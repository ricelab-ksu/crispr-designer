# ============================================================================
# tests/testthat/test-arm_resolution.R — homology arms come from the BUNDLE
# ============================================================================
# The defect this fixes: attL donor fragments were built from plasmid flanks
# instead of the genomic arms carried by the imported cas9_guide_bundle, which
# violates rules R1/R4.  These tests prove the bundle arms are used, the
# 250-1000 bp window is enforced, and the fallback still works.
# ============================================================================

library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "arm_resolution.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "plasmid_assembly.R"))
source(file.path(project_root, "R", "guide_bundle.R"))
source(file.path(project_root, "R", "gateway_attl.R"))
source(file.path(project_root, "R", "arm_resolution.R"))

mk <- function(n, ch = "A") strrep(ch, n)

make_guide_df <- function(a5, a3, gid = "g1") {
  data.frame(guide_id = gid, arm_5p_seq = a5, arm_3p_seq = a3,
             arm_5p_len = nchar(a5), arm_3p_len = nchar(a3),
             stringsAsFactors = FALSE)
}

# ---- trim helper ------------------------------------------------------------

test_that("arm_trim_to_window keeps the cut-site-adjacent end", {
  five <- paste0(strrep("A", 1400L), strrep("G", 100L))   # 3' end = cut-adjacent
  r5 <- arm_trim_to_window(five, "5p", 1000L)
  expect_true(r5$trimmed)
  expect_equal(r5$length_bp, 1000L)
  expect_equal(substr(r5$sequence, 901, 1000), strrep("G", 100L))  # kept

  three <- paste0(strrep("G", 100L), strrep("A", 1400L))  # 5' end = cut-adjacent
  r3 <- arm_trim_to_window(three, "3p", 1000L)
  expect_true(r3$trimmed)
  expect_equal(r3$length_bp, 1000L)
  expect_equal(substr(r3$sequence, 1, 100), strrep("G", 100L))     # kept
})

test_that("arm_trim_to_window leaves in-window arms untouched", {
  r <- arm_trim_to_window(mk(500), "5p", 1000L)
  expect_false(r$trimmed)
  expect_equal(r$length_bp, 500L)
})

# ---- bundle source (the fix) ------------------------------------------------

test_that("resolve_homology_arms uses the bundle arms, not plasmid flanks", {
  bundle <- make_guide_df(mk(600, "A"), mk(700, "C"))
  res <- resolve_homology_arms(bundle = bundle,
                               flank_5p = mk(500, "T"), flank_3p = mk(500, "T"))
  expect_equal(res$source, "bundle")
  expect_equal(res$left_arm, mk(600, "A"))
  expect_equal(res$right_arm, mk(700, "C"))
  expect_match(res$provenance, "cas9_guide_bundle")
  expect_equal(res$lengths[["arm_5p"]], 600L)
  expect_equal(res$lengths[["arm_3p"]], 700L)
})

test_that("resolve_homology_arms trims over-long bundle arms to 1000 bp", {
  bundle <- make_guide_df(mk(1500, "A"), mk(1300, "C"))
  res <- resolve_homology_arms(bundle = bundle)
  expect_true(res$trimmed)
  expect_equal(res$lengths[["arm_5p"]], 1000L)
  expect_equal(res$lengths[["arm_3p"]], 1000L)
  expect_equal(length(res$messages), 2L)
})

test_that("resolve_homology_arms rejects bundle arms below 250 bp", {
  bundle <- make_guide_df(mk(200), mk(600))
  expect_error(resolve_homology_arms(bundle = bundle), "validated minimum")
})

test_that("resolve_homology_arms selects the requested guide from a multi-row bundle", {
  g <- rbind(make_guide_df(mk(300, "A"), mk(300, "T"), "g1"),
             make_guide_df(mk(400, "C"), mk(400, "G"), "g2"))
  res <- resolve_homology_arms(bundle = g, guide_id = "g2")
  expect_equal(res$left_arm, mk(400, "C"))
  expect_equal(res$right_arm, mk(400, "G"))
  expect_match(res$provenance, "g2")
})

# ---- fallback source --------------------------------------------------------

test_that("resolve_homology_arms falls back to plasmid flanks when no bundle", {
  res <- resolve_homology_arms(flank_5p = mk(500, "T"), flank_3p = mk(500, "G"))
  expect_equal(res$source, "flanks")
  expect_match(res$provenance, "fallback")
  expect_equal(res$left_arm, mk(500, "T"))
})

test_that("resolve_homology_arms errors when neither source has arms", {
  expect_error(resolve_homology_arms(), "No homology arms available")
  expect_error(resolve_homology_arms(source = "bundle"), "no usable")
})

# ---- integration with fragment assembly ------------------------------------

test_that("bundle arms flow into pjat_assemble_fragments in the right order", {
  bundle <- make_guide_df(mk(600, "A"), mk(700, "C"))
  res <- resolve_homology_arms(bundle = bundle)
  attl <- pjat_attl_template_sequences()
  frags <- pjat_assemble_fragments(res$left_arm, res$right_arm,
                                   strrep("G", 20L),
                                   attl3 = attl$attl3, attl4 = attl$attl4,
                                   attl5 = attl$attl5, attl6 = attl$attl6,
                                   trna = attl$trna)
  expect_equal(frags$status, "READY")
  f1 <- frags$fragments$fragment1_attL3_arm_trna_attL4
  f2 <- frags$fragments$fragment2_attL5_arm_attL6
  # fragment 1 = attL3 + 5' arm + array + attL4
  expect_true(startsWith(f1, attl$attl3))
  expect_true(grepl(mk(600, "A"), f1, fixed = TRUE))
  # fragment 2 = attL5 + 3' arm + attL6
  expect_true(startsWith(f2, attl$attl5))
  expect_true(grepl(mk(700, "C"), f2, fixed = TRUE))
})

test_that("pjat_assemble_fragments now rejects arms above the 1000 bp ceiling", {
  attl <- pjat_attl_template_sequences()
  expect_error(
    pjat_assemble_fragments(mk(1200), mk(600), strrep("G", 20L),
                            attl3 = attl$attl3, attl4 = attl$attl4,
                            attl5 = attl$attl5, attl6 = attl$attl6,
                            trna = attl$trna),
    "validated maximum")
})
