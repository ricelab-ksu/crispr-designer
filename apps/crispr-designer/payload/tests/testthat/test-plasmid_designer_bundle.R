# ============================================================================
# tests/testthat/test-plasmid_designer_bundle.R — cas9_guide_bundle/1.1 producer
# ============================================================================

library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "plasmid_designer_bundle.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "plasmid_designer_bundle.R"))

# ---- fixtures ---------------------------------------------------------------

make_candidates <- function(gid = "ebony_g0001", n = 3L) {
  data.frame(
    guide_id          = c(gid, sprintf("ebony_g%04d", seq(2L, n))),
    chr               = "3R",
    accession         = "NC_004354.4",
    strand            = "+",
    protospacer_start = c(21231000L, 21232000L, 21233000L)[seq_len(n)],
    protospacer_end   = c(21231019L, 21232019L, 21233019L)[seq_len(n)],
    pam_start         = c(21231020L, 21232020L, 21233020L)[seq_len(n)],
    pam_end           = c(21231022L, 21232022L, 21233022L)[seq_len(n)],
    cut_site          = c(21231017L, 21232017L, 21233017L)[seq_len(n)],
    spacer_seq        = strrep("G", 20),
    pam_seq           = "NGG",
    on_target_score   = 72.5,
    cfd_specificity_score = 91.0,
    stringsAsFactors  = FALSE
  )
}

make_arm_design <- function(gid = "ebony_g0001",
                            a5 = strrep("A", 1000L), a3 = strrep("C", 1000L)) {
  list(
    guide_id      = gid,
    species       = "melanogaster",
    chr           = "3R",
    accession     = "NC_004354.4",
    strand        = "+",
    spacer_seq    = strrep("G", 20),
    pam_seq       = "NGG",
    arm_5p_seq    = a5,
    arm_3p_seq    = a3,
    arm_5p_start  = 21230000L,
    arm_5p_end    = 21230999L,
    arm_3p_start  = 21231023L,
    arm_3p_end    = 21232022L,
    arm_3p_anchor = "cut_site",
    pam_shielding_summary = "3' arm anchored at cut site"
  )
}

make_meta <- function() {
  list(species = "melanogaster", query = "ebony",
       assembly_accession = "GCF_000001215.4")
}

# ---- unit tests -------------------------------------------------------------

test_that("build_plasmid_designer_bundle emits all 27 required columns", {
  b <- build_plasmid_designer_bundle(make_candidates(), make_arm_design(), make_meta())
  expect_equal(nrow(b$guides), 1L)
  expect_setequal(names(b$guides), .gbu_cols)
  expect_true(all(.gbu_cols %in% names(b$guides)))
})

test_that("schema and sidecar keys match the consumer contract", {
  b <- build_plasmid_designer_bundle(make_candidates(), make_arm_design(), make_meta())
  expect_true(all(b$guides$guide_bundle_schema == "cas9_guide_bundle/1.1"))
  expect_true(all(c("schema", "app_version", "created_utc", "species", "assembly",
                    "gene_symbol", "n_guides", "arm_policy") %in% names(b$sidecar)))
  # v1.1 additions: the protocol-validated arm window and trim provenance.
  expect_equal(b$sidecar$arm_policy, "250-1000 bp validated homology arms")
  expect_equal(b$sidecar$arm_min_bp, 250L)
  expect_equal(b$sidecar$arm_max_bp, 1000L)
  expect_false(b$sidecar$arms_trimmed)
  expect_equal(b$sidecar$n_guides, 1L)
  expect_equal(b$sidecar$species, "melanogaster")
  expect_equal(b$sidecar$gene_symbol, "ebony")
})

test_that("cfd_specificity_score is renamed to cfd_specificity", {
  b <- build_plasmid_designer_bundle(make_candidates(), make_arm_design(), make_meta())
  expect_true("cfd_specificity" %in% names(b$guides))
  expect_false("cfd_specificity_score" %in% names(b$guides))
  expect_equal(b$guides$cfd_specificity, 91.0)
})

test_that("only the guide with arms is exported", {
  b <- build_plasmid_designer_bundle(make_candidates(n = 3L), make_arm_design(), make_meta())
  expect_equal(b$n_guides, 1L)
  expect_equal(b$guides$guide_id, "ebony_g0001")
})

test_that("rejects non-20-nt spacer and arms < 250 bp", {
  ad <- make_arm_design(); ad$spacer_seq <- strrep("A", 21L)
  expect_error(build_plasmid_designer_bundle(make_candidates(), ad, make_meta()),
               "20 nt")
  ad2 <- make_arm_design(a5 = strrep("A", 200L))
  expect_error(build_plasmid_designer_bundle(make_candidates(), ad2, make_meta()),
               "250 bp")
  ad3 <- make_arm_design(); ad3$arm_3p_seq <- ""
  expect_error(build_plasmid_designer_bundle(make_candidates(), ad3, make_meta()),
               "empty")
})

test_that("rejects missing/unknown guide_id and empty grna_df", {
  expect_error(build_plasmid_designer_bundle(make_candidates(), NULL, make_meta()),
               "arm_design is required")
  ad <- make_arm_design(); ad$guide_id <- "nope_g9999"
  expect_error(build_plasmid_designer_bundle(make_candidates(), ad, make_meta()),
               "not found")
  expect_error(build_plasmid_designer_bundle(make_candidates(n = 3L)[0, ], make_arm_design(), make_meta()),
               "non-empty")
})

test_that("write_plasmid_designer_bundle writes CSV + JSON sidecar", {
  td <- tempfile("pdb_"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  f <- file.path(td, "Cas9_ebony_melanogaster_guides_for_plasmid_designer.csv")
  out <- write_plasmid_designer_bundle(make_candidates(), make_arm_design(), make_meta(), file = f)
  expect_true(file.exists(out$csv))
  expect_true(file.exists(out$json))
  sc <- jsonlite::fromJSON(out$json, simplifyVector = FALSE)
  expect_equal(sc$schema, "cas9_guide_bundle/1.1")
  expect_equal(sc$arm_min_bp, 250L)
  expect_equal(sc$arm_max_bp, 1000L)
  expect_equal(sc$gene_symbol, "ebony")
})

test_that("arms longer than 1000 bp are trimmed to the R2 ceiling on export", {
  ad <- make_arm_design(a5 = strrep("A", 1500L), a3 = strrep("C", 1400L))
  b <- build_plasmid_designer_bundle(make_candidates(), ad, make_meta())
  expect_equal(b$guides$arm_5p_len, 1000L)
  expect_equal(b$guides$arm_3p_len, 1000L)
  expect_true(b$sidecar$arms_trimmed)
  # cut-site-adjacent bases must be preserved: 5' arm keeps its 3' end,
  # 3' arm keeps its 5' end (both abut the cut site).
  expect_equal(substr(b$guides$arm_5p_seq, 1, 1), "A")
  expect_equal(substr(b$guides$arm_3p_seq, 1, 1), "C")
})


# ---- round-trip through the real Plasmid Designer validator -----------------

test_that("exported bundle round-trips through read_guide_bundle(strict = TRUE)", {
  pd_root <- plasmid_designer_root
  skip_if_not(file.exists(file.path(pd_root, "R", "guide_bundle.R")),
              "Plasmid_Designer not present")
  source(file.path(pd_root, "R", "pjat_plasmids.R"))
  source(file.path(pd_root, "R", "guide_bundle.R"))

  td <- tempfile("pdb_rt_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  f <- file.path(td, "Cas9_ebony_melanogaster_guides_for_plasmid_designer.csv")
  write_plasmid_designer_bundle(make_candidates(), make_arm_design(), make_meta(), file = f)

  res <- read_guide_bundle(f, strict = TRUE)
  expect_equal(nrow(res$guides), 1L)
  expect_equal(res$guides$guide_id, "ebony_g0001")
  expect_equal(res$guides$spacer_seq, strrep("G", 20))
  expect_true(all(nchar(res$guides$arm_5p_seq) >= 250L))
  expect_true(all(nchar(res$guides$arm_3p_seq) >= 250L))
  expect_equal(res$bundle_meta$n_guides, 1L)
  expect_match(res$provenance, "cas9_designer")
})
