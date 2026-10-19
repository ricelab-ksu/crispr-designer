library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "pjat_plasmids.R"))) {
  project_root <- "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design"
}
source(file.path(project_root, "R", "pjat_plasmids.R"))

PJAT_TEST_DIR <- "/Volumes/Backup_Plus/projects/bioinformatics/reference_files/janelia_atalanta_plasmids"
HAS_DATA <- dir.exists(PJAT_TEST_DIR)

test_that("read_pjat_genbank parses all plasmids with size matching sequence", {
  skip_if_not(HAS_DATA, "pJAT collection not present")
  files <- list.files(PJAT_TEST_DIR, pattern = "^pJAT[0-9]+", full.names = TRUE)
  files <- files[!grepl("(^|/)\\._", files)]
  expect_gte(length(files), 24L)
  for (f in files) {
    gb <- suppressWarnings(read_pjat_genbank(f))
    expect_false(is.na(gb$size_bp), info = basename(f))
    expect_equal(nchar(gb$sequence), gb$size_bp, info = basename(f))
    expect_gt(nrow(gb$features), 0L)
  }
})

test_that("hazard H1: Geneious HTML-in-notes record parses as GenBank", {
  skip_if_not(HAS_DATA, "pJAT collection not present")
  f <- file.path(PJAT_TEST_DIR, "pJAT55.gbff")
  skip_if_not(file.exists(f), "pJAT55.gbff not present")
  gb <- read_pjat_genbank(f)
  expect_equal(gb$size_bp, 13951L)
  expect_equal(nchar(gb$sequence), 13951L)
  att <- find_pjat_recombination_sites(gb)
  expect_true(all(c("attR3","attR4","attR5","attR6") %in% att$att))
})

test_that("hazard H2: protein_bind+label and custom-key att annotations both found", {
  skip_if_not(HAS_DATA, "pJAT collection not present")
  a <- find_pjat_recombination_sites(read_pjat_genbank(
    file.path(PJAT_TEST_DIR, "pJAT74.gbk")))
  expect_equal(sort(a$att), c("attR3","attR4","attR5","attR6"))
  expect_true(all(a$key == "protein_bind"))
  expect_true(all(a$width == 124L))

  b <- find_pjat_recombination_sites(read_pjat_genbank(
    file.path(PJAT_TEST_DIR, "pJAT55.gbff")))
  expect_equal(sort(b$att), c("attR3","attR4","attR5","attR6"))
  expect_true(all(b$width == 124L))
})

test_that("manifest validates every plasmid and reports attR coordinates", {
  skip_if_not(HAS_DATA, "pJAT collection not present")
  man <- pjat_manifest(PJAT_TEST_DIR)
  expect_gte(nrow(man), 24L)
  expect_true(all(man$status == "OK"))
  p74 <- man[man$plasmid == "pJAT74", ]
  expect_equal(p74$attR3_start, 3980L)
  expect_equal(p74$attR4_end,   5967L)
  expect_equal(man$n_u6[man$plasmid == "pJAT77"], 0L)
  expect_true(all(man$n_u6[man$plasmid %in%
    c("pJAT30","pJAT31","pJAT32","pJAT39")] == 2L))
  expect_equal(man$n_u6[man$plasmid == "pJAT74"], 1L)
  expect_true(man$has_attp[man$plasmid == "pJAT11"])
  expect_true(man$has_attp[man$plasmid == "pJAT7"])
  expect_false(man$has_attp[man$plasmid == "pJAT74"])
})

test_that("validate_pjat_collection flags a missing att site", {
  m <- data.frame(plasmid = c("good","bad"),
                  attR3_start = c(10L,10L), attR4_start = c(20L,NA),
                  attR5_start = c(30L,30L), attR6_start = c(40L,40L))
  out <- validate_pjat_collection(m)
  expect_equal(out$status[1], "OK")
  expect_match(out$status[2], "MISSING:attR4")
})

test_that("synthesis QC flags >10 bp homopolymers and extreme GC", {
  set.seed(1)
  ok_seq <- paste(sample(c("A","C","G","T"), 400, replace = TRUE), collapse = "")
  qc <- pjat_synthesis_qc(ok_seq)
  expect_true(qc$ok)
  expect_lte(qc$max_homopolymer, 10L)

  bad <- paste0("ACGT", strrep("A", 12), "ACGT")
  qc2 <- pjat_synthesis_qc(bad)
  expect_false(qc2$ok)
  expect_equal(qc2$max_homopolymer, 13L)
  expect_match(qc2$warnings[1], "Homopolymer run of 13")
})

test_that("spacer guard enforces 20 nt and rejects PAM inclusion (rule R3)", {
  expect_true(pjat_validate_spacer(strrep("A", 20)))
  expect_error(pjat_validate_spacer(strrep("A", 21)), "exactly 20 nt")
  expect_error(pjat_validate_spacer(strrep("A", 19)), "exactly 20 nt")
  # spacer + PAM concatenated must be rejected
  expect_error(pjat_validate_spacer(paste0(strrep("A",20), "AGG")), "exactly 20 nt")
})

test_that("fragment assembly returns READY only when attL sites are supplied", {
  set.seed(2)
  L <- paste(sample(c("A","C","G","T"), 250, replace = TRUE), collapse = "")
  R <- paste(sample(c("A","C","G","T"), 260, replace = TRUE), collapse = "")

  res <- pjat_assemble_fragments(L, R, strrep("G", 20))
  expect_equal(res$status, "ARM_ONLY")
  expect_false(res$checks$attl_sites_supplied)
  expect_true(res$checks$spacer_is_20nt)
  expect_true(res$checks$left_arm_ge_250)
  expect_true(res$checks$right_arm_ge_250)

  res2 <- pjat_assemble_fragments(L, R, strrep("G", 20),
    attl3 = strrep("C",10), attl4 = strrep("T",10),
    attl5 = strrep("A",10), attl6 = strrep("G",10))
  expect_equal(res2$status, "READY")
  # fragment 1 must begin with attL3 and contain the left arm
  expect_true(startsWith(res2$fragments$fragment1_attL3_arm_trna_attL4,
                         strrep("C",10)))
})

test_that("catalog loads 24 plasmids and matches attR coordinates", {
  skip_if_not(file.exists(file.path(project_root, "donor_templates", "pjat_plasmids.yml")),
              "pjat_plasmids.yml not present")
  cat <- pjat_load_catalog()
  expect_equal(nrow(cat), 24L)
  p74 <- cat[cat$name == "pJAT74", ]
  expect_equal(p74$attR3_start, 3980L)
  expect_equal(p74$attR3_strand, "+")
  expect_equal(p74$attR4_strand, "-")
  expect_equal(cat$u6_promoters[cat$name == "pJAT77"], 0L)
  expect_true(cat$has_attp[cat$name == "pJAT11"])
})

test_that("recommendations depend only on gRNA count (U6 promoters)", {
  skip_if_not(file.exists(file.path(project_root, "donor_templates", "pjat_plasmids.yml")),
              "pjat_plasmids.yml not present")
  expect_equal(pjat_recommend_backbones(0)$name, "pJAT77")
  expect_equal(pjat_recommend_backbones(2)$name, "pJAT32")
  r1 <- pjat_recommend_backbones(1)
  expect_true(all(r1$u6_promoters >= 1L))
  expect_true("pJAT74" %in% r1$name)
})

test_that("compatibility guard errors when U6 count is insufficient", {
  skip_if_not(file.exists(file.path(project_root, "donor_templates", "pjat_plasmids.yml")),
              "pjat_plasmids.yml not present")
  expect_error(pjat_check_compatible("pJAT74", 2), "pJAT32")
  expect_error(pjat_check_compatible("pJAT77", 1), "Use one of")
  expect_true(pjat_check_compatible("pJAT32", 2))
  expect_true(pjat_check_compatible("pJAT74", 1))
  expect_error(pjat_check_compatible("pJAT99", 1), "Unknown")
  expect_error(pjat_recommend_backbones(3), "0, 1 or 2")
})

test_that("fragment assembly enforces the 250 bp validated minimum (rule R2)", {
  set.seed(3)
  short <- paste(sample(c("A","C","G","T"), 200, replace = TRUE), collapse = "")
  long  <- paste(sample(c("A","C","G","T"), 300, replace = TRUE), collapse = "")
  expect_error(pjat_assemble_fragments(short, long, strrep("G",20)), ">=250 bp")
  expect_error(pjat_assemble_fragments(long, short, strrep("G",20)), ">=250 bp")
})

test_that("reverse complement is self-inverse", {
  expect_equal(.pjat_revcomp("ACGT"), "ACGT")
  expect_equal(.pjat_revcomp("AACCGGTT"), "AACCGGTT")
  expect_equal(.pjat_revcomp("ATGC"), "GCAT")
})