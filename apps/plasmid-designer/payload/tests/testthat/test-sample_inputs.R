library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "sample_inputs.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "guide_bundle.R"))
source(file.path(project_root, "R", "sample_inputs.R"))

test_that("pjd_load_sample_guide_bundle loads w_g0351 with 20-nt spacer and 1000-bp arms", {
  csv_path <- pjd_sample_guide_bundle_path(project_root)
  expect_true(file.exists(csv_path))

  bundle <- pjd_load_sample_guide_bundle(project_root)
  expect_type(bundle, "list")
  expect_true("guides" %in% names(bundle))
  expect_s3_class(bundle$guides, "data.frame")
  expect_gte(nrow(bundle$guides), 1L)

  expect_equal(bundle$guides$guide_id[1], "w_g0351")
  expect_equal(nchar(bundle$guides$spacer_seq), rep(20L, nrow(bundle$guides)))
  expect_true(all(grepl("^[ACGT]{20}$", bundle$guides$spacer_seq)))

  expect_equal(bundle$guides$arm_5p_len, rep(1000L, nrow(bundle$guides)))
  expect_equal(bundle$guides$arm_3p_len, rep(1000L, nrow(bundle$guides)))
  expect_equal(nchar(bundle$guides$arm_5p_seq), rep(1000L, nrow(bundle$guides)))
  expect_equal(nchar(bundle$guides$arm_3p_seq), rep(1000L, nrow(bundle$guides)))
})

test_that("pjd_load_sample_donor_fasta parses 3xP3_dsRed payload and matches header coords", {
  fasta_path <- pjd_sample_donor_fasta_path(project_root)
  expect_true(file.exists(fasta_path))

  raw_lines <- readLines(fasta_path, warn = FALSE)
  expect_equal(raw_lines[1], ">3xP3_dsRed|source=pJAT30|addgene=204294|coords=4749..6609")
  seq_lines <- raw_lines[-1]
  expect_true(all(nchar(seq_lines) <= 80L))

  donor <- pjd_load_sample_donor_fasta(project_root)
  expect_type(donor, "list")
  expect_named(donor, c("header", "sequence", "length_bp", "provenance"), ignore.order = TRUE)
  expect_equal(donor$header, "3xP3_dsRed|source=pJAT30|addgene=204294|coords=4749..6609")
  expect_true(grepl("^[ACGTN]+$", donor$sequence))

  expected_len <- 6609L - 4749L + 1L
  expect_equal(donor$length_bp, expected_len)
  expect_equal(nchar(donor$sequence), expected_len)
  expect_equal(donor$provenance$id, "3xP3_dsRed")
  expect_equal(donor$provenance$source, "pJAT30")
  expect_equal(donor$provenance$addgene, "204294")
  expect_equal(donor$provenance$coords, "4749..6609")
  expect_equal(donor$provenance$start, 4749L)
  expect_equal(donor$provenance$end, 6609L)
})

test_that("both sample loaders error clearly when files or roots are missing or invalid", {
  missing_root <- file.path(tempdir(), "nonexistent_plasmid_designer_root_xyz")
  expect_error(
    pjd_load_sample_guide_bundle(missing_root),
    "Sample guide bundle CSV not found"
  )
  expect_error(
    pjd_load_sample_donor_fasta(missing_root),
    "Sample donor FASTA not found"
  )
  expect_error(
    pjd_sample_guide_bundle_path(""),
    "non-empty character scalar"
  )

  # Verify donor loader catches coordinate/length mismatch and non-ACGTN chars
  tmp_root <- tempfile("sample_root_")
  dir.create(file.path(tmp_root, "uploads"), recursive = TRUE)
  bad_fasta <- pjd_sample_donor_fasta_path(tmp_root)

  writeLines(c(">3xP3_dsRed|source=pJAT30|addgene=204294|coords=4749..6609", "ACGTACGT"), bad_fasta)
  expect_error(pjd_load_sample_donor_fasta(tmp_root), "does not match header coords")

  writeLines(c(">3xP3_dsRed|source=pJAT30|addgene=204294|coords=1..8", "ACGTXYZ!"), bad_fasta)
  expect_error(pjd_load_sample_donor_fasta(tmp_root), "invalid non-ACGTN")
})
