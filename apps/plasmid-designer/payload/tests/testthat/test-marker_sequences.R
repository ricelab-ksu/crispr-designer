library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "marker_sequences.R"))) {
  project_root <- project_root
}

source(file.path(project_root, "R", "marker_sequences.R"))

CAT_YML  <- file.path(project_root, "donor_templates", "marker_cassettes.yml")
SEQ_DIR  <- file.path(project_root, "donor_templates", "marker_sequences")

test_that("catalog parses and exposes the expected columns", {
  cat <- pjd_marker_catalog()
  expect_gte(nrow(cat), 8L)
  expect_true(all(c("id", "label", "sequence_status", "sequence_file",
                    "length_bp", "md5", "sequence_provenance") %in% names(cat)))
  expect_true("none" %in% cat$id)
  expect_equal(cat$sequence_status[cat$id == "none"], "none")
})

test_that("every bundled marker has a valid FASTA matching the catalog", {
  cat <- pjd_marker_catalog()
  bundled <- cat[cat$sequence_status == "bundled", , drop = FALSE]
  skip_if(nrow(bundled) == 0L, "no bundled marker sequences yet")

  for (i in seq_len(nrow(bundled))) {
    id  <- bundled$id[i]
    fp  <- file.path(project_root, bundled$sequence_file[i])
    expect_true(file.exists(fp), info = id)
    seq <- pjd_read_fasta_sequence(fp)$sequence
    expect_match(seq, "^[ACGTN]+$", info = id)
    expect_equal(nchar(seq), bundled$length_bp[i], info = id)
    expect_equal(unname(tools::md5sum(fp)), bundled$md5[i], info = id)
    expect_false(grepl("unset", bundled$sequence_provenance[i]), info = id)
  }
})

test_that("the local pJAT-derived cassettes are the expected size", {
  cat <- pjd_marker_catalog()
  skip_if(!all(c("3xP3_dsRed", "ie1_dsRed") %in% cat$id))
  expect_equal(cat$length_bp[cat$id == "3xP3_dsRed"], 1861L)
  expect_equal(cat$length_bp[cat$id == "ie1_dsRed"],  2306L)
  expect_match(cat$sequence_provenance[cat$id == "3xP3_dsRed"], "pJAT30")
  expect_match(cat$sequence_provenance[cat$id == "ie1_dsRed"],  "pJAT11")
})

test_that("pjd_get_marker_sequence returns sequence + provenance", {
  skip_if_not(pjd_marker_sequence_available("neo"), "neo not bundled")
  out <- pjd_get_marker_sequence("neo")
  expect_equal(out$id, "neo")
  expect_equal(nchar(out$sequence), 795L)
  expect_equal(out$length_bp, 795L)
  expect_match(out$provenance, "NCBI")
  expect_equal(out$status, "bundled")
})

test_that("unknown or unavailability markers fail loudly", {
  expect_error(pjd_get_marker_sequence("not_a_marker"), "marker|not found|Unknown",
               ignore.case = TRUE)
  expect_false(pjd_marker_sequence_available("not_a_marker"))
  expect_false(pjd_marker_sequence_available("none"))
})

test_that("missing catalog file errors with an actionable message", {
  expect_error(pjd_marker_catalog(root = tempdir()), "not found")
})
