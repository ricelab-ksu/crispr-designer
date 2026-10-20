library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "pjat_plasmids.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "guide_bundle.R"))

PJAT_TEST_DIR <- pjats_dir

make_valid_bundle <- function(tmpdir, n = 2, bad = FALSE) {
  l <- paste0(sample(c("A","C","G","T"), 250, replace = TRUE), collapse = "")
  r <- paste0(sample(c("A","C","G","T"), 260, replace = TRUE), collapse = "")
  spacer <- strrep("G", 20)
  rows <- lapply(seq_len(n), function(i) {
    data.frame(
      guide_id        = sprintf("ebony_g%04d", i),
      gene_symbol     = "ebony",
      species         = "melanogaster",
      assembly        = "Release_6_plus_ISO1_MT",
      chr             = "3R",
      accession       = "NC_004354.4",
      strand          = "+",
      protospacer_start = 21231000L + i*1000,
      protospacer_end   = 21231000L + i*1000 + 19L,
      pam_start       = 21231000L + i*1000 + 20L,
      pam_end         = 21231000L + i*1000 + 22L,
      cut_site        = 21231000L + i*1000 + 17L,
      spacer_seq      = if (bad && i == n) strrep("A", 21) else spacer,
      pam_seq         = "NGG",
      on_target_score = 72.5,
      cfd_specificity = 91.0,
      arm_5p_seq      = l,
      arm_5p_start    = 21230800L + i*1000,
      arm_5p_end      = 21231050L + i*1000,
      arm_5p_len      = nchar(l),
      arm_3p_seq      = r,
      arm_3p_start    = 21231100L + i*1000,
      arm_3p_end      = 21231360L + i*1000,
      arm_3p_len      = nchar(r),
      arm_3p_anchor   = "cut_site",
      pam_shielding_summary = "none",
      guide_bundle_schema   = GUIDE_BUNDLE_SCHEMA,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  csv <- file.path(tmpdir, "Cas9_ebony_melanogaster_guides_for_plasmid_designer.csv")
  utils::write.csv(df, csv, row.names = FALSE)
  json <- sub("\\.csv$", ".json", csv)
  sidecar <- list(
    schema = GUIDE_BUNDLE_SCHEMA,
    app_version = "1.0.0",
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
    species = "melanogaster",
    assembly = "Release_6_plus_ISO1_MT",
    gene_symbol = "ebony",
    n_guides = n,
    arm_policy = "250-1000 bp validated homology arms"
  )
  yaml::write_yaml(sidecar, json)
  list(csv = csv, json = json, n = n)
}

test_that("read_guide_bundle validates schema and required columns", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  out <- read_guide_bundle(b$csv, b$json, strict = TRUE)
  expect_equal(nrow(out$guides), b$n)
  expect_true(all(out$guides$spacer_seq == strrep("G", 20)))
  expect_equal(out$bundle_meta$n_guides, b$n)
  expect_false(is.null(out$sidecar))
  expect_match(out$provenance, "cas9_designer")
})

test_that("read_guide_bundle rejects missing required column", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  df <- readr::read_csv(b$csv, show_col_types = FALSE)
  df$guide_bundle_schema <- NULL
  csv2 <- file.path(td, "bad.csv")
  utils::write.csv(df, csv2, row.names = FALSE)
  expect_error(read_guide_bundle(csv2), "missing required columns")
})

test_that("read_guide_bundle rejects wrong schema version", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  df <- readr::read_csv(b$csv, show_col_types = FALSE)
  df$guide_bundle_schema <- "cas9_guide_bundle/0.9"
  csv2 <- file.path(td, "bad2.csv")
  utils::write.csv(df, csv2, row.names = FALSE)
  expect_error(read_guide_bundle(csv2), "Unsupported guide_bundle_schema")
})

test_that("read_guide_bundle rejects PAM-included spacer (rule R3)", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td, n = 1, bad = TRUE)
  expect_error(read_guide_bundle(b$csv, b$json, strict = TRUE), "valid 20-nt")
})

test_that("read_guide_bundle rejects arms below 250 bp (rule R2)", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  df <- readr::read_csv(b$csv, show_col_types = FALSE)
  df$arm_5p_seq <- strrep("A", 200)
  df$arm_5p_len <- 200L
  csv2 <- file.path(td, "bad3.csv")
  utils::write.csv(df, csv2, row.names = FALSE)
  expect_error(read_guide_bundle(csv2), "arms < 250 bp")
})

test_that("guide_bundle_check_compatible reports per-guide backbone compatibility", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td, n = 1)
  out <- read_guide_bundle(b$csv, b$json, strict = TRUE)
  chk <- guide_bundle_check_compatible(out)
  expect_equal(nrow(chk), 1L)
  expect_true(chk$compatible[1])
  expect_match(chk$compatible_backbones[1], "pJAT77|pJAT74|pJAT32")
})

test_that("read_guide_bundle works without a sidecar", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  file.remove(b$json)
  out <- read_guide_bundle(b$csv, json = NULL, strict = TRUE)
  expect_null(out$sidecar)
  expect_match(out$provenance, "user-supplied CSV")
})

test_that("read_guide_bundle imports an extension-less upload path", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  no_ext <- file.path(td, "uploaded_tmp")
  file.copy(b$csv, no_ext)
  out <- read_guide_bundle(no_ext, json = NULL, strict = TRUE)
  expect_equal(nrow(out$guides), b$n)
  expect_null(out$sidecar)
})

test_that("read_guide_bundle survives a malformed sidecar", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  writeLines("this is not a sidecar document", b$json)
  out <- suppressWarnings(read_guide_bundle(b$csv, b$json, strict = TRUE))
  expect_equal(nrow(out$guides), b$n)
  expect_null(out$sidecar)
  expect_match(out$provenance, "user-supplied CSV")
})

test_that("read_guide_bundle accepts a JSON sidecar (producer format)", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  td <- tempfile("bundle_"); dir.create(td)
  b <- make_valid_bundle(td)
  sc <- yaml::read_yaml(b$json)
  writeLines(jsonlite::toJSON(sc, pretty = TRUE, auto_unbox = TRUE), b$json)
  out <- read_guide_bundle(b$csv, b$json, strict = TRUE)
  expect_false(is.null(out$sidecar))
  expect_equal(out$sidecar$schema, GUIDE_BUNDLE_SCHEMA)
  expect_match(out$provenance, "cas9_designer")
})
