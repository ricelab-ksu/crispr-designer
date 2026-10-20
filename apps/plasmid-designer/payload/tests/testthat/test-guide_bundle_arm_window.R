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

# ---- v1.1 schema: 250-1000 bp arm window (rule R2) -------------------------

test_that("read_guide_bundle accepts both 1.0 (legacy) and 1.1 schemas", {
  expect_true("cas9_guide_bundle/1.0" %in% GUIDE_BUNDLE_ACCEPTED_SCHEMAS)
  expect_true("cas9_guide_bundle/1.1" %in% GUIDE_BUNDLE_ACCEPTED_SCHEMAS)
  expect_equal(GUIDE_BUNDLE_ARM_MIN_BP, 250L)
  expect_equal(GUIDE_BUNDLE_ARM_MAX_BP, 1000L)
})

test_that("read_guide_bundle rejects arms longer than 1000 bp (R2 ceiling)", {
  mk <- function(n) substr(paste(rep("ACGT", length.out = ceiling(n / 4) + 1),
                                 collapse = ""), 1, n)
  # Build a minimal valid bundle CSV with a 1500 bp arm.
  cols <- GUIDE_BUNDLE_REQUIRED_COLS
  row <- as.list(rep("x", length(cols))); names(row) <- cols
  row$guide_id <- "g1"; row$gene_symbol <- "w"; row$species <- "melanogaster"
  row$assembly <- "GCF_test"; row$chr <- "3R"; row$accession <- "NC_1"
  row$strand <- "+"; row$protospacer_start <- 1L; row$protospacer_end <- 20L
  row$pam_start <- 21L; row$pam_end <- 23L; row$cut_site <- 17L
  row$spacer_seq <- "GAGTCTGACCTGACTGACTG"; row$pam_seq <- "TGG"
  row$on_target_score <- 50; row$cfd_specificity <- 90
  row$arm_5p_seq <- mk(1500); row$arm_5p_start <- 1L; row$arm_5p_end <- 1500L
  row$arm_5p_len <- 1500L
  row$arm_3p_seq <- mk(1500); row$arm_3p_start <- 2000L; row$arm_3p_end <- 3500L
  row$arm_3p_len <- 1500L
  row$arm_3p_anchor <- "cut_site"; row$pam_shielding_summary <- "none"
  row$guide_bundle_schema <- GUIDE_BUNDLE_SCHEMA_V11
  df <- as.data.frame(row, stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".csv")
  utils::write.csv(df, f, row.names = FALSE)
  expect_error(read_guide_bundle(f), "validated maximum")
  unlink(f)
})

test_that("a 1.1 bundle with in-window arms imports and reports the policy", {
  mk <- function(n) substr(paste(rep("ACGT", length.out = ceiling(n / 4) + 1),
                                 collapse = ""), 1, n)
  cols <- GUIDE_BUNDLE_REQUIRED_COLS
  row <- as.list(rep("x", length(cols))); names(row) <- cols
  row$guide_id <- "g1"; row$gene_symbol <- "w"; row$species <- "melanogaster"
  row$assembly <- "GCF_test"; row$chr <- "3R"; row$accession <- "NC_1"
  row$strand <- "+"; row$protospacer_start <- 1L; row$protospacer_end <- 20L
  row$pam_start <- 21L; row$pam_end <- 23L; row$cut_site <- 17L
  row$spacer_seq <- "GAGTCTGACCTGACTGACTG"; row$pam_seq <- "TGG"
  row$on_target_score <- 50; row$cfd_specificity <- 90
  row$arm_5p_seq <- mk(600); row$arm_5p_start <- 1L; row$arm_5p_end <- 600L
  row$arm_5p_len <- 600L
  row$arm_3p_seq <- mk(700); row$arm_3p_start <- 2000L; row$arm_3p_end <- 2700L
  row$arm_3p_len <- 700L
  row$arm_3p_anchor <- "cut_site"; row$pam_shielding_summary <- "none"
  row$guide_bundle_schema <- GUIDE_BUNDLE_SCHEMA_V11
  df <- as.data.frame(row, stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".csv")
  utils::write.csv(df, f, row.names = FALSE)
  b <- read_guide_bundle(f)
  expect_equal(b$bundle_meta$schema, GUIDE_BUNDLE_SCHEMA_V11)
  expect_equal(b$bundle_meta$arm_min_bp, 250L)
  expect_equal(b$bundle_meta$arm_max_bp, 1000L)
  expect_equal(b$guides$arm_5p_len, 600L)
  unlink(f)
})
