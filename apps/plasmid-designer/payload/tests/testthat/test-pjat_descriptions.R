# ============================================================================
# tests/testthat/test-pjat_descriptions.R — curated plasmid descriptions
# ============================================================================
# Every pJAT plasmid must carry a curated 2-4 sentence description so users can
# tell what a plasmid is for without leaving the app.  Descriptions are
# additive: attR coordinates must be untouched.
# ============================================================================

library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "pjat_plasmids.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "pjat_plasmids.R"))

test_that("every pJAT plasmid has a substantial curated description", {
  yml <- file.path(project_root, "donor_templates", "pjat_plasmids.yml")
  skip_if_not(file.exists(yml), "pJAT catalog not available")
  cat_df <- pjat_load_catalog(yml = yml)
  expect_true("description" %in% names(cat_df))
  expect_equal(nrow(cat_df), 24L)
  for (i in seq_len(nrow(cat_df))) {
    d <- cat_df$description[i]
    expect_false(is.na(d), info = cat_df$name[i])
    expect_gte(nchar(d), 120L)          # a real explanation, not a fragment
  }
})

test_that("descriptions are distinct per plasmid (no copy-paste)", {
  yml <- file.path(project_root, "donor_templates", "pjat_plasmids.yml")
  skip_if_not(file.exists(yml), "pJAT catalog not available")
  cat_df <- pjat_load_catalog(yml = yml)
  expect_equal(anyDuplicated(cat_df$description), 0L)
})

test_that("the pJAT adapter carries descriptions through to the registry", {
  root <- project_root
  yml <- file.path(root, "donor_templates", "pjat_plasmids.yml")
  ref_dir <- pjats_dir
  skip_if_not(file.exists(yml) && dir.exists(ref_dir), "pJAT assets not mounted")
  source(file.path(root, "R", "plasmid_assembly.R"))
  source(file.path(root, "R", "backbone_registry.R"))
  reg <- backbone_pjat_adapter(yml = yml, ref_dir = ref_dir)
  e <- reg[["pJAT74"]]
  expect_true(nzchar(as.character(e$provenance$purpose)))
})
