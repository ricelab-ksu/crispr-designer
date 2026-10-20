# ============================================================================
# tests/testthat/test-onboarding.R — help content with worked examples
# ============================================================================

library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "onboarding.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "onboarding.R"))

test_that("the concatenation example shows literal arithmetic", {
  txt <- pjd_help_concat_example()
  expect_match(txt, "ATGC \\+ CGGG = ATGCCGGG", fixed = FALSE)
  expect_match(txt, "ATGCCGGGTTAA")
})

test_that("the orientation example shows the reverse complement", {
  txt <- pjd_help_orientation_example()
  expect_match(txt, "ATGC")
  expect_match(txt, "GCAT")
})

test_that("the window example states the length identity", {
  txt <- pjd_help_window_example()
  expect_match(txt, "final length = backbone length")
  expect_match(txt, "end = start - 1")
})

test_that("the arm example names both fragments and the 250-1000 window", {
  txt <- pjd_help_arm_example()
  expect_match(txt, "attL3")
  expect_match(txt, "attL6")
  expect_match(txt, "250-1000")
})

test_that("the glossary covers the core vocabulary", {
  gl <- pjd_help_glossary()
  expect_true(all(c("term", "meaning") %in% names(gl)))
  for (t in c("Spacer", "PAM", "Homology arm", "Payload", "Marker")) {
    expect_true(t %in% gl$term, info = t)
  }
})

test_that("the onboarding card renders for both apps", {
  for (app in c("plasmid", "crispr")) {
    html <- paste(as.character(pjd_onboarding_card(app)), collapse = "")
    expect_true(grepl("How this app works", html), info = app)
    expect_true(grepl("Glossary", html), info = app)
  }
})
