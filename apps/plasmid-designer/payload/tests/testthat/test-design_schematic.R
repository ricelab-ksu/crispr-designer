# ============================================================================
# tests/testthat/test-design_schematic.R — live construct-design schematic
# ============================================================================

library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "design_schematic.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "design_schematic.R"))

render_svg <- function(x) paste(as.character(x), collapse = "")

test_that("pjd_design_revcomp reverse-complements correctly", {
  expect_equal(.pjd_design_revcomp("ATGC"), "GCAT")
  expect_equal(.pjd_design_revcomp("AAAACCCCGGGGTTTT"), "AAAACCCCGGGGTTTT")
  expect_equal(.pjd_design_revcomp(""), "")
  expect_equal(.pjd_design_revcomp("AN"), "NT")
})

test_that("the schematic renders backbone, slot, attR and orientation blocks", {
  html <- render_svg(pjd_design_schematic(
    backbone_name = "pJAT74", spacer = strrep("G", 20L),
    payload_seq = strrep("A", 500L), payload_label = "3xP3-DsRed",
    orient = "+", slot_label = "linker window (attR4-attR5)",
    grna_site = "1 x U6-3 promoter"))
  expect_true(grepl("pJAT74", html))
  expect_true(grepl("3xP3-DsRed", html))
  expect_true(grepl("attR4", html))
  expect_true(grepl("attR5", html))
  expect_true(grepl("1 x U6-3 promoter", html))
  expect_true(grepl("500 bp payload", html))
})

test_that("plus orientation says no reverse complement; minus says reversed", {
  plus <- render_svg(pjd_design_schematic(payload_seq = "AAAACCCCGGGGTTTT",
                                          orient = "+"))
  expect_true(grepl("no reverse complement", plus))

  minus <- render_svg(pjd_design_schematic(payload_seq = "AAAACCCCGGGGTTTT",
                                           orient = "-"))
  expect_true(grepl("reverse complement", minus))
  # AAAA CCCC GGGG TTTT revcomp is itself, so use an asymmetric payload too.
  asym <- render_svg(pjd_design_schematic(payload_seq = "AAAAAAAAAACCCCCCCCCC",
                                          orient = "-"))
  expect_true(grepl("reverse complement", asym))
})

test_that("orientation flip changes the rendered readout text", {
  seq <- "ACGTACGTACGTACGTACGT"
  plus <- render_svg(pjd_design_schematic(payload_seq = seq, orient = "+"))
  minus <- render_svg(pjd_design_schematic(payload_seq = seq, orient = "-"))
  expect_false(identical(plus, minus))
})

test_that("assembled flag adds the assembled note", {
  html <- render_svg(pjd_design_schematic(payload_seq = "ACGT", assembled = TRUE))
  expect_true(grepl("Assembled", html))
})

test_that("pjd_grna_site_label describes the gRNA insertion site", {
  expect_match(pjd_grna_site_label(0L), "none")
  expect_match(pjd_grna_site_label(0L), "pJAT77")
  expect_match(pjd_grna_site_label(1L), "1 .* U6-3 promoter")
  expect_match(pjd_grna_site_label(2L), "2 .* U6-3 promoters")
  expect_match(pjd_grna_site_label(NA_integer_), "unknown")
  expect_match(pjd_grna_site_label(1L, has_trna = TRUE), "tRNA array")
})

test_that("the schematic tolerates a missing payload and spacer", {
  html <- render_svg(pjd_design_schematic())
  expect_true(nzchar(html))
  expect_true(grepl("0 bp payload", html))
})
