library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "gateway_attl.R"))) {
  project_root <- project_root
}

source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "gateway_attl.R"))

test_that("published attL templates load with the documented lengths", {
  t <- pjat_load_attl_templates()
  expect_equal(nchar(t$attl3), 96L)
  expect_equal(nchar(t$attl4), 97L)
  expect_equal(nchar(t$attl5), 97L)
  expect_equal(nchar(t$attl6), 97L)
  expect_equal(nchar(t$trna),  72L)
  expect_true(all(vapply(t[c("attl3","attl4","attl5","attl6","trna")],
                         function(s) grepl("^[ACGT]+$", s), logical(1))))
  expect_true(isTRUE(t$provenance$verified))
  expect_match(t$provenance$url, "protocols\\.io")
})

test_that("forward and inward orientations are reverse complements", {
  t <- pjat_load_attl_templates()
  rc <- function(x) paste(rev(strsplit(chartr("ACGT", "TGCA", x), "")[[1L]]),
                          collapse = "")
  expect_equal(rc(t$attl4), t$attl4_forward)
  expect_equal(rc(t$attl6), t$attl6_forward)
})

test_that("templates drive pjat_assemble_fragments to READY", {
  s <- pjat_attl_template_sequences()
  arm <- paste(rep("ACGT", 75), collapse = "")
  f <- pjat_assemble_fragments(arm, arm, "GAGTCTGACCTGACTGACTG",
                               attl3 = s$attl3, attl4 = s$attl4,
                               attl5 = s$attl5, attl6 = s$attl6, trna = s$trna)
  expect_equal(f$status, "READY")
  expect_true(isTRUE(f$checks$attl_sites_supplied))
  expect_false(grepl("[attL", paste(unlist(f$fragments), collapse = ""),
                     fixed = TRUE))
  expect_false(grepl("[tRNA]", paste(unlist(f$fragments), collapse = ""),
                     fixed = TRUE))
  expect_true(grepl("GAGTCTGACCTGACTGACTG", f$fragments$fragment1_attL3_arm_trna_attL4,
                    fixed = TRUE))
})

test_that("without templates the assembly degrades to ARM_ONLY", {
  arm <- paste(rep("ACGT", 75), collapse = "")
  f <- pjat_assemble_fragments(arm, arm, "GAGTCTGACCTGACTGACTG")
  expect_equal(f$status, "ARM_ONLY")
  expect_false(isTRUE(f$checks$attl_sites_supplied))
})

test_that("the pJAT backbone carries attR3-6 (124 bp), not attL", {
  gb = suppressWarnings(read_pjat_genbank(file.path(DEFAULT_PJAT_DIR, "pJAT74.gbk")))
  a <- pjat_extract_backbone_attr_sequences(gb)
  expect_setequal(a$att, c("attR3", "attR4", "attR5", "attR6"))
  expect_true(all(a$length_bp == 124L))
  expect_true(all(grepl("^[ACGT]+$", a$sequence)))
})

test_that("a corrupt template file is rejected with a clear message", {
  bad <- tempfile(fileext = ".yml")
  writeLines(c("attl3: \"ACGT\"", "attl4: \"ACGT\""), bad)
  expect_error(pjat_load_attl_templates(bad), "missing required key")
  expect_error(pjat_load_attl_templates(tempfile()), "not found")
})

test_that("the attL/attR explanation mentions both site families", {
  txt <- pjat_attl_explanation()
  expect_match(txt, "attL")
  expect_match(txt, "attR")
  expect_match(txt, "attB")
})
