library(testthat)
suppressPackageStartupMessages(library(shiny))

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "plasmid_map.R"))) {
  project_root <- project_root
}

source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "plasmid_map.R"))

PJAT_DIR <- pjats_dir

test_that("GenBank records normalise into a plottable map", {
  skip_if_not(file.exists(file.path(PJAT_DIR, "pJAT74.gbk")), "pJAT data missing")
  gb <- suppressWarnings(read_pjat_genbank(file.path(PJAT_DIR, "pJAT74.gbk")))
  rec <- pjd_map_record_from_genbank(gb)

  expect_equal(rec$size_bp, 7986L)
  expect_gt(nrow(rec$features), 10L)
  expect_true(all(c("type", "start", "end", "strand", "label", "class",
                    "color", "row") %in% names(rec$features)))
  expect_true(all(rec$features$start >= 1L))
  expect_true(all(rec$features$end <= rec$size_bp))
  expect_true(all(rec$features$end >= rec$features$start))
  expect_true(all(rec$features$color %in% PJD_MAP_COLORS))
  expect_true(any(grepl("attR", rec$features$label)))
  expect_true(all(rec$features$class[grepl("^attR", rec$features$label)] == "attR"))
})

test_that("origin-spanning features are split into two segments", {
  seq <- paste(rep("ACGT", 250), collapse = "")          # 1000 bp
  gb <- list(sequence = seq, size_bp = 1000L, name = "synthetic",
             features = data.frame(
               key = c("misc_feature", "CDS"),
               start = c(900L, 100L),
               end = c(100L, 200L),
               strand = c("+", "+"),
               label = c("wraps", "normal"),
               product = c("", ""),
               stringsAsFactors = FALSE))
  rec <- pjd_map_record_from_genbank(gb)
  wraps <- rec$features[rec$features$label == "wraps", , drop = FALSE]
  expect_equal(nrow(wraps), 2L)
  expect_true(any(wraps$start == 900L & wraps$end == 1000L))
  expect_true(any(wraps$start == 1L & wraps$end == 100L))
  expect_true(all(wraps$end >= wraps$start))
})

test_that("assembled constructs render as map records", {
  seq <- paste(rep("ACGT", 500), collapse = "")
  asm <- list(sequence = seq,
              insert_start = 101L, insert_end = 300L,
              orient = "+",
              features = data.frame(type = "insert", start = 101L, end = 300L,
                                    strand = "+", label = "user insert",
                                    stringsAsFactors = FALSE))
  rec <- pjd_map_record_from_assembly(asm, backbone = "pJAT74")
  expect_equal(rec$size_bp, 2000L)
  expect_equal(rec$title, "pJAT74 assembly")
  expect_true(any(rec$features$label == "user insert"))
  expect_equal(unname(rec$insert), c(101L, 300L))
  expect_equal(rec$features$class[rec$features$label == "user insert"], "insert")
})

test_that("circular and linear SVG output contains the expected structure", {
  skip_if_not(file.exists(file.path(PJAT_DIR, "pJAT74.gbk")), "pJAT data missing")
  rec <- pjd_map_record_from_genbank(
    suppressWarnings(read_pjat_genbank(file.path(PJAT_DIR, "pJAT74.gbk"))))

  circ <- as.character(.pjd_interactive_map_ui(rec))
  expect_match(circ, "<svg", fixed = TRUE)
  expect_match(circ, 'data-pjd-map="circular"', fixed = TRUE)
  expect_match(circ, 'data-pjd-map="feature"', fixed = TRUE)
  expect_match(circ, "attR3", fixed = TRUE)
  expect_match(circ, "pjd-map-tooltip", fixed = TRUE)
  expect_match(circ, "pjd-map-legend", fixed = TRUE)
  expect_match(circ, "<title>", fixed = TRUE)
  expect_match(circ, "kb", fixed = TRUE)

  lin <- as.character(.pjd_interactive_map_ui(rec, mode = "linear"))
  expect_match(lin, 'data-pjd-map="linear"', fixed = TRUE)
  expect_match(lin, "<rect", fixed = TRUE)
  expect_match(lin, "5' (+) strand", fixed = TRUE)
})

test_that("an empty record produces the empty-state card", {
  h <- as.character(.pjd_interactive_map_ui(NULL))
  expect_match(h, "No map data")
  expect_false(grepl("<svg", h, fixed = TRUE))
})

test_that("the interaction assets are present and non-trivial", {
  js  <- file.path(project_root, "www", "pjd_assets", "pjd_map.js")
  css <- file.path(project_root, "www", "pjd_assets", "pjd_map.css")
  expect_true(file.exists(js))
  expect_true(file.exists(css))
  expect_gt(file.info(js)$size, 1000)
  expect_gt(file.info(css)$size, 500)
  expect_match(paste(readLines(js, warn = FALSE), collapse = "\n"),
               "pjd_map_selected_feature", fixed = TRUE)
})

test_that("feature classes are assigned by biological role", {
  expect_equal(pjd_map_feature_class("protein_bind", "attR4"), "attR")
  expect_equal(pjd_map_feature_class("CDS", "AmpR"), "resistance")
  expect_equal(pjd_map_feature_class("promoter", "dU6-3 promoter"), "promoter")
  expect_equal(pjd_map_feature_class("CDS", "EGFP"), "cds")
  expect_equal(pjd_map_feature_class("rep_origin", "ori"), "origin")
  expect_equal(pjd_map_feature_class("CDS", "ccdB"), "negative")
})
