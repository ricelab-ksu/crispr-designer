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
source(file.path(project_root, "R", "addgene_maps.R"))
source(file.path(project_root, "R", "plasmid_map.R"))

PJAT_DIR <- pjats_dir

# Bring in the two new UI helpers by sourcing just their definitions from the
# main script (avoids launching the Shiny app / CLI at source time).
.load_map_blocks <- function() {
  src <- readLines(file.path(project_root, "Plasmid_Designer.R"), warn = FALSE)
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) > 0L && all(!is.na(a))) a else b
  }
  e <- new.env(parent = globalenv())
  e$`%||%` <- `%||%`
  for (nm in c("addgene_map_src", "addgene_page_url", "pjd_map_feature_class",
               "PJD_MAP_COLORS", ".pjd_addgene_link")) {
    if (exists(nm, mode = "function") || exists(nm)) {
      assign(nm, get(nm), envir = e)
    }
  }
  start <- grep("^\\.pjd_addgene_link <- function", src)[1L]
  end <- grep("^\\.pjd_backbone_panel <- function", src)[1L] - 1L
  eval(parse(text = paste(src[start:end], collapse = "\n")), envir = e)
  e
}
blk <- .load_map_blocks()

make_cat <- function() {
  data.frame(name = "pJAT74", addgene_id = 204312L, map_available = TRUE,
             map_sequence_id = 406186L, size_bp = 7986L,
             stringsAsFactors = FALSE)
}

test_that("enhanced PNG block renders container and zoom", {
  html <- as.character(blk$.pjd_addgene_map_block(make_cat(), zoom_pct = 150))
  expect_match(html, "pjd-png-view", fixed = TRUE)
  expect_match(html, "width:150%", fixed = TRUE)
  expect_match(html, "sequence ID 406186", fixed = TRUE)
  expect_false(grepl("pjd-png-overlay", html, fixed = TRUE))
  expect_false(grepl("pjd-layer-strip", html, fixed = TRUE))
})

test_that("PNG block renders at 100% by default", {
  html <- as.character(blk$.pjd_addgene_map_block(make_cat()))
  expect_match(html, "pjd-png-view", fixed = TRUE)
  expect_match(html, "width:100%", fixed = TRUE)
})

test_that("PNG block clamps out-of-range zoom values", {
  html <- as.character(blk$.pjd_addgene_map_block(make_cat(), zoom_pct = 9999))
  expect_match(html, "width:300%", fixed = TRUE)
  html2 <- as.character(blk$.pjd_addgene_map_block(make_cat(), zoom_pct = 1))
  expect_match(html2, "width:25%", fixed = TRUE)
})

test_that("plasmidplot parses the real pJAT74 GenBank file", {
  skip_if_not_installed("plasmidplot")
  p <- plasmidplot::read_genbank(file.path(PJAT_DIR, "pJAT74.gbk"))
  expect_s3_class(p, "plasmid")
  txt <- capture.output(print(p))
  expect_true(any(grepl("7986|7,986", txt)))
})

test_that("plasmidplot renders circular and linear maps without error", {
  skip_if_not_installed("plasmidplot")
  suppressMessages(library(plasmidplot))
  p <- plasmidplot::read_genbank(file.path(PJAT_DIR, "pJAT74.gbk"))
  f1 <- tempfile(fileext = ".png")
  f2 <- tempfile(fileext = ".png")
  on.exit({ unlink(f1); unlink(f2) }, add = TRUE)
  grDevices::png(f1, width = 800, height = 800)
  expect_no_error(plot(p))
  grDevices::dev.off()
  expect_gt(file.info(f1)$size, 10000)
  grDevices::png(f2, width = 1200, height = 500)
  expect_no_error(plot(p, style = plasmidplot::pp_style("minimal", layout = "linear")))
  grDevices::dev.off()
  expect_gt(file.info(f2)$size, 5000)
})

test_that("plasmidplot block helper returns the plot slot", {
  skip_if_not(file.exists(file.path(PJAT_DIR, "pJAT74.gbk")), "pJAT data missing")
  # When plasmidplot IS installed, the helper returns the plotOutput slot.
  if (requireNamespace("plasmidplot", quietly = TRUE)) {
    html <- as.character(blk$.pjd_plasmidplot_map_block(make_cat()))
    expect_match(html, "plasmidplot_map", fixed = TRUE)
  } else {
    html <- as.character(blk$.pjd_plasmidplot_map_block(make_cat()))
    expect_match(html, "not installed", fixed = TRUE)
  }
})
