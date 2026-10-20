library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "addgene_maps.R"))) {
  project_root <- project_root
}

source(file.path(project_root, "R", "addgene_maps.R"))

MAPS_YML <- file.path(project_root, "donor_templates", "addgene_maps.yml")
MAPS_DIR <- file.path(project_root, "www", "addgene_maps")

# ---- URL extraction ---------------------------------------------------------

test_that("map URL is extracted from a plain src attribute", {
  html <- paste0('<img src="https://sequences.addgene.org/snapgene-media/',
                 'v3.101.0/sequences/406181/8ebe0074-1de2-4f33-98ff-87bacda6d3ba/',
                 'addgene-plasmid-204289-sequence-406181-map.png" alt="406181_map">')
  url <- addgene_extract_map_url(html, 204289)
  expect_false(is.na(url))
  expect_match(url, "addgene-plasmid-204289-sequence-406181-map\\.png$")
  expect_equal(addgene_map_sequence_id(url), 406181L)
})

test_that("map URL is extracted when JSON-escaped", {
  html <- paste0('{"map":"https:\\/\\/sequences.addgene.org\\/snapgene-media\\/',
                 'v3.101.0\\/sequences\\/406181\\/8ebe0074-1de2-4f33-98ff-87bacda6d3ba\\/',
                 'addgene-plasmid-204289-sequence-406181-map.png"}')
  url <- addgene_extract_map_url(html, 204289)
  expect_false(is.na(url))
  expect_false(grepl("\\\\", url))
})

test_that("map URL for a different Addgene ID is rejected", {
  html <- paste0('<img src="https://sequences.addgene.org/snapgene-media/',
                 'v3.101.0/sequences/406181/8ebe0074/',
                 'addgene-plasmid-204289-sequence-406181-map.png">')
  expect_true(is.na(addgene_extract_map_url(html, 204290)))
  expect_true(is.na(addgene_extract_map_url("", 204289)))
  expect_true(is.na(addgene_extract_map_url(NULL, 204289)))
})

# ---- Paths ------------------------------------------------------------------

test_that("filename, page URL and Shiny src are consistent", {
  expect_equal(addgene_map_filename(204289), "addgene_204289.png")
  expect_equal(addgene_page_url(204289), "https://www.addgene.org/204289/")

  tmp <- tempfile("pjd_maps_")
  dir.create(file.path(tmp, "www", ADDGENE_MAPS_DIRNAME), recursive = TRUE)
  expect_null(addgene_map_src(204289, root = tmp))

  file.create(file.path(tmp, "www", ADDGENE_MAPS_DIRNAME,
                        addgene_map_filename(204289)))
  expect_equal(addgene_map_src(204289, root = tmp),
               "addgene_maps/addgene_204289.png")
  unlink(tmp, recursive = TRUE)
})

# ---- Manifest I/O -----------------------------------------------------------

test_that("missing manifest yields a zero-row frame with the documented columns", {
  m <- addgene_maps_read(tempfile())
  expect_equal(nrow(m), 0L)
  expect_true(all(c("addgene_id", "plasmid", "page_url", "map_url",
                    "sequence_id", "map_file", "width_px", "height_px",
                    "bytes", "md5", "fetched", "status") %in% names(m)))
})

test_that("manifest round-trips through YAML", {
  fixture <- data.frame(
    addgene_id = c(204289L, 204290L),
    plasmid = c("pJAT39", "pJAT11"),
    page_url = c(addgene_page_url(204289), addgene_page_url(204290)),
    map_url = c("https://sequences.addgene.org/a/addgene-plasmid-204289-sequence-1-map.png",
                "https://sequences.addgene.org/b/addgene-plasmid-204290-sequence-2-map.png"),
    sequence_id = c(1L, 2L),
    map_file = c("addgene_maps/addgene_204289.png",
                 "addgene_maps/addgene_204290.png"),
    width_px = c(907L, 900L), height_px = c(753L, 700L),
    bytes = c(123456, 234567), md5 = c("aaa", "bbb"),
    fetched = c("2026-10-17", "2026-10-17"),
    status = c("ok", "ok"), stringsAsFactors = FALSE
  )
  yml <- tempfile(fileext = ".yml")
  addgene_maps_write(fixture, yml)
  expect_true(file.exists(yml))
  back <- addgene_maps_read(yml)
  expect_equal(nrow(back), 2L)
  expect_equal(back$addgene_id, c(204289L, 204290L))
  expect_equal(back$sequence_id, c(1L, 2L))
  expect_equal(back$map_url, fixture$map_url)
  unlink(yml)

# ---- Catalog join -----------------------------------------------------------

test_that("join leaves the catalog untouched when there is no manifest", {
  cat <- data.frame(name = c("pJAT74", "pJAT77"),
                    addgene_id = c(204296L, 204295L),
                    stringsAsFactors = FALSE)
  out <- addgene_maps_join(cat, maps = addgene_maps_read(tempfile()),
                           root = tempdir())
  expect_equal(nrow(out), 2L)
  expect_equal(out$name, cat$name)
  expect_true(all(is.na(out$map_url)))
  expect_true(all(is.na(out$map_file)))
  expect_true(all(is.na(out$map_sequence_id)))
  expect_true(all(!out$map_available))
})

test_that("join matches by Addgene ID without reordering or dropping rows", {
  cat <- data.frame(name = c("pJAT77", "pJAT74"),
                    addgene_id = c(204295L, 204296L),
                    stringsAsFactors = FALSE)
  maps <- data.frame(
    addgene_id = c(204296L, 204295L),
    plasmid = c("pJAT74", "pJAT77"),
    page_url = c(addgene_page_url(204296), addgene_page_url(204295)),
    map_url = c("u74", "u77"),
    sequence_id = c(2L, 1L),
    map_file = c("addgene_maps/addgene_204296.png",
                 "addgene_maps/addgene_204295.png"),
    width_px = NA_integer_, height_px = NA_integer_,
    bytes = NA_real_, md5 = NA_character_, fetched = NA_character_,
    status = c("ok", "ok"), stringsAsFactors = FALSE
  )
  out <- addgene_maps_join(cat, maps = maps, root = tempdir())
  expect_equal(out$name, c("pJAT77", "pJAT74"))
  expect_equal(out$map_url, c("u77", "u74"))
  expect_equal(out$map_sequence_id, c(1L, 2L))
  expect_true(all(!out$map_available))   # no files under tempdir()
})

# ---- PNG helpers ------------------------------------------------------------

test_that("png_dimensions reads width/height from IHDR", {
  skip_if_not(requireNamespace("png", quietly = TRUE), "png package not installed")
  f <- tempfile(fileext = ".png")
  png::writePNG(array(255, dim = c(120, 240, 1)), f)   # rows = height
  d <- png_dimensions(f)
  expect_equal(d$width, 240L)
  expect_equal(d$height, 120L)
  expect_true(is_addgene_map_png(f, min_bytes = 1L))
  expect_false(is_addgene_map_png(f))                  # below 5 kB floor
  unlink(f)
})

test_that("non-PNG payloads are rejected", {
  f <- tempfile(fileext = ".png")
  writeLines("<html>not an image</html>", f)
  expect_true(is.na(png_dimensions(f)$width))
  expect_false(is_addgene_map_png(f, min_bytes = 1L))
  expect_false(is_addgene_map_png(tempfile()))         # missing file
  unlink(f)
})

# ---- Generated cache (integration) -----------------------------------------

test_that("generated manifest covers the whole pJAT collection", {
  skip_if_not(file.exists(MAPS_YML), "addgene_maps.yml not generated")
  m <- addgene_maps_read(MAPS_YML)
  expect_equal(nrow(m), 24L)
  expect_setequal(m$addgene_id, 204289:204312)
  expect_true(all(m$status %in% c("ok", "cached")))
  expect_true(all(vapply(seq_len(nrow(m)), function(i) {
    grepl(sprintf("addgene-plasmid-%d-sequence-[0-9]+-map\\.png$", m$addgene_id[i]),
          m$map_url[i])
  }, logical(1))))
  expect_true(all(m$width_px > 100L))
  expect_true(all(m$height_px > 100L))
  expect_false(any(duplicated(m$addgene_id)))
})

test_that("downloaded maps are valid PNGs on disk", {
  skip_if_not(dir.exists(MAPS_DIR), "map cache not generated")
  files <- list.files(MAPS_DIR, pattern = "^addgene_[0-9]+\\.png$",
                      full.names = TRUE)
  expect_equal(length(files), 24L)
  for (f in files) expect_true(is_addgene_map_png(f), info = basename(f))
})

test_that("pjat_load_catalog exposes map columns for every backbone", {
  src_pjat <- file.path(project_root, "R", "pjat_plasmids.R")
  skip_if_not(file.exists(src_pjat), "pjat_plasmids.R missing")
  source(src_pjat)
  cat <- pjat_load_catalog()
  expect_equal(nrow(cat), 24L)
  expect_true(all(c("map_url", "map_file", "map_sequence_id",
                    "map_available") %in% names(cat)))
  skip_if_not(dir.exists(MAPS_DIR), "map cache not generated")
  expect_true(all(cat$map_available))
  expect_false(any(is.na(cat$map_sequence_id)))
})

})
