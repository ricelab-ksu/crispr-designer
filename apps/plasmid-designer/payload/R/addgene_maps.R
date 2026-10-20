# ============================================================================
# R/addgene_maps.R — Addgene / SnapGene vector-map images for the pJAT catalog
# ============================================================================
# Addgene serves the interactive "vector map" of every plasmid as a PNG under
# sequences.addgene.org. The path is NOT derivable from the Addgene ID alone:
#
#   https://sequences.addgene.org/snapgene-media/v<version>/sequences/<seq_id>/
#       <uuid>/addgene-plasmid-<addgene_id>-sequence-<seq_id>-map.png
#
# Both <seq_id> and <uuid> are per-plasmid and only discoverable from the public
# plasmid page (https://www.addgene.org/<addgene_id>/). tools/fetch_addgene_maps.R
# scrapes that page once per plasmid, downloads the PNG into www/addgene_maps/
# and records the resolved URL in donor_templates/addgene_maps.yml.
#
# This module is deliberately pure/offline: it never touches the network. It only
# reads the generated YAML, joins it to the pJAT catalog, and answers "is there a
# local image for this backbone, and what is its Shiny-served URL?".
#
# Images are (c) Addgene and are redistributed here for local, non-commercial
# viewing with attribution; the app always links back to the source page.
# ============================================================================


# ---- Portable path helpers --------------------------------------------------
# app_root() and the reference-store resolvers live in R/app_paths.R, which
# sits next to this file. Source it here so this module can be loaded on its
# own (as the test suite does).
if (!exists("app_root", mode = "function")) {
  .self_dir <- local({
    for (i in seq_along(sys.frames())) {
      of <- tryCatch(sys.frames()[[i]]$ofile, error = function(e) NULL)
      if (!is.null(of) && nzchar(of)) {
        return(dirname(normalizePath(of, mustWork = FALSE)))
      }
    }
    NA_character_
  })
  .cands <- c(
    if (!is.na(.self_dir)) file.path(.self_dir, "app_paths.R"),
    file.path("R", "app_paths.R"),
    file.path(getwd(), "R", "app_paths.R")
  )
  for (p in unique(.cands)) {
    if (!is.na(p) && file.exists(p)) {
      tryCatch(source(p, local = FALSE), error = function(e) NULL)
      break
    }
  }
}
if (exists("load_app_paths", mode = "function")) load_app_paths()
ADDGENE_MAPS_VERSION <- "1.0.0"

ADDGENE_MAPS_DIRNAME <- "addgene_maps"

# `%||%` normally comes from Plasmid_Designer.R; define it when this module is
# sourced standalone (unit tests, tools). Only pass scalars through it: other
# modules in this app redefine `%||%` with a scalar-only test.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) == 1L && !is.na(a)) a else b
  }
}

#' Resolve the Plasmid_Designer project root (independent of the working dir).
#' @keywords internal
.addgene_maps_project_root <- function() {
  candidates <- unique(c(
    getwd(),
    app_root()
  ))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (file.exists(script_path)) {
      candidates <- unique(c(dirname(normalizePath(script_path, mustWork = FALSE)),
                             candidates))
    }
  }
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "R")) || dir.exists(file.path(cand, "donor_templates"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  candidates[1]
}

ADDGENE_MAPS_ROOT <- .addgene_maps_project_root()

#' Path to the generated Addgene map manifest.
#' @export
addgene_maps_yml <- function(root = ADDGENE_MAPS_ROOT) {
  file.path(root, "donor_templates", "addgene_maps.yml")
}

#' Directory holding the downloaded map PNGs (Shiny serves `www/` at `/`).
#' @export
addgene_maps_dir <- function(root = ADDGENE_MAPS_ROOT) {
  file.path(root, "www", ADDGENE_MAPS_DIRNAME)
}

#' Public Addgene plasmid page for an ID.
#' @export
addgene_page_url <- function(addgene_id) {
  sprintf("https://www.addgene.org/%d/", as.integer(addgene_id))
}

# ----------------------------------------------------------------------------
# HTML scraping helpers (pure: used by the fetch tool and by the unit tests)
# ----------------------------------------------------------------------------

#' Extract the SnapGene map PNG URL from an Addgene plasmid page.
#'
#' Handles the URL appearing raw in `src="..."`, JSON-escaped as `https:\/\/`,
#' or double-escaped inside a script payload. Returns NA when the page carries
#' no map (some records are sequence-only).
#'
#' @param html Character scalar: the page HTML (or any text containing the URL).
#' @param addgene_id Addgene plasmid ID the map must belong to.
#' @return Character scalar URL, or NA_character_.
#' @export
addgene_extract_map_url <- function(html, addgene_id) {
  if (is.null(html) || length(html) == 0L || is.na(html[1L]) || !nzchar(html[1L])) {
    return(NA_character_)
  }
  txt <- html[1L]
  # Undo JSON / HTML escaping so one pattern covers every embedding.
  txt <- gsub("\\\\+/", "/", txt)
  txt <- gsub("&amp;", "&", txt, fixed = TRUE)
  id <- as.integer(addgene_id)
  pat <- paste0(
    "https?://sequences\\.addgene\\.org/[^\"'<>[:space:]\\\\]*?",
    "addgene-plasmid-", id, "-sequence-[0-9]+-map\\.png"
  )
  m <- regmatches(txt, regexpr(pat, txt, perl = TRUE))
  if (length(m) == 0L || !nzchar(m[1L])) return(NA_character_)
  m[1L]
}

#' Extract the Addgene sequence ID from a map URL.
#' @export
addgene_map_sequence_id <- function(url) {
  m <- regmatches(url, regexpr("/sequences/[0-9]+/", url))
  if (length(m) == 0L) return(NA_integer_)
  as.integer(gsub("[^0-9]", "", m[1L]))
}

#' Local filename for a plasmid's map PNG.
#' @export
addgene_map_filename <- function(addgene_id) {
  sprintf("addgene_%d.png", as.integer(addgene_id))
}

#' Shiny-served relative URL for a plasmid's map (NULL when not downloaded).
#'
#' Shiny serves the app's `www/` directory at the site root, so
#' `www/addgene_maps/addgene_204289.png` is reachable as
#' `addgene_maps/addgene_204289.png`.
#'
#' @export
addgene_map_src <- function(addgene_id, root = ADDGENE_MAPS_ROOT) {
  p <- file.path(addgene_maps_dir(root), addgene_map_filename(addgene_id))
  if (!file.exists(p)) return(NULL)
  paste0(ADDGENE_MAPS_DIRNAME, "/", basename(p))
}


# ----------------------------------------------------------------------------
# Manifest I/O
# ----------------------------------------------------------------------------

.addgene_maps_empty <- function() {
  data.frame(
    addgene_id = integer(0),
    plasmid = character(0),
    page_url = character(0),
    map_url = character(0),
    sequence_id = integer(0),
    map_file = character(0),
    width_px = integer(0),
    height_px = integer(0),
    bytes = numeric(0),
    md5 = character(0),
    fetched = character(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
}

#' Read the generated Addgene map manifest.
#'
#' Returns a zero-row frame with the documented columns when the manifest has not
#' been generated yet, so the app degrades gracefully instead of erroring.
#'
#' @param yml Path to `addgene_maps.yml`.
#' @return data.frame, one row per attempted plasmid.
#' @export
addgene_maps_read <- function(yml = addgene_maps_yml()) {
  if (!file.exists(yml)) return(.addgene_maps_empty())
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read the Addgene map manifest.")
  }
  y <- yaml::read_yaml(yml)
  # NOTE: do not use `%||%` here — guide_bundle.R / plasmid_assembly.R redefine it
  # with a scalar-only test (`!is.na(a)`), which errors on a 24-element list.
  recs <- y$maps
  if (is.null(recs) || length(recs) == 0L) return(.addgene_maps_empty())

  chr <- function(v) if (is.null(v) || length(v) == 0L) NA_character_ else as.character(v)[1L]
  int <- function(v) {
    if (is.null(v) || length(v) == 0L) return(NA_integer_)
    suppressWarnings(as.integer(v[[1L]]))
  }
  num <- function(v) {
    if (is.null(v) || length(v) == 0L) return(NA_real_)
    suppressWarnings(as.numeric(v[[1L]]))
  }

  out <- data.frame(
    addgene_id  = vapply(recs, function(r) int(r$addgene_id), integer(1)),
    plasmid     = vapply(recs, function(r) chr(r$plasmid), character(1)),
    page_url    = vapply(recs, function(r) chr(r$page_url), character(1)),
    map_url     = vapply(recs, function(r) chr(r$map_url), character(1)),
    sequence_id = vapply(recs, function(r) int(r$sequence_id), integer(1)),
    map_file    = vapply(recs, function(r) chr(r$map_file), character(1)),
    width_px    = vapply(recs, function(r) int(r$width_px), integer(1)),
    height_px   = vapply(recs, function(r) int(r$height_px), integer(1)),
    bytes       = vapply(recs, function(r) num(r$bytes), numeric(1)),
    md5         = vapply(recs, function(r) chr(r$md5), character(1)),
    fetched     = vapply(recs, function(r) chr(r$fetched), character(1)),
    status      = vapply(recs, function(r) chr(r$status), character(1)),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$addgene_id), , drop = FALSE]
}

#' Write the Addgene map manifest (deterministic key order + provenance header).
#'
#' @param maps data.frame from `addgene_maps_read()` / the fetch tool.
#' @param yml Destination path.
#' @return Invisibly, the path written.
#' @export
addgene_maps_write <- function(maps, yml = addgene_maps_yml()) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to write the Addgene map manifest.")
  }
  dir.create(dirname(yml), recursive = TRUE, showWarnings = FALSE)
  recs <- lapply(seq_len(nrow(maps)), function(i) {
    r <- maps[i, , drop = FALSE]
    list(
      addgene_id  = as.integer(r$addgene_id),
      plasmid     = r$plasmid,
      page_url    = r$page_url,
      map_url     = r$map_url,
      sequence_id = if (is.na(r$sequence_id)) NULL else as.integer(r$sequence_id),
      map_file    = r$map_file,
      width_px    = if (is.na(r$width_px)) NULL else as.integer(r$width_px),
      height_px   = if (is.na(r$height_px)) NULL else as.integer(r$height_px),
      bytes       = if (is.na(r$bytes)) NULL else as.numeric(r$bytes),
      md5         = r$md5,
      fetched     = r$fetched,
      status      = r$status
    )
  })
  body <- yaml::as.yaml(
    list(version = ADDGENE_MAPS_VERSION,
         generator = "tools/fetch_addgene_maps.R",
         source = "https://www.addgene.org/<addgene_id>/ (vector maps: SnapGene via sequences.addgene.org)",
         maps = recs),
    indent = 2
  )
  header <- paste0(
    "# Addgene vector-map manifest for the Janelia Atalanta (pJAT) collection\n",
    "# Generated by tools/fetch_addgene_maps.R from the public Addgene plasmid pages.\n",
    "# Do not hand-edit: re-run the fetch tool instead.\n",
    "# Images are (c) Addgene, shown locally with attribution and a link to the\n",
    "# source page in the app.\n",
    "\n"
  )
  writeLines(paste0(header, body), con = yml, useBytes = TRUE)
  invisible(yml)
}

# ----------------------------------------------------------------------------
# Catalog join
# ----------------------------------------------------------------------------

#' Attach Addgene vector-map columns to the pJAT catalog.
#'
#' Adds `map_url`, `map_file`, `map_sequence_id` and `map_available` to a
#' `pjat_load_catalog()` data.frame. Rows are never dropped or reordered; a
#' plasmid without a manifest entry simply gets NA values and
#' `map_available = FALSE`.
#'
#' @param catalog data.frame from `pjat_load_catalog()`.
#' @param maps Optional manifest from `addgene_maps_read()`.
#' @param root Project root used to check that the PNG exists on disk.
#' @return `catalog` with the four map columns appended.
#' @export
addgene_maps_join <- function(catalog, maps = NULL, root = ADDGENE_MAPS_ROOT) {
  if (is.null(catalog) || nrow(catalog) == 0L) return(catalog)
  if (is.null(maps)) maps <- addgene_maps_read()

  pick <- function(col, default) {
    if (nrow(maps) == 0L) return(rep(default, nrow(catalog)))
    idx <- match(as.integer(catalog$addgene_id), as.integer(maps$addgene_id))
    maps[[col]][idx]
  }

  catalog$map_url         <- pick("map_url", NA_character_)
  catalog$map_file        <- pick("map_file", NA_character_)
  catalog$map_sequence_id <- pick("sequence_id", NA_integer_)

  catalog$map_available <- vapply(seq_len(nrow(catalog)), function(i) {
    id <- catalog$addgene_id[i]
    if (is.na(id)) return(FALSE)
    file.exists(file.path(addgene_maps_dir(root), addgene_map_filename(id)))
  }, logical(1))

  catalog
}

# ----------------------------------------------------------------------------
# PNG helpers (no external package: IHDR is bytes 16-23, big-endian)
# ----------------------------------------------------------------------------

#' Read the pixel dimensions of a PNG without loading the image.
#'
#' @param path Path to a PNG file.
#' @return list(width, height), with NA_integer_ values if not a readable PNG.
#' @export
png_dimensions <- function(path) {
  na <- list(width = NA_integer_, height = NA_integer_)
  if (!file.exists(path)) return(na)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  sig <- readBin(con, "raw", n = 8L)
  if (length(sig) < 8L) return(na)
  if (!identical(as.integer(sig),
                 c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))) return(na)
  hdr <- readBin(con, "raw", n = 16L)   # length(4) + "IHDR"(4) + width(4) + height(4)
  if (length(hdr) < 16L) return(na)
  be_int <- function(b) sum(as.integer(b) * c(256^3, 256^2, 256, 1))
  list(width  = as.integer(be_int(hdr[9:12])),
       height = as.integer(be_int(hdr[13:16])))
}

#' Is this file a plausible SnapGene map PNG (signature + sane dimensions)?
#' @export
is_addgene_map_png <- function(path, min_bytes = 5000) {
  if (!file.exists(path)) return(FALSE)
  if (file.info(path)$size < min_bytes) return(FALSE)
  d <- png_dimensions(path)
  !is.na(d$width) && !is.na(d$height) && d$width > 100L && d$height > 100L
}
