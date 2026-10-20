#!/usr/bin/env Rscript
# ============================================================================
# tools/fetch_addgene_maps.R — download Addgene/SnapGene vector maps for pJAT
# ============================================================================
# For each plasmid in donor_templates/pjat_plasmids.yml this tool:
#   1. GETs the public page  https://www.addgene.org/<addgene_id>/
#   2. extracts the SnapGene map PNG URL
#      (https://sequences.addgene.org/snapgene-media/.../
#       addgene-plasmid-<id>-sequence-<seq_id>-map.png)
#   3. downloads the PNG to  www/addgene_maps/addgene_<id>.png
#   4. records page URL, image URL, sequence ID, dimensions, size, md5 and
#      status in donor_templates/addgene_maps.yml
#
# The sequence ID and UUID are per-plasmid, so they cannot be templated from the
# Addgene ID — scraping the public page is the only reliable route.
#
# USAGE
#   Rscript tools/fetch_addgene_maps.R                 # fetch everything missing
#   Rscript tools/fetch_addgene_maps.R --force         # re-download everything
#   Rscript tools/fetch_addgene_maps.R --ids 204289,204290
#   Rscript tools/fetch_addgene_maps.R --dry-run       # resolve URLs, write nothing
#
# Options: --out <dir>  --yml <file>  --sleep <seconds>  --quiet
#
# Exit codes: 0 = every requested plasmid has an image on disk,
#             1 = at least one plasmid failed (details in the manifest/status).
#
# Images are (c) Addgene. They are cached locally for offline viewing inside the
# Plasmid Designer only, with attribution and a link back to the source page.
# Re-running is idempotent; cached files are reused unless --force is given.
# ============================================================================

suppressPackageStartupMessages({
  ok <- requireNamespace("httr", quietly = TRUE)
})

# ---- Resolve project root --------------------------------------------------

.script_path <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(sub("^--file=", "", a[1L]))
  NA_character_
}

.resolve_tool_root <- function() {
  sp <- .script_path()
  cands <- unique(c(
    if (!is.na(sp)) dirname(dirname(normalizePath(sp, mustWork = FALSE))) else NULL,
    getwd(),
    "/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer"
  ))
  for (c in cands) {
    if (file.exists(file.path(c, "donor_templates", "pjat_plasmids.yml"))) {
      return(normalizePath(c, mustWork = FALSE))
    }
  }
  stop("Could not locate the Plasmid_Designer project root (donor_templates/pjat_plasmids.yml).")
}

ROOT <- .resolve_tool_root()

# addgene_maps.R uses %||%, which normally comes from Plasmid_Designer.R.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) > 0L && all(!is.na(a))) a else b
  }
}
source(file.path(ROOT, "R", "addgene_maps.R"))

if (!ok) {
  stop("Package 'httr' is required by tools/fetch_addgene_maps.R ",
       "(install.packages('httr')).")
}

# ---- CLI arguments ---------------------------------------------------------

.parse_flag <- function(args, flag, default = NULL) {
  eq <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(eq)) return(sub(paste0("^", flag, "="), "", eq[1L]))
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1L])
  default
}

ARGS <- commandArgs(trailingOnly = TRUE)
if (any(grepl("^(--help|-h)$", ARGS))) {
  cat("Usage: Rscript tools/fetch_addgene_maps.R [--force] [--dry-run]",
      "[--ids 204289,204290] [--out <dir>] [--yml <file>]",
      "[--sleep <seconds>] [--quiet]\n")
  quit(status = 0L)
}

FORCE   <- "--force" %in% ARGS
DRY_RUN <- "--dry-run" %in% ARGS
QUIET   <- "--quiet" %in% ARGS
OUT_DIR <- .parse_flag(ARGS, "--out", addgene_maps_dir(ROOT))
YML     <- .parse_flag(ARGS, "--yml", addgene_maps_yml(ROOT))
SLEEP   <- as.numeric(.parse_flag(ARGS, "--sleep", "1.0"))
IDS_ARG <- .parse_flag(ARGS, "--ids", NULL)

if (!is.finite(SLEEP) || SLEEP < 0) SLEEP <- 1.0
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

say <- function(...) if (!QUIET) cat(..., "\n", sep = "")

# ---- Catalog (YAML only: no GenBank parsing needed) ------------------------

CAT_YML <- file.path(ROOT, "donor_templates", "pjat_plasmids.yml")
if (!file.exists(CAT_YML)) stop("Catalog not found: ", CAT_YML)
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required (install.packages('yaml')).")
}

catalog <- yaml::read_yaml(CAT_YML)$plasmids
if (is.null(catalog) || length(catalog) == 0L) stop("No plasmids in ", CAT_YML)

wanted <- vapply(catalog, function(p) as.integer(p$addgene_id), integer(1))
if (!is.null(IDS_ARG)) {
  keep <- wanted %in% as.integer(trimws(strsplit(IDS_ARG, ",")[[1L]]))
  catalog <- catalog[keep]
  wanted <- wanted[keep]
  if (length(catalog) == 0L) stop("--ids matched no plasmid in the catalog.")
}

say(sprintf("Addgene map fetch: %d plasmid(s) -> %s", length(catalog), OUT_DIR))

# Previous manifest: lets a cached run keep the resolved URL / sequence ID
# without re-scraping every page.
prev <- tryCatch(addgene_maps_read(YML), error = function(e) NULL)

# ---- HTTP helpers ----------------------------------------------------------

UA <- paste0("PlasmidDesigner/", ADDGENE_MAPS_VERSION,
             " (local research app; +https://www.addgene.org/)")

.http_get <- function(url, accept = "*/*", tries = 3L) {
  last <- NULL
  for (i in seq_len(tries)) {
    r <- tryCatch(
      httr::GET(url,
                httr::user_agent(UA),
                httr::accept(accept),
                httr::timeout(45)),
      error = function(e) e
    )
    if (!inherits(r, "error") && httr::status_code(r) == 200L) return(r)
    last <- if (inherits(r, "error")) conditionMessage(r) else
      paste0("HTTP ", httr::status_code(r))
    if (i < tries) Sys.sleep(i * 2)
  }
  stop(last %||% "request failed")
}

# ---- Main loop -------------------------------------------------------------

rows <- vector("list", length(catalog))
failures <- 0L

for (i in seq_along(catalog)) {
  p  <- catalog[[i]]
  id <- as.integer(p$addgene_id)
  nm <- as.character(p$name %||% NA_character_)
  png_path <- file.path(OUT_DIR, addgene_map_filename(id))
  page <- addgene_page_url(id)

  row <- list(
    addgene_id = id, plasmid = nm, page_url = page,
    map_url = NA_character_, sequence_id = NA_integer_,
    map_file = paste0(ADDGENE_MAPS_DIRNAME, "/", addgene_map_filename(id)),
    width_px = NA_integer_, height_px = NA_integer_,
    bytes = NA_real_, md5 = NA_character_,
    fetched = format(Sys.Date()), status = NA_character_
  )

  # Reuse a valid cached image unless --force.
  if (!FORCE && is_addgene_map_png(png_path)) {
    d <- png_dimensions(png_path)
    row$width_px  <- d$width
    row$height_px <- d$height
    row$bytes     <- as.numeric(file.info(png_path)$size)
    row$md5       <- unname(tools::md5sum(png_path))

    # Keep the previously resolved URL / sequence ID so a cached run never loses
    # provenance; scrape the page only when the manifest has no URL yet.
    pr <- if (!is.null(prev) && nrow(prev) > 0L) {
      prev[match(id, prev$addgene_id), , drop = FALSE]
    } else NULL

    if (!is.null(pr) && nrow(pr) == 1L && !is.na(pr$map_url)) {
      row$map_url     <- pr$map_url
      row$sequence_id <- pr$sequence_id
      if (!is.na(pr$fetched)) row$fetched <- pr$fetched
    } else {
      row <- tryCatch({
        html <- httr::content(.http_get(page, accept = "text/html"), "text",
                              encoding = "UTF-8")
        u <- addgene_extract_map_url(html, id)
        row$map_url     <- u
        row$sequence_id <- addgene_map_sequence_id(u)
        row
      }, error = function(e) row)
    }

    row$status <- "cached"
    rows[[i]] <- row
    say(sprintf("  [%2d/%2d] %-8s #%d  cached", i, length(catalog), nm, id))
    next
  }

  res <- tryCatch({
    html <- httr::content(.http_get(page, accept = "text/html"), "text",
                          encoding = "UTF-8")
    map_url <- addgene_extract_map_url(html, id)
    if (is.na(map_url)) stop("no SnapGene map URL found on the Addgene page")

    if (DRY_RUN) {
      row$map_url <- map_url
      row$sequence_id <- addgene_map_sequence_id(map_url)
      row$status <- "resolved (dry-run)"
      row
    } else {
      img <- .http_get(map_url, accept = "image/png,image/*")
      ctype <- httr::headers(img)[["content-type"]] %||% ""
      body <- httr::content(img, "raw")
      # Write the temp file next to the destination: the project volume is often a
      # different device from TMPDIR, where file.rename() fails with EXDEV.
      tmp <- file.path(OUT_DIR, paste0(".", addgene_map_filename(id), ".part"))
      writeBin(body, tmp)

      if (!is_addgene_map_png(tmp)) {
        unlink(tmp)
        stop(sprintf("downloaded payload is not a valid map PNG (content-type: %s, %d bytes)",
                     ctype, length(body)))
      }

      if (!file.rename(tmp, png_path)) {
        ok_copy <- file.copy(tmp, png_path, overwrite = TRUE)
        unlink(tmp)
        if (!ok_copy) stop("could not write ", png_path)
      }
      d <- png_dimensions(png_path)
      row$map_url     <- map_url
      row$sequence_id <- addgene_map_sequence_id(map_url)
      row$width_px    <- d$width
      row$height_px   <- d$height
      row$bytes       <- as.numeric(file.info(png_path)$size)
      row$md5         <- unname(tools::md5sum(png_path))
      row$status      <- "ok"
      row
    }
  }, error = function(e) {
    row$status <- paste0("error: ", conditionMessage(e))
    row
  })

  if (startsWith(res$status, "error")) failures <- failures + 1L
  rows[[i]] <- res
  say(sprintf("  [%2d/%2d] %-8s #%d  %s", i, length(catalog), nm, id, res$status))

  if (i < length(catalog) && SLEEP > 0) Sys.sleep(SLEEP)
}

maps <- do.call(rbind, lapply(rows, function(r) {
  as.data.frame(r, stringsAsFactors = FALSE)
}))

# ---- Write manifest --------------------------------------------------------

if (DRY_RUN) {
  say("\nDry run: manifest not written. Resolved URLs:")
  for (i in seq_len(nrow(maps))) {
    say(sprintf("  #%d %s -> %s", maps$addgene_id[i], maps$plasmid[i],
                ifelse(is.na(maps$map_url[i]), "<none>", maps$map_url[i])))
  }
} else {
  addgene_maps_write(maps, YML)
  say("\nWrote ", YML)
}

n_ok <- sum(maps$status %in% c("ok", "cached"))
say(sprintf("Summary: %d/%d with a local map, %d failed, %d resolved only%s",
            n_ok, nrow(maps), failures,
            sum(grepl("^resolved", maps$status)),
            if (DRY_RUN) " (dry-run)" else ""))

if (failures > 0L) {
  say("\nFailures:")
  bad <- maps[startsWith(maps$status, "error"), , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    say(sprintf("  #%d %s: %s", bad$addgene_id[i], bad$plasmid[i], bad$status[i]))
  }
  quit(status = 1L)
}
quit(status = 0L)
