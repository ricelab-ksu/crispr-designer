# ============================================================================
# R/sample_inputs.R — Bundled sample inputs for Plasmid Designer
# ============================================================================
# Provides offline loaders for the bundled sample inputs stored in `uploads/`:
#   1. Cas9_w_melanogaster_w_g0351_guides_for_plasmid_designer.csv
#   2. sample_donor_3xP3_dsRed.fasta (pJAT30 Addgene #204294, coords 4749..6609)
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
SAMPLE_GUIDE_BUNDLE_FILENAME <- "Cas9_w_melanogaster_w_g0351_guides_for_plasmid_designer.csv"
SAMPLE_DONOR_FASTA_FILENAME  <- "sample_donor_3xP3_dsRed.fasta"

.pjd_sample_resolve_root <- function(root = NULL) {
  if (!is.null(root)) {
    if (!is.character(root) || length(root) != 1L || is.na(root) || !nzchar(trimws(root))) {
      stop("`root` must be a non-empty character scalar path.", call. = FALSE)
    }
    return(root)
  }
  candidates <- unique(c(
    getwd(),
    app_root()
  ))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (file.exists(script_path)) {
      candidates <- unique(c(
        dirname(normalizePath(script_path, mustWork = FALSE)),
        candidates
      ))
    }
  }
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "R")) || dir.exists(file.path(cand, "uploads"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  candidates[1]
}

#' Return the path to the bundled sample Cas9 guide bundle CSV.
#' @export
pjd_sample_guide_bundle_path <- function(root = NULL) {
  resolved <- .pjd_sample_resolve_root(root)
  file.path(resolved, "uploads", SAMPLE_GUIDE_BUNDLE_FILENAME)
}

#' Load and validate the bundled sample Cas9 guide bundle.
#' @export
pjd_load_sample_guide_bundle <- function(root = NULL) {
  csv_path <- pjd_sample_guide_bundle_path(root)
  if (!file.exists(csv_path)) {
    stop(
      "Sample guide bundle CSV not found at '", csv_path, "'. ",
      "Ensure 'uploads/", SAMPLE_GUIDE_BUNDLE_FILENAME, "' exists under the project root.",
      call. = FALSE
    )
  }
  if (!exists("read_guide_bundle", mode = "function")) {
    resolved <- .pjd_sample_resolve_root(root)
    pjat_r <- file.path(resolved, "R", "pjat_plasmids.R")
    gb_r   <- file.path(resolved, "R", "guide_bundle.R")
    if (file.exists(pjat_r) && !exists("pjat_validate_spacer", mode = "function")) {
      source(pjat_r, local = FALSE)
    }
    if (file.exists(gb_r)) {
      source(gb_r, local = FALSE)
    }
  }
  if (!exists("read_guide_bundle", mode = "function")) {
    stop(
      "Function `read_guide_bundle()` is not available; source 'R/guide_bundle.R' first.",
      call. = FALSE
    )
  }
  read_guide_bundle(csv_path, strict = TRUE)
}

#' Return the path to the bundled sample donor FASTA file.
#' @export
pjd_sample_donor_fasta_path <- function(root = NULL) {
  resolved <- .pjd_sample_resolve_root(root)
  file.path(resolved, "uploads", SAMPLE_DONOR_FASTA_FILENAME)
}

#' Load and validate the bundled sample donor FASTA.
#' @export
pjd_load_sample_donor_fasta <- function(root = NULL) {
  fasta_path <- pjd_sample_donor_fasta_path(root)
  if (!file.exists(fasta_path)) {
    stop(
      "Sample donor FASTA not found at '", fasta_path, "'. ",
      "Ensure 'uploads/", SAMPLE_DONOR_FASTA_FILENAME, "' exists under the project root.",
      call. = FALSE
    )
  }

  lines <- readLines(fasta_path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) {
    stop("Sample donor FASTA is empty: ", fasta_path, call. = FALSE)
  }

  hdr_idx <- which(substr(lines, 1L, 1L) == ">")
  if (length(hdr_idx) == 0L || hdr_idx[1L] != 1L) {
    stop(
      "Sample donor FASTA must start with a '>' header line: ", fasta_path,
      call. = FALSE
    )
  }
  if (length(hdr_idx) > 1L) {
    stop(
      "Sample donor FASTA must contain a single FASTA record (found ",
      length(hdr_idx), " headers): ", fasta_path,
      call. = FALSE
    )
  }

  header <- sub("^>\\s*", "", lines[1L])
  if (!nzchar(header)) {
    stop("Sample donor FASTA has an empty '>' header: ", fasta_path, call. = FALSE)
  }

  seq_lines <- lines[-1L]
  if (length(seq_lines) == 0L) {
    stop("Sample donor FASTA contains no sequence lines: ", fasta_path, call. = FALSE)
  }

  raw_seq <- toupper(paste(gsub("\\s+", "", seq_lines), collapse = ""))
  if (!nzchar(raw_seq)) {
    stop("Sample donor FASTA sequence is empty: ", fasta_path, call. = FALSE)
  }
  if (grepl("[^ACGTN]", raw_seq)) {
    bad_chars <- unique(strsplit(gsub("[ACGTN]", "", raw_seq), "")[[1L]])
    stop(
      "Sample donor FASTA contains invalid non-ACGTN nucleotide character(s) [",
      paste(bad_chars, collapse = ", "), "] in: ", fasta_path,
      call. = FALSE
    )
  }

  len_bp <- as.integer(nchar(raw_seq))
  prov <- .pjd_parse_donor_header_provenance(header, fasta_path)

  if (!is.na(prov$start) && !is.na(prov$end)) {
    expected_len <- as.integer(prov$end - prov$start + 1L)
    if (len_bp != expected_len) {
      stop(
        "Sample donor FASTA length (", len_bp, " bp) does not match header coords ",
        prov$coords, " (expected ", expected_len, " bp) in: ", fasta_path,
        call. = FALSE
      )
    }
  }

  list(
    header     = header,
    sequence   = raw_seq,
    length_bp  = len_bp,
    provenance = prov
  )
}

.pjd_parse_donor_header_provenance <- function(header, fasta_path = NA_character_) {
  parts <- trimws(strsplit(header, "\\|", fixed = FALSE)[[1L]])
  id <- if (length(parts) >= 1L) parts[1L] else header
  kv <- list()
  if (length(parts) > 1L) {
    for (p in parts[-1L]) {
      eq <- regexpr("=", p, fixed = TRUE)
      if (eq > 1L) {
        kv[[trimws(substr(p, 1L, eq - 1L))]] <- trimws(substr(p, eq + 1L, nchar(p)))
      }
    }
  }
  coords <- if (!is.null(kv$coords)) kv$coords else NA_character_
  start_pos <- NA_integer_
  end_pos   <- NA_integer_
  if (!is.na(coords) && grepl("^[0-9]+\\.\\.[0-9]+$", coords)) {
    bounds <- as.integer(strsplit(coords, "\\.\\.")[[1L]])
    start_pos <- bounds[1L]
    end_pos   <- bounds[2L]
  }
  list(
    raw     = header,
    id      = id,
    source  = if (!is.null(kv$source)) kv$source else NA_character_,
    addgene = if (!is.null(kv$addgene)) kv$addgene else NA_character_,
    coords  = coords,
    start   = start_pos,
    end     = end_pos,
    path    = normalizePath(fasta_path, mustWork = FALSE)
  )
}
