# ============================================================================
# R/marker_sequences.R — Offline Selectable Marker Catalog & Sequence Loader
# ============================================================================
# Pure, offline module that reads donor_templates/marker_cassettes.yml and the
# bundled FASTA files under donor_templates/marker_sequences/.
#
# Never initiates network connections. All sequences on disk are produced by
# tools/fetch_marker_sequences.R with full provenance (local pJAT GenBank
# extraction coordinates or NCBI E-utilities accession + coordinates).
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
MARKER_SEQUENCES_VERSION <- "1.0.0"

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) > 0L && all(!is.na(a))) a else b
  }
}

.marker_scalar_or <- function(x, default) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) {
    return(default)
  }
  x[[1L]]
}

#' Resolve the Plasmid_Designer project root.
#' @keywords internal
.marker_sequences_project_root <- function() {
  candidates <- unique(c(
    getwd(),
    app_root()
  ))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1L])
    if (file.exists(script_path)) {
      candidates <- unique(c(
        dirname(normalizePath(script_path, mustWork = FALSE)),
        dirname(dirname(normalizePath(script_path, mustWork = FALSE))),
        candidates
      ))
    }
  }
  for (cand in candidates) {
    if (file.exists(file.path(cand, "donor_templates", "marker_cassettes.yml")) ||
        dir.exists(file.path(cand, "R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  candidates[1L]
}

PJD_MARKER_ROOT <- .marker_sequences_project_root()

#' Path to donor_templates/marker_cassettes.yml.
#' @param root Project root directory.
#' @return Character(1) file path.
#' @export
marker_cassettes_yml <- function(root = PJD_MARKER_ROOT) {
  file.path(root, "donor_templates", "marker_cassettes.yml")
}

#' Directory holding bundled marker FASTA files.
#' @param root Project root directory.
#' @return Character(1) directory path.
#' @export
marker_sequences_dir <- function(root = PJD_MARKER_ROOT) {
  file.path(root, "donor_templates", "marker_sequences")
}

#' Read and validate a single-record marker FASTA file.
#'
#' @param path Path to a FASTA file.
#' @return list(header, sequence, length_bp, md5).
#' @export
pjd_read_fasta_sequence <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Invalid FASTA file path.")
  }
  if (!file.exists(path)) {
    stop("Marker FASTA file not found: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0L) {
    stop("Marker FASTA file is empty: ", path)
  }
  hdr_idx <- which(startsWith(lines, ">"))
  if (length(hdr_idx) == 0L || hdr_idx[1L] != 1L) {
    stop("Marker FASTA file lacks a leading '>' header line: ", path)
  }
  if (length(hdr_idx) > 1L) {
    stop("Marker FASTA file must contain exactly one FASTA record: ", path)
  }
  header <- sub("^>", "", trimws(lines[1L]))
  raw_seq <- paste(gsub("\\s+", "", lines[-1L]), collapse = "")
  seq <- toupper(raw_seq)
  if (!nzchar(seq)) {
    stop("Marker FASTA file contains an empty sequence: ", path)
  }
  if (grepl("[^ACGTN]", seq)) {
    bad <- unique(strsplit(gsub("[ACGTN]", "", seq), "")[[1L]])
    stop("Marker FASTA contains non-ACGTN nucleotide characters (",
         paste(bad, collapse = ","), "): ", path)
  }
  list(
    header    = header,
    sequence  = seq,
    length_bp = as.integer(nchar(seq)),
    md5       = unname(tools::md5sum(path))
  )
}


#' Write a single-record marker FASTA file wrapped at 80 columns.
#'
#' Header format: `><id>|<source>|<accession>|<length>`
#'
#' @param id Marker ID.
#' @param sequence Nucleotide sequence (ACGTN).
#' @param source Source tag (e.g. `pJAT30.gbk` or `NCBI_nuccore`).
#' @param accession Accession / plasmid identifier.
#' @param path Destination FASTA path.
#' @param width Line wrap width (default 80).
#' @return Invisible list(path, length_bp, md5, header).
#' @export
pjd_write_marker_fasta <- function(id, sequence, source, accession, path,
                                   width = 80L) {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id),
            is.character(sequence), length(sequence) == 1L)
  seq <- toupper(gsub("\\s+", "", sequence))
  if (!nzchar(seq)) stop("Cannot write empty sequence for marker '", id, "'.")
  if (grepl("[^ACGTN]", seq)) {
    stop("Sequence for marker '", id, "' contains non-ACGTN characters.")
  }
  n <- nchar(seq)
  hdr <- sprintf(">%s|%s|%s|%d", id, source, accession, n)
  starts <- seq.int(1L, n, by = as.integer(width))
  ends   <- pmin(starts + as.integer(width) - 1L, n)
  chunks <- substring(seq, starts, ends)

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(hdr, chunks), con = path, useBytes = TRUE)
  list(
    path      = path,
    length_bp = as.integer(n),
    md5       = unname(tools::md5sum(path)),
    header    = sub("^>", "", hdr)
  )
}

#' Load the selectable marker catalog from donor_templates/marker_cassettes.yml.
#'
#' @param root Project root directory containing `donor_templates/marker_cassettes.yml`.
#' @return data.frame with one row per marker.
#' @export
pjd_marker_catalog <- function(root = PJD_MARKER_ROOT) {
  if (is.null(root) || length(root) != 1L || is.na(root) || !nzchar(root)) {
    stop("`root` must be a non-empty directory path.")
  }
  yml <- marker_cassettes_yml(root)
  if (!file.exists(yml)) {
    stop("Marker catalog YAML not found: ", yml)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read ", yml)
  }
  y <- yaml::read_yaml(yml)
  if (is.null(y$markers) || length(y$markers) == 0L) {
    stop("Invalid marker catalog YAML (missing 'markers' list): ", yml)
  }

  rows <- lapply(y$markers, function(m) {
    id   <- as.character(.marker_scalar_or(m$id, NA_character_))
    lbl  <- as.character(.marker_scalar_or(m$label, NA_character_))
    catg <- as.character(.marker_scalar_or(m$category, NA_character_))
    desc <- as.character(.marker_scalar_or(m$description, NA_character_))
    sel  <- as.character(.marker_scalar_or(m$selection, NA_character_))
    src  <- as.character(.marker_scalar_or(m$source, NA_character_))
    st   <- as.character(.marker_scalar_or(m$sequence_status, NA_character_))
    ff   <- as.character(.marker_scalar_or(m$fetch_from, NA_character_))
    sf   <- as.character(.marker_scalar_or(m$sequence_file, NA_character_))
    lbp  <- suppressWarnings(as.integer(.marker_scalar_or(m$length_bp, NA_integer_)))
    md5  <- as.character(.marker_scalar_or(m$md5, NA_character_))
    prov <- as.character(.marker_scalar_or(m$sequence_provenance, NA_character_))
    ef   <- if (is.null(m$expected_features)) character(0) else as.character(unlist(m$expected_features))

    data.frame(
      id                  = id,
      label               = lbl,
      category            = catg,
      description         = desc,
      selection           = sel,
      source              = src,
      sequence_status     = st,
      fetch_from          = ff,
      expected_features   = I(list(ef)),
      sequence_file       = sf,
      length_bp           = lbp,
      md5                 = md5,
      sequence_provenance = prov,
      stringsAsFactors    = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Check whether a selectable marker has a verified bundled FASTA on disk.
#'
#' @param marker_id Character(1) marker identifier (e.g. `"3xP3_dsRed"`).
#' @param root Project root directory.
#' @return Logical(1): TRUE if bundled and the FASTA file exists on disk, FALSE otherwise.
#' @export
pjd_marker_sequence_available <- function(marker_id, root = PJD_MARKER_ROOT) {
  if (is.null(marker_id) || length(marker_id) != 1L || is.na(marker_id) || !nzchar(marker_id)) {
    return(FALSE)
  }
  cat_df <- tryCatch(pjd_marker_catalog(root = root), error = function(e) NULL)
  if (is.null(cat_df) || nrow(cat_df) == 0L) return(FALSE)
  row <- cat_df[cat_df$id == marker_id, , drop = FALSE]
  if (nrow(row) == 0L) return(FALSE)
  row <- row[1L, , drop = FALSE]
  if (!identical(row$sequence_status, "bundled")) return(FALSE)
  if (is.na(row$sequence_file) || !nzchar(row$sequence_file)) return(FALSE)
  fasta_path <- file.path(root, row$sequence_file)
  file.exists(fasta_path) && isTRUE(file.info(fasta_path)$size > 0)
}

#' Retrieve a marker's nucleotide sequence and provenance from the local bundle.
#'
#' @param marker_id Character(1) marker identifier.
#' @param root Project root directory.
#' @return `list(id, label, sequence, length_bp, provenance, status)`.
#' @export
pjd_get_marker_sequence <- function(marker_id, root = PJD_MARKER_ROOT) {
  if (is.null(marker_id) || length(marker_id) != 1L || is.na(marker_id) || !nzchar(marker_id)) {
    stop("`marker_id` must be a non-empty character scalar.")
  }
  cat_df <- pjd_marker_catalog(root = root)
  row <- cat_df[cat_df$id == marker_id, , drop = FALSE]
  if (nrow(row) == 0L) {
    stop("Unknown marker_id '", marker_id, "'. Available IDs: ",
         paste(cat_df$id, collapse = ", "))
  }
  row <- row[1L, , drop = FALSE]

  if (!identical(row$sequence_status, "bundled")) {
    stop("Sequence for marker '", marker_id, "' is not bundled (status = '",
         row$sequence_status, "'; provenance = '", row$sequence_provenance, "').")
  }
  if (is.na(row$sequence_file) || !nzchar(row$sequence_file)) {
    stop("Marker '", marker_id, "' is marked bundled but has no sequence_file in marker_cassettes.yml.")
  }

  fasta_path <- file.path(root, row$sequence_file)
  if (!file.exists(fasta_path)) {
    stop("Bundled FASTA file for marker '", marker_id, "' does not exist: ", fasta_path)
  }

  fa <- pjd_read_fasta_sequence(fasta_path)

  if (!is.na(row$length_bp) && fa$length_bp != row$length_bp) {
    stop("Length mismatch for marker '", marker_id, "': FASTA has ",
         fa$length_bp, " bp but YAML records ", row$length_bp, " bp.")
  }
  if (!is.na(row$md5) && nzchar(row$md5) && !identical(tolower(fa$md5), tolower(row$md5))) {
    stop("MD5 checksum mismatch for marker '", marker_id, "': FASTA has ",
         fa$md5, " but YAML records ", row$md5, ".")
  }

  list(
    id         = as.character(row$id),
    label      = as.character(row$label),
    sequence   = fa$sequence,
    length_bp  = as.integer(fa$length_bp),
    provenance = as.character(row$sequence_provenance),
    status     = as.character(row$sequence_status)
  )
}
