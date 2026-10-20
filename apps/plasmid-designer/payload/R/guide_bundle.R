# ============================================================================
# R/guide_bundle.R — cas9_guide_bundle/1.0 reader & validator
# ============================================================================
# Contract exported by the Cas9 sgRNA Designer's
#   "Export guides for Plasmid Designer" button:
#
#   Cas9_<gene>_<species>_guides_for_plasmid_designer.csv  (+ .json sidecar)
#
# The bundle is the interop channel: the Plasmid Designer imports a set of
# validated SpCas9 guides (spacer_seq + flanking homology arms) and clones them
# into a pJAT backbone.  This module reads, validates, and normalises bundles
# without trusting them.
#
# Schema version: cas9_guide_bundle/1.0
# Required CSV columns:
#   guide_id, gene_symbol, species, assembly, chr, accession,
#   strand, protospacer_start, protospacer_end, pam_start, pam_end, cut_site,
#   spacer_seq, pam_seq, on_target_score, cfd_specificity,
#   arm_5p_seq, arm_5p_start, arm_5p_end, arm_5p_len,
#   arm_3p_seq, arm_3p_start, arm_3p_end, arm_3p_len,
#   arm_3p_anchor, pam_shielding_summary, guide_bundle_schema
#
# Required JSON sidecar keys:
#   schema, app_version, created_utc, species, assembly, gene_symbol,
#   n_guides, arm_policy
# ============================================================================

GUIDE_BUNDLE_SCHEMA <- "cas9_guide_bundle/1.0"
GUIDE_BUNDLE_SCHEMA_V11 <- "cas9_guide_bundle/1.1"
# Accepted on read (v1.0 legacy + v1.1 current). v1.1 adds the protocol-validated
# arm window as a hard contract: 250 bp minimum AND 1000 bp maximum.
GUIDE_BUNDLE_ACCEPTED_SCHEMAS <- c(GUIDE_BUNDLE_SCHEMA, GUIDE_BUNDLE_SCHEMA_V11)
# Protocol-validated arm window (Stern et al. 2024 / protocols.io v3, rule R2):
# "have not observed a difference in efficiency between arms containing 250 bp
#  and 1000 bp" -> 250 bp floor, 1000 bp ceiling.
GUIDE_BUNDLE_ARM_MIN_BP <- 250L
GUIDE_BUNDLE_ARM_MAX_BP <- 1000L
GUIDE_BUNDLE_ARM_POLICY <- "250-1000 bp validated homology arms"
GUIDE_BUNDLE_REQUIRED_COLS <- c(
  "guide_id", "gene_symbol", "species", "assembly", "chr", "accession",
  "strand", "protospacer_start", "protospacer_end", "pam_start", "pam_end",
  "cut_site", "spacer_seq", "pam_seq", "on_target_score", "cfd_specificity",
  "arm_5p_seq", "arm_5p_start", "arm_5p_end", "arm_5p_len",
  "arm_3p_seq", "arm_3p_start", "arm_3p_end", "arm_3p_len",
  "arm_3p_anchor", "pam_shielding_summary", "guide_bundle_schema"
)

GUIDE_BUNDLE_JSON_REQUIRED_KEYS <- c(
  "schema", "app_version", "created_utc", "species", "assembly",
  "gene_symbol", "n_guides", "arm_policy"
)

#' Read and validate a cas9_guide_bundle/1.0 export.
#'
#' Reads the CSV (and, if present, the matching .json sidecar), validates the
#' schema version, checks every required column, normalises sequences, and
#' recomputes spacer length so callers cannot import a stale or corrupt bundle.
#' Never fabricates sequence: missing spacer_seq or arm sequence is an error.
#'
#' @param csv Path to the exported guides CSV.
#' @param json Optional path to the .json sidecar (defaults to csv with the
#'   \`.csv\` suffix replaced by \`.json\`).
#' @param strict If TRUE (default), the spacer must pass `pjat_validate_spacer()`
#'   (exactly 20 nt, PAM excluded).
#' @return A list: `guides` (data.frame), `bundle_meta` (list),
#'   `sidecar` (list or NULL), `provenance` (character).
#' @export
read_guide_bundle <- function(csv, json = NULL, strict = TRUE) {
  if (!file.exists(csv)) stop("Guide bundle CSV not found: ", csv)
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required to read guide bundles.")
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read guide bundle sidecars.")
  }

  csv <- normalizePath(csv, mustWork = TRUE)
  df <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
  if (!is.data.frame(df) || nrow(df) == 0L) {
    stop("Guide bundle CSV is empty or unreadable: ", csv)
  }

  missing <- setdiff(GUIDE_BUNDLE_REQUIRED_COLS, names(df))
  if (length(missing)) {
    stop(
      "Guide bundle is missing required columns: ",
      paste(missing, collapse = ", "),
      ". Re-export from the Cas9 Designer."
    )
  }

  bundle_schema <- as.character(df$guide_bundle_schema[1L])
  if (!bundle_schema %in% GUIDE_BUNDLE_ACCEPTED_SCHEMAS) {
    stop(
      "Unsupported guide_bundle_schema: '", bundle_schema,
      "' (expected one of: ", paste(GUIDE_BUNDLE_ACCEPTED_SCHEMAS, collapse = ", "),
      "). Re-export from the Cas9 Designer."
    )
  }

  if (is.null(json)) {
    json <- sub("\\.csv$", ".json", csv)
  }
  # Never treat the CSV itself as a sidecar (e.g. extension-less upload paths),
  # and never let a malformed/optional sidecar abort the import.
  sidecar <- NULL
  json_is_csv <- !is.null(json) && identical(normalizePath(json, mustWork = FALSE),
                                             normalizePath(csv, mustWork = FALSE))
  if (!is.null(json) && !json_is_csv && file.exists(json)) {
    # Accept JSON (producer) or YAML (legacy/consumer tests); try JSON first.
    sidecar <- tryCatch(
      {
        sc <- if (requireNamespace("jsonlite", quietly = TRUE)) {
          txt <- paste(readLines(json, warn = FALSE), collapse = "\n")
          if (grepl("^\\s*[\\[{]", txt)) jsonlite::fromJSON(txt, simplifyVector = FALSE)
          else yaml::read_yaml(json)
        } else {
          yaml::read_yaml(json)
        }
        if (!is.list(sc)) {
          warning("Guide bundle sidecar is not a key/value document; ignoring it.",
                  call. = FALSE)
          NULL
        } else {
          sc
        }
      },
      error = function(e) {
        warning("Could not parse guide bundle sidecar '", json, "' (",
                conditionMessage(e), "); importing CSV without provenance.",
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(sidecar)) {
      sm <- as.character(sidecar$schema %||% NA_character_)
      if (!sm %in% GUIDE_BUNDLE_ACCEPTED_SCHEMAS) {
        warning(
          "Sidecar schema '", sm,
          "' does not match an accepted schema (",
          paste(GUIDE_BUNDLE_ACCEPTED_SCHEMAS, collapse = ", "), "); ",
          "proceeding with CSV validation only."
        )
      }
    }
  }

  # Normalise sequences (upper-case, strip non-ACGTN) and recompute lengths.
  df$spacer_seq <- .gbu_clean_dna(df$spacer_seq)
  df$arm_5p_seq <- .gbu_clean_dna(df$arm_5p_seq)
  df$arm_3p_seq <- .gbu_clean_dna(df$arm_3p_seq)

  # Recompute arm lengths from sequence so a hand-edited CSV cannot lie.
  df$arm_5p_len <- nchar(df$arm_5p_seq)
  df$arm_3p_len <- nchar(df$arm_3p_seq)

  bad_spacer <- !nzchar(df$spacer_seq)
  if (any(bad_spacer)) {
    stop(
      "guide bundle has empty spacer_seq for row(s): ",
      paste(df$guide_id[bad_spacer], collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(strict)) {
    for (i in seq_len(nrow(df))) {
      ok <- tryCatch(
        pjat_validate_spacer(df$spacer_seq[i]),
        error = function(e) FALSE
      )
      if (!isTRUE(ok)) {
        stop(
          "spacer_seq for guide '", df$guide_id[i],
          "' is not a valid 20-nt PAM-excluded spacer: ",
          df$spacer_seq[i],
          call. = FALSE
        )
      }
    }
  }

  for (nm in c("arm_5p_seq", "arm_3p_seq")) {
    bad <- !nzchar(df[[nm]])
    if (any(bad)) {
      stop(
        "guide bundle has empty ", nm, " for row(s): ",
        paste(df$guide_id[bad], collapse = ", "),
        call. = FALSE
      )
    }
  }

  bad_arm <- df$arm_5p_len < GUIDE_BUNDLE_ARM_MIN_BP |
             df$arm_3p_len < GUIDE_BUNDLE_ARM_MIN_BP
  if (any(bad_arm)) {
    stop(
      "guide bundle has arms < ", GUIDE_BUNDLE_ARM_MIN_BP,
      " bp (validated minimum) for row(s): ",
      paste(df$guide_id[bad_arm], collapse = ", "),
      call. = FALSE
    )
  }

  # Rule R2 ceiling: arms longer than 1000 bp are outside the validated window.
  long_arm <- df$arm_5p_len > GUIDE_BUNDLE_ARM_MAX_BP |
              df$arm_3p_len > GUIDE_BUNDLE_ARM_MAX_BP
  if (any(long_arm)) {
    stop(
      "guide bundle has arms > ", GUIDE_BUNDLE_ARM_MAX_BP,
      " bp (validated maximum) for row(s): ",
      paste(df$guide_id[long_arm], collapse = ", "),
      ". Trim the arms to <= ", GUIDE_BUNDLE_ARM_MAX_BP,
      " bp before importing (protocols.io v3, rule R2).",
      call. = FALSE
    )
  }

  provenance <- character(nrow(df))
  if (!is.null(sidecar)) {
    provenance <- paste0(
      "cas9_designer v", sidecar$app_version %||% "?",
      " (", sidecar$species %||% "?", ", ",
      sidecar$gene_symbol %||% "?", ")",
      " exported ", sidecar$created_utc %||% "unknown"
    )
  } else {
    provenance <- "user-supplied CSV (no sidecar)"
  }

  cols <- c("guide_id", setdiff(names(df), "guide_id"))
  df <- df[, cols, drop = FALSE]
  rownames(df) <- NULL

  bundle_meta <- list(
    n_guides = nrow(df),
    species = unique(df$species),
    gene_symbol = unique(df$gene_symbol),
    schema = bundle_schema,
    arm_policy = GUIDE_BUNDLE_ARM_POLICY,
    arm_min_bp = GUIDE_BUNDLE_ARM_MIN_BP,
    arm_max_bp = GUIDE_BUNDLE_ARM_MAX_BP,
    imported_from = csv,
    imported_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
  )

  list(
    guides = df,
    bundle_meta = bundle_meta,
    sidecar = sidecar,
    provenance = provenance
  )
}

#' Check whether every guide in a bundle is compatible with a chosen backbone.
#'
#' A guide requires n U6 promoters equal to the number of guides the user
#' wants to express from this plasmid. Advisory only: the user's own backbone
#' selection is never overridden.
#'
#' @param bundle Output of read_guide_bundle().
#' @param catalog Optional pJAT catalog from pjat_load_catalog().
#' @return data.frame: one row per guide with compatible_backbones and
#'   compatible flag.
#' @export
guide_bundle_check_compatible <- function(bundle, catalog = NULL) {
  if (is.null(catalog)) catalog <- pjat_load_catalog()
  df <- bundle$guides
  out <- data.frame(
    guide_id = df$guide_id,
    n_guides_needed = 1L,
    compatible_backbones = character(nrow(df)),
    compatible = logical(nrow(df)),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(df))) {
    rec <- tryCatch(
      pjat_recommend_backbones(1L, catalog),
      error = function(e) data.frame(name = character(0), u6_promoters = integer(0))
    )
    out$compatible_backbones[i] <- paste(rec$name, collapse = ", ")
    out$compatible[i] <- nrow(rec) > 0L
  }
  out
}

# ---- internal helpers -------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && length(a)) a else b

.gbu_clean_dna <- function(x) {
  x <- as.character(x)
  x <- toupper(gsub("[^ACGTN]", "", x))
  x[is.na(x)] <- ""
  x
}
