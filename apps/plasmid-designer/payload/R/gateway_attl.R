# ============================================================================
# R/gateway_attl.R — Gateway attL donor-fragment templates
# ============================================================================
# The Janelia Atalanta protocol builds two synthesized fragments per construct:
#   Fragment 1: attL3 + 5' homology arm + tRNA-sgRNA-tRNA + attL4
#   Fragment 2: attL5 + 3' homology arm + attL6
# and clones them into a pJAT destination vector by dual-Gateway LR Clonase:
#   attL x attR -> attB  (so pJAT plasmids carry attR3-6, never attL).
#
# This module supplies the published attL3-6 / tRNA sequences (see
# donor_templates/gateway_attl_templates.yml) so pjat_assemble_fragments() can
# produce READY, synthesizable fragments instead of ARM_ONLY placeholders.
# No network access: everything is read from the bundled YAML.
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
GATEWAY_ATTL_VERSION <- "1.0.0"

GATEWAY_ATTL_LENGTHS <- c(attl3 = 96L, attl4 = 97L, attl5 = 97L,
                          attl6 = 97L, trna = 72L)

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a) && length(a) == 1L && !is.na(a)) a else b
}

.gateway_attl_root <- function() {
  candidates <- unique(c(
    getwd(),
    app_root()
  ))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    sp <- sub("^--file=", "", file_arg[1L])
    if (file.exists(sp)) {
      d <- dirname(normalizePath(sp, mustWork = FALSE))
      candidates <- unique(c(d, dirname(d), candidates))
    }
  }
  for (cand in candidates) {
    if (file.exists(file.path(cand, "donor_templates",
                              "gateway_attl_templates.yml"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  candidates[1L]
}

#' Path to the bundled attL template YAML.
#' @export
gateway_attl_yml <- function(root = .gateway_attl_root()) {
  file.path(root, "donor_templates", "gateway_attl_templates.yml")
}

#' Load and validate the published attL3-6 + tRNA templates.
#'
#' @param yml Path to the template YAML.
#' @param strict Validate sequence alphabet and expected lengths (default TRUE).
#' @return list with attl3, attl4, attl5, attl6, trna, trna_3p_full_78bp,
#'   cassette_parts and provenance.
#' @export
pjat_load_attl_templates <- function(yml = gateway_attl_yml(), strict = TRUE) {
  if (!file.exists(yml)) {
    stop("Gateway attL template file not found: ", yml,
         "\nRestore donor_templates/gateway_attl_templates.yml.")
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read ", yml)
  }
  y <- yaml::read_yaml(yml)

  required <- c("attl3", "attl4", "attl5", "attl6", "trna")
  missing_keys <- required[vapply(required, function(k) {
    is.null(y[[k]]) || !nzchar(as.character(y[[k]]))
  }, logical(1))]
  if (length(missing_keys)) {
    stop("Gateway template file is missing required key(s): ",
         paste(missing_keys, collapse = ", "))
  }
  if (is.null(y$provenance) || is.null(y$provenance$url)) {
    stop("Gateway template file must document provenance (source + url).")
  }

  out <- list(
    version = as.character(y$version %||% NA_character_),
    label = as.character(y$label %||% "template from published protocol"),
    attl3 = toupper(gsub("[^ACGTN]", "", as.character(y$attl3))),
    attl4 = toupper(gsub("[^ACGTN]", "", as.character(y$attl4))),
    attl5 = toupper(gsub("[^ACGTN]", "", as.character(y$attl5))),
    attl6 = toupper(gsub("[^ACGTN]", "", as.character(y$attl6))),
    trna = toupper(gsub("[^ACGTN]", "", as.character(y$trna))),
    trna_3p_full_78bp = if (is.null(y$trna_3p_full_78bp)) NULL else
      toupper(gsub("[^ACGTN]", "", as.character(y$trna_3p_full_78bp))),
    attl4_forward = if (is.null(y$attl4_forward)) NULL else
      toupper(gsub("[^ACGTN]", "", as.character(y$attl4_forward))),
    attl6_forward = if (is.null(y$attl6_forward)) NULL else
      toupper(gsub("[^ACGTN]", "", as.character(y$attl6_forward))),
    cassette_parts = y$cassette_parts,
    provenance = y$provenance
  )

  if (isTRUE(strict)) {
    for (k in names(GATEWAY_ATTL_LENGTHS)) {
      expected <- GATEWAY_ATTL_LENGTHS[[k]]
      actual <- nchar(out[[k]])
      if (!identical(as.integer(actual), as.integer(expected))) {
        stop(sprintf("Gateway template '%s' is %d bp, expected %d bp (%s)",
                     k, actual, expected, yml))
      }
    }
    if (!isTRUE(out$provenance$verified)) {
      warning("Gateway attL templates are not marked verified in ", yml)
    }
  }
  out
}

#' Named list ready to pass into pjat_assemble_fragments().
#' @export
pjat_attl_template_sequences <- function(templates = pjat_load_attl_templates()) {
  list(
    attl3 = templates$attl3,
    attl4 = templates$attl4,
    attl5 = templates$attl5,
    attl6 = templates$attl6,
    trna  = templates$trna
  )
}

#' Short, UI-ready explanation of attL vs attR.
#' @export
pjat_attl_explanation <- function() {
  paste0("Donor fragments carry attL3/attL4 (fragment 1) and attL5/attL6 ",
         "(fragment 2); the pJAT destination backbone carries attR3-attR6. ",
         "LR Clonase converts attL x attR to attB, so the attL sites are ",
         "consumed by cloning and are absent from finished pJAT plasmids.")
}

#' Extract the four attR sequences (124 bp each) from a parsed pJAT record.
#'
#' The expression vector carries attR, not attL; showing both sides makes the
#' pairing explicit in the UI.
#'
#' @param gb Parsed GenBank record from read_pjat_genbank().
#' @return data.frame(att, start, end, strand, length_bp, sequence).
#' @export
pjat_extract_backbone_attr_sequences <- function(gb) {
  if (is.null(gb) || is.null(gb$sequence) || is.null(gb$features)) {
    stop("gb must be a parsed GenBank record (sequence + features).")
  }
  if (!exists("find_pjat_recombination_sites", mode = "function")) {
    stop("Source R/pjat_plasmids.R before calling this function.")
  }
  att <- find_pjat_recombination_sites(gb)
  if (nrow(att) == 0L) return(att)

  seqs <- vapply(seq_len(nrow(att)), function(i) {
    s <- substr(gb$sequence, att$start[i], att$end[i])
    if (identical(att$strand[i], "-")) {
      paste(rev(strsplit(chartr("ACGTN", "TGCAN", s), "")[[1L]]), collapse = "")
    } else {
      s
    }
  }, character(1))

  out <- att[, intersect(c("att", "start", "end", "strand"), names(att)),
             drop = FALSE]
  out$length_bp <- nchar(seqs)
  out$sequence <- seqs
  out
}

