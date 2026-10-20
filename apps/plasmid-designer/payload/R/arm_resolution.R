# ============================================================================
# R/arm_resolution.R — Homology-arm resolution for donor-fragment assembly
# ============================================================================
# THE FIX for the "dead bundle arms" defect.  Previously the attL donor
# fragments were built from plasmid-backbone flanks extracted around the
# insertion point, which violates Atalanta rules R1/R4: the arms must be
# GENOMIC sequence that flanks the Cas9 cut site, in genomic orientation.
#
# This module resolves arms from the correct source, in priority order:
#   1. "bundle"  — arm_5p_seq / arm_3p_seq carried by the imported
#                  cas9_guide_bundle (the CRISPR Designer's genomic arms)
#   2. "flanks"  — plasmid flanks around the insertion point (fallback, used
#                  only for pasted spacers where no genomic arms exist)
#
# Rule R2 window (250-1000 bp) is enforced on BOTH paths: arms shorter than
# 250 bp are rejected; arms longer than 1000 bp are trimmed from their OUTER
# end so the bases abutting the cut site are preserved exactly.
#
# Requires: R/guide_bundle.R (bundle shape), R/pjat_plasmids.R (R2 constants
# live in the assembly module). Sourced after plasmid_assembly.R.
# ============================================================================

ARM_MIN_BP <- 250L
ARM_MAX_BP <- 1000L
ARM_POLICY <- "250-1000 bp validated homology arms"

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) == 1L && !is.na(a)) a else b
  }
}

#' Trim an arm to the R2 ceiling, preserving the cut-site-adjacent end.
#'
#' @param seq Character(1) arm sequence.
#' @param side "5p" (keep the 3' end) or "3p" (keep the 5' end).
#' @param max_bp Integer ceiling (default 1000).
#' @return list(sequence, trimmed, original_len, length_bp).
#' @export
arm_trim_to_window <- function(seq, side = c("5p", "3p"), max_bp = ARM_MAX_BP) {
  side <- match.arg(side)
  seq <- toupper(gsub("[^ACGTN]", "", as.character(seq[1L])))
  n <- nchar(seq)
  if (n <= max_bp) {
    return(list(sequence = seq, trimmed = FALSE,
                original_len = n, length_bp = n))
  }
  # 5' arm abuts the cut site at its 3' end; 3' arm abuts it at its 5' end.
  kept <- if (identical(side, "5p")) {
    substr(seq, n - max_bp + 1L, n)
  } else {
    substr(seq, 1L, max_bp)
  }
  list(sequence = kept, trimmed = TRUE, original_len = n, length_bp = nchar(kept))
}

#' Resolve homology arms for donor-fragment assembly.
#'
#' @param bundle Optional cas9_guide_bundle from read_guide_bundle(), or a
#'   one-row guide data.frame carrying arm_5p_seq/arm_3p_seq.
#' @param guide_id Optional guide id to select when the bundle has >1 row.
#' @param flank_5p,flank_3p Plasmid flank sequences (fallback source).
#' @param source "auto", "bundle", or "flanks".
#' @param trim_to Integer ceiling for trimming (default 1000).
#' @return list(left_arm, right_arm, source, provenance, trimmed, lengths,
#'   ok, messages).
#' @export
resolve_homology_arms <- function(bundle = NULL, guide_id = NULL,
                                  flank_5p = NULL, flank_3p = NULL,
                                  source = c("auto", "bundle", "flanks"),
                                  trim_to = ARM_MAX_BP) {
  source <- match.arg(source)
  messages <- character(0)

  # ---- pull arms out of the bundle -----------------------------------------
  bundle_arms <- NULL
  if (!is.null(bundle)) {
    g <- if (is.data.frame(bundle)) bundle else bundle$guides
    if (is.data.frame(g) && nrow(g) > 0L &&
        all(c("arm_5p_seq", "arm_3p_seq") %in% names(g))) {
      row <- if (!is.null(guide_id) && guide_id %in% g$guide_id) {
        g[g$guide_id == guide_id, , drop = FALSE]
      } else {
        g[1L, , drop = FALSE]
      }
      a5 <- toupper(gsub("[^ACGTN]", "", as.character(row$arm_5p_seq[1])))
      a3 <- toupper(gsub("[^ACGTN]", "", as.character(row$arm_3p_seq[1])))
      if (nzchar(a5) && nzchar(a3)) {
        bundle_arms <- list(a5 = a5, a3 = a3,
                            guide_id = as.character(row$guide_id[1] %||% NA))
      }
    }
  }

  use_bundle <- switch(source,
    bundle = TRUE,
    flanks = FALSE,
    auto   = !is.null(bundle_arms)
  )

  if (use_bundle && is.null(bundle_arms)) {
    stop("source = 'bundle' but the bundle carries no usable arm_5p_seq/arm_3p_seq.")
  }

  if (use_bundle) {
    left <- bundle_arms$a5; right <- bundle_arms$a3
    prov <- paste0("imported cas9_guide_bundle arms (guide ",
                   bundle_arms$guide_id, ")")
    src <- "bundle"
  } else {
    left  <- toupper(gsub("[^ACGTN]", "", as.character(flank_5p %||% "")))
    right <- toupper(gsub("[^ACGTN]", "", as.character(flank_3p %||% "")))
    if (!nzchar(left) || !nzchar(right)) {
      stop("No homology arms available: supply a guide bundle with arms, or ",
           "plasmid flanks for the fallback path.")
    }
    prov <- "plasmid flanks around the insertion point (fallback; no genomic arms)"
    src <- "flanks"
  }

  l5 <- arm_trim_to_window(left, "5p", trim_to)
  l3 <- arm_trim_to_window(right, "3p", trim_to)
  if (l5$trimmed) messages <- c(messages, sprintf("5' arm trimmed %d -> %d bp",
                                                  l5$original_len, l5$length_bp))
  if (l3$trimmed) messages <- c(messages, sprintf("3' arm trimmed %d -> %d bp",
                                                  l3$original_len, l3$length_bp))

  if (l5$length_bp < ARM_MIN_BP || l3$length_bp < ARM_MIN_BP) {
    stop(sprintf("arms below the validated minimum (%d bp): 5'=%d, 3'=%d. ",
                 ARM_MIN_BP, l5$length_bp, l3$length_bp),
         "Increase arm length (protocols.io v3, rule R2).")
  }

  list(
    left_arm = l5$sequence,
    right_arm = l3$sequence,
    source = src,
    provenance = prov,
    trimmed = l5$trimmed || l3$trimmed,
    lengths = c(arm_5p = l5$length_bp, arm_3p = l3$length_bp),
    policy = ARM_POLICY,
    ok = TRUE,
    messages = messages
  )
}
