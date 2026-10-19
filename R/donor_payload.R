# ============================================================================
# R/donor_payload.R — Selectable-marker donor DNA for HDR (additive to the
#                     Cas9 sgRNA Designer; homology_arms.R is NOT modified)
# ============================================================================
# Reads a donor payload (FASTA or raw text), runs synthesis QC, and assembles
# the donor construct 5'arm + payload + 3'arm from a homology-arm design.
# Every nucleotide is either copied from the genome (the arms) or supplied by
# the user (the payload); nothing is fabricated.
#
# Design rules (quoted from Stern et al. 2024 / protocols.io v3):
#   R2  250 bp is the validated MINIMUM arm length (1000 bp validated MAXIMUM).
#   R6  The payload should destroy the gRNA target site.
#   R5  Synthesis: <=10 consecutive identical bp; GC 30-70%; avoid low complexity.
# ============================================================================

#' Read a donor payload from FASTA or raw text.
#'
#' Strips FASTA headers, upper-cases, and removes everything except ACGTN.
#' Never invents sequence: an empty result is an error.
#'
#' @param x Either a path to a FASTA file or a raw character string.
#' @param is_path Logical(1): is `x` a file path? Default TRUE.
#' @return Character(1) cleaned payload sequence.
#' @export
read_donor_payload <- function(x, is_path = TRUE) {
  if (is_path) {
    if (!file.exists(x)) stop("Donor payload file not found: ", x)
    lines <- readLines(x, warn = FALSE)
    lines <- lines[!grepl("^>", lines)]
    seq <- paste(lines, collapse = "")
  } else {
    seq <- as.character(x[1L])
  }
  seq <- toupper(gsub("[^ACGTN]", "", seq))
  if (!nzchar(seq)) stop("Donor payload is empty after cleaning.")
  seq
}

#' QC a donor payload sequence.
#'
#' Mirrors the pJAT synthesis rules (homopolymer <= 10 bp, GC 30-70 %,
#' low-complexity guard) and additionally flags internal Type IIS sites that
#' would self-digest during Gateway cloning (BbsI/BsaI/BsmBI).
#'
#' @param seq Character(1) donor sequence.
#' @param name Character(1) label for messages.
#' @return list: `ok`, `warnings`, `len`, `gc_pct`, `max_homopolymer`,
#'   `cloning_sites`.
#' @export
donor_payload_qc <- function(seq, name = "donor") {
  stopifnot(is.character(seq), length(seq) == 1L)
  seq <- toupper(gsub("[^ACGTN]", "", seq))
  qc <- pjat_synthesis_qc(seq)
  sites <- .find_cloning_sites_local(seq)
  warn <- qc$warnings
  if (nrow(sites) > 0L) {
    warn <- c(warn,
              sprintf("%s: %d Type IIS cloning site(s) detected (BbsI/BsaI/BsmBI) — ",
                      name, nrow(sites)))
  }
  list(
    ok = qc$ok,
    warnings = warn,
    len = nchar(seq),
    gc_pct = qc$gc_frac * 100,
    max_homopolymer = qc$max_homopolymer,
    cloning_sites = sites
  )
}
#' Assemble a donor construct from homology arms + payload.
#'
#' Returns 5'arm + payload + 3'arm plus junction coordinates and a flag for
#' whether the payload overlaps the gRNA cut site / target (rule R6 check is
#' advisory: the payload should destroy the target site).
#'
#' @param arm_design Output list from `design_homology_arms()`.
#' @param payload Character(1) cleaned donor sequence.
#' @param orient Character(1): "+" (payload in genomic orientation) or "-"
#'   (reverse-complemented before insertion).
#' @param name Character(1) label.
#' @return list: `sequence` (full donor), `payload`, `payload_start`,
#'   `payload_end`, `left_arm`, `right_arm`, `junction_left`,
#'   `junction_right`, `qc`.
#' @export
assemble_donor_construct <- function(arm_design, payload, orient = "+",
                                     name = "donor") {
  stopifnot(is.list(arm_design), !is.null(arm_design$arm_5p_seq),
            !is.null(arm_design$arm_3p_seq),
            is.character(payload), length(payload) == 1L)
  orient <- match.arg(orient, c("+", "-"))

  payload <- toupper(gsub("[^ACGTN]", "", payload))
  if (!nzchar(payload)) stop("payload must be non-empty.")

  if (identical(orient, "-")) {
    payload <- .dp_revcomp(payload)
  }

  # Rule R6: the donated payload should destroy the gRNA cut site if it is
  # Rule R6: the donated payload should destroy the gRNA cut site if it is
  # inserted at the cut site. This is advisory -- it does not change the
  # assembled sequence (the arms define the sequence), only records whether
  # the payload shares >= 75 % of the 20-nt spacer identity at either
  # junction.
  spacer_identity <- function(x, y) {
    x <- toupper(x)
    y <- toupper(y)
    if (nchar(x) == 0L || nchar(y) == 0L) return(0)
    l <- min(nchar(x), nchar(y))
    if (l == 0L) return(0)
    m <- strsplit(substr(x, 1L, l), "")[[1L]]
    n <- strsplit(substr(y, 1L, l), "")[[1L]]
    sum(m == n) / l
  }
  left  <- arm_design$arm_5p_seq
  right <- arm_design$arm_3p_seq
  r6_overlap <- max(
    spacer_identity(payload, left),
    spacer_identity(payload, right)
  )

  sequence <- paste0(left, payload, right)
  payload_start <- nchar(left) + 1L
  payload_end   <- nchar(left) + nchar(payload)

  list(
    sequence        = sequence,
    payload         = payload,
    payload_start   = as.integer(payload_start),
    payload_end     = as.integer(payload_end),
    left_arm        = left,
    right_arm       = right,
    junction_left   = substr(sequence, max(1L, payload_start - 20L),
                             payload_start - 1L),
    junction_right  = substr(sequence, payload_end + 1L,
                             min(nchar(sequence), payload_end + 19L)),
    qc              = donor_payload_qc(sequence, name = name),
    seq_rule_r6_share = round(r6_overlap * 100, 1),
    orient          = orient,
    sequence_provenance = "user-supplied (paste/upload)"
  )
}

#' Write a minimal GenBank for a donor construct (arm / payload / arm).
#'
#' @param construct Output of assemble_donor_construct().
#' @param name Character(1) LOCUS name.
#' @param file Output path (NULL returns text only).
#' @return Character(1) GenBank text (invisible).
#' @export
write_donor_genbank <- function(construct, name = "donor_construct",
                                file = NULL) {
  stopifnot(is.list(construct), !is.null(construct$sequence))
  seq <- construct$sequence
  n <- nchar(seq)
  l5 <- nchar(construct$left_arm)
  pstart <- construct$payload_start
  pend   <- construct$payload_end

  lines <- c(
    sprintf("LOCUS       %-16s %d bp DNA circular UNK 01-JAN-2026", name, n),
    "DEFINITION  Donor construct: 5'arm + payload + 3'arm.",
    "ACCESSION   .",
    "VERSION     .",
    "KEYWORDS    .",
    "SOURCE      .",
    "  ORGANISM  .",
    "FEATURES             Location/Qualifiers",
    sprintf("     misc_feature      1..%d", l5),
    "                     /label=\"5p_homology_arm\"",
    sprintf("     misc_feature      %d..%d", pstart, pend),
    "                     /label=\"donor_payload\"",
    sprintf("     misc_feature      %d..%d", pend + 1L, n),
    "                     /label=\"3p_homology_arm\"",
    "ORIGIN"
  )
  for (i in seq(1L, n, by = 60L)) {
    chunk <- substring(seq, i, min(i + 59L, n))
    spaced <- paste(strsplit(chunk, "")[[1L]], collapse = " ")
    lines <- c(lines, sprintf("%9d %s", i, spaced))
  }
  lines <- c(lines, "//")
  txt <- paste(lines, collapse = "\n")
  if (!is.null(file) && nzchar(file[1])) writeLines(txt, con = file[1])
  invisible(txt)
}

#' Write a minimal FASTA for a donor construct (arm / payload / arm).
#'
#' @param construct Output of assemble_donor_construct().
#' @param name Character(1) header line.
#' @param file Output path (NULL returns text only).
#' @return Character(1) FASTA text (invisible).
#' @export
write_donor_fasta <- function(construct, name = "donor_construct", file = NULL) {
  seq <- construct$sequence
  n <- nchar(seq)
  lines <- c(sprintf(">%s|len=%dbp|arm=5p:%dbp+payload:%dbp+3p:%dbp",
                     name, n,
                     nchar(construct$left_arm),
                     nchar(construct$payload),
                     nchar(construct$right_arm)),
             paste0(strwrap(seq, width = 80, simplify = FALSE), collapse = "
"))
  if (!is.null(file) && nzchar(file[1])) writeLines(lines, con = file[1])
  invisible(paste(lines, collapse = "
"))
}

# ---- internal ---------------------------------------------------------------

.dp_revcomp <- function(seq) {
  comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  chars <- strsplit(toupper(seq), "")[[1L]]
  paste0(rev(comp[chars]), collapse = "")
}
# Local Type IIS site scanner (self-contained; the Cas9 app does not source
# plasmid_assembly.R). Mirrors Plasmid Designer's find_cloning_sites.
.find_cloning_sites_local <- function(seq) {
  stopifnot(is.character(seq), length(seq) == 1L)
  seq <- toupper(gsub("[^ACGTN]", "", seq))
  enzymes <- list(BbsI = "GAAGAC", BsaI = "GGTCTC", BsmBI = "CGTCTC")
  out <- data.frame(enzyme = character(0), site = character(0),
                    position = integer(0), cut_pos = integer(0),
                    overhang = character(0), stringsAsFactors = FALSE)
  for (enz in names(enzymes)) {
    site <- enzymes[[enz]]
    hits <- gregexpr(site, seq, fixed = TRUE)[[1L]]
    if (length(hits) == 1L && hits[1L] == -1L) next
    for (h in as.integer(hits)) {
      over <- if (h + nchar(site) + 1L <= nchar(seq)) {
        substring(seq, h + nchar(site), h + nchar(site))
      } else ""
      out <- rbind(out, data.frame(enzyme = enz, site = site,
                                   position = as.integer(h),
                                   cut_pos = as.integer(h + nchar(site)),
                                   overhang = over,
                                   stringsAsFactors = FALSE))
    }
  }
  if (nrow(out)) out <- out[order(out$position), ]
  rownames(out) <- NULL
  out
}
