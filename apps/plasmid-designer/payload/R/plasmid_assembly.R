# ============================================================================
# R/plasmid_assembly.R — pJAT vector assembly (no sequences fabricated)
# ============================================================================
# Operations:
#   insert_cassette(backbone_seq, start, end, insert, strand) — replace a
#       window of the backbone with an insert strand-aware.
#   assemble_gRNA_cassette(spacer, scaffold, terminator) — build the
#       gRNA expression cassette from a validated 20-nt spacer.
#   assemble_expression_vector(backbone_gb, spacer, strategy, coords, orient)
#       — full plasmid assembly: cassette inserted at attR4_end+1..attR5_start-1
#       (default), a restriction site, or explicit coordinates; returns the
#       complete circular sequence plus annotation features.
#   find_cloning_sites(seq, enzymes) — enumerate BbsI/BsaI/BsmBI sites with
#       cut offsets so the user can pick one deliberately (never auto-cut).
#   write_plasmid_genbank(seq, features, ...) — minimal circular GenBank.
#   plasmid_qc(...) — length / restriction-site / synthesis QC.
#
# Design invariant: every nucleotide in the output is either copied verbatim
# from the backbone or derived from user-supplied sequence. No marker or
# scaffold sequence is ever fabricated.

GRRNA_SCAFFOLD <- "GTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
GRRNA_TERMINATOR <- "TTTTTT"
CLONING_ENZYMES <- list(
  BbsI  = list(site = "GAAGAC", cut_after = 2L),
  BsaI  = list(site = "GGTCTC", cut_after = 1L),
  BsmBI = list(site = "CGTCTC", cut_after = 1L)
)

#' Insert a cassette into a backbone window, strand-aware.
#'
#' The replacement window [start, end] (1-based inclusive) is excised from
#' `backbone_seq` and `insert` is spliced in.  When `strand == "-"`, the
#' insert is reverse-complemented before insertion so it reads 5'->3' in the
#' backbone's transcriptional orientation.
#'
#' @param backbone_seq Character(1), full plasmid sequence.
#' @param start,end Integer(1), 1-based inclusive window to replace.
#' @param insert Character(1), insert sequence.
#' @param strand Character(1), "+" or "-".
#' @return List: `sequence`, `insert_start`, `insert_end`, `insert_len`,
#'   `backbone_len`, `final_len`, `junction_left`, `junction_right`.
#' @export
insert_cassette <- function(backbone_seq, start, end, insert, strand = "+") {
  stopifnot(is.character(backbone_seq), length(backbone_seq) == 1L,
            is.character(insert), length(insert) == 1L,
            start > 0L, end <= nchar(backbone_seq))
  if (end < start - 1L) {
    stop("invalid insertion range: end (", end, ") < start (", start, ") - 1; ",
         "use end = start - 1 for a pure insertion")
  }
  strand <- match.arg(strand, c("+", "-"))

  backbone_seq <- toupper(gsub("[^ACGTN]", "", backbone_seq))
  insert <- toupper(gsub("[^ACGTN]", "", insert))

  ins_final <- insert
  if (identical(strand, "-")) {
    ins_final <- .paf_revcomp(insert)
  }

  before <- substring(backbone_seq, 1L, start - 1L)
  after  <- substring(backbone_seq, end + 1L, nchar(backbone_seq))
  final  <- paste0(before, ins_final, after)

  list(
    sequence         = final,
    insert_start     = as.integer(start),
    insert_end       = as.integer(start + nchar(ins_final) - 1L),
    insert_len       = as.integer(nchar(ins_final)),
    backbone_len     = as.integer(nchar(backbone_seq)),
    final_len        = as.integer(nchar(final)),
    strand           = strand,
    junction_left    = substr(final, max(1L, start - 20L), start - 1L),
    junction_right   = substr(final, start + nchar(ins_final),
                              min(nchar(final), start + nchar(ins_final) + 19L))
  )
}

#' Assemble a SpCas9 gRNA expression cassette from a validated spacer.
#'
#' @param spacer Character(1), exactly 20 nt, PAM excluded (rule R3).
#' @param scaffold Character(1), constant sgRNA handle.
#' @param terminator Character(1), Pol III terminator.
#' @return Character(1): spacer + scaffold + terminator.
#' @export
assemble_gRNA_cassette <- function(spacer,
                                   scaffold = GRRNA_SCAFFOLD,
                                   terminator = GRRNA_TERMINATOR) {
  stopifnot(is.character(spacer), length(spacer) == 1L)
  spacer <- toupper(gsub("[^ACGTN]", "", spacer))
  if (nchar(spacer) != 20L) {
    stop("spacer must be exactly 20 nt (PAM excluded, rule R3); got ", nchar(spacer), " nt")
  }
  paste0(spacer, toupper(scaffold), toupper(terminator))
}

#' Find Type IIS cloning sites in a sequence with cut offsets.
#'
#' Returns one row per hit: enzyme, match position (1-based), cut position
#' (1-based, the first base *after* cleavage on the + strand), and the
#' overhang sequence produced by cleavage.  Never selects a site automatically
#' — the user chooses.
#'
#' @param seq Character(1), double-stranded sequence to scan (both strands).
#' @param enzymes Named list overriding `CLONING_ENZYMES`.
#' @return data.frame.
#' @export
find_cloning_sites <- function(seq, enzymes = CLONING_ENZYMES) {
  stopifnot(is.character(seq), length(seq) == 1L)
  seq <- toupper(gsub("[^ACGTN]", "", seq))
  out <- data.frame(
    enzyme = character(0), site = character(0), position = integer(0),
    cut_pos = integer(0), overhang = character(0), stringsAsFactors = FALSE
  )
  rc <- .paf_revcomp(seq)
  for (enz in names(enzymes)) {
    e <- enzymes[[enz]]
    site <- e$site; cut_after <- as.integer(e$cut_after)
    hits <- .paf_find_all(seq, site)
    for (h in hits) {
      over <- if (h + nchar(site) + cut_after - 1L <= nchar(seq)) {
        substring(seq, h + nchar(site), h + nchar(site) + cut_after - 1L)
      } else ""
      out <- rbind(out, data.frame(
        enzyme = enz, site = site, position = as.integer(h),
        cut_pos = as.integer(h + nchar(site) + cut_after - 1L),
        overhang = over, stringsAsFactors = FALSE))
    }
    hits_rc <- .paf_find_all(rc, site)
    for (h in hits_rc) {
      out <- rbind(out, data.frame(
        enzyme = enz, site = paste0(site, " (rc)"), position = as.integer(h),
        cut_pos = as.integer(h + nchar(site) + cut_after - 1L),
        overhang = "", stringsAsFactors = FALSE))
    }
  }
  if (nrow(out)) out <- out[order(out$position), ]
  rownames(out) <- NULL
  out
}

.paf_find_all <- function(haystack, needle) {
  if (nchar(needle) == 0L || nchar(needle) > nchar(haystack)) return(integer(0))
  pos <- gregexpr(needle, haystack, fixed = TRUE)[[1L]]
  if (length(pos) == 1L && pos[1L] == -1L) return(integer(0))
  as.integer(pos)
}

.paf_revcomp <- function(seq) {
  stopifnot(is.character(seq), length(seq) == 1L)
  comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  chars <- strsplit(toupper(seq), "")[[1L]]
  paste0(rev(comp[chars]), collapse = "")
}

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && length(a)) a else b

#' Assemble a complete pJAT expression vector.
#'
#' Reads a parsed pJAT GenBank record, locates the four attR sites, and inserts
#' a user-supplied cassette (default: a gRNA expression cassette built from a
#' validated spacer) into a chosen window:
#'   * "linker"   — attR4_end+1 .. attR5_start-1 (the guide-cassette window,
#'                  default; 158 bp in pJAT74/pJAT32)
#'   * "enzyme"   — a Type IIS site returned by find_cloning_sites(); the user
#'                  must confirm the exact position; never auto-selected
#'   * "custom"   — explicit start/end coordinates
#'
#' @param backbone_gb Parsed pJAT GenBank from read_pjat_genbank().
#' @param insert Character(1), insert sequence (or NULL to build a gRNA
#'   cassette from `spacer`).
#' @param spacer Character(1), 20-nt PAM-excluded spacer used only when
#'   insert is NULL.
#' @param strategy Character(1): "linker", "enzyme", or "custom".
#' @param enzyme Character(1), enzyme name when strategy == "enzyme".
#' @param site Integer(1), 1-based position of the enzyme cut when strategy
#'   == "enzyme" (defaults to the first hit).
#' @param start,end Integer(1), explicit window when strategy == "custom".
#' @param orient Character(1): "+" or "-".
#' @return List: `sequence` (full circular plasmid), `insert` (the final
#'   insert as inserted), `features` (data.frame of annotation features),
#'   `assembly` (summary), `qc` (synthesis QC of the insert).
#' @export
assemble_expression_vector <- function(backbone_gb,
                                       insert = NULL,
                                       spacer = NULL,
                                       strategy = c("linker", "enzyme", "custom"),
                                       enzyme = NULL,
                                       site = NULL,
                                       start = NULL,
                                       end = NULL,
                                       orient = "+") {
  stopifnot(is.list(backbone_gb), !is.null(backbone_gb$sequence))
  strategy <- match.arg(strategy)
  orient <- match.arg(orient, c("+", "-"))

  seq <- backbone_gb$sequence
  att <- find_pjat_recombination_sites(backbone_gb)
  # Build lookup: att -> list(start, end, strand)
  att_map <- lapply(seq_len(nrow(att)), function(i) {
    list(start = att$start[i], end = att$end[i], strand = att$strand[i])
  })
  names(att_map) <- att$att

  resolve_window <- function() {
    if (identical(strategy, "linker")) {
      # window between attR4 end and attR5 start (the cassette window)
      a4 <- att_map[["attR4"]]; a5 <- att_map[["attR5"]]
      if (is.null(a4) || is.null(a5)) {
        stop("backbone lacks attR4/attR5 annotation; cannot use linker strategy")
      }
      s <- a4$end + 1L
      e <- a5$start - 1L
      if (s > e) stop("attR4 and attR5 overlap; linker window is empty")
      list(start = s, end = e, label = "attR4_end+1 .. attR5_start-1 (linker)")
    } else if (identical(strategy, "enzyme")) {
      if (is.null(enzyme) || !(enzyme %in% names(CLONING_ENZYMES))) {
        stop("enzyme must be one of: ", paste(names(CLONING_ENZYMES), collapse = ", "))
      }
      sites <- find_cloning_sites(seq, enzymes = CLONING_ENZYMES[enzyme])
      hits <- sites[!grepl(" \\(rc\\)", sites$site), ]
      if (nrow(hits) == 0L) stop("no ", enzyme, " site found in backbone")
      idx <- if (!is.null(site) && site %in% hits$position) site else hits$position[1L]
      # cut window: 1 bp around the cut site (replace 0 bp = insert only)
      list(start = idx, end = idx - 1L,
           label = sprintf("%s cut at %d", enzyme, idx))
    } else {
      # custom
      if (is.null(start) || is.null(end)) {
        stop("start and end are required for strategy = 'custom'")
      }
      list(start = as.integer(start), end = as.integer(end),
           label = sprintf("custom %d..%d", start, end))
    }
  }

  win <- resolve_window()
  if (win$start > win$end && !identical(strategy, "enzyme")) {
    stop("empty insertion window: ", win$label)
  }

  if (is.null(insert)) {
    if (is.null(spacer)) stop("provide either `insert` or `spacer`")
    insert <- assemble_gRNA_cassette(spacer)
  }

  asm <- insert_cassette(seq, win$start, win$end, insert, strand = orient)

  # Annotation features for the insert.
  feats <- data.frame(
    type = c("insert", "spacer"),
    start = c(asm$insert_start, asm$insert_start),
    end = c(asm$insert_end, asm$insert_start + 19L),
    strand = orient,
    label = c("user_insert", "spacer_20nt"),
    stringsAsFactors = FALSE
  )

  qc <- pjat_synthesis_qc(insert)

  list(
    sequence   = asm$sequence,
    insert     = insert,
    insert_start = asm$insert_start,
    insert_end   = asm$insert_end,
    insert_len = asm$insert_len,
    backbone_len = asm$backbone_len,
    final_len  = asm$final_len,
    strategy   = strategy,
    window     = win,
    orient     = orient,
    features   = feats,
    junction_left  = asm$junction_left,
    junction_right = asm$junction_right,
    qc         = qc,
    assembly_summary = sprintf(
      "%s | %s | insert %d bp at %d..%d | final %d bp",
      backbone_gb$plasmid %||% backbone_gb$locus %||% "pJAT",
      win$label, asm$insert_len, asm$insert_start, asm$insert_end, asm$final_len
    )
  )
}

#' Write a minimal circular plasmid GenBank file.
#'
#' @param seq Character(1), full plasmid sequence.
#' @param features data.frame with columns type, start, end, strand, label.
#' @param name Character(1), plasmid name / LOCUS name.
#' @param topology Character(1): "circular" or "linear".
#' @param file Output path (NULL returns text only).
#' @return Character(1) GenBank text.
#' @export
write_plasmid_genbank <- function(seq, features = NULL,
                                  name = "pJAT_assembly",
                                  topology = c("circular", "linear"),
                                  file = NULL) {
  topology <- match.arg(topology)
  stopifnot(is.character(seq), length(seq) == 1L)
  seq <- toupper(gsub("[^ACGTN]", "", seq))
  n <- nchar(seq)

  lines <- c(
    sprintf("LOCUS       %-16s %d bp DNA %s   UNK 01-JAN-2026",
            name, n, topology),
    "DEFINITION  Plasmid Designer assembly.",
    "ACCESSION   .",
    "VERSION     .",
    "KEYWORDS    .",
    "SOURCE      .",
    "  ORGANISM  .",
    "FEATURES             Location/Qualifiers"
  )
  if (!is.null(features) && nrow(features) > 0L) {
    for (i in seq_len(nrow(features))) {
      f <- features[i, ]
      loc <- if (identical(f$strand, "-")) {
        sprintf("complement(%d..%d)", f$start, f$end)
      } else {
        sprintf("%d..%d", f$start, f$end)
      }
      lines <- c(lines, sprintf("     %-16s%s", f$type, loc))
      lines <- c(lines, sprintf("                     /label=\"%s\"", f$label))
    }
  }
  lines <- c(lines, "ORIGIN")
  for (i in seq(1L, n, by = 60L)) {
    chunk <- substring(seq, i, min(i + 59L, n))
    spaced <- paste(strsplit(chunk, "")[[1L]], collapse = " ")
    lines <- c(lines, sprintf("%9d %s", i, spaced))
  }
  lines <- c(lines, "//")
  txt <- paste(lines, collapse = "\n")
  if (!is.null(file) && nzchar(file[1])) {
    writeLines(txt, con = file[1])
  }
  invisible(txt)
}

#' QC a proposed insert / full plasmid.
#'
#' @param seq Character(1), sequence to QC.
#' @param name Character(1), label for messages.
#' @return list with `ok`, `warnings`, `len`, `gc_pct`, `max_homopolymer`,
#'   `cloning_sites` (data.frame from find_cloning_sites).
#' @export
plasmid_qc <- function(seq, name = "insert") {
  stopifnot(is.character(seq), length(seq) == 1L)
  seq <- toupper(gsub("[^ACGTN]", "", seq))
  qc <- pjat_synthesis_qc(seq)
  sites <- find_cloning_sites(seq)
  warn <- qc$warnings
  if (nrow(sites) > 0L) {
    warn <- c(warn,
              sprintf("%s: %d Type IIS cloning site(s) detected (BbsI/BsaI/BsmBI) — confirm before Gateway cloning.",
                      name, nrow(sites)))
  }
  list(
    ok = qc$ok && length(warn) == 0L,
    warnings = warn,
    len = nchar(seq),
    gc_pct = qc$gc_pct,
    max_homopolymer = qc$max_homopolymer,
    cloning_sites = sites
  )
}
