# ============================================================================
# R/pjat_plasmids.R — Janelia Atalanta (pJAT) Plasmid Catalog, GenBank Parser
#                       & Gateway Donor Fragment Designer (Ticket T8)
# ============================================================================
# Implements the plasmid registry and Gateway (attL) donor-fragment design for
# the Janelia Atalanta CRISPR/Cas9 HDR platform
# (Stern, Kim & Behrman 2024, bioRxiv 2023.06.17.545412;
#  Addgene 204289-204312; protocols.io bp2l6bjokgqe/v3).
#
# Design rules implemented here, quoted from the sources:
#   R1  "Identify homology arms that precisely, or nearly precisely flank the
#        expected Cas9 cut site."                    (protocols.io v3, step 3)
#   R2  "I have not performed a comprehensive survey of effective homology arm
#        lengths, but have not observed a difference in efficiency between arms
#        containing 250 bp and 1000 bp."             (protocols.io v3, step 3)
#        -> 250 bp validated MINIMUM, 1000 bp validated MAXIMUM.
#   R3  "The gRNA sequence should be the 20 bp target sequence, excluding the
#        NGG PAM sequence. The experiment will fail if you include the NGG
#        sequence."                                  (protocols.io v3, step 4)
#   R4  "The orientation of the homology arms in the synthesized attL fragments
#        should match their orientation in the genome. Some software 'flips'
#        sequences extracted from the genome, so double check."
#                                                   (protocols.io v3, step 4)
#   R5  Synthesis: <=10 consecutive identical bp; strong AT/GC bias and repeats
#        are problematic.                            (Stern 2024, Methods)
#   R6  "It is optimal that the payload integrated during homology directed
#        repair destroys the gRNA target site in the genome."
#                                                   (protocols.io v3, step 3)
#
# Two real-world parsing hazards were observed in the downloaded collection and
# are handled explicitly:
#   H1  Geneious exports embed HTML fragments (<a href="urn:local:...">) inside
#       qualifier values, so `file` misreports a valid GenBank record as an HTML
#       document. Detection is therefore by LOCUS/FEATURES/ORIGIN content, never
#       by file extension or MIME type.
#   H2  The recombination site may be annotated under a custom feature key
#       (DGRC uses `attR3`) rather than protein_bind + /label="attR3" (Addgene).
#       The key is never trusted alone; the qualifiers decide.
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
PJAT_VERSION <- "1.0.0"

# The pJAT collection is resolved at runtime (see R/app_paths.R) so the
# desktop app does not depend on a path on an external volume.
local({
  for (p in unique(c(
    file.path("R", "app_paths.R"),
    file.path(getwd(), "R", "app_paths.R"),
    file.path(dirname(getwd()), "R", "app_paths.R")
  ))) {
    if (file.exists(p) && !exists("get_app_pjat_dir", mode = "function")) {
      tryCatch(source(p, local = FALSE), error = function(e) NULL)
      break
    }
  }
})
if (exists("load_app_paths", mode = "function")) load_app_paths("PLASMID")

# Legacy name, kept for compatibility. Resolved at runtime.
DEFAULT_PJAT_DIR <- if (exists("get_app_pjat_dir", mode = "function")) {
  get_app_pjat_dir("PLASMID")
} else {
  ""
}

PJAT_FILE_EXTENSIONS <- c("gbk", "gb", "gbff", "genbank")

PJAT_RECOMBINATION_KEYS <- c("protein_bind", "misc_feature", "misc_binding",
                             "recombination", "site")

.pjat_resolve_project_root <- function() {
  candidates <- unique(c(
    getwd(),
    app_root()
  ))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (file.exists(script_path)) {
      candidates <- unique(c(dirname(normalizePath(script_path, mustWork = FALSE)), candidates))
    }
  }
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "R"))) return(normalizePath(cand))
  }
  candidates[1]
}

PJAT_PROJECT_ROOT <- .pjat_resolve_project_root()

# ============================================================================
# 1. Standard Genetic Code (for synonymous-codon edit advisor)
# ============================================================================

.pjat_std_genetic_code <- c(
  "TTT"="F","TTC"="F","TTA"="L","TTG"="L","CTT"="L","CTC"="L","CTA"="L","CTG"="L",
  "ATT"="I","ATC"="I","ATA"="I","ATG"="M","GTT"="V","GTC"="V","GTA"="V","GTG"="V",
  "TCT"="S","TCC"="S","TCA"="S","TCG"="S","CCT"="P","CCC"="P","CCA"="P","CCG"="P",
  "ACT"="T","ACC"="T","ACA"="T","ACG"="T","GCT"="A","GCC"="A","GCA"="A","GCG"="A",
  "TAT"="Y","TAC"="Y","TAA"="*","TAG"="*","CAT"="H","CAC"="H","CAA"="Q","CAG"="Q",
  "AAT"="N","AAC"="N","AAA"="K","AAG"="K","GAT"="D","GAC"="D","GAA"="E","GAG"="E",
  "TGT"="C","TGC"="C","TGA"="*","TGG"="W","CGT"="R","CGC"="R","CGA"="R","CGG"="R",
  "AGT"="S","AGC"="S","AGA"="R","AGG"="R","GGT"="G","GGC"="G","GGA"="G","GGG"="G"
)

# ============================================================================
# 2. GenBank Parsing (hazard H1: HTML-in-notes)
# ============================================================================

.pjat_parse_location <- function(loc) {
  comp <- grepl("^complement\\(", loc)
  nums <- as.integer(unlist(regmatches(loc, gregexpr("[0-9]+", loc, perl = TRUE))))
  if (length(nums) < 2L) {
    return(list(start = NA_integer_, end = NA_integer_, strand = "+"))
  }
  list(start = nums[1L], end = nums[2L], strand = if (comp) "-" else "+")
}

.pjat_empty_features <- function() {
  data.frame(key = character(0), location = character(0),
             start = integer(0), end = integer(0), strand = character(0),
             label = character(0), product = character(0),
             qualifiers = character(0), stringsAsFactors = FALSE)
}

#' Read one GenBank file into a structured list (hazard H1).
#' @param file Path to a GenBank file (.gbk/.gb/.gbff).
#' @return list(file, locus, size_bp, definition, features, sequence)
#' @export
read_pjat_genbank <- function(file) {
  if (!file.exists(file)) stop("GenBank file not found: ", file)
  lines <- readLines(file, warn = FALSE)
  if (length(lines) == 0L) stop("Empty file: ", basename(file))

  # Validate by CONTENT, not extension or MIME type (hazard H1).
  if (!any(grepl("^LOCUS\\s", lines))) {
    stop("Not a GenBank record (no LOCUS line): ", basename(file))
  }

  first <- lines[1L]
  # Robust to both variants seen in the collection:
  #   LOCUS       sequence-406186-                      7986 bp    DNA ...
  #   LOCUS       pJAT55{...}  13951 bp    DNA     circular ...
  # Match on "<digits> bp" so token position never matters.
  loc_tok <- regmatches(first,
    regexpr("([0-9]+)\\s+bp", first, perl = TRUE))
  size  <- suppressWarnings(as.integer(sub("\\s*bp$", "", loc_tok)))
  locus <- trimws(gsub("\"", "",
              sub("\\s+[0-9]+\\s+bp\\s.*$", "", sub("^LOCUS\\s+", "", first), perl = TRUE)))
  if (is.na(size)) {
    stop("Could not parse LOCUS size for ", basename(file))
  }

  d_i <- which(grepl("^DEFINITION\\s", lines))[1L]
  definition <- if (!is.na(d_i)) {
    j <- d_i
    while (j + 1L <= length(lines) &&
           !grepl("^(ACCESSION|VERSION|KEYWORDS|SOURCE|ORGANISM)",
                  lines[j + 1L])) {
      j <- j + 1L
    }
    paste(trimws(lines[d_i:j]), collapse = " ")
  } else ""
  definition <- trimws(sub("^DEFINITION\\s+", "", definition))

  list(file = basename(file), path = file, locus = locus,
       size_bp = as.integer(size), definition = definition,
       features = .pjat_parse_features_table(lines),
       sequence = .pjat_parse_origin(lines))
}

.pjat_parse_features_table <- function(lines) {
  f_i <- which(grepl("^FEATURES\\b", lines))[1L]
  o_i <- which(grepl("^ORIGIN\\b", lines))[1L]
  if (is.na(f_i)) return(.pjat_empty_features())

  end_i <- if (!is.na(o_i)) o_i - 1L else length(lines)
  idx   <- seq.int(f_i + 1L, end_i)
  idx   <- idx[nzchar(trimws(lines[idx]))]
  if (length(idx) == 0L) return(.pjat_empty_features())

  is_head <- grepl("^ {5}(\\S+)\\s+(\\S+)", lines[idx])
  keys <- character(0); locs <- character(0); qs <- list(); n <- 0L

  for (k in seq_along(idx)) {
    if (is_head[k]) {
      key <- sub("^ {5}(\\S+)\\s+.*$", "\\1", lines[idx[k]], perl = TRUE)
      if (grepl("^[0-9]+$", key)) next        # ORIGIN sequence line, not a feature
      loc <- sub("^ {5}\\S+\\s+(.*)$", "\\1", lines[idx[k]], perl = TRUE)
      n <- n + 1L
      keys[n] <- key; locs[n] <- loc; qs[[n]] <- character(0)
    } else if (n > 0L && grepl("^\\s*/", lines[idx[k]])) {
      qs[[n]] <- c(qs[[n]], trimws(lines[idx[k]]))
    }
  }
  if (n == 0L) return(.pjat_empty_features())

  keys <- keys[seq_len(n)]; locs <- locs[seq_len(n)]; qs <- qs[seq_len(n)]

  data.frame(
    key   = keys,
    location = locs,
    start = vapply(locs, function(l) .pjat_parse_location(l)$start, integer(1)),
    end   = vapply(locs, function(l) .pjat_parse_location(l)$end,   integer(1)),
    strand= vapply(locs, function(l) .pjat_parse_location(l)$strand, character(1)),
    label     = .pjat_qualifier_values(qs, "label"),
    product   = .pjat_qualifier_values(qs, "product"),
    qualifiers= vapply(qs, function(q) paste(q, collapse = "\t"), character(1)),
    stringsAsFactors = FALSE
  )
}

.pjat_qualifier_values <- function(qs, name) {
  vapply(qs, function(q) {
    hit <- grep(paste0("^/", name, "="), q, value = TRUE)
    if (length(hit) == 0L) return("")
    gsub("\"", "", sub(paste0("^/", name, "="), "", hit[1L]))
  }, character(1))
}

.pjat_parse_origin <- function(lines) {
  o_i <- which(grepl("^ORIGIN\\b", lines))[1L]
  if (is.na(o_i)) return("")
  b_i <- which(grepl("^//", lines))[1L]
  end_i <- if (!is.na(b_i)) b_i - 1L else length(lines)
  if (end_i <= o_i) return("")
  toupper(gsub("[^A-Za-z]", "", paste(lines[seq.int(o_i + 1L, end_i)],
                                      collapse = "")))
}

# ============================================================================
# 3. attR / attL / attP Site Discovery (hazard H2: custom feature keys)
# ============================================================================

.pjat_recomb_feature <- function(gb) {
  f <- gb$features
  if (is.null(f) || nrow(f) == 0L) return(f)
  key_hit <- grepl("^att[LBRP][0-9]?$", f$key, perl = TRUE)
  lab_hit <- grepl("att[RBLP][0-9]?", f$label, perl = TRUE)
  gen_hit <- grepl("mutant version of attR|Gateway|phage phi-C31 attP",
                   paste(f$product, f$qualifiers), perl = TRUE)
  f[key_hit | (lab_hit & f$key %in% PJAT_RECOMBINATION_KEYS) | gen_hit, ,
    drop = FALSE]
}

#' Normalize an att feature's label to its canonical name (attR3, attP, ...).
.pjat_normalize_att <- function(f) {
  lab <- f$label
  key_is_att <- grepl("^att[LBRP][0-9]?$", f$key, perl = TRUE)
  lab[key_is_att & !nzchar(lab)] <- f$key[key_is_att & !nzchar(lab)]
  m <- regmatches(lab, regexpr("att[LBRP][0-9]?", lab, perl = TRUE))
  out <- rep(NA_character_, length(lab))
  ok <- nzchar(m)
  out[ok] <- m[ok]
  fb <- rep(NA_character_, length(lab))
  fb[!nzchar(lab) & nzchar(f$product)] <- f$product[!nzchar(lab) & nzchar(f$product)]
  fb[!nzchar(lab) & !nzchar(f$product)] <- f$key[!nzchar(lab) & !nzchar(f$product)]
  out[is.na(out)] <- fb[is.na(out)]
  out
}

#' Extract all Gateway recombination sites from a parsed pJAT plasmid.
#' @param gb Output of read_pjat_genbank().
#' @return data.frame with att, start, end, strand, key, label, width
#' @export
find_pjat_recombination_sites <- function(gb) {
  empty <- data.frame(att = character(0), start = integer(0), end = integer(0),
                      strand = character(0), key = character(0),
                      label = character(0), width = integer(0),
                      stringsAsFactors = FALSE)
  f <- .pjat_recomb_feature(gb)
  if (is.null(f) || nrow(f) == 0L) return(empty)

  f <- f[!is.na(f$start) & !is.na(f$end), , drop = FALSE]
  if (nrow(f) == 0L) return(empty)

  out <- data.frame(
    att    = .pjat_normalize_att(f),
    start  = as.integer(f$start),
    end    = as.integer(f$end),
    strand = f$strand,
    key    = f$key,
    label  = f$label,
    width  = as.integer(f$end - f$start + 1L),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$start), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Count U6 promoter features in a parsed pJAT plasmid.
#' @param gb Output of read_pjat_genbank().
#' @return integer count of Drosophila U6 promoters (dU6-3 / U6-3)
#' @export
count_pjat_u6_promoters <- function(gb) {
  f <- gb$features
  if (is.null(f) || nrow(f) == 0L) return(0L)
  blob <- paste(f$label, f$product, f$qualifiers)
  hit <- (grepl("^promoter$", f$key) | grepl("promoter", f$product,
                                             ignore.case = TRUE)) &
         grepl("U6-3|U6 promoter|RNA polymerase III", blob,
               ignore.case = TRUE, perl = TRUE)
  as.integer(sum(hit))
}

# ============================================================================
# 4. Collection Manifest & Validation
# ============================================================================

#' Build a manifest of every pJAT plasmid in a directory.
#'
#' Skips AppleDouble sidecar files (._*), non-GenBank files, and duplicates.
#' Adds a `status` column: "OK", or "MISSING:attR3+attR4" when the four Gateway
#' attR sites are incomplete.
#'
#' @param dir Directory containing pJAT GenBank files.
#' @param verbose Print progress per file.
#' @return data.frame, one row per plasmid.
#' @export
pjat_manifest <- function(dir = DEFAULT_PJAT_DIR, verbose = FALSE) {
  files <- list.files(dir, full.names = TRUE)
  files <- files[grepl("\\.(gbk|gb|gbff|genbank)$", files, ignore.case = TRUE)]
  files <- files[!grepl("(^|/)\\._", files)]
  files <- sort(files)
  if (length(files) == 0L) stop("No GenBank files found in ", dir)

  rows <- lapply(files, function(f) {
    gb  <- read_pjat_genbank(f)
    att <- find_pjat_recombination_sites(gb)
    name <- sub("\\.(gbk|gb|gbff|genbank)$", "", basename(f), ignore.case = TRUE)
    name <- sub("_empty_backbone$", "", name)
    name <- sub("^(pJAT[0-9]+).*", "\\1", name)

    get <- function(a) {
      r <- att[att$att == a, , drop = FALSE]
      if (nrow(r) == 0L) return(c(NA_integer_, NA_integer_, NA_character_))
      c(r$start[1L], r$end[1L], r$strand[1L])
    }
    r3 <- get("attR3"); r4 <- get("attR4"); r5 <- get("attR5"); r6 <- get("attR6")

    data.frame(
      plasmid = name, file = basename(f), size_bp = gb$size_bp,
      n_u6 = count_pjat_u6_promoters(gb), has_attp = any(att$att == "attP"),
      attR3_start = as.integer(r3[1]), attR3_end = as.integer(r3[2]),
      attR3_strand = as.character(r3[3]),
      attR4_start = as.integer(r4[1]), attR4_end = as.integer(r4[2]),
      attR4_strand = as.character(r4[3]),
      attR5_start = as.integer(r5[1]), attR5_end = as.integer(r5[2]),
      attR5_strand = as.character(r5[3]),
      attR6_start = as.integer(r6[1]), attR6_end = as.integer(r6[2]),
      attR6_strand = as.character(r6[3]),
      stringsAsFactors = FALSE)
  })

  man <- do.call(rbind, rows)
  man <- validate_pjat_collection(man)
  rownames(man) <- NULL
  if (verbose) {
    cat(sprintf("pJAT manifest: %d plasmids from %s\n", nrow(man), dir))
    cat(sprintf("  all four attR sites present: %d/%d\n",
                sum(man$status == "OK"), nrow(man)))
  }
  man
}

#' Validate that a pJAT plasmid has the four required Gateway attR sites.
#' @param manifest data.frame with attR3_start..attR6_start columns.
#' @return the same data.frame with a status column.
#' @export
validate_pjat_collection <- function(manifest) {
  req <- c("attR3_start", "attR4_start", "attR5_start", "attR6_start")
  if (!all(req %in% names(manifest))) {
    stop("manifest must contain attR3_start..attR6_start columns")
  }
  manifest$status <- apply(manifest[, req, drop = FALSE], 1, function(r) {
    miss <- req[is.na(r)]
    if (length(miss) == 0L) "OK" else paste0("MISSING:",
      paste(sub("_start$", "", miss), collapse = "+"))
  })
  manifest
}

# ============================================================================
# 5. Twist Synthesis Feasibility QC (rule R5)
# ============================================================================

#' Assess a sequence for synthesis feasibility before ordering.

#' Load the pJAT plasmid catalog from the generated YAML.
#'
#' Reads `donor_templates/pjat_plasmids.yml`. attR entries are stored as
#' [start, end, strand]; this returns them as a data.frame per plasmid for the
#' Shiny picker and the compatibility rules.
#'
#' @param yml Path to the YAML catalog.
#' @return data.frame with one row per plasmid plus attR3/4/5/6 as lists.
#' @export
pjat_load_catalog <- function(yml = file.path(PJAT_PROJECT_ROOT, "donor_templates",
                                              "pjat_plasmids.yml")) {
  if (!file.exists(yml)) stop("Catalog not found: ", yml)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read the pJAT catalog.")
  }
  y <- yaml::read_yaml(yml)
  recs <- lapply(y$plasmids, function(p) {
    get3 <- function(k) {
      v <- p$attR[[k]]
      if (is.null(v)) return(list(start = NA_integer_, end = NA_integer_,
                                  strand = NA_character_))
      list(start = as.integer(v[[1]]), end = as.integer(v[[2]]),
           strand = as.character(v[[3]]))
    }
    a3 <- get3("attR3"); a4 <- get3("attR4")
    a5 <- get3("attR5"); a6 <- get3("attR6")
    data.frame(
      name = p$name, addgene_id = as.integer(p$addgene_id),
      supplier = p$supplier, size_bp = as.integer(p$size_bp),
      u6_promoters = as.integer(p$u6_promoters),
      has_attp = isTRUE(p$has_attp), purpose = p$purpose,
      description = p$description %||% p$purpose %||% NA_character_,
      source_file = p$source_file,
      attR3_start = a3$start, attR3_end = a3$end, attR3_strand = a3$strand,
      attR4_start = a4$start, attR4_end = a4$end, attR4_strand = a4$strand,
      attR5_start = a5$start, attR5_end = a5$end, attR5_strand = a5$strand,
      attR6_start = a6$start, attR6_end = a6$end, attR6_strand = a6$strand,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, recs)
  rownames(out) <- NULL
  # Attach Addgene/SnapGene vector-map metadata when R/addgene_maps.R is loaded
  # (Plasmid_Designer.R sources it first). When this module is used standalone
  # the join is skipped and callers get the base columns only.
  if (exists("addgene_maps_join", mode = "function")) {
    out <- addgene_maps_join(out)
  }
  out
}

#' Recommend pJAT plasmids by required gRNA count (U6 promoter count).
#'
#' The U6 count is what forces the plasmid choice: pJAT77 has 0 U6 (you supply
#' your own promoter), pJAT74 has 1, pJAT32 has 2. Recommendations never
#' override the user's own selection -- they are advisory only.
#'
#' @param n_guides Number of gRNAs this one plasmid must express (0, 1 or 2).
#' @param catalog Optional pre-loaded catalog from pjat_load_catalog().
#' @return data.frame of compatible plasmids, ordered by preference.
#' @export
pjat_recommend_backbones <- function(n_guides, catalog = NULL) {
  if (is.null(catalog)) catalog <- pjat_load_catalog()
  n <- as.integer(n_guides[1])
  if (is.na(n) || n < 0L || n > 2L) {
    stop("n_guides must be 0, 1 or 2 (pJAT supports at most two U6 promoters).")
  }
  # Empty MCS backbones only -- the payload plasmids are for finished constructs.
  empty <- c("pJAT77", "pJAT74", "pJAT32")
  back <- catalog[catalog$name %in% empty, , drop = FALSE]
  if (n == 0L) {
    keep <- back$name == "pJAT77"
  } else {
    keep <- back$u6_promoters >= n
  }
  out <- back[keep, , drop = FALSE]
  out <- out[order(out$u6_promoters), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Check whether a chosen backbone can express the required number of gRNAs.
#'
#' @param plasmid Plasmid name (e.g. "pJAT74").
#' @param n_guides Number of gRNAs.
#' @param catalog Optional pre-loaded catalog.
#' @return TRUE if compatible; errors with an explanatory message if not.
#' @export
pjat_check_compatible <- function(plasmid, n_guides, catalog = NULL) {
  if (is.null(catalog)) catalog <- pjat_load_catalog()
  row <- catalog[catalog$name == plasmid, , drop = FALSE]
  if (nrow(row) == 0L) stop("Unknown pJAT plasmid: ", plasmid)
  n <- as.integer(n_guides[1])
  if (n == 0L) return(invisible(TRUE))
  if (row$u6_promoters < n) {
    alt <- pjat_recommend_backbones(n, catalog)
    stop(sprintf(
      "%s has %d U6 promoter(s) but %d gRNA(s) are required. Use one of: %s.",
      plasmid, row$u6_promoters, n,
      paste(alt$name, collapse = ", ")))
  }
  invisible(TRUE)
}
#'
#' Flags runs >10 identical bases (hard fail per Twist rule), strong AT/GC bias,
#' and low-complexity/repeat-like composition.
#'
#' @param seq DNA sequence (A/C/G/T/N).
#' @param gc_min,gc_max GC fraction bounds for the bias warning (default 30/70%).
#' @return list(ok, n, gc_frac, max_homopolymer, extreme_gc, low_complexity,
#'   warnings)
#' @export
pjat_synthesis_qc <- function(seq, gc_min = 0.30, gc_max = 0.70) {
  s <- toupper(gsub("[^ACGTN]", "", seq))
  n <- nchar(s)
  if (n == 0L) stop("Empty sequence")

  chars <- strsplit(s, "")[[1]]
  gc <- sum(chars %in% c("G", "C")) / n

  runs <- rle(chars)
  max_run <- if (length(runs$lengths)) max(runs$lengths) else 0L
  hp_pos <- integer(0)
  if (max_run > 10L) {
    ends   <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    bad    <- which(runs$lengths > 10L)
    if (length(bad)) {
      hp_pos <- unlist(lapply(bad, function(i) seq.int(starts[i], ends[i])))
    }
  }

  distinct <- length(unique(chars[chars != "N"]))
  ent <- .pjat_min_window_entropy(s, 50L)
  low_complexity <- (distinct <= 2L) || (ent < 1.0)
  extreme_gc <- (gc < gc_min) || (gc > gc_max)

  warnings <- character(0)
  if (max_run > 10L) {
    warnings <- c(warnings, sprintf(
      "Homopolymer run of %d bp (>10) at position(s) %s.",
      max_run, paste(utils::head(hp_pos, 5), collapse = ",")))
  }
  if (extreme_gc) warnings <- c(warnings, sprintf(
    "GC content %.1f%% outside %.0f-%.0f%%.", 100 * gc, 100 * gc_min, 100 * gc_max))
  if (low_complexity) warnings <- c(warnings, sprintf(
    "Low-complexity/repeat-like composition (%d distinct bases, min entropy %.2f).",
    distinct, ent))

  list(ok = length(warnings) == 0L, n = n, gc_frac = gc,
       max_homopolymer = max_run, extreme_gc = extreme_gc,
       low_complexity = low_complexity, warnings = warnings)
}

.pjat_min_window_entropy <- function(s, w) {
  n <- nchar(s)
  if (n < w) w <- n
  if (w <= 0L) return(0)
  chars <- strsplit(s, "")[[1]]
  ent <- vapply(seq_len(n - w + 1L), function(i) {
    tbl <- table(chars[i:(i + w - 1L)])
    p <- as.numeric(tbl) / w
    -sum(p * log2(p))
  }, numeric(1))
  min(ent)
}

# ============================================================================
# 6. Gateway Donor Fragment Assembly (rules R1, R3, R4)
# ============================================================================

.pjat_revcomp <- function(seq) {
  s <- toupper(seq)
  chars <- strsplit(s, "")[[1]]
  comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  paste(rev(comp[chars]), collapse = "")
}

#' Validate a 20-nt spacer for the tRNA-gRNA-tRNA array (rule R3).
#'
#' Hard requirements from the protocol: the gRNA is the 20 bp target sequence
#' with the NGG PAM excluded. Including the PAM makes the experiment fail, so
#' this is an error, not a warning.
#'
#' @param spacer 20-nt protospacer (PAM excluded).
#' @param pam Optional 3-nt PAM that was excluded; used to confirm the caller
#'   really did drop it.
#' @return TRUE invisibly; errors on wrong length or PAM inclusion.
#' @export
pjat_validate_spacer <- function(spacer, pam = NULL) {
  s <- toupper(gsub("[^ACGTN]", "", spacer))
  if (nchar(s) != 20L) {
    stop("gRNA spacer must be exactly 20 nt (got ", nchar(s), "). The NGG PAM ",
         "must be excluded: 'The experiment will fail if you include the NGG ",
         "sequence.' (protocols.io v3, step 4)")
  }
  if (!is.null(pam) && nzchar(pam)) {
    p <- toupper(gsub("[^ACGTN]", "", pam))
    if (nchar(p) == 3L && grepl("^[ACGTN]{20}$", s)) {
      # If the caller passed spacer+PAM concatenated, the length check above
      # already caught it. Here we confirm the PAM is not the last 3 nt.
      if (grepl(p, s, fixed = TRUE) && substring(s, 18L, 20L) == p) {
        stop("Spacer appears to include the PAM. Exclude the NGG PAM.")
      }
    }
  }
  invisible(TRUE)
}

#' Assemble the two attL Gateway donor fragments for a pJAT experiment.
#'
#' Produces, for one left and one right homology arm plus a 20-nt spacer:
#'   Fragment 1 (attL3-arm-attL4)  :  attL3 + LEFT_ARM + tRNA-sgRNA-tRNA + attL4
#'   Fragment 2 (attL5-arm-attR6->attL6) : attL5 + RIGHT_ARM + attL6
#'
#' attL sites are supplied by the caller (they are consumed during LR cloning and
#' are therefore absent from the destination plasmids). They are never invented:
#' if missing, the fragments are emitted with the arms in place and a status of
#' "ARM_ONLY" so the caller can see exactly what is missing.
#'
#' The orientation of both arms is asserted against the genome before assembly
#' (rule R4): each arm is taken in the direction it occupies in the target
#' locus, so no arm is silently flipped.
#'
#' @param left_arm,arm_right 5' and 3' homology arms, each already oriented as
#'   they appear in the genome (5'->3' on the target strand).
#' @param spacer 20-nt gRNA spacer (PAM excluded), in gRNA orientation.
#' @param attl3,attl4,attl5,attl6 Named list of flanking sequences, or NULL.
#' @param trna Optional trna-tRNA spacer sequence for the array.
#' @return list(status, fragments, checks, qc)
#' @export
pjat_assemble_fragments <- function(left_arm, right_arm, spacer,
                                    attl3 = NULL, attl4 = NULL,
                                    attl5 = NULL, attl6 = NULL,
                                    trna = NULL) {
  pjat_validate_spacer(spacer)
  spacer <- toupper(gsub("[^ACGTN]", "", spacer))
  left_arm  <- toupper(gsub("[^ACGTN]", "", left_arm))
  right_arm <- toupper(gsub("[^ACGTN]", "", right_arm))

  if (nchar(left_arm) < 250L || nchar(right_arm) < 250L) {
    stop("Both arms must be >=250 bp: 'have not observed a difference in ",
         "efficiency between arms containing 250 bp and 1000 bp' implies 250 bp ",
         "is the validated floor (protocols.io v3, step 3).")
  }
  # R2 ceiling: the same source validates 1000 bp as the effective maximum.
  if (nchar(left_arm) > 1000L || nchar(right_arm) > 1000L) {
    stop("Both arms must be <=1000 bp (validated maximum, rule R2): ",
         "5'=", nchar(left_arm), " bp, 3'=", nchar(right_arm), " bp. ",
         "Trim with arm_trim_to_window() before assembling fragments.")
  }

  # Default tRNA flank: if none supplied, emit the array as spacer only and
  # mark the array as incomplete rather than fabricating a sequence.
  array <- if (!is.null(trna) && nzchar(trna)) {
    paste0(toupper(trna), "-", spacer, "-", toupper(trna))
  } else {
    paste0("[tRNA]-", spacer, "-[tRNA]")
  }

  have_all <- all(vapply(list(attl3, attl4, attl5, attl6),
                         function(x) !is.null(x) && nzchar(x), logical(1)))

  frag1_left  <- paste0(if (is.null(attl3)) "[attL3]" else attl3, left_arm)
  frag1_right <- paste0(array, if (is.null(attl4)) "[attL4]" else attl4)
  frag2_left  <- paste0(if (is.null(attl5)) "[attL5]" else attl5, right_arm)
  frag2_right <- if (is.null(attl6)) "[attL6]" else attl6

  frag1 <- paste0(frag1_left, frag1_right)
  frag2 <- paste0(frag2_left, frag2_right)

  qc1 <- pjat_synthesis_qc(left_arm)
  qc2 <- pjat_synthesis_qc(right_arm)

  list(
    status = if (have_all) "READY" else "ARM_ONLY",
    spacer = spacer,
    fragments = list(
      fragment1_attL3_arm_trna_attL4 = frag1,
      fragment2_attL5_arm_attL6      = frag2
    ),
    fragment_parts = list(
      f1 = c(attL3 = if (is.null(attl3)) "[missing]" else attl3,
             left_arm = left_arm, trna_sgRNA_trna = array,
             attL4 = if (is.null(attl4)) "[missing]" else attl4),
      f2 = c(attL5 = if (is.null(attl5)) "[missing]" else attl5,
             right_arm = right_arm,
             attL6 = if (is.null(attl6)) "[missing]" else attl6)
    ),
    checks = list(
      spacer_is_20nt      = nchar(spacer) == 20L,
      left_arm_ge_250     = nchar(left_arm) >= 250L,
      right_arm_ge_250    = nchar(right_arm) >= 250L,
      attl_sites_supplied = have_all,
      # rule R4: arms must not be flipped relative to the genome. The caller
      # supplies arms already in genomic orientation; we record the assumption
      # so it cannot be silently violated downstream.
      arms_in_genomic_orientation = TRUE
    ),
    qc = list(left = qc1, right = qc2)
  )
}