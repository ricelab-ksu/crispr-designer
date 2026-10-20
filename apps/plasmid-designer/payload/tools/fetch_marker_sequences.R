#!/usr/bin/env Rscript
# ============================================================================
# tools/fetch_marker_sequences.R — bundle real marker sequences for the app
# ============================================================================
# marker_cassettes.yml ships metadata only. This tool turns every marker into a
# provenance-tracked FASTA under donor_templates/marker_sequences/ and records
# sequence_file / length_bp / md5 / sequence_provenance in the catalog.
#
# Two sources, in this order:
#   1. LOCAL  — Drosophila cassettes already present in the pJAT GenBank files
#      (reference_files/janelia_atalanta_plasmids/). Extracted between attR4 and
#      attR5, with circular wrap-around handled for origin-spanning payloads.
#   2. NCBI   — standard bacterial/yeast markers fetched from NuCCORE with
#      E-utilities (efetch GenBank, then the CDS feature matching the gene).
#
# Honesty rules (enforced in code):
#   * No sequence is ever invented or hand-typed.
#   * A marker whose sequence cannot be obtained is written with
#     sequence_status: unavailable and a reason — never a placeholder sequence.
#   * Every stored FASTA records its source (file + coordinates, or accession).
#
# USAGE
#   Rscript tools/fetch_marker_sequences.R
#   Rscript tools/fetch_marker_sequences.R --force      # re-extract everything
#   Rscript tools/fetch_marker_sequences.R --offline    # local pJAT only
#   Rscript tools/fetch_marker_sequences.R --sleep 0.5  # NCBI politeness delay
#
# Exit code 0 = catalog updated (some markers may still be unavailable).
# ============================================================================

suppressPackageStartupMessages({
  HAS_YAML <- requireNamespace("yaml", quietly = TRUE)
  HAS_HTTR <- requireNamespace("httr", quietly = TRUE)
})

.script_path <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) sub("^--file=", "", a[1L]) else NA_character_
}

.resolve_root <- function() {
  sp <- .script_path()
  cands <- unique(c(
    if (!is.na(sp)) dirname(dirname(normalizePath(sp, mustWork = FALSE))) else NULL,
    getwd(),
    "/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer"
  ))
  for (c in cands) {
    if (file.exists(file.path(c, "donor_templates", "marker_cassettes.yml"))) {
      return(normalizePath(c, mustWork = FALSE))
    }
  }
  stop("Could not locate the Plasmid_Designer project root.")
}

ROOT <- .resolve_root()
source(file.path(ROOT, "R", "pjat_plasmids.R"))

if (!HAS_YAML) stop("Package 'yaml' is required.")

ARGS <- commandArgs(trailingOnly = TRUE)
if (any(grepl("^(--help|-h)$", ARGS))) {
  cat("Usage: Rscript tools/fetch_marker_sequences.R",
      "[--force] [--offline] [--sleep <sec>]\n")
  quit(status = 0L)
}

parse_flag <- function(flag, default = NULL) {
  eq <- grep(paste0("^", flag, "="), ARGS, value = TRUE)
  if (length(eq)) return(sub(paste0("^", flag, "="), "", eq[1L]))
  i <- match(flag, ARGS)
  if (!is.na(i) && i < length(ARGS)) return(ARGS[i + 1L])
  default
}

FORCE   <- "--force" %in% ARGS
OFFLINE <- "--offline" %in% ARGS
SLEEP   <- as.numeric(parse_flag("--sleep", "0.4"))
if (!is.finite(SLEEP) || SLEEP < 0) SLEEP <- 0.4

CAT_YML <- file.path(ROOT, "donor_templates", "marker_cassettes.yml")
SEQ_DIR <- file.path(ROOT, "donor_templates", "marker_sequences")
dir.create(SEQ_DIR, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(..., "\n", sep = "")

# ---- helpers ---------------------------------------------------------------

revcomp <- function(x) {
  paste(rev(strsplit(chartr("ACGTNacgtn", "TGCANtgcan", x), "")[[1L]]), collapse = "")
}

clean_seq <- function(x) toupper(gsub("[^ACGTN]", "", x))

wrap80 <- function(x) paste0(strwrap(x, width = 80, simplify = FALSE), collapse = "\n")

#' FASTA writer with a provenance-bearing header.
write_marker_fasta <- function(id, sequence, provenance, path) {
  header <- sprintf(">%s|length=%dbp|%s", id, nchar(sequence), provenance)
  writeLines(c(header, wrap80(sequence)), con = path, useBytes = TRUE)
  invisible(path)
}

#' Extract [start,end] from a circular sequence, wrapping past the origin.
circular_slice <- function(sequence, start, end) {
  n <- nchar(sequence)
  if (start <= end) return(substr(sequence, start, end))
  paste0(substr(sequence, start, n), substr(sequence, 1L, end))
}


# ---- source specifications -------------------------------------------------

PJAT_DIR <- DEFAULT_PJAT_DIR

# Drosophila cassettes: payload between attR4 and attR5 of the named plasmid.
LOCAL_SPECS <- list(
  list(id = "3xP3_dsRed", file = "pJAT30.gbk"),
  list(id = "ie1_dsRed",  file = "pJAT11.gbk")
)

# Standard markers: direct accession, or gene lookup (esearch), then the CDS
# whose annotation matches `match`.
NCBI_SPECS <- list(
  list(id = "neo",   accession = "U00004.1",
       match = "aph|aminoglycoside|kanamycin|neomycin"),
  list(id = "puro",  queries = c("pac[Gene] AND Streptomyces alboniger[Organism]",
                                 "puromycin N-acetyltransferase"),
       match = "pac|puromycin|acetyltransferase"),
  list(id = "hygro", queries = c("hph[Gene]", "hygromycin B phosphotransferase"),
       match = "hph|hygromycin"),
  list(id = "blast", queries = c("bsd[Gene] AND Aspergillus terreus[Organism]",
                                 "blasticidin S deaminase"),
       match = "bsd|blasticidin|deaminase"),
  list(id = "URA3",  queries = c("URA3[Gene] AND Saccharomyces cerevisiae[Organism]"),
       match = "URA3|orotidine"),
  list(id = "LEU2",  queries = c("LEU2[Gene] AND Saccharomyces cerevisiae[Organism]"),
       match = "LEU2|isopropylmalate"),
  list(id = "KanMX", accession = "11344909",
       match = "KanMX|kanamycin|aminoglycoside|aph")
)

MIN_BP <- 300L
MAX_BP <- 4000L

# ---- local extraction ------------------------------------------------------

#' Payload between attR4 and attR5 of a pJAT plasmid.
local_cassette <- function(id, file) {
  path <- file.path(PJAT_DIR, file)
  if (!file.exists(path)) stop("missing local GenBank file: ", path)
  gb  <- suppressWarnings(read_pjat_genbank(path))
  att <- find_pjat_recombination_sites(gb)
  a4  <- att[att$att == "attR4", , drop = FALSE]
  a5  <- att[att$att == "attR5", , drop = FALSE]
  if (nrow(a4) == 0L || nrow(a5) == 0L) stop("attR4/attR5 not found in ", file)
  start <- a4$end[1L] + 1L
  end   <- a5$start[1L] - 1L
  seq   <- clean_seq(circular_slice(gb$sequence, start, end))
  if (nchar(seq) < MIN_BP || nchar(seq) > MAX_BP) {
    stop(sprintf("extracted region is %d bp (expected %d..%d); refusing to bundle it",
                 nchar(seq), MIN_BP, MAX_BP))
  }
  list(
    sequence = seq,
    provenance = sprintf("pJAT local GenBank %s bp=%d..%d (attR4..attR5 payload)",
                         file, start, end)
  )
}

# ---- NCBI ------------------------------------------------------------------

NCBI_BASE <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA <- "PlasmidDesigner/1.0 (local research app; marker sequence cache)"

http_text <- function(url, tries = 3L) {
  if (!HAS_HTTR) stop("Package 'httr' is required for NCBI downloads")
  last <- NULL
  for (i in seq_len(tries)) {
    r <- tryCatch(httr::GET(url, httr::user_agent(UA), httr::timeout(60)),
                  error = function(e) e)
    if (!inherits(r, "error") && httr::status_code(r) == 200L) {
      return(httr::content(r, "text", encoding = "UTF-8"))
    }
    last <- if (inherits(r, "error")) conditionMessage(r) else
      paste0("HTTP ", httr::status_code(r))
    if (i < tries) Sys.sleep(i * 2)
  }
  stop(last)
}

#' Search NuCCORE and return up to `retmax` record IDs.
ncbi_search <- function(term, retmax = 6L) {
  url <- sprintf("%s/esearch.fcgi?db=nuccore&term=%s&retmax=%d&retmode=json",
                 NCBI_BASE, utils::URLencode(term, reserved = TRUE), retmax)
  js <- jsonlite::fromJSON(http_text(url), simplifyVector = FALSE)
  unlist(js$esearchresult$idlist, use.names = FALSE)
}

.gb_slice_one <- function(sequence, part) {
  part <- gsub("[<>]", "", part)
  b <- as.integer(strsplit(part, "\\.\\.")[[1L]])
  if (length(b) == 1L) b <- c(b, b)
  substr(sequence, min(b), max(b))
}

#' Slice a GenBank location: "123..456", "complement(...)", "join(a..b,c..d)".
gb_slice <- function(sequence, location, strand = "+") {
  comp <- grepl("^complement", location)
  loc <- sub("^complement\\(", "", location)
  if (grepl("\\)$", loc)) loc <- sub("\\)$", "", loc)
  if (grepl("^join", loc)) {
    inner <- sub("^join\\(", "", loc)
    inner <- sub("\\)$", "", inner)
    parts <- strsplit(inner, ",")[[1L]]
    s <- paste0(vapply(parts, function(p) .gb_slice_one(sequence, p),
                       character(1)), collapse = "")
  } else {
    s <- .gb_slice_one(sequence, loc)
  }
  if (comp || identical(strand, "-")) revcomp(s) else s
}

#' Fetch a GenBank record and return the CDS best matching `pattern`.
ncbi_cds <- function(accession, pattern) {
  url <- sprintf("%s/efetch.fcgi?db=nuccore&id=%s&rettype=gb&retmode=text",
                 NCBI_BASE, utils::URLencode(accession, reserved = TRUE))
  gb_text <- http_text(url)
  tmp <- tempfile(fileext = ".gb")
  writeLines(gb_text, tmp, useBytes = TRUE)
  on.exit(unlink(tmp), add = TRUE)
  gb <- suppressWarnings(read_pjat_genbank(tmp))

  feats <- gb$features
  cds <- feats[feats$key == "CDS", , drop = FALSE]
  if (nrow(cds) == 0L) stop("no CDS feature in ", accession)

  haystack <- paste(cds$label, cds$product,
                    gsub("[\t\n]", " ", cds$qualifiers), sep = " ")
  hit <- grepl(pattern, haystack, ignore.case = TRUE)
  if (!any(hit)) stop(sprintf("no CDS in %s matched /%s/", accession, pattern))

  cand <- cds[hit, , drop = FALSE]
  cand$len <- abs(cand$end - cand$start) + 1L
  best <- cand[which.max(cand$len), , drop = FALSE]

  seq <- clean_seq(gb_slice(gb$sequence, best$location[1L], best$strand[1L]))
  lab <- best$label[1L]
  list(
    sequence = seq,
    provenance = sprintf("NCBI NuCCORE %s CDS %s%s",
                         accession, best$location[1L],
                         if (!is.na(lab) && nzchar(lab)) paste0(" ", lab) else "")
  )
}

# ---- catalog update --------------------------------------------------------

update_catalog <- function(results, yml) {
  y <- yaml::read_yaml(yml)
  for (i in seq_along(y$markers)) {
    id <- as.character(y$markers[[i]]$id)
    r <- results[[id]]
    if (is.null(r)) next
    y$markers[[i]]$sequence_status     <- r$status
    y$markers[[i]]$sequence_file       <- r$sequence_file
    y$markers[[i]]$length_bp           <- r$length_bp
    y$markers[[i]]$md5                 <- r$md5
    y$markers[[i]]$sequence_provenance <- r$provenance
  }
  header <- paste0(
    "# ============================================================================\n",
    "# Selectable marker library for Plasmid Designer.\n",
    "#\n",
    "# Sequence fields (sequence_file, length_bp, md5, sequence_provenance,\n",
    "# sequence_status) are GENERATED by tools/fetch_marker_sequences.R.\n",
    "# Re-run that tool instead of hand-editing them.\n",
    "#\n",
    "# Drosophila cassettes are extracted from the local pJAT GenBank files;\n",
    "# standard markers come from NCBI NuCCORE (accession recorded per marker).\n",
    "# sequence_status: bundled = FASTA present | unavailable = could not obtain\n",
    "# (never a placeholder sequence) | none = no marker by design.\n",
    "# ============================================================================\n\n"
  )
  writeLines(paste0(header, yaml::as.yaml(y, indent = 2)), yml, useBytes = TRUE)
  invisible(yml)
}

# ---- main ------------------------------------------------------------------

say(sprintf("Marker sequences -> %s", SEQ_DIR))
results <- list()

record <- function(id, status, path = NA_character_, sequence = NA_character_,
                   provenance, reason = NULL) {
  if (identical(status, "bundled")) {
    results[[id]] <<- list(
      status = "bundled",
      sequence_file = sub(paste0("^", ROOT, "/?"), "", path),
      length_bp = nchar(sequence),
      md5 = unname(tools::md5sum(path)),
      provenance = provenance
    )
    say(sprintf("  %-10s bundled  %5d bp", id, nchar(sequence)))
  } else {
    results[[id]] <<- list(
      status = "unavailable",
      sequence_file = NA_character_,
      length_bp = NA_integer_,
      md5 = NA_character_,
      provenance = paste0("unavailable: ", reason)
    )
    say(sprintf("  %-10s UNAVAILABLE (%s)", id, reason))
  }
}

for (spec in LOCAL_SPECS) {
  path <- file.path(SEQ_DIR, paste0(spec$id, ".fasta"))
  if (!FORCE && file.exists(path) && file.info(path)$size > 0) {
    seq <- clean_seq(paste(readLines(path, warn = FALSE)[-1L], collapse = ""))
    record(spec$id, "bundled", path, seq,
           provenance = "cached (use --force to re-extract)")
    next
  }
  res <- tryCatch(local_cassette(spec$id, spec$file), error = function(e) e)
  if (inherits(res, "error")) {
    record(spec$id, "unavailable", reason = conditionMessage(res))
  } else {
    write_marker_fasta(spec$id, res$sequence, res$provenance, path)
    record(spec$id, "bundled", path, res$sequence, provenance = res$provenance)
  }
}

if (OFFLINE) {
  say("Offline mode: skipping NCBI markers.")
} else {
  for (spec in NCBI_SPECS) {
    path <- file.path(SEQ_DIR, paste0(spec$id, ".fasta"))
    if (!FORCE && file.exists(path) && file.info(path)$size > 0) {
      seq <- clean_seq(paste(readLines(path, warn = FALSE)[-1L], collapse = ""))
      record(spec$id, "bundled", path, seq,
             provenance = "cached (use --force to re-fetch)")
      next
    }
    res <- tryCatch({
      ids <- if (!is.null(spec$accession)) {
        spec$accession
      } else {
        unique(unlist(lapply(spec$queries, function(q) {
          tryCatch(ncbi_search(q), error = function(e) character(0))
        }), use.names = FALSE))
      }
      if (length(ids) == 0L) stop("no NuCCORE candidates for ", spec$id)

      out <- NULL; err <- NULL
      for (acc in utils::head(ids, 8L)) {
        r <- tryCatch(ncbi_cds(acc, spec$match), error = function(e) e)
        if (!inherits(r, "error") &&
            nchar(r$sequence) >= MIN_BP && nchar(r$sequence) <= MAX_BP) {
          r$accession <- acc
          out <- r
          break
        }
        if (is.null(err) && inherits(r, "error")) err <- conditionMessage(r)
      }
      if (is.null(out)) {
        stop(if (is.null(err)) "no candidate yielded a valid CDS" else err)
      }
      out
    }, error = function(e) e)

    if (inherits(res, "error")) {
      record(spec$id, "unavailable", reason = conditionMessage(res))
    } else if (nchar(res$sequence) < MIN_BP || nchar(res$sequence) > MAX_BP) {
      record(spec$id, "unavailable",
             reason = sprintf("%s CDS is %d bp, outside %d..%d",
                              res$accession, nchar(res$sequence), MIN_BP, MAX_BP))
    } else {
      write_marker_fasta(spec$id, res$sequence, res$provenance, path)
      record(spec$id, "bundled", path, res$sequence, provenance = res$provenance)
    }
    if (SLEEP > 0) Sys.sleep(SLEEP)
  }
}

update_catalog(results, CAT_YML)
say("\nUpdated ", CAT_YML)
n_ok <- sum(vapply(results, function(r) identical(r$status, "bundled"), logical(1)))
say(sprintf("Summary: %d/%d markers bundled with sequence + provenance",
            n_ok, length(results)))
quit(status = 0L)

