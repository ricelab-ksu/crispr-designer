# ============================================================================
# R/backbone_registry.R — Pluggable backbone registry (any plasmid, not just pJAT)
# ============================================================================
# A registry entry describes ONE plasmid backbone as a container of named
# "cloning slots" (gRNA array, donor payload with homology arms, selectable
# marker swap, restriction site).  The registry is the abstraction that lets
# Plasmid Designer accept *any* GenBank/FASTA backbone: pJAT is just the first
# (default) adapter.
#
# Design invariant (shared with plasmid_assembly.R): every nucleotide in an
# assembly is either copied verbatim from the backbone/genome or supplied by
# the user.  Nothing is ever fabricated.
#
# Data structures (all exported constructors validate):
#   backbone_entry()   list(id, name, source_path, format, topology, sequence,
#                        size_bp, features, slots, provenance)
#   backbone_slot()    one row: slot_id, slot_type, label, start, end, strand,
#                        enzyme, notes   (slot_type in {grna_array, donor_payload,
#                        marker_swap, restriction_site})
#   backbone_registry() named list of backbone_entry objects
#
# Contract: backbone_registry/1.0
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
BACKBONE_REGISTRY_SCHEMA <- "backbone_registry/1.0"

VALID_SLOT_TYPES <- c("grna_array", "donor_payload", "marker_swap",
                      "restriction_site", "recombination_site")

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) == 1L && !is.na(a)) a else b
  }
}

# ---- Project-root resolution (same pattern as other modules) ----------------

.backbone_registry_root <- function() {
  candidates <- unique(c(
    getwd(),
    app_root()
  ))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1L])
    if (file.exists(script_path)) {
      candidates <- unique(c(dirname(normalizePath(script_path,
                                             mustWork = FALSE)), candidates))
    }
  }
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  candidates[1L]
}

BACKBONE_REGISTRY_ROOT <- .backbone_registry_root()

# ---- Constructors -----------------------------------------------------------

# Canonical empty slots data.frame (keeps column contract everywhere).
backbone_empty_slots <- function() {
  data.frame(slot_id = character(0), slot_type = character(0),
             label = character(0), start = integer(0), end = integer(0),
             strand = character(0), enzyme = character(0),
             notes = character(0), stringsAsFactors = FALSE)
}

#' Build one validated cloning slot row.
#'
#' @param slot_id Character scalar, unique within one backbone.
#' @param slot_type One of `VALID_SLOT_TYPES`.
#' @param label Human-readable label.
#' @param start,end Integer 1-based inclusive bounds of the slot (NA if unknown).
#' @param strand "+" or "-" or NA.
#' @param enzyme For restriction_site slots: the enzyme name, else NA.
#' @param notes Free-text notes.
#' @return One-row data.frame.
#' @export
backbone_slot <- function(slot_id, slot_type, label = NA_character_,
                          start = NA_integer_, end = NA_integer_,
                          strand = NA_character_, enzyme = NA_character_,
                          notes = NA_character_) {
  stopifnot(is.character(slot_id), length(slot_id) == 1L, nzchar(slot_id),
            is.character(slot_type), length(slot_type) == 1L,
            slot_type %in% VALID_SLOT_TYPES)
  data.frame(
    slot_id   = slot_id,
    slot_type = slot_type,
    label     = as.character(label %||% NA_character_),
    start     = as.integer(start),
    end       = as.integer(end),
    strand    = as.character(strand %||% NA_character_),
    enzyme    = as.character(enzyme %||% NA_character_),
    notes     = as.character(notes %||% NA_character_),
    stringsAsFactors = FALSE
  )
}
#' Build one validated backbone registry entry.
#'
#' @param id,character scalar unique id.
#' @param name Human-readable plasmid name.
#' @param source_path Path to the source file (GenBank/FASTA).
#' @param format "genbank" or "fasta".
#' @param topology "circular" or "linear".
#' @param sequence Full nucleotide sequence (ACGTN only).
#' @param size_bp Integer length.
#' @param features data.frame from read_pjat_genbank() (may be empty for FASTA).
#' @param slots data.frame of backbone_slot rows (may be empty).
#' @param provenance Named list of provenance metadata.
#' @return list (backbone entry).
#' @export
backbone_entry <- function(id, name, source_path = NA_character_,
                           format = c("genbank", "fasta"),
                           topology = c("circular", "linear"),
                           sequence, size_bp = NA_integer_,
                           features = NULL, slots = NULL,
                           provenance = list()) {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id),
            is.character(name), length(name) == 1L,
            is.character(sequence), length(sequence) == 1L,
            nzchar(sequence))
  format <- match.arg(format)
  topology <- match.arg(topology)
  sequence <- toupper(gsub("[^ACGTN]", "", toupper(sequence)))
  if (!nzchar(sequence)) stop("backbone sequence is empty after cleaning.")
  if (is.na(size_bp)) size_bp <- nchar(sequence)
  size_bp <- as.integer(size_bp)
  stopifnot(length(size_bp) == 1L, size_bp == nchar(sequence))
  if (is.null(features)) {
    features <- data.frame(key = character(0), start = integer(0),
                           end = integer(0), strand = character(0),
                           label = character(0), product = character(0),
                           qualifiers = character(0),
                           stringsAsFactors = FALSE)
  }
  if (is.null(slots)) {
    slots <- backbone_empty_slots()
  } else {
    stopifnot(is.data.frame(slots),
              all(c("slot_id", "slot_type") %in% names(slots)))
  }
  list(
    schema        = BACKBONE_REGISTRY_SCHEMA,
    id            = id,
    name          = name,
    source_path   = as.character(source_path),
    format        = format,
    topology      = topology,
    sequence      = sequence,
    size_bp       = size_bp,
    features      = features,
    slots         = slots,
    provenance    = provenance
  )
}

#' Create a named backbone registry (list of backbone_entry).
#'
#' @param ... Named backbone_entry objects.
#' @return named list (class "backbone_registry").
#' @export
backbone_registry <- function(...) {
  entries <- list(...)
  if (length(entries) && is.null(names(entries))) {
    stop("backbone_registry() entries must be named (id = entry).")
  }
  ids <- vapply(entries, function(e) e$id, character(1))
  if (anyDuplicated(ids) || anyDuplicated(names(entries))) {
    dup <- unique(c(ids[duplicated(ids)], names(entries)[duplicated(names(entries))]))
    stop("duplicate backbone id(s): ", paste(dup, collapse = ", "))
  }
  structure(entries, class = "backbone_registry")
}

#' Print method for backbone_registry.
#' @export
print.backbone_registry <- function(x, ...) {
  cat(sprintf("backbone registry (%d backbone%s):\n",
              length(x), if (length(x) != 1L) "s" else ""))
  for (e in x) {
    cat(sprintf("  %-10s %-20s %s %5d bp  (%d slot%s)\n",
                e$id, e$name, e$format, e$size_bp,
                nrow(e$slots), if (nrow(e$slots) != 1L) "s" else ""))
  }
  invisible(x)
}

# ---- Accessors --------------------------------------------------------------

#' Get one entry by id (or name).
#' @param registry backbone_registry object.
#' @param id Backbone id or name.
#' @return backbone_entry or NULL.
#' @export
backbone_get <- function(registry, id) {
  stopifnot(inherits(registry, "backbone_registry"))
  if (is.null(id) || length(id) != 1L || !nzchar(id)) return(NULL)
  if (id %in% names(registry)) return(registry[[id]])
  hit <- vapply(registry, function(e) identical(e$name, id), logical(1))
  if (any(hit)) return(registry[hit][[1L]])
  NULL
}

#' Find slots of a given type in one backbone.
#' @param registry backbone_registry object.
#' @param id Backbone id or name.
#' @param slot_type One of VALID_SLOT_TYPES (NULL = all).
#' @return data.frame of matching slot rows (possibly empty).
#' @export
backbone_slots <- function(registry, id, slot_type = NULL) {
  e <- backbone_get(registry, id)
  if (is.null(e)) return(data.frame())
  s <- e$slots
  if (nrow(s) == 0L) return(s)
  if (!is.null(slot_type)) {
    s <- s[s$slot_type == slot_type, , drop = FALSE]
  }
  s
}

# ---- Feature -> slot discovery ----------------------------------------------

# Is a feature key/label plausibly a U6 / dU6 Pol III promoter (gRNA array)?
#
# Deliberately narrow: only Pol III U6-class promoters drive a gRNA array.
# Broad "promoter" matching would sweep in T7/SP6/lac/AmpR promoters, which are
# NOT gRNA-array slots (a real precision bug worth avoiding in the UI).
.is_u6_promoter <- function(key, label, product) {
  txt <- tolower(paste(key, label, product, sep = " "))
  grepl("\\bu6\\b|du6|u6[-_]?[0-9]|u6[-_]?promoter|pol\\s?iii|snrna:u6", txt)
}

# Is a feature plausibly a selectable / resistance marker (marker swap)?
#
# Matches marker CDS/promoter features only. A bare "promoter" is not enough:
# the feature must name a marker or resistance determinant.
.is_marker_cds <- function(key, label, product) {
  txt <- tolower(paste(key, label, product, sep = " "))
  marker_tok <- "ampr|amp[-_]?r|ampr|chlor|cmr|kanr|kan[-_]?r|hygr|hyg|puro|bsd|blast|aad|aph|resistance|ccdb|sacb|dsred|egfp|ecfp|eyfp|gfp|marker|neo|g418|kanmx|ura3|leu2"
  if (!grepl(marker_tok, txt)) return(FALSE)
  # Exclude promoter-only features that merely sit next to a marker gene.
  key_l <- tolower(as.character(key[1]))
  !grepl("^promoter$|terminator|primer_bind|^misc", key_l)
}

# Is a feature key/label an attR recombination site?
.is_att_site <- function(key, label) {
  txt <- tolower(paste(key, label, sep = " "))
  grepl("attr[3-6]|attp|attl", txt)
}

#' Auto-discover cloning slots from a parsed GenBank feature table.
#'
#' Heuristics (never fabricates; only annotates what is already annotated):
#'   * dU6/U6/promoter features            -> grna_array slot
#'   * attR3/4/5/6 / attP / attL features  -> recombination_site slot
#'   * the window attR4_end+1 .. attR5_start-1  -> donor_payload slot
#'   * marker / resistance CDS features    -> marker_swap slot
#'   * Type IIS sites (BbsI/BsaI/BsmBI)    -> restriction_site slot
#'
#' @param features data.frame (key, start, end, strand, label, product).
#' @param sequence Character(1) backbone sequence (for restriction scan).
#' @param topology "circular" or "linear".
#' @return data.frame of backbone_slot rows.
#' @export
backbone_discover_slots <- function(features, sequence, topology = "circular") {
  topology <- match.arg(topology)
  slots <- list()
  n <- 1L

  add <- function(slot_id, slot_type, label, start, end, strand,
                  enzyme = NA_character_, notes = NA_character_) {
    slots[[n]] <<- backbone_slot(slot_id, slot_type, label,
                                 as.integer(start), as.integer(end),
                                 as.character(strand), enzyme,
                                 as.character(notes))
    n <<- n + 1L
  }

  if (!is.null(features) && nrow(features) > 0L) {
    feats <- features
    if (!"product" %in% names(feats)) feats$product <- ""
    if (!"qualifiers" %in% names(feats)) feats$qualifiers <- ""

    for (i in seq_len(nrow(feats))) {
      f <- feats[i, , drop = FALSE]
      key <- as.character(f$key[1])
      label <- as.character(f$label[1] %||% "")
      product <- as.character(f$product[1] %||% "")
      s <- as.integer(f$start[1]); e <- as.integer(f$end[1])
      strd <- as.character(f$strand[1] %||% "+")

      if (.is_att_site(key, label)) {
        add(paste0("att_", tolower(label)),
            "recombination_site",
            paste0("recombination site ", label),
            if (is.na(s)) NA_integer_ else s,
            if (is.na(e)) NA_integer_ else e,
            strd, NA_character_,
            "Gateway / recombination site from backbone annotation")
      } else if (.is_u6_promoter(key, label, product)) {
        add(paste0("grna_", n),
            "grna_array",
            paste0("gRNA array / promoter: ", label),
            s, e, strd, NA_character_,
            "Pol III promoter for gRNA expression")
      } else if (.is_marker_cds(key, label, product)) {
        add(paste0("marker_", n),
            "marker_swap",
            paste0("selectable marker: ", label),
            s, e, strd, NA_character_,
            "Selectable/resistance marker cassette")
      }
    }

    # Donor-payload window between attR4 and attR5 (the classic pJAT linker).
    #
    # On inverted-cassette plasmids (pJAT39/43/77/18/21...) attR4 and attR5 sit on
    # opposite strands, so the naive attR4_end+1..attR5_start-1 window is empty
    # or nonsensical. In that case the payload window is the segment BETWEEN the
    # two sites along the strand of attR4: from attR4's end (or start, if it runs
    # the other way) to attR5's near edge. We take the interval between the two
    # sites on the attR4 strand and validate it is non-empty and in-range.
    att <- tryCatch(find_pjat_recombination_sites(
                        list(features = feats, sequence = sequence)),
                    error = function(e) NULL)
    if (is.data.frame(att) && nrow(att) > 0L &&
        all(c("attR4", "attR5") %in% att$att)) {
      a4 <- att[att$att == "attR4", , drop = FALSE]
      a5 <- att[att$att == "attR5", , drop = FALSE]
      if (nrow(a4) && nrow(a5)) {
        same_strand <- identical(as.character(a4$strand[1]),
                                 as.character(a5$strand[1]))
        start <- a4$end[1] + 1L
        end   <- a5$start[1] - 1L
        if (!same_strand || is.na(start) || is.na(end) || start > end) {
          # Fall back to the interval strictly between the two sites, taken on
          # whichever side leaves a positive-length window.
          if (a4$start[1] < a5$start[1]) {
            start <- a4$end[1] + 1L; end <- a5$start[1] - 1L
          } else {
            start <- a5$end[1] + 1L; end <- a4$start[1] - 1L
          }
          note <- "interval between attR4 and attR5 (inverted cassette)"
        } else {
          note <- "attR4_end+1 .. attR5_start-1 linker window"
        }
        if (!is.na(start) && !is.na(end) && start <= end &&
            start >= 1L && end <= nchar(sequence)) {
          add("donor_payload_linker", "donor_payload",
              "donor payload / homology-arm insertion window",
              start, end, as.character(a4$strand[1] %||% "+"),
              NA_character_, note)
        }
      }
    }
  }

  # Restriction-site slots from the sequence (Type IIS, never auto-selected).
  if (exists("find_cloning_sites", mode = "function") && nzchar(sequence)) {
    sites <- tryCatch(find_cloning_sites(sequence), error = function(e) NULL)
    if (is.data.frame(sites) && nrow(sites) > 0L) {
      for (i in seq_len(nrow(sites))) {
        site <- sites[i, , drop = FALSE]
        add(paste0("restriction_", site$enzyme[1], "_", n),
            "restriction_site",
            paste0(site$enzyme[1], " site at ", site$position[1]),
            as.integer(site$position[1]),
            as.integer(site$position[1]) + nchar(as.character(site$site[1])) - 1L,
            "+", as.character(site$enzyme[1]),
            "Type IIS cloning site (user must confirm cut position)")
      }
    }
  }

  if (length(slots) == 0L) {
    return(backbone_empty_slots())
  }
  do.call(rbind, slots)
}

# ---- Registration from files ------------------------------------------------

#' Register a backbone from a GenBank file (reuses read_pjat_genbank()).
#'
#' @param file Path to a .gbk/.gb/.gbff file.
#' @param id Backbone id (defaults to file stem).
#' @param name Name (defaults to id).
#' @param topology "circular" or "linear".
#' @param slots "auto" (discover) or a pre-built slots data.frame or NULL.
#' @param provenance Named list of extra provenance metadata.
#' @return backbone_entry.
#' @export
backbone_register_gbk <- function(file, id = NULL, name = NULL,
                                  topology = c("circular", "linear"),
                                  slots = "auto",
                                  provenance = list()) {
  topology <- match.arg(topology)
  if (!file.exists(file)) stop("GenBank file not found: ", file)
  if (!exists("read_pjat_genbank", mode = "function")) {
    stop("read_pjat_genbank() is required; source R/pjat_plasmids.R first.")
  }
  gb <- read_pjat_genbank(file)
  if (is.null(id))   id   <- sub("\\.(gbk|gb|gbff|genbank)$", "", basename(file),
                                 ignore.case = TRUE)
  if (is.null(name)) name <- id
  if (isTRUE(slots == "auto")) {
    slots <- backbone_discover_slots(gb$features, gb$sequence, topology)
  } else if (is.null(slots)) {
    slots <- backbone_empty_slots()
  }
  backbone_entry(id = id, name = name, source_path = normalizePath(file),
                 format = "genbank", topology = topology,
                 sequence = gb$sequence, size_bp = gb$size_bp %||%
                 nchar(gb$sequence), features = gb$features,
                 slots = slots,
                 provenance = utils::modifyList(provenance,
                           list(locus = gb$locus %||% NA_character_,
                                definition = gb$definition %||% NA_character_)))
}

#' Register a backbone from a FASTA file (sequence-only; no features).
#'
#' @param file Path to a .fa/.fasta file.
#' @param id Backbone id (defaults to file stem).
#' @param name Name (defaults to id).
#' @param topology "circular" or "linear".
#' @param provenance Named list of extra provenance metadata.
#' @return backbone_entry.
#' @export
backbone_register_fasta <- function(file, id = NULL, name = NULL,
                                    topology = c("circular", "linear"),
                                    provenance = list()) {
  topology <- match.arg(topology)
  if (!file.exists(file)) stop("FASTA file not found: ", file)
  lines <- readLines(file, warn = FALSE)
  lines <- lines[!grepl("^>", lines)]
  sequence <- toupper(gsub("[^ACGTN]", "", paste(lines, collapse = "")))
  if (!nzchar(sequence)) stop("FASTA contains no sequence: ", file)
  if (is.null(id))   id   <- sub("\\.(fa|fasta|fa\\.gz|fasta\\.gz)$", "",
                                 basename(file), ignore.case = TRUE)
  if (is.null(name)) name <- id
  backbone_entry(id = id, name = name, source_path = normalizePath(file),
                 format = "fasta", topology = topology,
                 sequence = sequence, size_bp = nchar(sequence),
                 features = data.frame(), slots = backbone_empty_slots(),
                 provenance = provenance)
}

# ---- pJAT adapter (pJAT becomes the first registry entry) -------------------

#' Build a backbone registry from the installed pJAT collection.
#'
#' Reads `donor_templates/pjat_plasmids.yml`, then loads each referenced
#' GenBank file via `backbone_register_gbk()`.  This makes pJAT the default
#' registry entry with behavior identical to the current `pjat_load_catalog()`,
#' while adding auto-discovered cloning slots.
#'
#' @param yml Path to the pJAT catalog YAML.
#' @param ref_dir Directory containing the referenced GenBank files.
#' @return backbone_registry (named by plasmid id, e.g. "pJAT74").
#' @export
backbone_pjat_adapter <- function(
  yml = file.path(BACKBONE_REGISTRY_ROOT, "donor_templates", "pjat_plasmids.yml"),
  ref_dir = NULL) {
  if (!file.exists(yml)) stop("pJAT catalog not found: ", yml)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read the pJAT catalog.")
  }
  if (is.null(ref_dir)) {
    y <- yaml::read_yaml(yml)
    ref_dir <- y$source_dir %||%
      get_app_pjat_dir("PLASMID")
  }
  cat <- pjat_load_catalog(yml = yml)
  entries <- list()
  for (i in seq_len(nrow(cat))) {
    r <- cat[i, , drop = FALSE]
    gb_file <- file.path(ref_dir, r$source_file)
    if (!file.exists(gb_file)) {
      warning("skipping ", r$name, ": GenBank not found at ", gb_file)
      next
    }
    e <- tryCatch(
      backbone_register_gbk(gb_file, id = r$name, name = r$name,
                            topology = "circular",
                            provenance = list(
                              addgene_id = r$addgene_id %||% NA_integer_,
                              supplier   = r$supplier %||% NA_character_,
                              purpose    = r$purpose %||% NA_character_,
                              u6_promoters = r$u6_promoters %||% NA_integer_,
                              has_attp   = r$has_attp %||% NA)),
      error = function(err) {
        warning("skipping ", r$name, ": ", conditionMessage(err))
        NULL
      })
    if (!is.null(e)) entries[[r$name]] <- e
  }
  if (length(entries) == 0L) stop("no pJAT backbones could be registered.")
  do.call(backbone_registry, entries)
}

#' Default backbone registry (pJAT adapter unless a user registry exists).
#'
#' @return backbone_registry.
#' @export
backbone_registry_default <- function() {
  backbone_pjat_adapter()
}

