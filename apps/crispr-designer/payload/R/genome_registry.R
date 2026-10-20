# ============================================================================
# R/genome_registry.R — Drosophila Genome Registry, FAI Indexing & Fast Subseq
# ============================================================================
# Supports all 7 Drosophila species in reference_files/drosophila_species_genomes:
#   - d_melanogaster : GCF_000001215.4_Release_6_plus_ISO1_MT
#   - d_simulans     : GCF_016746395.2_Prin_Dsim_3.1
#   - d_eugracilis   : GCF_018153835.1_ASM1815383v1
#   - d_erecta       : GCF_003286155.1_DereRS2
#   - d_teissieri    : GCF_016746235.2_Prin_Dtei_1.1
#   - d_ananassae    : GCF_017639315.1_ASM1763931v2
#   - d_virilis      : GCF_003285735.1_DvirRS2
# ============================================================================

# The reference store is resolved at runtime (see R/app_paths.R) rather than
# hard-coded to a path on an external volume, so the desktop app works on any
# machine. Load the shared resolver when it is available.
local({
  for (p in unique(c(
    file.path("R", "app_paths.R"),
    file.path(getwd(), "R", "app_paths.R"),
    file.path(dirname(getwd()), "R", "app_paths.R")
  ))) {
    if (file.exists(p) && !exists("get_app_reference_dir", mode = "function")) {
      tryCatch(source(p, local = FALSE), error = function(e) NULL)
      break
    }
  }
})
if (exists("load_app_paths", mode = "function")) load_app_paths("CRISPR")

# Legacy name, kept for compatibility. Resolved at runtime.
DEFAULT_REFERENCE_DIR <- if (exists("get_app_reference_dir", mode = "function")) {
  get_app_reference_dir("CRISPR")
} else {
  ""
}

SPECIES_GCF <- c(
  melanogaster = "GCF_000001215.4_Release_6_plus_ISO1_MT",
  simulans     = "GCF_016746395.2_Prin_Dsim_3.1",
  eugracilis   = "GCF_018153835.1_ASM1815383v1",
  erecta       = "GCF_003286155.1_DereRS2",
  teissieri    = "GCF_016746235.2_Prin_Dtei_1.1",
  ananassae    = "GCF_017639315.1_ASM1763931v2",
  virilis      = "GCF_003285735.1_DvirRS2"
)

SPECIES_ORGANISMS <- c(
  melanogaster = "Drosophila melanogaster",
  simulans     = "Drosophila simulans",
  eugracilis   = "Drosophila eugracilis",
  erecta       = "Drosophila erecta",
  teissieri    = "Drosophila teissieri",
  tesseri      = "Drosophila teissieri",
  ananassae    = "Drosophila ananassae",
  virilis      = "Drosophila virilis"
)

SPECIES_DL_MB <- c(
  melanogaster = 73,
  simulans     = 61,
  eugracilis   = 66,
  erecta       = 62,
  teissieri    = 61,
  ananassae    = 86,
  virilis      = 75
)

KNOWN_CHR_ALIASES <- list(
  melanogaster = c(
    "X"   = "NC_004354.4",
    "2L"  = "NT_033779.5",
    "2R"  = "NT_033778.4",
    "3L"  = "NT_037436.4",
    "3R"  = "NT_033777.3",
    "4"   = "NC_004353.4",
    "Y"   = "NC_024512.1",
    "M"   = "NC_024511.2",
    "MT"  = "NC_024511.2",
    "mt"  = "NC_024511.2",
    "mitochondrion" = "NC_024511.2",
    "dmel_mitochondrion_genome" = "NC_024511.2"
  ),
  simulans = c(
    "2L"  = "NC_052520.2",
    "2R"  = "NC_052521.2",
    "3L"  = "NC_052522.2",
    "3R"  = "NC_052523.2",
    "4"   = "NC_052524.2",
    "X"   = "NC_052525.2"
  ),
  teissieri = c(
    "2L"  = "NC_053029.1",
    "2R"  = "NC_053030.1",
    "3L"  = "NC_053031.1",
    "3R"  = "NC_053032.1",
    "4"   = "NC_053033.1",
    "X"   = "NC_053034.1"
  )
)

.genome_registry_cache <- new.env(parent = emptyenv())

get_reference_dir <- function() {
  # An explicit override always wins, exactly as it did before the reference
  # store became runtime-resolved. The Electron main process sets
  # CRISPR_REFERENCE_DIR for the desktop builds; HCR_REFERENCE_DIR is honoured
  # so the two applications can be pointed at one store.
  for (ev in c("CRISPR_REFERENCE_DIR", "CAS9_REFERENCE_DIR", "HCR_REFERENCE_DIR")) {
    val <- Sys.getenv(ev, unset = "")
    if (nzchar(val)) return(val)
  }
  if (exists("get_app_reference_dir", mode = "function")) {
    return(get_app_reference_dir("CRISPR"))
  }
  DEFAULT_REFERENCE_DIR
}

canonical_species <- function(sp) {
  if (is.null(sp) || length(sp) == 0 || is.na(sp[1])) return(NA_character_)
  s <- tolower(trimws(as.character(sp[1])))
  s <- sub("^drosophila[ _\\.]+", "", s)
  s <- sub("^d[ _\\.]+", "", s)
  alias_map <- c(
    tesseri    = "teissieri",
    tessieri   = "teissieri",
    teissier   = "teissieri",
    teissierii = "teissieri",
    dmel       = "melanogaster",
    dsim       = "simulans",
    deug       = "eugracilis",
    dere       = "erecta",
    dtei       = "teissieri",
    dana       = "ananassae",
    dvir       = "virilis"
  )
  if (s %in% names(alias_map)) return(unname(alias_map[[s]]))
  s
}

organism_for <- function(species) {
  sp <- canonical_species(species)
  if (is.na(sp) || !nzchar(sp)) return(NA_character_)
  if (sp %in% names(SPECIES_ORGANISMS)) return(unname(SPECIES_ORGANISMS[[sp]]))
  if (grepl(" ", trimws(species))) return(trimws(species))
  paste("Drosophila", sp)
}

species_ref_dir <- function(species, ref_dir = get_reference_dir()) {
  sp <- canonical_species(species)
  file.path(ref_dir, paste0("d_", sp))
}

.list_valid_files <- function(sd, pattern) {
  if (!dir.exists(sd)) return(character(0))
  hits <- list.files(sd, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (!length(hits)) return(character(0))
  base <- basename(hits)
  hits <- hits[!grepl("^\\.", base)]
  hits
}

.find_species_dirs <- function(species, ref_dir = get_reference_dir()) {
  if (!dir.exists(ref_dir)) return(character(0))
  key <- canonical_species(species)
  if (is.na(key) || !nzchar(key)) return(character(0))
  exact_dir <- file.path(ref_dir, paste0("d_", key))
  subs <- list.dirs(ref_dir, recursive = FALSE, full.names = TRUE)
  subs <- subs[!grepl("^\\.", basename(subs))]
  matched <- subs[grepl(key, tolower(basename(subs)), fixed = TRUE)]
  if (dir.exists(exact_dir)) {
    matched <- unique(c(exact_dir, matched))
  }
  matched
}

find_genome_file <- function(species, ref_dir = get_reference_dir()) {
  subs <- .find_species_dirs(species, ref_dir = ref_dir)
  for (sd in subs) {
    hits <- .list_valid_files(sd, "(^|_)genomic\\.fna(\\.gz)?$")
    hits <- hits[!grepl("cds_from|rna_from", basename(hits), ignore.case = TRUE)]
    if (length(hits)) {
      is_gz <- grepl("\\.gz$", hits, ignore.case = TRUE)
      ord <- order(is_gz, nchar(basename(hits)))
      return(hits[ord[1]])
    }
  }
  NULL
}

find_gtf_file <- function(species, ref_dir = get_reference_dir()) {
  subs <- .find_species_dirs(species, ref_dir = ref_dir)
  for (sd in subs) {
    hits <- .list_valid_files(sd, "(^|_)genomic\\.(gtf|gff)(\\.gz)?$")
    if (length(hits)) {
      b <- tolower(basename(hits))
      is_gz <- grepl("\\.gz$", b)
      is_gff <- grepl("\\.gff", b)
      is_exact_gtf <- b %in% c("genomic.gtf", "genomic.gtf.gz")
      ord <- order(is_gz, is_gff, !is_exact_gtf, nchar(b))
      return(hits[ord[1]])
    }
  }
  NULL
}

find_rna_file <- function(species, ref_dir = get_reference_dir()) {
  subs <- .find_species_dirs(species, ref_dir = ref_dir)
  for (sd in subs) {
    hits <- .list_valid_files(sd, "(^|_)rna\\.fna(\\.gz)?$")
    if (length(hits)) {
      is_gz <- grepl("\\.gz$", hits, ignore.case = TRUE)
      ord <- order(is_gz, nchar(basename(hits)))
      return(hits[ord[1]])
    }
  }
  NULL
}

find_gene_ref_file <- function(species, ref_dir = get_reference_dir()) {
  subs <- .find_species_dirs(species, ref_dir = ref_dir)
  for (sd in subs) {
    hits <- .list_valid_files(sd, "(^|_)gene_reference\\.csv$")
    if (length(hits)) {
      ord <- order(nchar(basename(hits)))
      return(hits[ord[1]])
    }
  }
  NULL
}

get_species_registry <- function(ref_dir = get_reference_dir()) {
  known_keys <- names(SPECIES_GCF)
  extra_keys <- character(0)
  if (dir.exists(ref_dir)) {
    dirs <- list.dirs(ref_dir, recursive = FALSE, full.names = FALSE)
    dirs <- dirs[grepl("^d_[a-z0-9_]+$", dirs, ignore.case = TRUE)]
    detected <- vapply(dirs, canonical_species, character(1))
    extra_keys <- setdiff(unique(detected), known_keys)
  }
  all_keys <- c(known_keys, extra_keys)

  rows <- lapply(all_keys, function(sp) {
    fna <- find_genome_file(sp, ref_dir = ref_dir)
    gtf <- find_gtf_file(sp, ref_dir = ref_dir)
    rna <- find_rna_file(sp, ref_dir = ref_dir)
    gref <- find_gene_ref_file(sp, ref_dir = ref_dir)

    fai <- NA_character_
    if (!is.null(fna) && !grepl("\\.gz$", fna, ignore.case = TRUE)) {
      cand_fai <- paste0(fna, ".fai")
      if (file.exists(cand_fai) && file.size(cand_fai) > 0) {
        fai <- cand_fai
      }
    }

    status <- if (is.null(fna) || is.null(gtf)) {
      "needs_download"
    } else if (grepl("\\.gz$", fna, ignore.case = TRUE) ||
               grepl("\\.gz$", gtf, ignore.case = TRUE)) {
      "compressed"
    } else {
      "ready"
    }

    acc <- if (sp %in% names(SPECIES_GCF)) unname(SPECIES_GCF[[sp]]) else NA_character_

    data.frame(
      species_key        = sp,
      folder_name        = paste0("d_", sp),
      organism_name      = organism_for(sp),
      assembly_accession = acc,
      genomic_fna_path   = if (!is.null(fna)) fna else NA_character_,
      genomic_gtf_path   = if (!is.null(gtf)) gtf else NA_character_,
      rna_fna_path       = if (!is.null(rna)) rna else NA_character_,
      fai_path           = fai,
      gene_ref_path      = if (!is.null(gref)) gref else NA_character_,
      status             = status,
      stringsAsFactors   = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- out$species_key
  out
}

download_species_ref <- function(species,
                                 with_genome = TRUE,
                                 ref_dir = get_reference_dir(),
                                 progress = message,
                                 download_fn = utils::download.file) {
  old_to <- getOption("timeout")
  on.exit(options(timeout = old_to), add = TRUE)
  options(timeout = max(900, as.numeric(old_to)))

  sp <- canonical_species(species)
  gcf <- if (!is.na(sp) && sp %in% names(SPECIES_GCF)) SPECIES_GCF[[sp]] else NULL
  if (is.null(gcf)) {
    return(paste0("No known assembly for '", species, "'."))
  }

  acc <- sub("^(GCF_[0-9]+)\\..*$", "\\1", gcf)
  nums <- sub("^GCF_", "", acc)
  segs <- substring(nums, c(1, 4, 7), c(3, 6, 9))
  base <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/",
    segs[1], "/", segs[2], "/", segs[3], "/", gcf, "/"
  )
  dest <- species_ref_dir(sp, ref_dir = ref_dir)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  components <- list(
    rna            = paste0(gcf, "_rna.fna.gz"),
    genome         = paste0(gcf, "_genomic.fna.gz"),
    annotation_gtf = paste0(gcf, "_genomic.gtf.gz"),
    annotation_gff = paste0(gcf, "_genomic.gff.gz")
  )
  if (!with_genome) {
    components <- components["rna"]
  }

  saved <- character(0)
  for (nm in names(components)) {
    f <- switch(
      nm,
      rna            = "rna.fna.gz",
      genome         = "genomic.fna.gz",
      annotation_gtf = "genomic.gtf.gz",
      annotation_gff = "genomic.gff.gz"
    )
    out <- file.path(dest, f)

    existing <- switch(
      nm,
      rna            = find_rna_file(sp, ref_dir = ref_dir),
      genome         = find_genome_file(sp, ref_dir = ref_dir),
      annotation_gtf = {
        g <- find_gtf_file(sp, ref_dir = ref_dir)
        if (!is.null(g) && grepl("\\.gtf(\\.gz)?$", g, ignore.case = TRUE)) g else NULL
      },
      annotation_gff = {
        g <- find_gtf_file(sp, ref_dir = ref_dir)
        if (!is.null(g) && grepl("\\.gff(\\.gz)?$", g, ignore.case = TRUE)) g else NULL
      }
    )
    if ((file.exists(out) && file.size(out) > 0) ||
        (!is.null(existing) && file.exists(existing) && file.size(existing) > 0)) {
      saved <- c(saved, if (file.exists(out)) f else basename(existing))
      next
    }

    part <- paste0(out, ".part")
    on.exit(if (file.exists(part)) unlink(part), add = TRUE)
    url <- paste0(base, components[[nm]])
    if (is.function(progress)) {
      progress(paste0("Downloading ", sp, " ", components[[nm]], " ..."))
    }
    ok <- tryCatch({
      download_fn(url, part, mode = "wb", quiet = TRUE)
      file.exists(part) && file.size(part) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)

    if (!ok) {
      if (file.exists(part)) unlink(part)
      if (nm == "rna") {
        return(paste0(
          "Download failed for ", sp, " (rna.fna.gz). Check your ",
          "internet connection, or download manually from ", base
        ))
      }
      if (is.function(progress)) {
        progress(paste0("  (", components[[nm]], " unavailable for ", sp, "; continuing)"))
      }
      next
    }
    file.rename(part, out)
    saved <- c(saved, f)
  }

  paste0(sp, " reference saved to ", dest, " (", paste(saved, collapse = ", "), ")")
}

decompress_gzip_file <- function(gz_path,
                                 dest_path = sub("\\.gz$", "", gz_path, ignore.case = TRUE),
                                 overwrite = FALSE,
                                 remove_gz = FALSE,
                                 chunk_size = 16777216L) {
  if (!file.exists(gz_path)) {
    stop("Gzip file does not exist: ", gz_path)
  }
  if (!overwrite && file.exists(dest_path) && file.size(dest_path) > 0) {
    return(dest_path)
  }
  part_path <- paste0(dest_path, ".part")
  on.exit(if (file.exists(part_path)) unlink(part_path), add = TRUE)

  in_con <- gzfile(gz_path, open = "rb")
  on.exit(try(close(in_con), silent = TRUE), add = TRUE)
  out_con <- file(part_path, open = "wb")
  on.exit(try(close(out_con), silent = TRUE), add = TRUE)

  repeat {
    buf <- readBin(in_con, what = "raw", n = chunk_size)
    if (length(buf) == 0L) break
    writeBin(buf, out_con)
  }
  close(out_con)
  on.exit(NULL)
  try(close(in_con), silent = TRUE)

  if (!file.rename(part_path, dest_path)) {
    file.copy(part_path, dest_path, overwrite = TRUE)
    unlink(part_path)
  }
  if (remove_gz && file.exists(dest_path) && file.size(dest_path) > 0) {
    unlink(gz_path)
  }
  dest_path
}

ensure_fai_index <- function(fna_path,
                             fai_path = paste0(fna_path, ".fai"),
                             force = FALSE) {
  if (is.null(fna_path) || is.na(fna_path) || !file.exists(fna_path)) {
    stop("FASTA file does not exist: ", fna_path)
  }
  if (grepl("\\.gz$", fna_path, ignore.case = TRUE)) {
    stop("Cannot build .fai index directly on a gzip-compressed FASTA: ", fna_path)
  }
  if (!force && file.exists(fai_path) && file.size(fai_path) > 0) {
    return(fai_path)
  }

  con <- file(fna_path, open = "rb")
  on.exit(close(con), add = TRUE)

  chunk_size <- 33554432L
  file_pos <- 0
  leftover <- raw(0)

  names_vec <- character(0)
  lengths_vec <- numeric(0)
  offsets_vec <- numeric(0)
  linebases_vec <- integer(0)
  linewidth_vec <- integer(0)

  cur_name <- NULL
  cur_len <- 0
  cur_offset <- 0
  cur_linebases <- 0L
  cur_linewidth <- 0L

  nl_raw <- as.raw(0x0a)
  cr_raw <- as.raw(0x0d)
  gt_raw <- as.raw(0x3e)

  process_lines <- function(raw_vec, base_offset, is_eof = FALSE) {
    nl_idx <- which(raw_vec == nl_raw)
    if (!length(nl_idx) && !is_eof) {
      return(raw_vec)
    }

    line_starts <- c(1L, nl_idx + 1L)
    line_ends <- c(nl_idx, if (is_eof && (length(nl_idx) == 0L || nl_idx[length(nl_idx)] < length(raw_vec))) length(raw_vec) else integer(0))

    n_complete <- length(line_ends)
    if (n_complete > 0L) {
      ls_vec <- line_starts[seq_len(n_complete)]
      le_vec <- line_ends[seq_len(n_complete)]
      first_bytes <- raw_vec[ls_vec]
      has_nl <- (raw_vec[le_vec] == nl_raw)
      has_cr <- logical(n_complete)
      gt1 <- which(le_vec > ls_vec)
      if (length(gt1)) {
        has_cr[gt1] <- (raw_vec[le_vec[gt1] - 1L] == cr_raw)
      }
      term_len <- as.integer(has_nl) + as.integer(has_cr)
      total_w <- le_vec - ls_vec + 1L
      content_w <- total_w - term_len

      hdr_pos <- which(first_bytes == gt_raw)
      if (!length(hdr_pos)) {
        if (!is.null(cur_name)) {
          if (cur_linebases == 0L) {
            non_empty <- which(content_w > 0L)
            if (length(non_empty)) {
              idx1 <- non_empty[1L]
              cur_offset <<- base_offset + (ls_vec[idx1] - 1)
              cur_linebases <<- content_w[idx1]
              cur_linewidth <<- total_w[idx1]
            }
          }
          cur_len <<- cur_len + sum(as.numeric(content_w))
        }
      } else {
        boundaries <- c(0L, hdr_pos, n_complete + 1L)
        for (b_idx in seq_len(length(boundaries) - 1L)) {
          seg_start <- boundaries[b_idx]
          seg_end <- boundaries[b_idx + 1L] - 1L

          if (seg_start > 0L) {
            if (!is.null(cur_name)) {
              names_vec <<- c(names_vec, cur_name)
              lengths_vec <<- c(lengths_vec, cur_len)
              offsets_vec <<- c(offsets_vec, cur_offset)
              linebases_vec <<- c(linebases_vec, cur_linebases)
              linewidth_vec <<- c(linewidth_vec, cur_linewidth)
            }
            h_start <- ls_vec[seg_start] + 1L
            h_end <- ls_vec[seg_start] + content_w[seg_start] - 1L
            hdr_str <- if (h_end >= h_start) rawToChar(raw_vec[h_start:h_end]) else ""
            hdr_str <- trimws(hdr_str)
            cur_name <<- sub("\\s.*$", "", hdr_str)
            cur_len <<- 0
            cur_offset <<- base_offset + le_vec[seg_start]
            cur_linebases <<- 0L
            cur_linewidth <<- 0L
          }

          seq_lo <- seg_start + 1L
          if (seq_lo <= seg_end && !is.null(cur_name)) {
            rng <- seq_lo:seg_end
            cw <- content_w[rng]
            tw <- total_w[rng]
            if (cur_linebases == 0L) {
              non_empty <- which(cw > 0L)
              if (length(non_empty)) {
                idx1 <- rng[non_empty[1L]]
                cur_offset <<- base_offset + (ls_vec[idx1] - 1)
                cur_linebases <<- content_w[idx1]
                cur_linewidth <<- total_w[idx1]
              }
            }
            cur_len <<- cur_len + sum(as.numeric(cw))
          }
        }
      }
    }

    if (is_eof) {
      return(raw(0))
    }
    last_end <- line_ends[n_complete]
    if (last_end < length(raw_vec)) {
      raw_vec[(last_end + 1L):length(raw_vec)]
    } else {
      raw(0)
    }
  }

  repeat {
    buf <- readBin(con, what = "raw", n = chunk_size)
    if (length(buf) == 0L) {
      if (length(leftover) > 0L) {
        process_lines(leftover, file_pos - length(leftover), is_eof = TRUE)
      }
      break
    }
    combined <- if (length(leftover) > 0L) c(leftover, buf) else buf
    base_off <- file_pos - length(leftover)
    file_pos <- file_pos + length(buf)
    leftover <- process_lines(combined, base_off, is_eof = FALSE)
  }

  if (!is.null(cur_name)) {
    names_vec <- c(names_vec, cur_name)
    lengths_vec <- c(lengths_vec, cur_len)
    offsets_vec <- c(offsets_vec, cur_offset)
    linebases_vec <- c(linebases_vec, cur_linebases)
    linewidth_vec <- c(linewidth_vec, cur_linewidth)
  }

  fai_lines <- sprintf(
    "%s\t%.0f\t%.0f\t%d\t%d",
    names_vec, lengths_vec, offsets_vec, linebases_vec, linewidth_vec
  )

  target_fai <- tryCatch({
    part <- paste0(fai_path, ".part")
    writeLines(fai_lines, part)
    if (!file.rename(part, fai_path)) {
      file.copy(part, fai_path, overwrite = TRUE)
      unlink(part)
    }
    fai_path
  }, error = function(e) {
    fallback <- file.path(tempdir(), paste0(basename(fna_path), ".fai"))
    writeLines(fai_lines, fallback)
    fallback
  })

  cache_key <- paste0(".fai:", fna_path)
  if (exists(cache_key, envir = .genome_registry_cache, inherits = FALSE)) {
    rm(list = cache_key, envir = .genome_registry_cache)
  }

  target_fai
}

read_fai_index <- function(fai_path) {
  if (!file.exists(fai_path)) {
    stop("FAI index file does not exist: ", fai_path)
  }
  df <- utils::read.table(
    fai_path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("name", "length", "offset", "linebases", "linewidth"),
    colClasses = c("character", "numeric", "numeric", "integer", "integer"),
    quote = "",
    comment.char = ""
  )
  rownames(df) <- df$name
  df
}

ensure_uncompressed_or_indexed <- function(species,
                                           ref_dir = get_reference_dir(),
                                           decompress_gz = FALSE) {
  sp <- canonical_species(species)
  fna <- find_genome_file(sp, ref_dir = ref_dir)
  if (is.null(fna)) {
    return(list(
      species = sp,
      fna_path = NA_character_,
      fai_path = NA_character_,
      is_compressed = NA,
      status = "needs_download"
    ))
  }

  if (grepl("\\.gz$", fna, ignore.case = TRUE)) {
    if (decompress_gz) {
      fna <- decompress_gzip_file(fna)
    } else {
      return(list(
        species = sp,
        fna_path = fna,
        fai_path = NA_character_,
        is_compressed = TRUE,
        status = "compressed"
      ))
    }
  }

  fai <- ensure_fai_index(fna)
  list(
    species = sp,
    fna_path = fna,
    fai_path = fai,
    is_compressed = FALSE,
    status = "ready"
  )
}

.get_fai_table <- function(fna_path, cache = NULL) {
  env <- if (is.environment(cache)) cache else .genome_registry_cache
  key <- paste0(".fai:", fna_path)
  if (exists(key, envir = env, inherits = FALSE)) {
    return(get(key, envir = env, inherits = FALSE))
  }
  fai_path <- ensure_fai_index(fna_path)
  df <- read_fai_index(fai_path)
  assign(key, df, envir = env)
  df
}

.get_gz_genome <- function(fna_path, cache = NULL) {
  env <- if (is.environment(cache)) cache else .genome_registry_cache
  key <- paste0(".gz_genome:", fna_path)
  if (exists(key, envir = env, inherits = FALSE)) {
    return(get(key, envir = env, inherits = FALSE))
  }

  if (requireNamespace("Biostrings", quietly = TRUE)) {
    dna <- Biostrings::readDNAStringSet(fna_path)
    full_headers <- names(dna)
    ids <- sub("\\s.*$", "", full_headers)
    seqs <- as.character(dna)
    names(seqs) <- ids
    names(full_headers) <- ids
  } else {
    con <- if (grepl("\\.gz$", fna_path, ignore.case = TRUE)) gzfile(fna_path, "rt") else file(fna_path, "rt")
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    hdr_idx <- grep("^>", lines)
    ends <- c(hdr_idx[-1] - 1L, length(lines))
    full_headers <- sub("^>", "", lines[hdr_idx])
    ids <- sub("\\s.*$", "", trimws(full_headers))
    seqs <- vapply(seq_along(hdr_idx), function(i) {
      if (hdr_idx[i] + 1L > ends[i]) return("")
      toupper(gsub("[^A-Za-z]", "", paste0(lines[(hdr_idx[i] + 1L):ends[i]], collapse = "")))
    }, character(1))
    names(seqs) <- ids
    names(full_headers) <- ids
  }

  obj <- list(seqs = seqs, headers = full_headers, lengths = nchar(seqs))
  assign(key, obj, envir = env)
  obj
}

.read_fasta_headers <- function(fna_path, fai_df = NULL) {
  if (!file.exists(fna_path)) return(character(0))
  if (grepl("\\.gz$", fna_path, ignore.case = TRUE)) {
    gz_obj <- .get_gz_genome(fna_path)
    return(gz_obj$headers)
  }
  if (!is.null(fai_df) && nrow(fai_df) > 0) {
    idx <- which(grepl("^N[CT]_", fai_df$name))
    if (!length(idx)) idx <- seq_len(min(25L, nrow(fai_df)))
    con <- file(fna_path, open = "rb")
    on.exit(close(con), add = TRUE)
    hdrs <- character(length(idx))
    names(hdrs) <- fai_df$name[idx]
    for (k in seq_along(idx)) {
      i <- idx[k]
      seq_off <- fai_df$offset[i]
      read_start <- max(0, seq_off - 256)
      n_read <- as.integer(seq_off - read_start)
      seek(con, where = read_start, origin = "start")
      raw_bytes <- readBin(con, what = "raw", n = n_read)
      txt <- rawToChar(raw_bytes[raw_bytes != as.raw(0)])
      lines <- strsplit(txt, "[\r\n]+")[[1]]
      gt_lines <- grep("^>", lines, value = TRUE)
      if (length(gt_lines)) {
        hdrs[k] <- sub("^>", "", gt_lines[length(gt_lines)])
      } else {
        hdrs[k] <- fai_df$name[i]
      }
    }
    return(hdrs)
  }
  character(0)
}

get_chromosome_map <- function(species, cache = NULL, ref_dir = get_reference_dir()) {
  sp <- canonical_species(species)
  env <- if (is.environment(cache)) cache else .genome_registry_cache
  cache_key <- paste0(".chrmap:", ref_dir, ":", sp)
  if (exists(cache_key, envir = env, inherits = FALSE)) {
    return(get(cache_key, envir = env, inherits = FALSE))
  }

  alias_to_acc <- character(0)
  if (!is.na(sp) && sp %in% names(KNOWN_CHR_ALIASES)) {
    alias_to_acc <- KNOWN_CHR_ALIASES[[sp]]
  }

  fna <- find_genome_file(sp, ref_dir = ref_dir)
  all_accessions <- character(0)
  if (!is.null(fna) && file.exists(fna)) {
    if (!grepl("\\.gz$", fna, ignore.case = TRUE)) {
      fai_df <- .get_fai_table(fna, cache = env)
      all_accessions <- fai_df$name
      hdrs <- .read_fasta_headers(fna, fai_df = fai_df)
    } else {
      hdrs <- character(0)
    }

    if (length(hdrs)) {
      for (acc in names(hdrs)) {
        h <- hdrs[[acc]]
        if (grepl("\\bchromosome\\s+([A-Za-z0-9]+)\\b", h, ignore.case = TRUE)) {
          chr_name <- sub(".*\\bchromosome\\s+([A-Za-z0-9]+)\\b.*", "\\1", h, ignore.case = TRUE)
          if (!chr_name %in% names(alias_to_acc)) {
            alias_to_acc[[chr_name]] <- acc
          }
        } else if (grepl("mitochondrion", h, ignore.case = TRUE)) {
          for (m_alias in c("M", "MT", "mt", "mitochondrion")) {
            if (!m_alias %in% names(alias_to_acc)) {
              alias_to_acc[[m_alias]] <- acc
            }
          }
        }
      }
    }
  }

  base_aliases <- names(alias_to_acc)
  for (al in base_aliases) {
    if (!grepl("^chr", al, ignore.case = TRUE) && nchar(al) <= 4) {
      chr_al <- paste0("chr", al)
      if (!chr_al %in% names(alias_to_acc)) {
        alias_to_acc[[chr_al]] <- alias_to_acc[[al]]
      }
    }
  }

  acc_to_alias <- character(0)
  for (al in base_aliases) {
    acc <- alias_to_acc[[al]]
    if (!acc %in% names(acc_to_alias)) {
      acc_to_alias[[acc]] <- al
    }
  }

  unique_accs <- unique(c(unname(alias_to_acc), all_accessions))
  for (acc in unique_accs) {
    if (!acc %in% names(alias_to_acc)) {
      alias_to_acc[[acc]] <- acc
    }
    no_ver <- sub("\\.[0-9]+$", "", acc)
    if (!no_ver %in% names(alias_to_acc)) {
      alias_to_acc[[no_ver]] <- acc
    }
  }

  tbl <- data.frame(
    alias = names(alias_to_acc),
    accession = unname(alias_to_acc),
    stringsAsFactors = FALSE
  )

  res <- list(
    species = sp,
    alias_to_accession = alias_to_acc,
    accession_to_alias = acc_to_alias,
    table = tbl
  )
  assign(cache_key, res, envir = env)
  res
}

resolve_chromosome <- function(species, chr, cache = NULL, ref_dir = get_reference_dir()) {
  if (is.null(chr) || length(chr) == 0 || is.na(chr[1]) || !nzchar(trimws(chr[1]))) {
    return(NA_character_)
  }
  c_in <- trimws(as.character(chr[1]))
  cmap <- get_chromosome_map(species, cache = cache, ref_dir = ref_dir)
  a2a <- cmap$alias_to_accession

  if (c_in %in% names(a2a)) {
    return(unname(a2a[[c_in]]))
  }
  idx <- match(tolower(c_in), tolower(names(a2a)))
  if (!is.na(idx)) {
    return(unname(a2a[[idx]]))
  }
  c_stripped <- sub("^chr", "", c_in, ignore.case = TRUE)
  idx2 <- match(tolower(c_stripped), tolower(names(a2a)))
  if (!is.na(idx2)) {
    return(unname(a2a[[idx2]]))
  }
  fna <- find_genome_file(species, ref_dir = ref_dir)
  if (!is.null(fna) && grepl("\\.gz$", fna, ignore.case = TRUE)) {
    gz_obj <- .get_gz_genome(fna, cache = cache)
    if (c_in %in% names(gz_obj$seqs)) return(c_in)
    no_ver <- sub("\\.[0-9]+$", "", names(gz_obj$seqs))
    m <- match(tolower(c_in), tolower(no_ver))
    if (!is.na(m)) return(names(gz_obj$seqs)[m])
  }
  NA_character_
}

revcomp_dna <- function(seq) {
  if (is.null(seq) || !length(seq)) return(character(0))
  comp <- chartr(
    "ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
    "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn",
    seq
  )
  vapply(comp, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    rawToChar(rev(charToRaw(s)))
  }, character(1), USE.NAMES = FALSE)
}

get_chr_length <- function(species, chr, cache = NULL, ref_dir = get_reference_dir()) {
  fna <- find_genome_file(species, ref_dir = ref_dir)
  if (is.null(fna) || !file.exists(fna)) return(NA_integer_)
  acc <- resolve_chromosome(species, chr, cache = cache, ref_dir = ref_dir)
  if (is.na(acc)) return(NA_integer_)

  if (!grepl("\\.gz$", fna, ignore.case = TRUE)) {
    fai_df <- .get_fai_table(fna, cache = cache)
    if (acc %in% rownames(fai_df)) {
      return(as.integer(fai_df[acc, "length"]))
    }
    return(NA_integer_)
  } else {
    gz_obj <- .get_gz_genome(fna, cache = cache)
    if (acc %in% names(gz_obj$lengths)) {
      return(as.integer(gz_obj$lengths[[acc]]))
    }
    return(NA_integer_)
  }
}

extract_genomic_subseq <- function(species,
                                   chr,
                                   start,
                                   end,
                                   strand = "+",
                                   cache = NULL,
                                   ref_dir = get_reference_dir()) {
  fna <- find_genome_file(species, ref_dir = ref_dir)
  if (is.null(fna) || !file.exists(fna)) {
    stop("Genomic FASTA not found for species: ", species)
  }
  acc <- resolve_chromosome(species, chr, cache = cache, ref_dir = ref_dir)
  if (is.na(acc)) {
    stop("Chromosome '", chr, "' not found in species '", species, "'.")
  }
  if (!strand %in% c("+", "-")) {
    stop("strand must be '+' or '-'.")
  }

  s <- as.numeric(start[1])
  e <- as.numeric(end[1])
  if (is.na(s) || is.na(e)) {
    stop("start and end must be numeric coordinates.")
  }

  if (!grepl("\\.gz$", fna, ignore.case = TRUE)) {
    fai_df <- .get_fai_table(fna, cache = cache)
    if (!acc %in% rownames(fai_df)) {
      stop("Chromosome accession '", acc, "' not found in FAI index for ", fna)
    }
    row <- fai_df[acc, , drop = FALSE]
    chr_len <- row$length
    s <- max(1, floor(s))
    e <- min(chr_len, floor(e))
    if (s > e || chr_len == 0) return("")

    offset <- row$offset
    linebases <- as.numeric(row$linebases)
    linewidth <- as.numeric(row$linewidth)

    pos0 <- s - 1
    pos1 <- e - 1
    byte_start <- offset + floor(pos0 / linebases) * linewidth + (pos0 %% linebases)
    byte_end   <- offset + floor(pos1 / linebases) * linewidth + (pos1 %% linebases)
    n_bytes    <- as.integer(byte_end - byte_start + 1)

    con <- file(fna, open = "rb")
    on.exit(close(con), add = TRUE)
    seek(con, where = byte_start, origin = "start")
    raw_bytes <- readBin(con, what = "raw", n = n_bytes)
    keep <- raw_bytes != as.raw(0x0a) & raw_bytes != as.raw(0x0d)
    seq_str <- toupper(rawToChar(raw_bytes[keep]))
  } else {
    gz_obj <- .get_gz_genome(fna, cache = cache)
    if (!acc %in% names(gz_obj$seqs)) {
      stop("Chromosome accession '", acc, "' not found in compressed FASTA ", fna)
    }
    full_seq <- gz_obj$seqs[[acc]]
    chr_len <- nchar(full_seq)
    s <- max(1L, as.integer(floor(s)))
    e <- min(chr_len, as.integer(floor(e)))
    if (s > e || chr_len == 0L) return("")
    seq_str <- toupper(substr(full_seq, s, e))
  }

  if (identical(strand, "-")) {
    seq_str <- revcomp_dna(seq_str)
  }
  seq_str
}
