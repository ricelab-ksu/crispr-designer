# ============================================================================
# Multi-Species Drosophila Gene Lookup & Isoform Feature Model API (Ticket T3)
# ============================================================================
# Provides fast, case-insensitive gene lookup across all 7 Drosophila species
# by gene symbol, FlyBase ID (FBgn), CG number, NCBI Gene ID / LOC ID,
# synonym/alias, protein product name, and genomic coordinates, and extracts
# complete isoform & feature models (transcripts, exons_by_tx, cds_by_tx,
# shared_cds, five_prime_utr, three_prime_utr, promoter_regions,
# first_coding_exons) from genomic.gtf(.gz) in <1 second per query.
# ============================================================================

# Source R/genome_registry.R if present in the project
local({
  candidates <- c(
    "R/genome_registry.R",
    file.path(getwd(), "R", "genome_registry.R"),
    "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design/R/genome_registry.R"
  )
  for (p in unique(candidates)) {
    if (file.exists(p) && !exists(".GENOME_REGISTRY_LOADED", envir = .GlobalEnv)) {
      tryCatch(source(p, local = FALSE), error = function(e) NULL)
      break
    }
  }
})

# Module-level in-memory cache environment
if (!exists(".GENE_LOOKUP_CACHE", envir = .GlobalEnv)) {
  assign(".GENE_LOOKUP_CACHE", new.env(hash = TRUE, parent = emptyenv()), envir = .GlobalEnv)
}

# ---- Fallback helpers (compatible with R/genome_registry.R) ----------------

if (!exists("get_reference_dir", mode = "function")) {
  get_reference_dir <- function() {
    for (ev in c("CRISPR_REFERENCE_DIR", "CAS9_REFERENCE_DIR", "HCR_REFERENCE_DIR")) {
      val <- Sys.getenv(ev, unset = "")
      if (nzchar(val) && dir.exists(val)) return(normalizePath(val, mustWork = FALSE))
    }
    if (exists("get_app_reference_dir", mode = "function")) {
      return(get_app_reference_dir("CRISPR"))
    }
    candidates <- c(
      file.path(getwd(), "reference_files"),
      file.path(dirname(getwd()), "reference_files"),
      file.path(dirname(dirname(getwd())), "reference_files")
    )
    hit <- candidates[vapply(candidates, dir.exists, logical(1))]
    if (length(hit)) normalizePath(hit[1], mustWork = FALSE) else candidates[1]
  }
}

if (!exists("canonical_species", mode = "function")) {
  canonical_species <- function(species) {
    if (is.null(species) || !length(species) || is.na(species[1])) return("melanogaster")
    s <- tolower(trimws(as.character(species[1])))
    s <- sub("^drosophila[_\\s]+", "", s)
    s <- sub("^d[_\\.\\s]+", "", s)
    if (s %in% c("mel", "dmel")) return("melanogaster")
    if (s %in% c("sim", "dsim")) return("simulans")
    if (s %in% c("eug", "deug")) return("eugracilis")
    if (s %in% c("ere", "dere")) return("erecta")
    if (s %in% c("ana", "dana")) return("ananassae")
    if (s %in% c("vir", "dvir")) return("virilis")
    if (s %in% c("tei", "dtei", "tesseri", "tessieri", "teissier", "teissierii")) return("teissieri")
    s
  }
}

species_folder_path <- function(species, ref_dir = get_reference_dir()) {
  sp <- canonical_species(species)
  direct <- file.path(ref_dir, paste0("d_", sp))
  if (dir.exists(direct)) return(direct)
  if (dir.exists(ref_dir)) {
    subs <- list.dirs(ref_dir, recursive = FALSE, full.names = TRUE)
    hit <- subs[grepl(sp, tolower(basename(subs)), fixed = TRUE)]
    if (length(hit)) return(hit[1])
  }
  direct
}

if (!exists("find_gtf_file", mode = "function")) {
  find_gtf_file <- function(species, ref_dir = get_reference_dir()) {
    sd <- species_folder_path(species, ref_dir = ref_dir)
    if (!dir.exists(sd)) return(NULL)
    preferred <- c(
      file.path(sd, "genomic.gtf"),
      file.path(sd, "genomic.gtf.gz")
    )
    for (p in preferred) if (file.exists(p) && !startsWith(basename(p), "._")) return(p)
    hits <- list.files(sd, pattern = "(^|_)genomic\\.(gtf|gff)(\\.gz)?$",
                       ignore.case = TRUE, full.names = TRUE)
    hits <- hits[!startsWith(basename(hits), "._")]
    if (length(hits)) return(hits[1])
    NULL
  }
}

if (!exists("find_rna_file", mode = "function")) {
  find_rna_file <- function(species, ref_dir = get_reference_dir()) {
    sd <- species_folder_path(species, ref_dir = ref_dir)
    if (!dir.exists(sd)) return(NULL)
    preferred <- c(
      file.path(sd, "rna.fna"),
      file.path(sd, "rna.fna.gz")
    )
    for (p in preferred) if (file.exists(p) && !startsWith(basename(p), "._")) return(p)
    hits <- list.files(sd, pattern = "(^|_)rna\\.fna(\\.gz)?$",
                       ignore.case = TRUE, full.names = TRUE)
    hits <- hits[!startsWith(basename(hits), "._")]
    if (length(hits)) return(hits[1])
    NULL
  }
}

# ---- Chromosome Name & Accession Resolution --------------------------------

.SPECIES_CHR_MAP <- list(
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
    "MITOCHONDRION_GENOME" = "NC_024511.2"
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

.is_species_name_arg <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[1])) return(FALSE)
  s <- tolower(trimws(as.character(x[1])))
  s <- sub("^drosophila[_\\.\\s]+", "", s)
  s <- sub("^d[_\\.\\s]+", "", s)
  s %in% c(
    "melanogaster", "simulans", "eugracilis", "erecta",
    "teissieri", "tesseri", "tessieri", "teissier", "teissierii",
    "ananassae", "virilis",
    "dmel", "dsim", "deug", "dere", "dtei", "dana", "dvir",
    "mel", "sim", "eug", "ere", "tei", "ana", "vir"
  )
}

resolve_chromosome <- function(chr, species = "melanogaster", available_chrs = NULL,
                               cache = NULL, ref_dir = get_reference_dir()) {
  # Support both (chr, species) and (species, chr) positional calling conventions
  if (.is_species_name_arg(chr) && !.is_species_name_arg(species)) {
    tmp <- chr
    chr <- species
    species <- tmp
  }
  if (is.null(chr) || !length(chr) || is.na(chr[1])) return(NA_character_)
  raw_chr <- trimws(as.character(chr[1]))
  if (!nzchar(raw_chr)) return(NA_character_)
  sp <- canonical_species(species)

  if (!is.null(available_chrs)) {
    if (raw_chr %in% available_chrs) return(raw_chr)
    idx_ci <- match(tolower(raw_chr), tolower(available_chrs))
    if (!is.na(idx_ci)) return(available_chrs[idx_ci])
  }

  clean_arm <- toupper(sub("^(chr|chromosome[_\\s]*)", "", raw_chr, ignore.case = TRUE))
  cmap <- .SPECIES_CHR_MAP[[sp]]
  if (!is.null(cmap)) {
    if (clean_arm %in% names(cmap)) {
      mapped <- unname(cmap[[clean_arm]])
      if (is.null(available_chrs) || mapped %in% available_chrs) return(mapped)
    }
    no_ver <- sub("\\.[0-9]+$", "", raw_chr)
    cmap_no_ver <- sub("\\.[0-9]+$", "", unname(cmap))
    m_rev <- match(toupper(no_ver), toupper(cmap_no_ver))
    if (!is.na(m_rev)) {
      mapped <- unname(cmap[m_rev])
      if (is.null(available_chrs) || mapped %in% available_chrs) return(mapped)
    }
  }

  if (!is.null(available_chrs)) {
    no_ver <- sub("\\.[0-9]+$", "", raw_chr)
    avail_no_ver <- sub("\\.[0-9]+$", "", available_chrs)
    idx_nv <- match(tolower(no_ver), tolower(avail_no_ver))
    if (!is.na(idx_nv)) return(available_chrs[idx_nv])
  }

  if (exists("get_genome_headers", mode = "function")) {
    headers <- tryCatch(get_genome_headers(sp, cache = cache, ref_dir = ref_dir), error = function(e) NULL)
    if (!is.null(headers) && nrow(headers) > 0) {
      if (raw_chr %in% headers$accession) return(raw_chr)
      idx_alias <- which(toupper(headers$chr_alias) == clean_arm)
      if (length(idx_alias) > 0) return(headers$accession[idx_alias[1]])
    }
  }

  raw_chr
}

chromosome_arm_label <- function(chr, species = "melanogaster") {
  sp <- canonical_species(species)
  cmap <- .SPECIES_CHR_MAP[[sp]]
  if (!is.null(cmap)) {
    m <- match(chr, unname(cmap))
    if (!is.na(m)) return(names(cmap)[m])
  }
  chr
}

# ---- Classic Drosophila Gene Full-Name & Ortholog Product Dictionary -------

.CLASSIC_DMEL_PRODUCTS <- c(
  "y"      = "yellow; protein yellow",
  "w"      = "white; protein white",
  "e"      = "ebony; protein ebony; beta-alanyl-dopamine synthase; N-beta-alanyldopamine synthase",
  "Abd-B"  = "Abdominal B; homeobox protein abdominal-B; abdominal-B",
  "abd-A"  = "abdominal A; homeobox protein abdominal-A; abdominal-A",
  "Ubx"    = "Ultrabithorax; homeobox protein ultrabithorax",
  "Antp"   = "Antennapedia; homeobox protein antennapedia",
  "PlexB"  = "Plexin B; plexin-B; plexin B",
  "PlexA"  = "Plexin A; plexin-A; plexin A",
  "dsx"    = "doublesex; protein doublesex",
  "fru"    = "fruitless; sex determination protein fruitless",
  "tra"    = "transformer; female-specific protein transformer",
  "Sxl"    = "Sex-lethal; protein sex-lethal",
  "f"      = "forked; protein forked",
  "sn"     = "singed; fascin; protein singed",
  "v"      = "vermilion; tryptophan 2,3-dioxygenase",
  "cn"     = "cinnabar; kynurenine 3-monooxygenase",
  "bw"     = "brown; protein brown",
  "st"     = "scarlet; protein scarlet",
  "N"      = "Notch; neurogenic locus notch protein",
  "Dl"     = "Delta; neurogenic locus protein delta",
  "hh"     = "hedgehog; protein hedgehog",
  "wg"     = "wingless; protein wingless",
  "dpp"    = "decapentaplegic; protein decapentaplegic",
  "en"     = "engrailed; segmentation polarity homeobox protein engrailed",
  "ci"     = "cubitus interruptus; transcriptional activator cubitus interruptus",
  "pan"    = "pangolin; protein pangolin",
  "vas"    = "vasa; ATP-dependent RNA helicase vasa",
  "nos"    = "nanos; protein nanos",
  "bcd"    = "bicoid; homeotic protein bicoid"
)

# Parse coordinate query ('chr:start-end', 'chr:start..end', 'chr:start-end:strand')
parse_coordinate_query <- function(query, species = "melanogaster", available_chrs = NULL) {
  if (is.null(query) || !length(query) || is.na(query[1])) return(NULL)
  q <- trimws(as.character(query[1]))
  pat <- "^([A-Za-z0-9_\\.]+):([0-9][0-9,]*)(?:\\.\\.|-|\u2013)([0-9][0-9,]*)(?::([\\+\\-]))?$"
  if (!grepl(pat, q, perl = TRUE)) return(NULL)
  raw_chr <- sub(pat, "\\1", q, perl = TRUE)
  s_str   <- gsub(",", "", sub(pat, "\\2", q, perl = TRUE), fixed = TRUE)
  e_str   <- gsub(",", "", sub(pat, "\\3", q, perl = TRUE), fixed = TRUE)
  strd    <- sub(pat, "\\4", q, perl = TRUE)
  s <- suppressWarnings(as.integer(s_str))
  e <- suppressWarnings(as.integer(e_str))
  if (is.na(s) || is.na(e)) return(NULL)
  if (s > e) { tmp <- s; s <- e; e <- tmp }
  if (!nzchar(strd)) strd <- "+"
  resolved_chr <- resolve_chromosome(raw_chr, species = species, available_chrs = available_chrs)
  list(
    raw_chr = raw_chr,
    chr     = resolved_chr,
    start   = max(1L, s),
    end     = e,
    strand  = strd
  )
}
# ---- Cache Directory & Attribute Extraction Helpers ------------------------

get_gene_index_cache_dir <- function() {
  # The cache lives beside the app payload when there is one, and falls back
  # to a per-user location. It must never point at an external volume: the
  # packaged app is read-only and the volume may not be mounted.
  proj_cache <- file.path(
    if (exists("PJAT_PROJECT_ROOT", inherits = TRUE)) PJAT_PROJECT_ROOT else getwd(),
    ".gene_index_cache")
  if (dir.exists(dirname(proj_cache))) {
    dir.create(proj_cache, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(proj_cache)) return(proj_cache)
  }
  td <- file.path(tempdir(), "cas9_gene_index_cache")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  td
}

extract_gtf_attr_vec <- function(attr_vec, key) {
  n <- length(attr_vec)
  if (!n) return(character(0))
  out <- rep("", n)
  needle <- paste0(key, ' "')
  has_gtf <- grepl(needle, attr_vec, fixed = TRUE)
  if (any(has_gtf)) {
    sub_v <- attr_vec[has_gtf]
    if (key == "gene_id") {
      out[has_gtf] <- sub('^\\s*gene_id "([^"]*)".*$', "\\1", sub_v, perl = TRUE)
    } else {
      pat <- paste0('^(?:.*;\\s*)?', key, ' "([^"]*)".*$')
      val <- sub(pat, "\\1", sub_v, perl = TRUE)
      val[val == sub_v] <- ""
      out[has_gtf] <- val
    }
  }
  rem <- !nzchar(out)
  if (any(rem)) {
    needle_gff <- paste0(key, "=")
    has_gff <- rem & grepl(needle_gff, attr_vec, fixed = TRUE)
    if (any(has_gff)) {
      sub_v2 <- attr_vec[has_gff]
      pat_gff <- paste0("^(?:.*;\\s*)?", key, "=([^;]+).*$")
      val2 <- sub(pat_gff, "\\1", sub_v2, perl = TRUE)
      val2[val2 == sub_v2] <- ""
      out[has_gff] <- val2
    }
  }
  out
}

extract_dbxref_field <- function(attr_vec, prefix) {
  n <- length(attr_vec)
  if (!n) return(character(0))
  out <- rep("", n)
  has <- grepl(prefix, attr_vec, fixed = TRUE)
  if (any(has)) {
    pat <- paste0(".*", prefix, "([A-Za-z0-9_]+).*")
    out[has] <- sub(pat, "\\1", attr_vec[has], perl = TRUE)
  }
  out
}

# ---- D. melanogaster Gene Index Loader -------------------------------------

load_dmel_gene_index <- function(sd) {
  ref_csv <- file.path(sd, "dmel_gene_reference.csv")
  if (!file.exists(ref_csv)) return(NULL)
  dt <- if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(ref_csv, showProgress = FALSE, na.strings = c("", "NA")),
                  stringsAsFactors = FALSE)
  } else {
    utils::read.csv(ref_csv, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  }

  clean_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    trimws(x)
  }

  gene_symbol    <- clean_chr(dt$gene_symbol)
  flybase_id     <- clean_chr(dt$flybase_id)
  cg_number      <- clean_chr(dt$cg_number)
  locus_id       <- clean_chr(dt$ncbi_gene_id)   # e.g. Dmel_CG17245
  ncbi_geneid    <- clean_chr(dt$ncbi_geneid)    # e.g. 43766
  ncbi_combined  <- ifelse(nzchar(ncbi_geneid), ncbi_geneid, locus_id)
  chromosome     <- clean_chr(dt$chromosome)
  start_pos      <- as.integer(dt$start)
  end_pos        <- as.integer(dt$end)
  strand         <- clean_chr(dt$strand)
  biotype        <- clean_chr(dt$biotype)
  synonyms       <- clean_chr(dt$synonyms)
  transcript_ids <- clean_chr(dt$transcript_ids)
  product        <- rep("", nrow(dt))

  fbgn_csv <- file.path(sd, "fbgn_map.csv")
  if (file.exists(fbgn_csv)) {
    fb_df <- tryCatch(utils::read.csv(fbgn_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(fb_df) && all(c("FBgn", "symbol") %in% names(fb_df))) {
      m_fb <- match(flybase_id, fb_df$FBgn)
      miss_sym <- !nzchar(gene_symbol) & !is.na(m_fb)
      if (any(miss_sym)) gene_symbol[miss_sym] <- fb_df$symbol[m_fb[miss_sym]]
    }
  }

  ish_csv <- file.path(sd, "ISH_markers_pitt.csv")
  if (file.exists(ish_csv)) {
    ish_df <- tryCatch(utils::read.csv(ish_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(ish_df) && all(c("symbol", "name") %in% names(ish_df))) {
      m_ish <- match(gene_symbol, ish_df$symbol)
      has_ish <- !is.na(m_ish) & nzchar(ish_df$name[m_ish])
      product[has_ish] <- ish_df$name[m_ish[has_ish]]
    }
  }

  for (sym in names(.CLASSIC_DMEL_PRODUCTS)) {
    hit <- which(gene_symbol == sym)
    if (length(hit)) {
      full_prod <- .CLASSIC_DMEL_PRODUCTS[[sym]]
      primary_prod <- trimws(strsplit(full_prod, ";", fixed = TRUE)[[1]][1])
      for (h in hit) {
        if (!nzchar(product[h])) product[h] <- primary_prod
        existing_syns <- trimws(strsplit(synonyms[h], ";", fixed = TRUE)[[1]])
        new_syns <- trimws(strsplit(full_prod, ";", fixed = TRUE)[[1]])
        add_syns <- new_syns[!(tolower(new_syns) %in% tolower(existing_syns))]
        if (length(add_syns)) {
          synonyms[h] <- if (nzchar(synonyms[h])) {
            paste(c(synonyms[h], add_syns), collapse = "; ")
          } else {
            paste(add_syns, collapse = "; ")
          }
        }
      }
    }
  }

  has_locus <- nzchar(locus_id) & !grepl("Dmel_", synonyms, fixed = TRUE)
  synonyms[has_locus] <- ifelse(
    nzchar(synonyms[has_locus]),
    paste0(synonyms[has_locus], "; ", locus_id[has_locus]),
    locus_id[has_locus]
  )

  data.frame(
    gene_symbol    = gene_symbol,
    flybase_id     = flybase_id,
    cg_number      = cg_number,
    ncbi_gene_id   = ncbi_combined,
    chromosome     = chromosome,
    start          = start_pos,
    end            = end_pos,
    strand         = strand,
    biotype        = biotype,
    synonyms       = synonyms,
    product        = product,
    transcript_ids = transcript_ids,
    locus_id       = ifelse(nzchar(locus_id), locus_id, gene_symbol),
    stringsAsFactors = FALSE
  )
}
# ---- Non-Melanogaster GTF & RNA Fallback Gene Index Builders ----------------

clean_product_name <- function(prod) {
  if (!length(prod)) return(character(0))
  p <- sub(",\\s*transcript variant.*$", "", prod, ignore.case = TRUE, perl = TRUE)
  p <- sub(",?\\s*isoform\\s+[A-Za-z0-9]+$", "", p, ignore.case = TRUE, perl = TRUE)
  trimws(p)
}

build_gtf_gene_index <- function(gtf_path, species = "simulans", cache = NULL) {
  if (!file.exists(gtf_path)) return(NULL)
  is_gz <- grepl("\\.gz$", gtf_path, ignore.case = TRUE)

  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- tryCatch({
      if (is_gz) {
        cmd <- sprintf("gzip -dc %s | grep -v '^#'", shQuote(gtf_path))
        suppressWarnings(data.table::fread(
          cmd = cmd, sep = "\t", header = FALSE, quote = "",
          select = c(1L, 3L, 4L, 5L, 7L, 8L, 9L), showProgress = FALSE
        ))
      } else {
        suppressWarnings(data.table::fread(
          gtf_path, sep = "\t", header = FALSE, quote = "", skip = 4L,
          select = c(1L, 3L, 4L, 5L, 7L, 8L, 9L), showProgress = FALSE
        ))
      }
    }, error = function(e) NULL)
  } else {
    dt <- NULL
  }

  if (is.null(dt) || !nrow(dt)) {
    con <- if (is_gz) gzfile(gtf_path, "rt") else file(gtf_path, "rt")
    on.exit(close(con), add = TRUE)
    lns <- readLines(con, warn = FALSE)
    lns <- lns[!grepl("^#", lns) & nzchar(lns)]
    sp_list <- strsplit(lns, "\t", fixed = TRUE)
    sp_list <- Filter(function(x) length(x) >= 9, sp_list)
    dt <- data.frame(
      V1 = vapply(sp_list, `[`, "", 1),
      V3 = vapply(sp_list, `[`, "", 3),
      V4 = as.integer(vapply(sp_list, `[`, "", 4)),
      V5 = as.integer(vapply(sp_list, `[`, "", 5)),
      V7 = vapply(sp_list, `[`, "", 7),
      V8 = vapply(sp_list, `[`, "", 8),
      V9 = vapply(sp_list, `[`, "", 9),
      stringsAsFactors = FALSE
    )
    if (requireNamespace("data.table", quietly = TRUE)) {
      dt <- data.table::as.data.table(dt)
    }
  }

  names(dt) <- c("chr", "feature", "start", "end", "strand", "phase", "attr")

  keep_feat <- dt$feature %in% c("gene", "transcript", "mRNA", "exon", "CDS", "start_codon", "stop_codon")
  feat_dt <- dt[keep_feat, ]

  if (!is.null(cache)) {
    cache[[paste0(".gtf_dt:", normalizePath(gtf_path, mustWork = FALSE))]] <- feat_dt
  }
  gene_rows <- feat_dt[feat_dt$feature == "gene", ]
  if (!nrow(gene_rows)) {
    ex_rows <- feat_dt[feat_dt$feature == "exon", ]
    ex_rows$gene_id <- extract_gtf_attr_vec(ex_rows$attr, "gene_id")
    gene_rows <- do.call(rbind, lapply(split(ex_rows, ex_rows$gene_id), function(d) {
      d0 <- d[1, , drop = FALSE]
      d0$start <- min(d$start, na.rm = TRUE)
      d0$end   <- max(d$end, na.rm = TRUE)
      d0
    }))
  } else {
    gene_rows$gene_id <- extract_gtf_attr_vec(gene_rows$attr, "gene_id")
  }

  tx_rows <- feat_dt[feat_dt$feature %in% c("transcript", "mRNA"), ]
  if (!nrow(tx_rows)) {
    tx_rows <- feat_dt[feat_dt$feature == "exon" & grepl('exon_number "1"', feat_dt$attr, fixed = TRUE), ]
  }

  tx_gid  <- extract_gtf_attr_vec(tx_rows$attr, "gene_id")
  tx_id   <- extract_gtf_attr_vec(tx_rows$attr, "transcript_id")
  tx_prod <- clean_product_name(extract_gtf_attr_vec(tx_rows$attr, "product"))

  if (requireNamespace("data.table", quietly = TRUE)) {
    tmp_tx <- data.table::data.table(gid = tx_gid, tid = tx_id, prod = tx_prod)
    agg_tx <- tmp_tx[, .(
      product = {
        u <- unique(prod[nzchar(prod)])
        nu <- u[!grepl("^uncharacterized", u, ignore.case = TRUE)]
        if (length(nu)) paste(nu, collapse = "; ") else if (length(u)) paste(u, collapse = "; ") else ""
      },
      transcript_ids = {
        u <- unique(tid[nzchar(tid)])
        if (length(u)) paste(u, collapse = "; ") else ""
      }
    ), by = gid]
    prod_by_gene <- stats::setNames(agg_tx$product, agg_tx$gid)
    txid_by_gene <- stats::setNames(agg_tx$transcript_ids, agg_tx$gid)
  } else {
    prod_by_gene <- tapply(tx_prod, tx_gid, function(v) paste(unique(v[nzchar(v)]), collapse = "; "))
    txid_by_gene <- tapply(tx_id, tx_gid, function(v) paste(unique(v[nzchar(v)]), collapse = "; "))
  }

  gid        <- gene_rows$gene_id
  gsym       <- extract_gtf_attr_vec(gene_rows$attr, "gene")
  gsym[!nzchar(gsym)] <- gid[!nzchar(gsym)]
  locus_tag  <- extract_gtf_attr_vec(gene_rows$attr, "locus_tag")
  biotype    <- extract_gtf_attr_vec(gene_rows$attr, "gene_biotype")
  fb_id      <- extract_dbxref_field(gene_rows$attr, "FLYBASE:")
  ncbi_id    <- extract_dbxref_field(gene_rows$attr, "GeneID:")
  desc       <- extract_gtf_attr_vec(gene_rows$attr, "description")
  syn        <- extract_gtf_attr_vec(gene_rows$attr, "gene_synonym")

  prod_mapped <- unname(prod_by_gene[gid])
  prod_mapped[is.na(prod_mapped)] <- ""
  prod_mapped[!nzchar(prod_mapped) & nzchar(desc)] <- desc[!nzchar(prod_mapped) & nzchar(desc)]

  tx_mapped <- unname(txid_by_gene[gid])
  tx_mapped[is.na(tx_mapped)] <- ""

  cg_num <- ifelse(grepl("^CG[0-9]+$", gsym), gsym,
            ifelse(grepl("CG[0-9]+", prod_mapped),
                   sub(".*\\b(CG[0-9]+)\\b.*", "\\1", prod_mapped, perl = TRUE), ""))

  pl <- tolower(prod_mapped)
  extra_alias <- rep("", length(gid))
  extra_alias[grepl("^protein yellow(;|$)", pl)] <- "yellow; y"
  extra_alias[grepl("^protein white(;|$)", pl)]  <- "white; w"
  extra_alias[grepl("^plexin-b(;|$)|^plexin b(;|$)", pl)] <- "PlexB; plexin B; plex B"
  extra_alias[grepl("homeobox protein abdominal-b", pl)]  <- "Abd-B; Abdominal B; abdominal-B"
  extra_alias[grepl("homeobox protein abdominal-a", pl)]  <- "abd-A; abdominal A; abdominal-A"
  extra_alias[grepl("^protein doublesex(;|$)", pl)]       <- "dsx; doublesex"

  syn_full <- prod_mapped
  has_lt <- nzchar(locus_tag) & (locus_tag != gsym)
  syn_full[has_lt] <- ifelse(nzchar(syn_full[has_lt]), paste0(syn_full[has_lt], "; ", locus_tag[has_lt]), locus_tag[has_lt])
  has_s <- nzchar(syn)
  syn_full[has_s] <- ifelse(nzchar(syn_full[has_s]), paste0(syn_full[has_s], "; ", syn[has_s]), syn[has_s])
  has_ea <- nzchar(extra_alias)
  syn_full[has_ea] <- ifelse(nzchar(syn_full[has_ea]), paste0(syn_full[has_ea], "; ", extra_alias[has_ea]), extra_alias[has_ea])

  data.frame(
    gene_symbol    = as.character(gsym),
    flybase_id     = as.character(fb_id),
    cg_number      = as.character(cg_num),
    ncbi_gene_id   = ifelse(nzchar(ncbi_id), as.character(ncbi_id), as.character(gid)),
    chromosome     = as.character(gene_rows$chr),
    start          = as.integer(gene_rows$start),
    end            = as.integer(gene_rows$end),
    strand         = as.character(gene_rows$strand),
    biotype        = as.character(biotype),
    synonyms       = as.character(syn_full),
    product        = as.character(prod_mapped),
    transcript_ids = as.character(tx_mapped),
    locus_id       = as.character(gid),
    stringsAsFactors = FALSE
  )
}
# Fallback index builder from rna.fna(.gz) headers (for RNA-only species like
# d_ananassae and d_virilis before genomic.gtf is downloaded)
build_rna_fallback_index <- function(rna_path, species = "ananassae") {
  if (is.null(rna_path) || !file.exists(rna_path)) return(NULL)
  is_gz <- grepl("\\.gz$", rna_path, ignore.case = TRUE)
  hdrs <- tryCatch({
    cmd <- if (is_gz) {
      sprintf("gzip -dc %s | LC_ALL=C grep '^>'", shQuote(rna_path))
    } else {
      sprintf("LC_ALL=C grep '^>' %s", shQuote(rna_path))
    }
    system(cmd, intern = TRUE)
  }, error = function(e) character(0))

  if (!length(hdrs)) {
    con <- if (is_gz) gzfile(rna_path, "rt") else file(rna_path, "rt")
    on.exit(close(con), add = TRUE)
    lns <- readLines(con, warn = FALSE)
    hdrs <- lns[startsWith(lns, ">")]
  }
  if (!length(hdrs)) return(NULL)

  tx_id <- sub("^>(\\S+).*$", "\\1", hdrs, perl = TRUE)
  has_sym <- grepl("\\([A-Za-z0-9_\\-]+\\)", hdrs, perl = TRUE)
  sym <- rep("", length(hdrs))
  sym[has_sym] <- sub("^.*\\(([A-Za-z0-9_\\-]+)\\)[^\\(]*$", "\\1", hdrs[has_sym], perl = TRUE)
  sym[!nzchar(sym)] <- tx_id[!nzchar(sym)]

  prod <- sub("^>\\S+\\s+(?:PREDICTED:\\s+)?(?:Drosophila\\s+\\S+\\s+)?", "", hdrs, perl = TRUE)
  prod <- sub("\\s*\\([A-Za-z0-9_\\-]+\\).*$", "", prod, perl = TRUE)
  prod <- clean_product_name(prod)

  grp <- split(seq_along(hdrs), sym)
  g_syms <- names(grp)
  g_prod <- vapply(grp, function(idx) {
    u <- unique(prod[idx][nzchar(prod[idx])])
    paste(u, collapse = "; ")
  }, character(1))
  g_tx <- vapply(grp, function(idx) paste(unique(tx_id[idx]), collapse = "; "), character(1))
  ncbi_num <- sub("^LOC([0-9]+)$", "\\1", g_syms)

  data.frame(
    gene_symbol    = g_syms,
    flybase_id     = rep("", length(g_syms)),
    cg_number      = ifelse(grepl("CG[0-9]+", g_prod),
                            sub(".*\\b(CG[0-9]+)\\b.*", "\\1", g_prod, perl = TRUE), ""),
    ncbi_gene_id   = ncbi_num,
    chromosome     = rep("rna_transcriptome", length(g_syms)),
    start          = rep(1L, length(g_syms)),
    end            = rep(1000L, length(g_syms)),
    strand         = rep("+", length(g_syms)),
    biotype        = ifelse(grepl("^X[MR]_|^N[MR]_", g_tx), "protein_coding", "transcript"),
    synonyms       = g_prod,
    product        = g_prod,
    transcript_ids = g_tx,
    locus_id       = g_syms,
    stringsAsFactors = FALSE
  )
}

#' Load (and cache) the standardized gene lookup index for a Drosophila species
#'
#' @param species Species name (e.g., "melanogaster", "d_simulans", "eugracilis")
#' @param cache Optional environment for in-memory caching (defaults to .GENE_LOOKUP_CACHE)
#' @param ref_dir Base directory containing d_<species> reference genome folders
#' @return A standardized data.frame of genes for the requested species
load_gene_index <- function(species = "melanogaster", cache = NULL, ref_dir = get_reference_dir()) {
  sp <- canonical_species(species)
  cache_env <- if (is.environment(cache)) cache else .GENE_LOOKUP_CACHE
  mem_key <- paste0("gene_index:", sp, ":", normalizePath(ref_dir, mustWork = FALSE))

  if (exists(mem_key, envir = cache_env, inherits = FALSE)) {
    return(get(mem_key, envir = cache_env, inherits = FALSE))
  }

  sd <- species_folder_path(sp, ref_dir = ref_dir)
  gtf_path <- find_gtf_file(sp, ref_dir = ref_dir)
  dmel_csv <- file.path(sd, "dmel_gene_reference.csv")
  rna_path <- find_rna_file(sp, ref_dir = ref_dir)

  src_file <- if (sp == "melanogaster" && file.exists(dmel_csv)) {
    dmel_csv
  } else if (!is.null(gtf_path) && file.exists(gtf_path)) {
    gtf_path
  } else if (!is.null(rna_path) && file.exists(rna_path)) {
    rna_path
  } else {
    stop(sprintf("No gene reference CSV, genomic.gtf(.gz), or rna.fna(.gz) found for species '%s' in %s",
                 species, sd))
  }

  finfo <- file.info(src_file)
  disk_cache_dir <- get_gene_index_cache_dir()
  disk_rds <- file.path(
    disk_cache_dir,
    sprintf("%s_gene_index_%s.rds", sp, as.character(round(finfo$size)))
  )

  if (file.exists(disk_rds)) {
    cached_df <- tryCatch(readRDS(disk_rds), error = function(e) NULL)
    if (is.data.frame(cached_df) && nrow(cached_df) > 0) {
      assign(mem_key, cached_df, envir = cache_env)
      return(cached_df)
    }
  }

  idx_df <- if (sp == "melanogaster" && file.exists(dmel_csv)) {
    load_dmel_gene_index(sd)
  } else if (!is.null(gtf_path) && file.exists(gtf_path)) {
    build_gtf_gene_index(gtf_path, species = sp, cache = cache_env)
  } else {
    build_rna_fallback_index(rna_path, species = sp)
  }

  if (is.null(idx_df) || !nrow(idx_df)) {
    stop(sprintf("Failed to build gene index for species '%s' from %s", species, src_file))
  }

  tryCatch(saveRDS(idx_df, disk_rds), error = function(e) NULL)
  assign(mem_key, idx_df, envir = cache_env)
  idx_df
}
# ---- Standardized Gene Search API ------------------------------------------

.STANDARD_SEARCH_COLS <- c(
  "gene_symbol", "flybase_id", "cg_number", "ncbi_gene_id",
  "chromosome", "start", "end", "strand", "biotype",
  "synonyms", "product", "match_type"
)

empty_search_df <- function() {
  data.frame(
    gene_symbol  = character(0),
    flybase_id   = character(0),
    cg_number    = character(0),
    ncbi_gene_id = character(0),
    chromosome   = character(0),
    start        = integer(0),
    end          = integer(0),
    strand       = character(0),
    biotype      = character(0),
    synonyms     = character(0),
    product      = character(0),
    match_type   = character(0),
    stringsAsFactors = FALSE
  )
}

#' Search genes across any of the 7 Drosophila species
#'
#' @param query Search query (symbol, FBgn, CG, NCBI/LOC ID, synonym, product, or chr:start-end)
#' @param species Species name (default "melanogaster")
#' @param max_results Maximum number of rows to return (default 25)
#' @param cache Optional cache environment
#' @param ref_dir Reference genomes base directory
#' @return A standardized data.frame of matching genes with match_type column
search_genes <- function(query, species = "melanogaster", max_results = 25,
                         cache = NULL, ref_dir = get_reference_dir()) {
  if (is.null(query) || !length(query) || is.na(query[1])) return(empty_search_df())
  q <- trimws(as.character(query[1]))
  if (!nzchar(q)) return(empty_search_df())

  sp <- canonical_species(species)
  idx <- load_gene_index(sp, cache = cache, ref_dir = ref_dir)
  avail_chrs <- unique(idx$chromosome)

  coord <- parse_coordinate_query(q, species = sp, available_chrs = avail_chrs)
  if (!is.null(coord)) {
    ov <- which(idx$chromosome == coord$chr & idx$start <= coord$end & idx$end >= coord$start)
    if (length(ov)) {
      has_explicit_strand <- grepl(":[\\+\\-]$", q)
      strand_match <- if (has_explicit_strand) ifelse(idx$strand[ov] == coord$strand, 1L, 0L) else rep(1L, length(ov))
      ov_len <- pmin(idx$end[ov], coord$end) - pmax(idx$start[ov], coord$start) + 1L
      gene_len <- pmax(1L, idx$end[ov] - idx$start[ov] + 1L)
      union_len <- pmax(idx$end[ov], coord$end) - pmin(idx$start[ov], coord$start) + 1L
      jaccard <- ov_len / union_len
      gene_cov <- ov_len / gene_len
      is_coding <- ifelse(idx$biotype[ov] == "protein_coding", 1L, 0L)
      ov <- ov[order(-strand_match, -round(gene_cov, 2), -is_coding, -jaccard, -ov_len, idx$start[ov])]
      res <- idx[ov, intersect(names(idx), .STANDARD_SEARCH_COLS), drop = FALSE]
      res$match_type <- "coordinates"
      res <- res[, .STANDARD_SEARCH_COLS, drop = FALSE]
      rownames(res) <- NULL
      return(utils::head(res, max_results))
    } else {
      syn_row <- data.frame(
        gene_symbol  = sprintf("%s:%d-%d", coord$chr, coord$start, coord$end),
        flybase_id   = "",
        cg_number    = "",
        ncbi_gene_id = "",
        chromosome   = coord$chr,
        start        = coord$start,
        end          = coord$end,
        strand       = coord$strand,
        biotype      = "genomic_region",
        synonyms     = q,
        product      = sprintf("Genomic region %s:%d-%d (%s)", coord$chr, coord$start, coord$end, coord$strand),
        match_type   = "coordinates",
        stringsAsFactors = FALSE
      )
      return(syn_row)
    }
  }
  ql <- tolower(q)
  q_loc_stripped <- sub("^loc", "", ql)
  sym_l   <- tolower(idx$gene_symbol)
  fb_l    <- tolower(idx$flybase_id)
  cg_l    <- tolower(idx$cg_number)
  ncbi_l  <- tolower(idx$ncbi_gene_id)
  locus_l <- if ("locus_id" %in% names(idx)) tolower(idx$locus_id) else rep("", nrow(idx))
  prod_l  <- tolower(idx$product)
  syn_l   <- tolower(idx$synonyms)

  has_exact_token <- function(vec_l, target_l) {
    if (!nzchar(target_l)) return(rep(FALSE, length(vec_l)))
    esc <- gsub("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", target_l)
    pat <- paste0("(^|;\\s*)", esc, "(\\s*;|$)")
    grepl(pat, vec_l, perl = TRUE)
  }

  tier <- rep(Inf, nrow(idx))
  mtype <- rep("", nrow(idx))

  is_sym_exact <- idx$gene_symbol == q
  is_sym_ci    <- !is_sym_exact & (sym_l == ql)
  tier[is_sym_exact] <- 1.0; mtype[is_sym_exact] <- "gene_symbol"
  tier[is_sym_ci]    <- 1.1; mtype[is_sym_ci]    <- "gene_symbol"

  is_fb <- (tier == Inf) & nzchar(fb_l) & (fb_l == ql)
  tier[is_fb] <- 2.0; mtype[is_fb] <- "flybase_id"

  is_cg <- (tier == Inf) & nzchar(cg_l) & (cg_l == ql)
  tier[is_cg] <- 3.0; mtype[is_cg] <- "cg_number"

  is_ncbi <- (tier == Inf) & (
    (nzchar(ncbi_l) & (ncbi_l == ql | paste0("loc", ncbi_l) == ql | ncbi_l == q_loc_stripped)) |
    (nzchar(locus_l) & (locus_l == ql))
  )
  tier[is_ncbi] <- 4.0; mtype[is_ncbi] <- "ncbi_gene_id"

  is_syn_tok <- (tier == Inf) & has_exact_token(syn_l, ql)
  tier[is_syn_tok] <- 5.0; mtype[is_syn_tok] <- "synonym"

  is_prod_tok <- (tier == Inf) & (prod_l == ql | has_exact_token(prod_l, ql))
  tier[is_prod_tok] <- 6.0; mtype[is_prod_tok] <- "product"

  is_sym_pref <- (tier == Inf) & startsWith(sym_l, ql)
  tier[is_sym_pref] <- 7.0; mtype[is_sym_pref] <- "gene_symbol_prefix"

  escaped_q <- gsub("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", ql)
  wb_pat <- paste0("(^|[^a-z0-9])", escaped_q, "([^a-z0-9]|$)")
  is_prod_wb <- (tier == Inf) & nzchar(prod_l) & grepl(wb_pat, prod_l, perl = TRUE)
  tier[is_prod_wb] <- 8.0; mtype[is_prod_wb] <- "product"

  is_syn_wb <- (tier == Inf) & nzchar(syn_l) & grepl(wb_pat, syn_l, perl = TRUE)
  tier[is_syn_wb] <- 8.5; mtype[is_syn_wb] <- "synonym"

  if (nchar(ql) >= 2L) {
    is_sym_sub <- (tier == Inf) & grepl(ql, sym_l, fixed = TRUE)
    tier[is_sym_sub] <- 9.0; mtype[is_sym_sub] <- "gene_symbol_partial"

    is_prod_sub <- (tier == Inf) & grepl(ql, prod_l, fixed = TRUE)
    tier[is_prod_sub] <- 9.5; mtype[is_prod_sub] <- "product"

    is_syn_sub <- (tier == Inf) & grepl(ql, syn_l, fixed = TRUE)
    tier[is_syn_sub] <- 10.0; mtype[is_syn_sub] <- "synonym"
  }

  keep <- which(is.finite(tier))
  if (!length(keep)) return(empty_search_df())

  arm_labels <- vapply(idx$chromosome[keep], chromosome_arm_label, character(1), species = sp)
  is_x_yellow <- ifelse(ql %in% c("yellow", "protein yellow", "y") & arm_labels == "X", 0L, 1L)
  is_coding <- ifelse(idx$biotype[keep] == "protein_coding", 0L, 1L)
  n_tx_chars <- if ("transcript_ids" %in% names(idx)) nchar(idx$transcript_ids[keep]) else rep(0L, length(keep))

  ord <- order(
    tier[keep],
    is_x_yellow,
    is_coding,
    nchar(idx$product[keep]) == 0L,
    nchar(idx$gene_symbol[keep]),
    -n_tx_chars,
    idx$gene_symbol[keep]
  )
  sel <- keep[ord]
  if (length(sel) > max_results) sel <- sel[seq_len(max_results)]

  out <- idx[sel, intersect(names(idx), .STANDARD_SEARCH_COLS), drop = FALSE]
  out$match_type <- mtype[sel]
  out <- out[, .STANDARD_SEARCH_COLS, drop = FALSE]
  rownames(out) <- NULL
  out
}
# ---- Fast Sub-Second GTF Locus Extractor (Binary Search + Cache) -----------

.PRECOMPUTED_GTF_CHR_BOUNDS <- list(
  melanogaster_185199259 = data.frame(
    chr = c(
      "NC_004354.4", "NW_001845284.1", "NT_033779.5", "NT_033778.4",
      "NT_037436.4", "NT_033777.3", "NC_004353.4", "NC_024512.1",
      "NW_001845431.1", "NW_001845819.1", "NW_001846712.1", "NW_001846812.1",
      "NW_001847227.1", "NW_001846198.1", "NW_001845031.1", "NW_001845220.1",
      "NW_001845164.1", "NW_007931104.1", "NW_001844935.1", "NW_001846187.1",
      "NW_001845870.1", "NW_001845199.1", "NW_007931083.1", "NW_007931121.1",
      "NC_024511.2"
    ),
    start_byte = c(
      151, 33311782, 33315613, 66203872, 104664543, 139017189, 182168063,
      184905235, 185086306, 185087401, 185087887, 185090621, 185091710,
      185093032, 185094142, 185095258, 185096353, 185098534, 185098831,
      185103311, 185103602, 185104918, 185105212, 185120463, 185151290
    ),
    end_byte = c(
      33311782, 33315613, 66203872, 104664543, 139017189, 182168063,
      184905235, 185086306, 185087401, 185087887, 185090621, 185091710,
      185093032, 185094142, 185095258, 185096353, 185098534, 185098831,
      185103311, 185103602, 185104918, 185105212, 185120463, 185151290,
      185199259
    ),
    stringsAsFactors = FALSE
  ),
  simulans_167574976 = data.frame(
    chr = c(
      "NC_052520.2", "NC_052521.2", "NC_052522.2",
      "NC_052523.2", "NC_052524.2", "NC_052525.2"
    ),
    start_byte = c(169, 27315573, 64173220, 94539806, 133390859, 137116157),
    end_byte   = c(27315573, 64173220, 94539806, 133390859, 137116157, 166048284),
    stringsAsFactors = FALSE
  )
)

parse_gtf_lines_to_df <- function(lines) {
  if (!length(lines)) {
    return(data.frame(
      chr = character(0), feature = character(0), start = integer(0),
      end = integer(0), strand = character(0), phase = character(0),
      attr = character(0), gene_id = character(0), stringsAsFactors = FALSE
    ))
  }
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  sp <- strsplit(lines, "\t", fixed = TRUE)
  sp <- Filter(function(x) length(x) >= 9, sp)
  if (!length(sp)) {
    return(data.frame(
      chr = character(0), feature = character(0), start = integer(0),
      end = integer(0), strand = character(0), phase = character(0),
      attr = character(0), gene_id = character(0), stringsAsFactors = FALSE
    ))
  }
  attr_v <- vapply(sp, `[`, "", 9)
  data.frame(
    chr     = vapply(sp, `[`, "", 1),
    feature = vapply(sp, `[`, "", 3),
    start   = suppressWarnings(as.integer(vapply(sp, `[`, "", 4))),
    end     = suppressWarnings(as.integer(vapply(sp, `[`, "", 5))),
    strand  = vapply(sp, `[`, "", 7),
    phase   = vapply(sp, `[`, "", 8),
    attr    = attr_v,
    gene_id = extract_gtf_attr_vec(attr_v, "gene_id"),
    stringsAsFactors = FALSE
  )
}

extract_gtf_locus_rows <- function(gtf_path, species, chr, start_pos, end_pos,
                                   target_gene_id = NULL, cache = NULL) {
  if (is.null(gtf_path) || !file.exists(gtf_path)) return(parse_gtf_lines_to_df(character(0)))
  cache_env <- if (is.environment(cache)) cache else .GENE_LOOKUP_CACHE
  norm_gtf <- normalizePath(gtf_path, mustWork = FALSE)
  dt_key <- paste0(".gtf_dt:", norm_gtf)

  if (exists(dt_key, envir = cache_env, inherits = FALSE)) {
    feat_dt <- get(dt_key, envir = cache_env, inherits = FALSE)
    sub_dt <- as.data.frame(feat_dt[feat_dt$chr == chr & feat_dt$start <= end_pos & feat_dt$end >= start_pos, ],
                            stringsAsFactors = FALSE)
    if (nrow(sub_dt)) {
      sub_dt$gene_id <- extract_gtf_attr_vec(sub_dt$attr, "gene_id")
      if (!is.null(target_gene_id) && nzchar(target_gene_id)) {
        by_id <- sub_dt[sub_dt$gene_id == target_gene_id, , drop = FALSE]
        if (nrow(by_id)) return(by_id)
      }
      return(sub_dt)
    }
  }

  sp <- canonical_species(species)
  is_gz <- grepl("\\.gz$", gtf_path, ignore.case = TRUE)

  if (!is_gz) {
    fsize <- file.info(gtf_path)$size
    bkey <- paste0(sp, "_", as.character(round(fsize)))
    bounds_df <- .PRECOMPUTED_GTF_CHR_BOUNDS[[bkey]]
    has_pre <- !is.null(bounds_df) && (chr %in% bounds_df$chr)
    lo <- if (has_pre) as.numeric(bounds_df$start_byte[bounds_df$chr == chr][1]) else 0
    hi <- if (has_pre) as.numeric(bounds_df$end_byte[bounds_df$chr == chr][1])   else as.numeric(fsize)
    max_hi <- hi
    search_start <- max(1L, as.integer(start_pos) - 75000L)

    con <- file(gtf_path, "rb")
    on.exit(close(con), add = TRUE)

    can_bsearch <- has_pre || grepl("^NW_|^NC_|^NT_", chr)
    if (can_bsearch) {
      while ((hi - lo) > 65536) {
        mid <- floor((lo + hi) / 2)
        seek(con, mid)
        raw_buf <- readBin(con, "raw", n = 2048L)
        if (!length(raw_buf)) break
        txt <- rawToChar(raw_buf)
        lns <- strsplit(txt, "\n", fixed = TRUE)[[1]]
        if (mid > 0 && length(lns) > 1) lns <- lns[-1]
        lns <- lns[nzchar(lns) & substr(lns, 1, 1) != "#"]
        if (!length(lns)) { lo <- mid; next }
        parts <- strsplit(lns[1], "\t", fixed = TRUE)[[1]]
        c_mid <- parts[1]
        pos <- suppressWarnings(as.integer(parts[4]))
        if (has_pre) {
          if (!is.na(pos) && pos < search_start) lo <- mid else hi <- mid
        } else {
          if (c_mid < chr || (c_mid == chr && !is.na(pos) && pos < search_start)) lo <- mid else hi <- mid
        }
      }

      span_bp <- max(100000, as.numeric(end_pos - start_pos) + 100000)
      read_bytes <- min(max_hi - lo, max(393216, span_bp * 12))
      seek(con, lo)
      chunk_raw <- readBin(con, "raw", n = as.integer(read_bytes))
      chunk_txt <- rawToChar(chunk_raw)
      lns <- strsplit(chunk_txt, "\n", fixed = TRUE)[[1]]
      if (lo > 0 && length(lns) > 1) lns <- lns[-1]
      if (length(lns) > 1) lns <- lns[-length(lns)]

      if (!is.null(target_gene_id) && nzchar(target_gene_id)) {
        matched <- lns[grepl(paste0('"', target_gene_id, '"'), lns, fixed = TRUE)]
        if (length(matched)) return(parse_gtf_lines_to_df(matched))
      }
      df_loc <- parse_gtf_lines_to_df(lns)
      df_loc <- df_loc[df_loc$chr == chr & df_loc$start <= end_pos & df_loc$end >= start_pos, , drop = FALSE]
      if (nrow(df_loc)) return(df_loc)
    }
  }

  build_gtf_gene_index(gtf_path, species = sp, cache = cache_env)
  if (exists(dt_key, envir = cache_env, inherits = FALSE)) {
    feat_dt <- get(dt_key, envir = cache_env, inherits = FALSE)
    sub_dt <- as.data.frame(feat_dt[feat_dt$chr == chr & feat_dt$start <= end_pos & feat_dt$end >= start_pos, ],
                            stringsAsFactors = FALSE)
    if (nrow(sub_dt)) {
      sub_dt$gene_id <- extract_gtf_attr_vec(sub_dt$attr, "gene_id")
      if (!is.null(target_gene_id) && nzchar(target_gene_id)) {
        by_id <- sub_dt[sub_dt$gene_id == target_gene_id, , drop = FALSE]
        if (nrow(by_id)) return(by_id)
      }
      return(sub_dt)
    }
  }

  parse_gtf_lines_to_df(character(0))
}
# ---- Interval Arithmetic & Isoform Feature Geometry ------------------------

empty_interval_df <- function(with_tx = FALSE, with_phase = FALSE) {
  df <- data.frame(
    chr    = character(0),
    start  = integer(0),
    end    = integer(0),
    strand = character(0),
    stringsAsFactors = FALSE
  )
  if (with_phase) df$phase <- integer(0)
  if (with_tx)    df$transcript_id <- character(0)
  df
}

merge_intervals_df <- function(df, max_gap = 0L) {
  if (is.null(df) || !nrow(df)) return(empty_interval_df())
  df <- df[!is.na(df$start) & !is.na(df$end) & df$start <= df$end, , drop = FALSE]
  if (!nrow(df)) return(empty_interval_df())
  df <- df[order(df$chr, df$start, df$end), , drop = FALSE]
  if (nrow(df) == 1L) {
    rownames(df) <- NULL
    return(df[, intersect(c("chr", "start", "end", "strand"), names(df)), drop = FALSE])
  }
  out_chr <- character(0)
  out_s   <- integer(0)
  out_e   <- integer(0)
  out_st  <- character(0)

  cur_chr <- df$chr[1]
  cur_s   <- as.integer(df$start[1])
  cur_e   <- as.integer(df$end[1])
  cur_st  <- df$strand[1]

  for (i in 2:nrow(df)) {
    if (df$chr[i] == cur_chr && df$start[i] <= (cur_e + 1L + as.integer(max_gap))) {
      if (df$end[i] > cur_e) cur_e <- as.integer(df$end[i])
    } else {
      out_chr <- c(out_chr, cur_chr)
      out_s   <- c(out_s, cur_s)
      out_e   <- c(out_e, cur_e)
      out_st  <- c(out_st, cur_st)
      cur_chr <- df$chr[i]
      cur_s   <- as.integer(df$start[i])
      cur_e   <- as.integer(df$end[i])
      cur_st  <- df$strand[i]
    }
  }
  out_chr <- c(out_chr, cur_chr)
  out_s   <- c(out_s, cur_s)
  out_e   <- c(out_e, cur_e)
  out_st  <- c(out_st, cur_st)

  data.frame(
    chr    = out_chr,
    start  = out_s,
    end    = out_e,
    strand = out_st,
    stringsAsFactors = FALSE
  )
}

intersect_two_interval_dfs <- function(a, b) {
  if (is.null(a) || is.null(b) || !nrow(a) || !nrow(b)) return(empty_interval_df())
  out_chr <- character(0)
  out_s   <- integer(0)
  out_e   <- integer(0)
  out_st  <- character(0)

  for (i in seq_len(nrow(a))) {
    for (j in seq_len(nrow(b))) {
      if (a$chr[i] != b$chr[j]) next
      s <- max(as.integer(a$start[i]), as.integer(b$start[j]))
      e <- min(as.integer(a$end[i]), as.integer(b$end[j]))
      if (s <= e) {
        out_chr <- c(out_chr, a$chr[i])
        out_s   <- c(out_s, s)
        out_e   <- c(out_e, e)
        out_st  <- c(out_st, a$strand[i])
      }
    }
  }
  if (!length(out_s)) return(empty_interval_df())
  merge_intervals_df(data.frame(
    chr = out_chr, start = out_s, end = out_e, strand = out_st,
    stringsAsFactors = FALSE
  ))
}

# Compute shared CDS across coding isoforms (with fallback to largest overlapping
# subset or CDS union if isoforms are non-overlapping or only 1 coding isoform)
compute_shared_cds <- function(cds_by_tx) {
  coding_txs <- Filter(function(df) !is.null(df) && nrow(df) > 0, cds_by_tx)
  if (!length(coding_txs)) {
    res <- empty_interval_df()
    attr(res, "mode") <- "none"
    return(res)
  }
  merged_per_tx <- lapply(coding_txs, merge_intervals_df)
  if (length(merged_per_tx) == 1L) {
    res <- merged_per_tx[[1]]
    attr(res, "mode") <- "single_coding_isoform"
    return(res)
  }

  shared <- merged_per_tx[[1]]
  for (i in 2:length(merged_per_tx)) {
    shared <- intersect_two_interval_dfs(shared, merged_per_tx[[i]])
    if (!nrow(shared)) break
  }

  if (nrow(shared) > 0) {
    # Order in 5' -> 3' coding direction on the gene's strand
    strd <- shared$strand[1]
    shared <- if (identical(strd, "-")) {
      shared[order(shared$chr, -shared$end, -shared$start), , drop = FALSE]
    } else {
      shared[order(shared$chr, shared$start, shared$end), , drop = FALSE]
    }
    rownames(shared) <- NULL
    attr(shared, "mode") <- "intersection"
    return(shared)
  }

  # Fallback to CDS union if coding isoforms have mutually exclusive CDS regions
  all_cds <- do.call(rbind, lapply(merged_per_tx, function(d) d[, c("chr", "start", "end", "strand"), drop = FALSE]))
  union_df <- merge_intervals_df(all_cds)
  strd <- if (nrow(union_df)) union_df$strand[1] else "+"
  union_df <- if (identical(strd, "-")) {
    union_df[order(union_df$chr, -union_df$end, -union_df$start), , drop = FALSE]
  } else {
    union_df[order(union_df$chr, union_df$start, union_df$end), , drop = FALSE]
  }
  rownames(union_df) <- NULL
  attr(union_df, "mode") <- "fallback_union"
  union_df
}
# Extract 5' UTR, 3' UTR, promoter intervals, and 5'-most coding exons across
# all transcripts of a gene, strictly accounting for + vs - strand geometry.
compute_utr_promoter_and_first_cds <- function(exons_by_tx, cds_by_tx, promoter_bp = 500L) {
  utr5_list <- list()
  utr3_list <- list()
  prom_list <- list()
  fcds_list <- list()
  promoter_bp <- max(1L, as.integer(promoter_bp))

  tx_ids <- names(exons_by_tx)
  for (tid in tx_ids) {
    ex <- exons_by_tx[[tid]]
    if (is.null(ex) || !nrow(ex)) next
    strd <- ex$strand[1]
    chr  <- ex$chr[1]
    tx_s <- min(ex$start, na.rm = TRUE)
    tx_e <- max(ex$end, na.rm = TRUE)

    # Promoter: promoter_bp immediately upstream (5') of the TSS on the gene's strand
    if (identical(strd, "-")) {
      tss <- tx_e
      p_s <- tss + 1L
      p_e <- tss + promoter_bp
    } else {
      tss <- tx_s
      p_s <- max(1L, tss - promoter_bp)
      p_e <- max(1L, tss - 1L)
    }
    prom_list[[tid]] <- data.frame(
      chr           = chr,
      start         = as.integer(p_s),
      end           = as.integer(p_e),
      strand        = strd,
      transcript_id = tid,
      tss           = as.integer(tss),
      stringsAsFactors = FALSE
    )

    cds <- cds_by_tx[[tid]]
    if (is.null(cds) || !nrow(cds)) next

    cds_min <- min(cds$start, na.rm = TRUE)
    cds_max <- max(cds$end, na.rm = TRUE)

    # 5'-most coding exon (first CDS interval in 5' -> 3' coding order)
    if (identical(strd, "-")) {
      cds_ord <- cds[order(-cds$end, -cds$start), , drop = FALSE]
    } else {
      cds_ord <- cds[order(cds$start, cds$end), , drop = FALSE]
    }
    fcds_list[[tid]] <- cds_ord[1, c("chr", "start", "end", "strand", "phase", "transcript_id"), drop = FALSE]

    # 5' UTR and 3' UTR exonic intervals per exon
    for (r in seq_len(nrow(ex))) {
      es <- as.integer(ex$start[r])
      ee <- as.integer(ex$end[r])
      if (identical(strd, "-")) {
        # On minus strand: 5' is higher coordinates (> cds_max), 3' is lower coordinates (< cds_min)
        if (ee > cds_max) {
          u5_s <- max(es, cds_max + 1L)
          u5_e <- ee
          if (u5_s <= u5_e) {
            utr5_list[[length(utr5_list) + 1L]] <- data.frame(
              chr = chr, start = u5_s, end = u5_e, strand = strd,
              transcript_id = tid, stringsAsFactors = FALSE
            )
          }
        }
        if (es < cds_min) {
          u3_s <- es
          u3_e <- min(ee, cds_min - 1L)
          if (u3_s <= u3_e) {
            utr3_list[[length(utr3_list) + 1L]] <- data.frame(
              chr = chr, start = u3_s, end = u3_e, strand = strd,
              transcript_id = tid, stringsAsFactors = FALSE
            )
          }
        }
      } else {
        # On plus strand: 5' is lower coordinates (< cds_min), 3' is higher coordinates (> cds_max)
        if (es < cds_min) {
          u5_s <- es
          u5_e <- min(ee, cds_min - 1L)
          if (u5_s <= u5_e) {
            utr5_list[[length(utr5_list) + 1L]] <- data.frame(
              chr = chr, start = u5_s, end = u5_e, strand = strd,
              transcript_id = tid, stringsAsFactors = FALSE
            )
          }
        }
        if (ee > cds_max) {
          u3_s <- max(es, cds_max + 1L)
          u3_e <- ee
          if (u3_s <= u3_e) {
            utr3_list[[length(utr3_list) + 1L]] <- data.frame(
              chr = chr, start = u3_s, end = u3_e, strand = strd,
              transcript_id = tid, stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  five_prime_utr <- if (length(utr5_list)) do.call(rbind, utr5_list) else empty_interval_df(with_tx = TRUE)
  three_prime_utr <- if (length(utr3_list)) do.call(rbind, utr3_list) else empty_interval_df(with_tx = TRUE)
  promoter_regions <- if (length(prom_list)) do.call(rbind, prom_list) else {
    df <- empty_interval_df(with_tx = TRUE)
    df$tss <- integer(0)
    df
  }
  first_coding_exons <- if (length(fcds_list)) do.call(rbind, fcds_list) else empty_interval_df(with_tx = TRUE, with_phase = TRUE)

  rownames(five_prime_utr) <- NULL
  rownames(three_prime_utr) <- NULL
  rownames(promoter_regions) <- NULL
  rownames(first_coding_exons) <- NULL

  list(
    five_prime_utr     = five_prime_utr,
    three_prime_utr    = three_prime_utr,
    promoter_regions   = promoter_regions,
    first_coding_exons = first_coding_exons
  )
}

#' Convert exons_by_tx / cds_by_tx between data.frame-per-transcript and
#' HCR_probe_design_v41.R list-of-lists format for plot_isoform_tracks()
as_hcr_by_tx <- function(by_tx) {
  if (is.null(by_tx) || !length(by_tx)) return(NULL)
  lapply(by_tx, function(df) {
    if (is.data.frame(df)) {
      if (!nrow(df)) return(list())
      lapply(seq_len(nrow(df)), function(i) {
        list(
          chr           = as.character(df$chr[i]),
          start         = as.integer(df$start[i]),
          end           = as.integer(df$end[i]),
          strand        = as.character(df$strand[i]),
          transcript_id = as.character(df$transcript_id[i])
        )
      })
    } else if (is.list(df)) {
      df
    } else {
      list()
    }
  })
}
# ---- Complete Gene & Isoform Feature Model Lookup --------------------------

parse_variant_label <- function(attr_str, tx_id = "") {
  if (!length(attr_str) || !nzchar(attr_str)) return("")
  if (grepl("-R([A-Z0-9]+)", attr_str, perl = TRUE)) {
    return(sub(".*-R([A-Z0-9]+).*", "R\\1", attr_str, perl = TRUE))
  }
  if (grepl("transcript variant\\s+([A-Za-z0-9]+)", attr_str, ignore.case = TRUE, perl = TRUE)) {
    return(sub(".*transcript variant\\s+([A-Za-z0-9]+).*", "\\1", attr_str, ignore.case = TRUE, perl = TRUE))
  }
  if (grepl("isoform\\s+([A-Za-z0-9]+)", attr_str, ignore.case = TRUE, perl = TRUE)) {
    return(sub(".*isoform\\s+([A-Za-z0-9]+).*", "\\1", attr_str, ignore.case = TRUE, perl = TRUE))
  }
  tx_id
}

#' Resolve a single gene or genomic coordinate range and extract its complete
#' isoform & feature model from genomic.gtf(.gz)
#'
#' @param query Gene symbol, FBgn, CG, NCBI/LOC ID, synonym, product, or chr:start-end[:strand]
#' @param species Drosophila species name (default "melanogaster")
#' @param promoter_bp Upstream promoter window in bp (default 500L)
#' @param cache Optional environment for caching
#' @param ref_dir Reference genomes directory
#' @return A list containing gene metadata, transcripts, exons_by_tx, cds_by_tx,
#'   shared_cds, five_prime_utr, three_prime_utr, promoter_regions, and first_coding_exons
lookup_gene <- function(query, species = "melanogaster", promoter_bp = 500L,
                        cache = NULL, ref_dir = get_reference_dir()) {
  if (is.null(query) || !length(query) || is.na(query[1]) || !nzchar(trimws(query[1]))) {
    stop("Query must be a non-empty string.")
  }
  q <- trimws(as.character(query[1]))
  sp <- canonical_species(species)
  cache_env <- if (is.environment(cache)) cache else .GENE_LOOKUP_CACHE
  idx <- load_gene_index(sp, cache = cache_env, ref_dir = ref_dir)
  avail_chrs <- unique(idx$chromosome)

  coord <- parse_coordinate_query(q, species = sp, available_chrs = avail_chrs)
  hits <- search_genes(q, species = sp, max_results = 25L, cache = cache_env, ref_dir = ref_dir)
  if (!nrow(hits)) {
    stop(sprintf("No matching gene found for query '%s' in species '%s'.", q, sp))
  }

  best <- hits[1, , drop = FALSE]
  idx_row_i <- which(idx$gene_symbol == best$gene_symbol &
                     idx$chromosome  == best$chromosome &
                     idx$start       == best$start)[1]
  target_locus_id <- if (!is.na(idx_row_i) && "locus_id" %in% names(idx)) {
    idx$locus_id[idx_row_i]
  } else {
    best$gene_symbol
  }

  is_coord_query <- !is.null(coord)
  q_chr    <- if (is_coord_query) coord$chr else best$chromosome
  q_start  <- if (is_coord_query) min(coord$start, as.integer(best$start)) else as.integer(best$start)
  q_end    <- if (is_coord_query) max(coord$end,   as.integer(best$end))   else as.integer(best$end)
  q_strand <- if (is_coord_query && grepl(":[\\+\\-]$", q)) coord$strand else best$strand

  gtf_path <- find_gtf_file(sp, ref_dir = ref_dir)
  gtf_rows <- extract_gtf_locus_rows(
    gtf_path       = gtf_path,
    species        = sp,
    chr            = q_chr,
    start_pos      = q_start,
    end_pos        = q_end,
    target_gene_id = if (is_coord_query && !nzchar(best$flybase_id) && best$biotype == "genomic_region") NULL else target_locus_id,
    cache          = cache_env
  )

  if (nrow(gtf_rows) && nzchar(target_locus_id) && best$biotype != "genomic_region") {
    by_gid <- gtf_rows[gtf_rows$gene_id == target_locus_id, , drop = FALSE]
    if (!nrow(by_gid)) {
      g_attr <- extract_gtf_attr_vec(gtf_rows$attr, "gene")
      by_gid <- gtf_rows[tolower(g_attr) == tolower(best$gene_symbol), , drop = FALSE]
    }
    if (nrow(by_gid)) gtf_rows <- by_gid
  }

  exon_rows <- gtf_rows[gtf_rows$feature == "exon", , drop = FALSE]
  cds_rows  <- gtf_rows[gtf_rows$feature == "CDS", , drop = FALSE]
  tx_rows   <- gtf_rows[gtf_rows$feature %in% c("transcript", "mRNA"), , drop = FALSE]

  if (nrow(exon_rows)) {
    exon_rows$transcript_id <- extract_gtf_attr_vec(exon_rows$attr, "transcript_id")
    exon_rows$transcript_id[!nzchar(exon_rows$transcript_id)] <- best$gene_symbol
  }
  if (nrow(cds_rows)) {
    cds_rows$transcript_id <- extract_gtf_attr_vec(cds_rows$attr, "transcript_id")
    cds_rows$transcript_id[!nzchar(cds_rows$transcript_id)] <- best$gene_symbol
  }
  if (nrow(tx_rows)) {
    tx_rows$transcript_id <- extract_gtf_attr_vec(tx_rows$attr, "transcript_id")
  }

  tx_ids <- unique(c(
    if (nrow(tx_rows)) tx_rows$transcript_id else character(0),
    if (nrow(exon_rows)) exon_rows$transcript_id else character(0),
    if (nrow(cds_rows)) cds_rows$transcript_id else character(0)
  ))
  tx_ids <- tx_ids[nzchar(tx_ids)]

  exons_by_tx  <- list()
  cds_by_tx    <- list()
  tx_meta_list <- list()

  for (tid in tx_ids) {
    ex_t <- exon_rows[exon_rows$transcript_id == tid, , drop = FALSE]
    cd_t <- cds_rows[cds_rows$transcript_id == tid, , drop = FALSE]
    tx_t <- tx_rows[tx_rows$transcript_id == tid, , drop = FALSE]

    strd_t <- if (nrow(ex_t)) ex_t$strand[1] else if (nrow(tx_t)) tx_t$strand[1] else q_strand
    chr_t  <- if (nrow(ex_t)) ex_t$chr[1]    else if (nrow(tx_t)) tx_t$chr[1]    else q_chr

    if (nrow(ex_t)) {
      ex_t <- ex_t[order(ex_t$start, ex_t$end), , drop = FALSE]
      ex_df <- data.frame(
        chr           = as.character(ex_t$chr),
        start         = as.integer(ex_t$start),
        end           = as.integer(ex_t$end),
        strand        = as.character(ex_t$strand),
        transcript_id = rep(tid, nrow(ex_t)),
        stringsAsFactors = FALSE
      )
      rownames(ex_df) <- NULL
      exons_by_tx[[tid]] <- ex_df
    }

    if (nrow(cd_t)) {
      cd_t <- cd_t[order(cd_t$start, cd_t$end), , drop = FALSE]
      ph <- suppressWarnings(as.integer(cd_t$phase))
      ph[is.na(ph)] <- 0L
      cd_df <- data.frame(
        chr           = as.character(cd_t$chr),
        start         = as.integer(cd_t$start),
        end           = as.integer(cd_t$end),
        strand        = as.character(cd_t$strand),
        phase         = ph,
        transcript_id = rep(tid, nrow(cd_t)),
        stringsAsFactors = FALSE
      )
      rownames(cd_df) <- NULL
      cds_by_tx[[tid]] <- cd_df
    } else {
      cds_by_tx[[tid]] <- empty_interval_df(with_tx = TRUE, with_phase = TRUE)
    }

    attr_sample <- if (nrow(tx_t)) tx_t$attr[1] else if (nrow(ex_t)) ex_t$attr[1] else ""
    var_lbl <- parse_variant_label(attr_sample, tid)
    t_s <- if (nrow(ex_t)) min(ex_t$start) else if (nrow(tx_t)) min(tx_t$start) else q_start
    t_e <- if (nrow(ex_t)) max(ex_t$end)   else if (nrow(tx_t)) max(tx_t$end)   else q_end
    c_s <- if (nrow(cd_t)) min(cd_t$start) else NA_integer_
    c_e <- if (nrow(cd_t)) max(cd_t$end)   else NA_integer_

    tx_meta_list[[tid]] <- data.frame(
      transcript_id = tid,
      variant       = var_lbl,
      start         = as.integer(t_s),
      end           = as.integer(t_e),
      strand        = as.character(strd_t),
      n_exons       = as.integer(nrow(ex_t)),
      cds_start     = as.integer(c_s),
      cds_end       = as.integer(c_e),
      stringsAsFactors = FALSE
    )
  }

  transcripts_df <- if (length(tx_meta_list)) {
    df <- do.call(rbind, tx_meta_list)
    rownames(df) <- NULL
    df
  } else {
    data.frame(
      transcript_id = character(0), variant = character(0),
      start = integer(0), end = integer(0), strand = character(0),
      n_exons = integer(0), cds_start = integer(0), cds_end = integer(0),
      stringsAsFactors = FALSE
    )
  }

  shared_cds <- compute_shared_cds(cds_by_tx)
  utr_prom   <- compute_utr_promoter_and_first_cds(exons_by_tx, cds_by_tx, promoter_bp = promoter_bp)

  gene_meta <- list(
    symbol       = best$gene_symbol,
    gene_symbol  = best$gene_symbol,
    flybase_id   = best$flybase_id,
    cg_number    = best$cg_number,
    ncbi_gene_id = best$ncbi_gene_id,
    chr          = best$chromosome,
    chromosome   = best$chromosome,
    start        = as.integer(best$start),
    end          = as.integer(best$end),
    strand       = if (is_coord_query && grepl(":[\\+\\-]$", q)) coord$strand else best$strand,
    biotype      = best$biotype,
    synonyms     = best$synonyms,
    product      = best$product
  )

  list(
    gene               = gene_meta,
    symbol             = gene_meta$symbol,
    gene_symbol        = gene_meta$gene_symbol,
    flybase_id         = gene_meta$flybase_id,
    cg_number          = gene_meta$cg_number,
    ncbi_gene_id       = gene_meta$ncbi_gene_id,
    chr                = gene_meta$chr,
    chromosome         = gene_meta$chromosome,
    start              = gene_meta$start,
    end                = gene_meta$end,
    strand             = gene_meta$strand,
    biotype            = gene_meta$biotype,
    synonyms           = gene_meta$synonyms,
    product            = gene_meta$product,
    species            = sp,
    query              = q,
    match_type         = best$match_type,
    transcripts        = transcripts_df,
    exons_by_tx        = exons_by_tx,
    cds_by_tx          = cds_by_tx,
    shared_cds         = shared_cds,
    five_prime_utr     = utr_prom$five_prime_utr,
    three_prime_utr    = utr_prom$three_prime_utr,
    promoter_regions   = utr_prom$promoter_regions,
    first_coding_exons = utr_prom$first_coding_exons,
    hcr_exons_by_tx    = as_hcr_by_tx(exons_by_tx),
    hcr_cds_by_tx      = as_hcr_by_tx(cds_by_tx)
  )
}

