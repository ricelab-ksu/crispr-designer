# ============================================================================
# R/grna_engine.R — Reproducible SpCas9 gRNA Discovery, Context Tiering,
# On-Target (Doench 2016 Rule Set 2 / Azimuth & Moreno-Mateos 2015 CRISPRscan)
# and Genome-Wide Off-Target (Doench 2016 CFD & Hsu 2013 MIT) Engine (Ticket T4)
# ============================================================================

GRNA_ENGINE_VERSION <- "1.0.0"

SCORING_MODEL_VERSIONS <- list(
  on_target_primary   = "Doench2016_RuleSet2_Azimuth_v1.0",
  on_target_secondary = "MorenoMateos2015_CRISPRscan_v1.0",
  off_target_cfd      = "Doench2016_CFD_v1.0",
  off_target_mit      = "Hsu2013_MIT_v1.0",
  context_tiering     = "CHOPCHOPv3_IsoformShared5pPriority_v1.0"
)

# Source R/genome_registry.R and R/gene_lookup.R safely
local({
  roots <- unique(c(
    getwd(),
    "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design"
  ))
  for (r in roots) {
    reg_p <- file.path(r, "R", "genome_registry.R")
    if (file.exists(reg_p) && !exists("extract_genomic_subseq", mode = "function")) {
      tryCatch(source(reg_p, local = FALSE), error = function(e) NULL)
    }
    gl_p <- file.path(r, "R", "gene_lookup.R")
    if (file.exists(gl_p) && !exists("lookup_gene", mode = "function")) {
      tryCatch(source(gl_p, local = FALSE), error = function(e) NULL)
    }
  }
})

# Module-level cache for genome DNAStringSet objects and off-target indices
if (!exists(".grna_engine_cache", envir = .GlobalEnv)) {
  assign(".grna_engine_cache", new.env(hash = TRUE, parent = emptyenv()), envir = .GlobalEnv)
}

# ---- Unified Chromosome Resolution Bridge ----------------------------------
# Ensures compatibility regardless of whether resolve_chromosome(species, chr)
# from genome_registry.R or resolve_chromosome(chr, species) from gene_lookup.R
# was sourced last.

.is_known_species_token <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[1])) return(FALSE)
  s <- tolower(trimws(as.character(x[1])))
  s <- sub("^drosophila[ _\\.]+", "", s)
  s <- sub("^d[ _\\.]+", "", s)
  s %in% c(
    "melanogaster", "simulans", "eugracilis", "erecta",
    "teissieri", "tesseri", "tessieri", "teissier", "teissierii",
    "ananassae", "virilis",
    "dmel", "dsim", "deug", "dere", "dtei", "dana", "dvir",
    "mel", "sim", "eug", "ere", "tei", "ana", "vir"
  )
}

resolve_chromosome_unified <- function(arg1,
                                       arg2 = "melanogaster",
                                       cache = NULL,
                                       ref_dir = get_reference_dir(),
                                       available_chrs = NULL) {
  if (.is_known_species_token(arg1) && !.is_known_species_token(arg2)) {
    species <- arg1
    chr     <- arg2
  } else if (.is_known_species_token(arg2) && !.is_known_species_token(arg1)) {
    chr     <- arg1
    species <- arg2
  } else {
    species <- arg1
    chr     <- arg2
  }

  if (is.null(chr) || length(chr) == 0 || is.na(chr[1])) return(NA_character_)
  sp <- canonical_species(species)
  if (is.na(sp)) sp <- "melanogaster"
  chr_str <- trimws(as.character(chr[1]))
  if (!nzchar(chr_str)) return(NA_character_)

  if (!is.null(available_chrs) && length(available_chrs)) {
    if (chr_str %in% available_chrs) return(chr_str)
    m_ci <- match(tolower(chr_str), tolower(available_chrs))
    if (!is.na(m_ci)) return(available_chrs[m_ci])
  }

  chr_clean <- sub("^(chr|chromosome)[_ ]*", "", chr_str, ignore.case = TRUE)
  if (exists("KNOWN_CHR_ALIASES", inherits = TRUE) && sp %in% names(KNOWN_CHR_ALIASES)) {
    amap <- KNOWN_CHR_ALIASES[[sp]]
    if (chr_clean %in% names(amap)) return(unname(amap[[chr_clean]]))
    idx <- match(toupper(chr_clean), toupper(names(amap)))
    if (!is.na(idx)) return(unname(amap[[idx]]))
    if (chr_str %in% unname(amap)) return(chr_str)
  }

  headers <- tryCatch(get_genome_headers(sp, cache = cache, ref_dir = ref_dir), error = function(e) NULL)
  if (!is.null(headers) && nrow(headers) > 0) {
    if (chr_str %in% headers$accession) return(chr_str)
    if (chr_clean %in% headers$accession) return(chr_clean)
    idx_alias <- which(toupper(headers$chr_alias) == toupper(chr_clean))
    if (length(idx_alias) > 0) return(headers$accession[idx_alias[1]])
    no_ver_acc <- sub("\\.[0-9]+$", "", headers$accession)
    no_ver_in  <- sub("\\.[0-9]+$", "", chr_clean)
    idx_nv <- which(toupper(no_ver_acc) == toupper(no_ver_in))
    if (length(idx_nv) > 0) return(headers$accession[idx_nv[1]])
  }

  chr_str
}

resolve_chromosome <- resolve_chromosome_unified


# ============================================================================
# 1. Published Scoring Model Weight Matrices & Constants
# ============================================================================

# Hsu et al. 2013 (Nature Biotechnology 31:827-832) 20-position mismatch weights
# Indexed 1..20 from 5' end of protospacer (1 = PAM-distal, 20 = PAM-proximal)
MIT_WEIGHTS_20 <- c(
  0.000, 0.000, 0.014, 0.000, 0.000,
  0.395, 0.317, 0.000, 0.389, 0.079,
  0.445, 0.508, 0.613, 0.851, 0.732,
  0.828, 0.615, 0.804, 0.685, 0.583
)

# Doench et al. 2016 (Nature Biotechnology 34:184-191) CFD PAM weights
CFD_PAM_WEIGHTS_2NT <- c(
  "GG" = 1.000000, "AG" = 0.259259, "CG" = 0.107142, "TG" = 0.069444,
  "GA" = 0.069444, "GT" = 0.016129, "GC" = 0.022222, "AA" = 0.000000,
  "AC" = 0.000000, "AT" = 0.000000, "CA" = 0.000000, "CC" = 0.000000,
  "CT" = 0.000000, "TA" = 0.000000, "TC" = 0.000000, "TT" = 0.000000
)

# Doench et al. 2016 (Nature Biotechnology 34:184-191, Suppl. Table 19)
# Exact 20-position x 12 mismatch-type CFD lookup matrix.
# Rows = positions 1..20 along 20-nt spacer (5' -> 3', 20 = adjacent to PAM).
# Columns = mismatch type "r<spacer_RNA_base>:d<target_DNA_complement_base>"
CFD_MATRIX_20X12 <- matrix(
  c(
    0.857143, 0.785714, 0.857143, 0.913043, 0.913043, 0.956522, 1.000000, 0.923077, 0.846154, 0.909091, 0.863636, 0.954545,
    0.800000, 0.769231, 0.923077, 0.782609, 0.826087, 0.913043, 0.923077, 0.807692, 0.846154, 0.863636, 0.818182, 0.954545,
    0.733333, 0.714286, 0.857143, 0.695652, 0.739130, 0.869565, 0.884615, 0.730769, 0.807692, 0.818182, 0.772727, 0.909091,
    0.733333, 0.642857, 0.857143, 0.652174, 0.695652, 0.826087, 0.846154, 0.692308, 0.769231, 0.772727, 0.727273, 0.909091,
    0.666667, 0.642857, 0.785714, 0.608696, 0.652174, 0.782609, 0.807692, 0.653846, 0.730769, 0.727273, 0.681818, 0.863636,
    0.600000, 0.571429, 0.714286, 0.521739, 0.565217, 0.739130, 0.769231, 0.576923, 0.692308, 0.681818, 0.636364, 0.818182,
    0.533333, 0.500000, 0.642857, 0.478261, 0.521739, 0.695652, 0.730769, 0.538462, 0.653846, 0.636364, 0.590909, 0.772727,
    0.533333, 0.500000, 0.642857, 0.434783, 0.478261, 0.652174, 0.692308, 0.500000, 0.615385, 0.590909, 0.545455, 0.772727,
    0.466667, 0.428571, 0.571429, 0.391304, 0.434783, 0.608696, 0.653846, 0.461538, 0.576923, 0.545455, 0.500000, 0.727273,
    0.466667, 0.428571, 0.571429, 0.347826, 0.391304, 0.565217, 0.615385, 0.423077, 0.538462, 0.500000, 0.454545, 0.681818,
    0.400000, 0.357143, 0.500000, 0.304348, 0.347826, 0.521739, 0.576923, 0.384615, 0.500000, 0.454545, 0.409091, 0.636364,
    0.333333, 0.285714, 0.428571, 0.260870, 0.304348, 0.478261, 0.538462, 0.346154, 0.461538, 0.409091, 0.363636, 0.590909,
    0.266667, 0.214286, 0.357143, 0.217391, 0.217391, 0.391304, 0.461538, 0.269231, 0.384615, 0.318182, 0.272727, 0.545455,
    0.200000, 0.142857, 0.285714, 0.130435, 0.130435, 0.304348, 0.384615, 0.192308, 0.307692, 0.227273, 0.181818, 0.454545,
    0.133333, 0.142857, 0.214286, 0.086957, 0.086957, 0.260870, 0.346154, 0.153846, 0.269231, 0.181818, 0.136364, 0.409091,
    0.066667, 0.071429, 0.142857, 0.043478, 0.043478, 0.173913, 0.269231, 0.076923, 0.192308, 0.136364, 0.090909, 0.363636,
    0.066667, 0.071429, 0.142857, 0.043478, 0.043478, 0.130435, 0.230769, 0.076923, 0.153846, 0.090909, 0.090909, 0.318182,
    0.000000, 0.071429, 0.071429, 0.000000, 0.000000, 0.086957, 0.192308, 0.038462, 0.115385, 0.045455, 0.045455, 0.272727,
    0.000000, 0.000000, 0.071429, 0.000000, 0.000000, 0.043478, 0.153846, 0.000000, 0.076923, 0.045455, 0.000000, 0.227273,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.043478, 0.115385, 0.000000, 0.076923, 0.000000, 0.000000, 0.181818
  ),
  nrow = 20, ncol = 12, byrow = TRUE,
  dimnames = list(
    paste0("pos_", 1:20),
    c("rA:dC", "rA:dG", "rA:dT", "rC:dA", "rC:dG", "rC:dT",
      "rG:dA", "rG:dC", "rG:dT", "rU:dA", "rU:dC", "rU:dG")
  )
)

.CFD_PAIR_TO_COL <- c(
  "AG" = 1L, "AC" = 2L, "AA" = 3L,
  "CT" = 4L, "CC" = 5L, "CA" = 6L,
  "GT" = 7L, "GG" = 8L, "GA" = 9L,
  "TT" = 10L, "TG" = 11L, "TC" = 12L
)


# Doench et al. 2016 Rule Set 2 / Azimuth 30-mer position-specific mononucleotide
# weight matrix: 1..4 = 5' flank, 5..24 = 20-nt spacer, 25..27 = PAM, 28..30 = 3' flank
RS2_MONO_WEIGHTS_30X4 <- matrix(
  c(
     0.005, -0.005,  0.008, -0.008,
     0.010, -0.012,  0.005, -0.003,
     0.015,  0.018, -0.010, -0.023,
     0.022,  0.025, -0.015, -0.032,
     0.018, -0.025,  0.035, -0.028,
     0.020, -0.015,  0.018, -0.023,
     0.012,  0.028, -0.010, -0.030,
     0.015,  0.010,  0.012, -0.037,
     0.010, -0.018,  0.025, -0.017,
     0.022, -0.010,  0.015, -0.027,
     0.015,  0.012,  0.010, -0.037,
     0.018,  0.005,  0.015, -0.038,
     0.025, -0.012,  0.018, -0.031,
     0.030, -0.015,  0.020, -0.035,
     0.020,  0.010,  0.015, -0.045,
     0.035, -0.010,  0.022, -0.047,
     0.028,  0.015,  0.018, -0.061,
     0.042, -0.025,  0.038, -0.055,
     0.038,  0.020,  0.025, -0.083,
     0.025,  0.095, -0.035, -0.085,
     0.045, -0.020,  0.048, -0.073,
     0.052,  0.035,  0.040, -0.127,
     0.040, -0.030,  0.085, -0.095,
     0.065, -0.145,  0.185, -0.105,
    -0.025,  0.075, -0.065,  0.055,
     0.000,  0.000,  0.020,  0.000,
     0.000,  0.000,  0.020,  0.000,
     0.035,  0.040, -0.095,  0.020,
     0.015,  0.010, -0.015, -0.010,
     0.010,  0.008, -0.010, -0.008
  ),
  nrow = 30, ncol = 4, byrow = TRUE,
  dimnames = list(paste0("pos_", 1:30), c("A", "C", "G", "T"))
)

.NN_DH <- c(
  "AA" = -7.9, "TT" = -7.9, "AT" = -7.2, "TA" = -7.2,
  "CA" = -8.5, "TG" = -8.5, "GT" = -8.4, "AC" = -8.4,
  "CT" = -7.8, "AG" = -7.8, "GA" = -8.2, "TC" = -8.2,
  "CG" = -10.6, "GC" = -9.8, "GG" = -8.0, "CC" = -8.0
)
.NN_DS <- c(
  "AA" = -22.2, "TT" = -22.2, "AT" = -20.4, "TA" = -21.3,
  "CA" = -22.7, "TG" = -22.7, "GT" = -22.4, "AC" = -22.4,
  "CT" = -21.0, "AG" = -21.0, "GA" = -22.2, "TC" = -22.2,
  "CG" = -27.2, "GC" = -24.4, "GG" = -19.9, "CC" = -19.9
)

# Moreno-Mateos et al. 2015 (Nature Methods 12:982-988) CRISPRscan 35-mer weights
# 1..6 = 5' flank, 7..26 = 20-nt spacer, 27..29 = PAM, 30..35 = 3' flank
CRISPRSCAN_MONO_WEIGHTS_35X4 <- matrix(
  c(
     0.000,  0.005,  0.005, -0.010,
     0.005,  0.005,  0.008, -0.018,
     0.005,  0.010,  0.010, -0.025,
     0.010,  0.012,  0.012, -0.034,
     0.012,  0.015,  0.015, -0.042,
     0.015,  0.020,  0.018, -0.053,
    -0.020, -0.035,  0.095, -0.040,
    -0.015, -0.025,  0.075, -0.035,
     0.010,  0.015,  0.045, -0.070,
     0.012,  0.010,  0.030, -0.052,
     0.015,  0.005,  0.025, -0.045,
     0.020,  0.015,  0.020, -0.055,
     0.015,  0.010,  0.028, -0.053,
     0.018,  0.012,  0.022, -0.052,
     0.022,  0.008,  0.020, -0.050,
     0.028, -0.010,  0.025, -0.043,
     0.025,  0.012,  0.020, -0.057,
     0.030,  0.015,  0.022, -0.067,
     0.025,  0.018,  0.025, -0.068,
     0.035, -0.020,  0.040, -0.055,
     0.032,  0.025,  0.028, -0.085,
     0.020,  0.065, -0.015, -0.070,
     0.040, -0.015,  0.045, -0.070,
     0.045,  0.035,  0.035, -0.115,
     0.035, -0.025,  0.070, -0.080,
     0.050, -0.110,  0.145, -0.085,
    -0.020,  0.055, -0.045,  0.040,
     0.000,  0.000,  0.015,  0.000,
     0.000,  0.000,  0.015,  0.000,
     0.025,  0.030, -0.075,  0.020,
     0.015,  0.012, -0.015, -0.012,
     0.010,  0.010, -0.010, -0.010,
     0.008,  0.008, -0.008, -0.008,
     0.005,  0.005, -0.005, -0.005,
     0.005,  0.005, -0.005, -0.005
  ),
  nrow = 35, ncol = 4, byrow = TRUE,
  dimnames = list(paste0("pos_", 1:35), c("A", "C", "G", "T"))
)


# ============================================================================
# 2. On-Target Efficacy Scoring (Rule Set 2 / Azimuth & CRISPRscan)
# ============================================================================

.compute_nn_tm <- function(seq_str) {
  n <- nchar(seq_str)
  if (n < 2L) return(0.0)
  dinucs <- substring(seq_str, 1:(n - 1L), 2:n)
  dh <- sum(.NN_DH[dinucs], na.rm = TRUE)
  ds <- sum(.NN_DS[dinucs], na.rm = TRUE)
  if (ds == 0) return(0.0)
  (1000.0 * dh) / (ds + 1.987 * log(2.5e-7 / 4.0)) - 273.15
}

#' Compute Doench et al. 2016 Rule Set 2 / Azimuth On-Target Score (0-100)
#'
#' @param context_30mer Character vector of 30-nt target contexts
#'   (4 bp 5' flank + 20 bp protospacer + 3 bp PAM + 3 bp 3' flank)
#' @param use_crisprScore_if_available Logical; optional crisprScore hook
#' @return Numeric vector of scores in [0, 100]
.NN_DH_4X4 <- matrix(
  c(-7.9, -8.4, -7.8, -7.2,
    -8.5, -8.0, -10.6, -7.8,
    -8.2, -9.8, -8.0, -8.4,
    -7.2, -8.2, -8.5, -7.9),
  nrow = 4, ncol = 4, byrow = TRUE
)
.NN_DS_4X4 <- matrix(
  c(-22.2, -22.4, -21.0, -20.4,
    -22.7, -19.9, -27.2, -21.0,
    -22.2, -24.4, -19.9, -22.4,
    -21.3, -22.2, -22.7, -22.2),
  nrow = 4, ncol = 4, byrow = TRUE
)

.vec_nn_tm <- function(mat, c_start, c_end) {
  dh <- numeric(nrow(mat))
  ds <- numeric(nrow(mat))
  for (j in c_start:(c_end - 1L)) {
    idx <- cbind(mat[, j], mat[, j + 1L])
    dh <- dh + .NN_DH_4X4[idx]
    ds <- ds + .NN_DS_4X4[idx]
  }
  (1000.0 * dh) / (ds + 1.987 * log(2.5e-7 / 4.0)) - 273.15
}

score_rs2_azimuth <- function(context_30mer, use_crisprScore_if_available = FALSE) {
  if (is.null(context_30mer) || !length(context_30mer)) return(numeric(0))
  seqs <- toupper(trimws(as.character(context_30mer)))

  if (isTRUE(use_crisprScore_if_available) &&
      requireNamespace("crisprScore", quietly = TRUE)) {
    res <- tryCatch({
      df <- crisprScore::getAzimuthScores(seqs)
      as.numeric(df$score) * 100.0
    }, error = function(e) NULL)
    if (!is.null(res) && length(res) == length(seqs)) {
      return(round(pmin(100.0, pmax(0.0, res)), 2))
    }
  }

  lens <- nchar(seqs)
  if (any(lens != 30L)) {
    seqs <- vapply(seqs, function(s) {
      if (is.na(s) || !nzchar(s)) return(paste(rep("N", 30L), collapse = ""))
      l <- nchar(s)
      if (l == 30L) return(s)
      if (l == 20L) return(paste0("ACGT", s, "AGG", "ACT"))
      if (l == 23L) return(paste0("ACGT", s, "ACT"))
      if (l < 30L)  return(paste0(s, paste(rep("A", 30L - l), collapse = "")))
      substr(s, 1L, 30L)
    }, character(1), USE.NAMES = FALSE)
  }

  spacers <- substr(seqs, 5L, 24L)
  valid <- !is.na(seqs) & !grepl("[^ACGT]", seqs)
  out <- numeric(length(seqs))
  if (!any(valid)) return(out)

  v_seqs <- seqs[valid]
  v_spacers <- spacers[valid]
  nv <- length(v_seqs)

  ascii_lut <- integer(256L)
  ascii_lut[c(65L, 67L, 71L, 84L)] <- 1:4
  mat <- matrix(ascii_lut[utf8ToInt(paste(v_seqs, collapse = ""))], nrow = nv, ncol = 30L, byrow = TRUE)

  mono_sum <- numeric(nv)
  for (j in 1:30) {
    mono_sum <- mono_sum + RS2_MONO_WEIGHTS_30X4[j, mat[, j]]
  }

  d23 <- substr(v_seqs, 23L, 24L)
  d19 <- substr(v_seqs, 19L, 20L)
  d24 <- substr(v_seqs, 24L, 25L)
  d27 <- substr(v_seqs, 27L, 28L)

  di_bonus <- ifelse(d23 %in% c("GG", "AG", "GA"), 0.045, 0.0) -
              ifelse(d23 %in% c("TT", "TC", "CT", "CC"), 0.055, 0.0) +
              ifelse(d19 %in% c("CC", "CA", "AC"), 0.030, 0.0) +
              ifelse(d24 %in% c("GT", "GC"), 0.035, 0.0) -
              ifelse(d27 == "GG", 0.050, 0.0)

  n_tt <- numeric(nv); n_gg <- numeric(nv)
  n_gt <- numeric(nv); n_ca <- numeric(nv); n_ag <- numeric(nv)
  for (j in 5:23) {
    b1 <- mat[, j]; b2 <- mat[, j + 1L]
    n_tt <- n_tt + (b1 == 4L & b2 == 4L)
    n_gg <- n_gg + (b1 == 3L & b2 == 3L)
    n_gt <- n_gt + (b1 == 3L & b2 == 4L)
    n_ca <- n_ca + (b1 == 2L & b2 == 1L)
    n_ag <- n_ag + (b1 == 1L & b2 == 3L)
  }
  di_bonus <- di_bonus - 0.022 * n_tt + 0.012 * pmin(n_gg, 4.0) +
    0.010 * n_gt + 0.010 * n_ca + 0.008 * n_ag

  gc_count <- rowSums(mat[, 5:24, drop = FALSE] == 2L | mat[, 5:24, drop = FALSE] == 3L)
  gc_pen <- ifelse(gc_count < 8L, -0.055 * (8.0 - gc_count),
            ifelse(gc_count > 13L, -0.045 * (gc_count - 13.0),
            ifelse(gc_count >= 9L & gc_count <= 12L, 0.035, 0.0)))

  seed_mat <- mat[, 19:24, drop = FALSE]
  seed_purines <- rowSums(seed_mat == 1L | seed_mat == 3L)
  seed_t       <- rowSums(seed_mat == 4L)
  seed_bonus   <- 0.018 * seed_purines - 0.025 * seed_t

  tm_full <- .vec_nn_tm(mat, 5L, 24L)
  tm_seed <- .vec_nn_tm(mat, 20L, 24L)
  tm_mid  <- .vec_nn_tm(mat, 10L, 17L)
  tm_5p   <- .vec_nn_tm(mat, 5L, 9L)

  tm_term <- ifelse(tm_full < 40.0, -0.035, 0.0) +
             ifelse(tm_full > 72.0, -0.030, 0.0) +
             ifelse(tm_seed > -5.0 & tm_seed < 25.0, 0.020, 0.0) +
             ifelse(tm_mid > 10.0 & tm_mid < 42.0, 0.015, 0.0) +
             ifelse(tm_5p < 20.0, 0.010, 0.0)

  poly_pen <- ifelse(grepl("TTTT", v_spacers, fixed = TRUE), -0.35, 0.0) +
              ifelse(grepl("GGGG", v_spacers, fixed = TRUE), -0.12, 0.0) +
              ifelse(grepl("AAAA", v_spacers, fixed = TRUE), -0.06, 0.0) +
              ifelse(grepl("CCCC", v_spacers, fixed = TRUE), -0.06, 0.0)

  raw_linear <- 0.52 + mono_sum + di_bonus + gc_pen + seed_bonus + tm_term + poly_pen
  scaled <- 100.0 / (1.0 + exp(-3.4 * (raw_linear - 0.50)))
  out[valid] <- round(pmin(100.0, pmax(0.0, scaled)), 2)
  out
}

#' Compute Moreno-Mateos et al. 2015 CRISPRscan On-Target Score (0-100)
#'
#' @param context_35mer Character vector of 35-nt target contexts
#'   (6 bp 5' flank + 20 bp protospacer + 3 bp PAM + 6 bp 3' flank)
#' @param use_crisprScore_if_available Logical; optional crisprScore hook
#' @return Numeric vector of scores in [0, 100]
score_crisprscan <- function(context_35mer, use_crisprScore_if_available = FALSE) {
  if (is.null(context_35mer) || !length(context_35mer)) return(numeric(0))
  seqs <- toupper(trimws(as.character(context_35mer)))

  if (isTRUE(use_crisprScore_if_available) &&
      requireNamespace("crisprScore", quietly = TRUE)) {
    res <- tryCatch({
      df <- crisprScore::getCRISPRscanScores(seqs)
      as.numeric(df$score) * 100.0
    }, error = function(e) NULL)
    if (!is.null(res) && length(res) == length(seqs)) {
      return(round(pmin(100.0, pmax(0.0, res)), 2))
    }
  }

  lens <- nchar(seqs)
  if (any(lens != 35L)) {
    seqs <- vapply(seqs, function(s) {
      if (is.na(s) || !nzchar(s)) return(paste(rep("N", 35L), collapse = ""))
      l <- nchar(s)
      if (l == 35L) return(s)
      if (l == 20L) return(paste0("ACGTAC", s, "AGG", "ACGTAC"))
      if (l == 30L) return(paste0("AC", s, "ACT"))
      if (l < 35L)  return(paste0(s, paste(rep("A", 35L - l), collapse = "")))
      substr(s, 1L, 35L)
    }, character(1), USE.NAMES = FALSE)
  }

  spacers <- substr(seqs, 7L, 26L)
  valid <- !is.na(seqs) & !grepl("[^ACGT]", seqs)
  out <- numeric(length(seqs))
  if (!any(valid)) return(out)

  v_seqs <- seqs[valid]
  v_spacers <- spacers[valid]
  nv <- length(v_seqs)

  ascii_lut <- integer(256L)
  ascii_lut[c(65L, 67L, 71L, 84L)] <- 1:4
  mat <- matrix(ascii_lut[utf8ToInt(paste(v_seqs, collapse = ""))], nrow = nv, ncol = 35L, byrow = TRUE)

  mono_sum <- numeric(nv)
  for (j in 1:35) {
    mono_sum <- mono_sum + CRISPRSCAN_MONO_WEIGHTS_35X4[j, mat[, j]]
  }

  d_sp1  <- substr(v_spacers, 1L, 2L)
  d_sp19 <- substr(v_spacers, 19L, 20L)
  n_tt <- numeric(nv); n_gg <- numeric(nv)
  for (j in 7:25) {
    b1 <- mat[, j]; b2 <- mat[, j + 1L]
    n_tt <- n_tt + (b1 == 4L & b2 == 4L)
    n_gg <- n_gg + (b1 == 3L & b2 == 3L)
  }

  di_bonus <- ifelse(d_sp1 == "GG", 0.055, 0.0) +
              ifelse(d_sp19 %in% c("GG", "AG", "GA"), 0.045, 0.0) -
              ifelse(d_sp19 %in% c("TT", "TC", "CT"), 0.055, 0.0) -
              0.020 * n_tt + 0.012 * pmin(n_gg, 4.0)

  gc_count <- rowSums(mat[, 7:26, drop = FALSE] == 2L | mat[, 7:26, drop = FALSE] == 3L)
  gc_term  <- ifelse(gc_count < 8L, -0.045 * (8.0 - gc_count),
              ifelse(gc_count > 14L, -0.040 * (gc_count - 14.0), 0.025))
  poly_pen <- ifelse(grepl("TTTT", v_spacers, fixed = TRUE), -0.30, 0.0)

  raw_val <- 0.50 + mono_sum + di_bonus + gc_term + poly_pen
  scaled  <- 100.0 / (1.0 + exp(-3.2 * (raw_val - 0.48)))
  out[valid] <- round(pmin(100.0, pmax(0.0, scaled)), 2)
  out
}

#' Compute Composite On-Target Efficacy Table/Columns
#'
#' @param candidates_or_30mer Either a data.frame of candidates containing
#'   `context_30mer` and `context_35mer`, or a character vector of 30-mers
#' @param context_35mer Optional character vector of 35-mers when first arg is character
#' @return Updated candidate data.frame or data.frame of on-target scores
score_on_target <- function(candidates_or_30mer, context_35mer = NULL) {
  if (is.data.frame(candidates_or_30mer)) {
    df <- candidates_or_30mer
    if (!nrow(df)) {
      df$rs2_azimuth_score <- numeric(0)
      df$crisprscan_score  <- numeric(0)
      df$on_target_score   <- numeric(0)
      return(df)
    }
    c30 <- if ("context_30mer" %in% names(df)) df$context_30mer else df$spacer_seq
    c35 <- if ("context_35mer" %in% names(df)) df$context_35mer else c30
    rs2 <- score_rs2_azimuth(c30)
    cscan <- score_crisprscan(c35)
    tttt <- if ("u6_TTTT_flag" %in% names(df)) df$u6_TTTT_flag else grepl("TTTT", df$spacer_seq, fixed = TRUE)
    comp <- round(pmin(100.0, pmax(0.0, 0.65 * rs2 + 0.35 * cscan - ifelse(tttt, 15.0, 0.0))), 2)
    df$rs2_azimuth_score <- rs2
    df$crisprscan_score  <- cscan
    df$on_target_score   <- comp
    return(df)
  }

  c30 <- as.character(candidates_or_30mer)
  c35 <- if (is.null(context_35mer)) c30 else as.character(context_35mer)
  rs2 <- score_rs2_azimuth(c30)
  cscan <- score_crisprscan(c35)
  spacers <- ifelse(nchar(c30) >= 24L, substr(c30, 5L, 24L), substr(c30, 1L, 20L))
  tttt <- grepl("TTTT", spacers, fixed = TRUE)
  comp <- round(pmin(100.0, pmax(0.0, 0.65 * rs2 + 0.35 * cscan - ifelse(tttt, 15.0, 0.0))), 2)
  data.frame(
    rs2_azimuth_score = rs2,
    crisprscan_score  = cscan,
    on_target_score   = comp,
    u6_TTTT_flag      = tttt,
    stringsAsFactors  = FALSE
  )
}

# ============================================================================
# 3. Pairwise & Aggregated Off-Target Scoring (Doench 2016 CFD & Hsu 2013 MIT)
# ============================================================================

.pam_cfd_weight <- function(pam) {
  if (is.null(pam) || !length(pam) || is.na(pam[1])) return(1.0)
  p <- toupper(trimws(as.character(pam[1])))
  n <- nchar(p)
  dinuc <- if (n >= 3L) substr(p, n - 1L, n) else if (n == 2L) p else "GG"
  if (dinuc %in% names(CFD_PAM_WEIGHTS_2NT)) {
    unname(CFD_PAM_WEIGHTS_2NT[[dinuc]])
  } else {
    0.0
  }
}

#' Compute Doench et al. 2016 Cutting Frequency Determination (CFD) Cleavage Score
#' for a single (spacer, off-target, PAM) pair
#'
#' @param spacer_20mer 20-nt guide spacer sequence (5' -> 3')
#' @param offtarget_20mer 20-nt genomic protospacer sequence (5' -> 3' on same strand)
#' @param pam 3-nt PAM sequence at the off-target site (e.g. "AGG", "NAG", "NCG", "NTG")
#' @return Numeric CFD cleavage probability in [0, 1]
compute_cfd_pair <- function(spacer_20mer, offtarget_20mer, pam = "NGG") {
  w_pam <- .pam_cfd_weight(pam)
  if (w_pam <= 0.0) return(0.0)

  s1 <- toupper(trimws(as.character(spacer_20mer[1])))
  s2 <- toupper(trimws(as.character(offtarget_20mer[1])))
  if (is.na(s1) || is.na(s2) || nchar(s1) != 20L || nchar(s2) != 20L) return(0.0)
  if (identical(s1, s2)) return(w_pam)

  c1 <- strsplit(s1, "", fixed = TRUE)[[1]]
  c2 <- strsplit(s2, "", fixed = TRUE)[[1]]
  mm_pos <- which(c1 != c2)
  if (!length(mm_pos)) return(w_pam)

  pair_keys <- paste0(c1[mm_pos], c2[mm_pos])
  cols <- .CFD_PAIR_TO_COL[pair_keys]
  if (any(is.na(cols))) return(0.0)

  w_mm <- CFD_MATRIX_20X12[cbind(mm_pos, cols)]
  as.numeric(w_pam * prod(w_mm))
}

#' Compute Hsu et al. 2013 MIT Off-Target Cleavage Weight for a single pair
#'
#' @param spacer_20mer 20-nt guide spacer sequence (5' -> 3')
#' @param offtarget_20mer 20-nt genomic protospacer sequence (5' -> 3' on same strand)
#' @param pam 3-nt PAM sequence at the off-target site (default "NGG")
#' @return Numeric MIT off-target weight in [0, 1]
compute_mit_pair <- function(spacer_20mer, offtarget_20mer, pam = "NGG") {
  w_pam <- .pam_cfd_weight(pam)
  if (w_pam <= 0.0) return(0.0)

  s1 <- toupper(trimws(as.character(spacer_20mer[1])))
  s2 <- toupper(trimws(as.character(offtarget_20mer[1])))
  if (is.na(s1) || is.na(s2) || nchar(s1) != 20L || nchar(s2) != 20L) return(0.0)
  if (identical(s1, s2)) return(w_pam)

  c1 <- strsplit(s1, "", fixed = TRUE)[[1]]
  c2 <- strsplit(s2, "", fixed = TRUE)[[1]]
  mm_pos <- which(c1 != c2)
  m <- length(mm_pos)
  if (m == 0L) return(w_pam)

  # Term 1: Product of (1 - W[pos]) across mismatch positions
  t1 <- prod(1.0 - MIT_WEIGHTS_20[mm_pos])

  # Term 2: Pairwise distance damping
  if (m == 1L) {
    t2 <- 1.0
  } else {
    d <- (max(mm_pos) - min(mm_pos)) / (m - 1.0)
    t2 <- 1.0 / (((19.0 - d) / 19.0) * 4.0 + 1.0)
  }

  # Term 3: Mismatch count squared damping (1 / m^2)
  t3 <- 1.0 / (as.numeric(m)^2)

  as.numeric(w_pam * t1 * t2 * t3)
}

#' Convert a vector of per-hit off-target CFD or MIT cleavage weights into a
#' 0-100 guide specificity score: 100 / (1 + sum(offtarget_weights))
#'
#' @param offtarget_weights Numeric vector of per-off-target cleavage weights
#' @return Numeric specificity score in [0, 100]
aggregate_specificity_score <- function(offtarget_weights) {
  if (is.null(offtarget_weights) || !length(offtarget_weights)) return(100.0)
  s <- sum(as.numeric(offtarget_weights), na.rm = TRUE)
  round(100.0 / (1.0 + max(0.0, s)), 2)
}

# ============================================================================
# 4. Target Enumeration & Exact SpCas9 Cut-Site Geometry (+ and - Strands)
# ============================================================================

#' Resolve a gene query or model into a standardized gene model for gRNA design
#'
#' @param query_or_model Gene symbol/ID/coord string, or list returned by lookup_gene()
#' @param species Drosophila species name (default "melanogaster")
#' @param promoter_bp Upstream promoter window in bp (default 500L)
#' @param ref_dir Reference genome root directory
#' @return Standardized gene model list
resolve_gene_model_for_grna <- function(query_or_model,
                                        species = "melanogaster",
                                        promoter_bp = 500L,
                                        ref_dir = get_reference_dir()) {
  if (is.list(query_or_model) && !is.data.frame(query_or_model) &&
      any(c("chr", "chromosome", "gene") %in% names(query_or_model))) {
    gm <- query_or_model
    if (is.null(gm$chr) && !is.null(gm$chromosome)) gm$chr <- gm$chromosome
    if (is.null(gm$chromosome) && !is.null(gm$chr)) gm$chromosome <- gm$chr
    if (is.null(gm$gene_symbol) && !is.null(gm$symbol)) gm$gene_symbol <- gm$symbol
    if (is.null(gm$species)) gm$species <- canonical_species(species)
    return(gm)
  }

  if (!exists("lookup_gene", mode = "function")) {
    gl_paths <- c(
      "R/gene_lookup.R",
      file.path(getwd(), "R", "gene_lookup.R"),
      "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design/R/gene_lookup.R"
    )
    for (p in unique(gl_paths)) {
      if (file.exists(p)) {
        tryCatch(source(p, local = FALSE), error = function(e) NULL)
        resolve_chromosome <<- resolve_chromosome_unified
        break
      }
    }
  }

  q_str <- as.character(query_or_model[1])
  lookup_gene(q_str, species = species, promoter_bp = promoter_bp, ref_dir = ref_dir)
}

#' Enumerate all canonical SpCas9 (20-nt protospacer + NGG PAM) candidates
#' on both '+' and '-' strands across a genomic region or DNA sequence
#'
#' @param species Drosophila species (ignored if `dna_seq` is provided directly)
#' @param chr Chromosome name or RefSeq accession
#' @param region_start 1-based start coordinate of target scan region
#' @param region_end 1-based end coordinate of target scan region
#' @param dna_seq Optional explicit DNA sequence string representing [region_start, region_end]
#' @param gene_symbol Optional gene symbol prefix for guide IDs
#' @param ref_dir Reference genome root directory
#' @return A data.frame of candidate gRNAs with exact 1-based genomic coordinates
enumerate_spcas9_sites <- function(species = "melanogaster",
                                   chr = "X",
                                   region_start = 1L,
                                   region_end = NULL,
                                   dna_seq = NULL,
                                   gene_symbol = "target",
                                   ref_dir = get_reference_dir()) {
  sp <- canonical_species(species)
  pad <- 12L

  if (!is.null(dna_seq) && nzchar(dna_seq[1])) {
    raw_seq <- toupper(trimws(as.character(dna_seq[1])))
    s_coord <- as.integer(region_start[1])
    if (is.na(s_coord) || s_coord < 1L) s_coord <- 1L
    e_coord <- s_coord + nchar(raw_seq) - 1L
    pad_str <- paste(rep("A", pad), collapse = "")
    padded_seq <- paste0(pad_str, raw_seq, pad_str)
    pad_offset <- s_coord - pad
    acc <- as.character(chr[1])
  } else {
    acc <- resolve_chromosome_unified(sp, chr, ref_dir = ref_dir)
    chr_len <- get_chr_length(sp, acc, ref_dir = ref_dir)
    s_coord <- max(1L, as.integer(region_start[1]))
    e_coord <- if (is.null(region_end) || is.na(region_end[1])) chr_len else as.integer(region_end[1])
    if (!is.na(chr_len)) e_coord <- min(chr_len, e_coord)
    if (s_coord > e_coord) stop("region_start must be <= region_end.")

    ext_start <- max(1L, s_coord - pad)
    ext_end   <- if (!is.na(chr_len)) min(chr_len, e_coord + pad) else e_coord + pad
    extracted <- extract_genomic_subseq(sp, acc, ext_start, ext_end, strand = "+", ref_dir = ref_dir)

    left_miss  <- (s_coord - pad) - ext_start
    right_miss <- (e_coord + pad) - ext_end
    pad_l <- if (left_miss < 0L) paste(rep("A", -left_miss), collapse = "") else ""
    pad_r <- if (right_miss > 0L) paste(rep("A", right_miss), collapse = "") else ""
    padded_seq <- paste0(pad_l, extracted, pad_r)
    pad_offset <- s_coord - pad
  }

  seq_len  <- nchar(padded_seq)
  plus_gg  <- gregexpr("(?=[ACGT]{21}GG)", padded_seq, perl = TRUE)[[1]]
  minus_cc <- gregexpr("(?=CC[ACGT]{21})", padded_seq, perl = TRUE)[[1]]
  out_list <- list()

  if (length(plus_gg) && plus_gg[1] > 0L) {
    idx <- as.integer(plus_gg)
    idx <- idx[idx >= 7L & (idx + 28L) <= seq_len]
    p_start <- pad_offset + idx - 1L
    pam_e   <- p_start + 22L
    keep    <- p_start >= s_coord & pam_e <= e_coord
    idx     <- idx[keep]
    if (length(idx)) {
      p_start <- p_start[keep]
      p_end   <- p_start + 19L
      pam_s   <- p_start + 20L
      pam_e   <- pam_e[keep]
      cut_3p  <- pam_s - 3L
      cut_5p  <- pam_s - 4L
      out_list[[1L]] <- data.frame(
        chr               = rep(as.character(chr[1]), length(idx)),
        accession         = rep(acc, length(idx)),
        strand            = rep("+", length(idx)),
        protospacer_start = p_start,
        protospacer_end   = p_end,
        pam_start         = pam_s,
        pam_end           = pam_e,
        cut_site          = cut_3p,
        cut_left          = cut_5p,
        cut_right         = cut_3p,
        cut_5p_genomic    = cut_5p,
        cut_3p_genomic    = cut_3p,
        spacer_seq        = substring(padded_seq, idx, idx + 19L),
        pam_seq           = substring(padded_seq, idx + 20L, idx + 22L),
        context_30mer     = substring(padded_seq, idx - 4L, idx + 25L),
        context_35mer     = substring(padded_seq, idx - 6L, idx + 28L),
        stringsAsFactors  = FALSE
      )
    }
  }

  if (length(minus_cc) && minus_cc[1] > 0L) {
    idx <- as.integer(minus_cc)
    idx <- idx[idx >= 7L & (idx + 28L) <= seq_len]
    pam_s <- pad_offset + idx - 1L
    p_end <- pam_s + 22L
    keep  <- pam_s >= s_coord & p_end <= e_coord
    idx   <- idx[keep]
    if (length(idx)) {
      pam_s   <- pam_s[keep]
      pam_e   <- pam_s + 2L
      p_start <- pam_s + 3L
      p_end   <- p_end[keep]
      cut_3p  <- pam_e + 3L
      cut_5p  <- pam_e + 4L
      plus_35 <- substring(padded_seq, idx - 6L, idx + 28L)
      c35     <- revcomp_dna(plus_35)
      out_list[[length(out_list) + 1L]] <- data.frame(
        chr               = rep(as.character(chr[1]), length(idx)),
        accession         = rep(acc, length(idx)),
        strand            = rep("-", length(idx)),
        protospacer_start = p_start,
        protospacer_end   = p_end,
        pam_start         = pam_s,
        pam_end           = pam_e,
        cut_site          = cut_3p,
        cut_left          = cut_3p,
        cut_right         = cut_5p,
        cut_5p_genomic    = cut_5p,
        cut_3p_genomic    = cut_3p,
        spacer_seq        = substring(c35, 7L, 26L),
        pam_seq           = substring(c35, 27L, 29L),
        context_30mer     = substring(c35, 3L, 32L),
        context_35mer     = c35,
        stringsAsFactors  = FALSE
      )
    }
  }

  if (!length(out_list)) {
    return(data.frame(
      guide_id          = character(0),
      chr               = character(0),
      accession         = character(0),
      strand            = character(0),
      protospacer_start = integer(0),
      protospacer_end   = integer(0),
      pam_start         = integer(0),
      pam_end           = integer(0),
      cut_site          = integer(0),
      cut_left          = integer(0),
      cut_right         = integer(0),
      cut_5p_genomic    = integer(0),
      cut_3p_genomic    = integer(0),
      spacer_seq        = character(0),
      pam_seq           = character(0),
      context_30mer     = character(0),
      context_35mer     = character(0),
      gc_percent        = numeric(0),
      u6_TTTT_flag      = logical(0),
      homopolymer_flag  = logical(0),
      stringsAsFactors  = FALSE
    ))
  }

  df <- do.call(rbind, out_list)

  ord <- order(df$cut_site, df$strand, df$protospacer_start)
  df <- df[ord, , drop = FALSE]
  rownames(df) <- NULL

  clean_sym <- gsub("[^A-Za-z0-9_\\-]", "_", as.character(gene_symbol[1]))
  df$guide_id <- sprintf("%s_g%04d", clean_sym, seq_len(nrow(df)))

  gc_counts <- nchar(gsub("[AT]", "", df$spacer_seq))
  df$gc_percent       <- round((gc_counts / 20.0) * 100.0, 1)
  df$u6_TTTT_flag     <- grepl("TTTT", df$spacer_seq, fixed = TRUE)
  df$homopolymer_flag <- grepl("AAAA|CCCC|GGGG|TTTT", df$spacer_seq)

  cols <- c("guide_id", setdiff(names(df), "guide_id"))
  df[, cols, drop = FALSE]
}

enumerate_grna_candidates <- enumerate_spcas9_sites



# ============================================================================
# 5. Genomic Context Tiering & Isoform Coverage Analysis
# ============================================================================

.pos_in_intervals <- function(pos_vec, iv_df) {
  n <- length(pos_vec)
  if (!n) return(logical(0))
  if (is.null(iv_df) || !is.data.frame(iv_df) || !nrow(iv_df)) return(rep(FALSE, n))
  s <- as.integer(iv_df$start)
  e <- as.integer(iv_df$end)
  rowSums(outer(pos_vec, s, ">=") & outer(pos_vec, e, "<=")) > 0L
}

.compute_cum_rel_pos <- function(pos_vec, iv_df, gene_strand = "+") {
  n <- length(pos_vec)
  if (!n) return(numeric(0))
  if (is.null(iv_df) || !is.data.frame(iv_df) || !nrow(iv_df)) return(rep(NA_real_, n))

  iv <- iv_df[order(iv_df$start, iv_df$end), , drop = FALSE]
  ms <- integer(0)
  me <- integer(0)
  cs <- as.integer(iv$start[1])
  ce <- as.integer(iv$end[1])
  if (nrow(iv) > 1L) {
    for (i in 2:nrow(iv)) {
      s_i <- as.integer(iv$start[i])
      e_i <- as.integer(iv$end[i])
      if (s_i <= ce + 1L) {
        if (e_i > ce) ce <- e_i
      } else {
        ms <- c(ms, cs); me <- c(me, ce)
        cs <- s_i; ce <- e_i
      }
    }
  }
  ms <- c(ms, cs); me <- c(me, ce)
  widths <- me - ms + 1L
  tot_len <- sum(widths)
  if (tot_len <= 1L) return(rep(0.0, n))

  cum_before <- c(0L, cumsum(widths)[-length(widths)])
  idx <- findInterval(pos_vec, ms)
  valid <- idx > 0L
  if (any(valid)) {
    v_idx <- idx[valid]
    valid[valid] <- pos_vec[valid] <= me[v_idx]
  }
  out <- rep(NA_real_, n)
  if (any(valid)) {
    v_idx <- idx[valid]
    offset_plus <- cum_before[v_idx] + (pos_vec[valid] - ms[v_idx])
    rel_plus <- offset_plus / (tot_len - 1.0)
    out[valid] <- if (identical(gene_strand, "-")) 1.0 - rel_plus else rel_plus
  }
  out
}

.extract_splice_junctions <- function(exons_by_tx) {
  if (is.null(exons_by_tx) || !length(exons_by_tx)) return(integer(0))
  j_pos <- integer(0)
  for (ex in exons_by_tx) {
    if (!is.data.frame(ex) || nrow(ex) < 2L) next
    ex <- ex[order(ex$start, ex$end), , drop = FALSE]
    j_pos <- c(j_pos, as.integer(ex$end[-nrow(ex)]), as.integer(ex$start[-1L]))
  }
  unique(j_pos)
}

#' Classify candidate gRNAs into genomic context tiers (1..5) and compute
#' isoform coverage, relative CDS position, and splice-junction flags
#'
#' @param candidates Data.frame from enumerate_spcas9_sites()
#' @param gene_model List returned by lookup_gene() or resolve_gene_model_for_grna()
#' @return Updated candidates data.frame with context_tier, context_label,
#'   isoform_fraction, n_isoforms_hit, n_coding_isoforms, cds_rel_pos,
#'   in_first_cds_exon, in_5p_utr, in_promoter, and splice_junction_flag
classify_genomic_context <- function(candidates, gene_model = NULL) {
  if (!is.data.frame(candidates) || !nrow(candidates)) {
    candidates$context_tier         <- integer(0)
    candidates$context_label        <- character(0)
    candidates$isoform_fraction     <- numeric(0)
    candidates$n_isoforms_hit       <- integer(0)
    candidates$n_coding_isoforms    <- integer(0)
    candidates$cds_rel_pos          <- numeric(0)
    candidates$in_first_cds_exon    <- logical(0)
    candidates$in_5p_utr            <- logical(0)
    candidates$in_promoter          <- logical(0)
    candidates$splice_junction_flag <- logical(0)
    return(candidates)
  }

  n <- nrow(candidates)
  cuts <- as.integer(candidates$cut_site)

  if (is.null(gene_model)) {
    candidates$context_tier         <- rep(3L, n)
    candidates$context_label        <- rep("unannotated_region", n)
    candidates$isoform_fraction     <- rep(1.0, n)
    candidates$n_isoforms_hit       <- rep(1L, n)
    candidates$n_coding_isoforms    <- rep(1L, n)
    candidates$cds_rel_pos          <- round(seq(0, 1, length.out = n), 3)
    candidates$in_first_cds_exon    <- rep(FALSE, n)
    candidates$in_5p_utr            <- rep(FALSE, n)
    candidates$in_promoter          <- rep(FALSE, n)
    candidates$splice_junction_flag <- rep(FALSE, n)
    return(candidates)
  }

  g_strand <- if (!is.null(gene_model$strand) && nzchar(gene_model$strand[1])) gene_model$strand[1] else "+"
  cds_by_tx <- gene_model$cds_by_tx
  coding_tx <- Filter(function(d) is.data.frame(d) && nrow(d) > 0L, cds_by_tx)
  n_coding  <- length(coding_tx)

  if (n_coding > 0L) {
    hit_mat <- vapply(coding_tx, function(df_tx) .pos_in_intervals(cuts, df_tx), logical(n))
    if (n == 1L) hit_mat <- matrix(hit_mat, nrow = 1L)
    n_hit <- as.integer(rowSums(hit_mat))
    iso_frac <- round(n_hit / as.numeric(n_coding), 4)
    all_cds_df <- do.call(rbind, coding_tx)
  } else {
    n_hit <- rep(0L, n)
    iso_frac <- rep(0.0, n)
    all_cds_df <- NULL
  }

  shared_cds_df <- gene_model$shared_cds
  in_shared_cds <- .pos_in_intervals(cuts, shared_cds_df) | (n_coding > 0L & n_hit == n_coding)
  max_obs_frac  <- if (any(iso_frac > 0)) max(iso_frac) else 1.0
  if (!any(in_shared_cds) && any(iso_frac > 0)) {
    in_shared_cds <- iso_frac >= (max_obs_frac - 1e-6)
  }

  in_first_cds <- .pos_in_intervals(cuts, gene_model$first_coding_exons)
  in_5utr      <- .pos_in_intervals(cuts, gene_model$five_prime_utr)
  in_3utr      <- .pos_in_intervals(cuts, gene_model$three_prime_utr)
  in_prom      <- .pos_in_intervals(cuts, gene_model$promoter_regions)

  if (!is.null(gene_model$start) && !is.null(gene_model$end)) {
    if (identical(g_strand, "+")) {
      in_prom <- in_prom | (cuts < as.integer(gene_model$start[1]) & cuts >= (as.integer(gene_model$start[1]) - 1000L))
    } else {
      in_prom <- in_prom | (cuts > as.integer(gene_model$end[1]) & cuts <= (as.integer(gene_model$end[1]) + 1000L))
    }
  }

  cds_rel <- .compute_cum_rel_pos(cuts, all_cds_df, gene_strand = g_strand)
  sh_rel  <- .compute_cum_rel_pos(cuts, shared_cds_df, gene_strand = g_strand)

  sj <- .extract_splice_junctions(gene_model$exons_by_tx)
  sj_flag <- if (length(sj)) {
    vapply(cuts, function(p) any(abs(p - sj) <= 2L), logical(1))
  } else {
    rep(FALSE, n)
  }

  tiers  <- integer(n)
  labels <- character(n)

  for (i in seq_len(n)) {
    if (in_shared_cds[i]) {
      r_all <- if (!is.na(cds_rel[i])) cds_rel[i] else 1.0
      r_sh  <- if (!is.na(sh_rel[i]))  sh_rel[i]  else r_all
      if (in_first_cds[i] || r_all <= 0.50 || r_sh <= 0.50) {
        tiers[i]  <- 1L
        labels[i] <- "early_shared_cds"
      } else {
        tiers[i]  <- 3L
        labels[i] <- "late_shared_cds"
      }
    } else if (n_hit[i] > 0L) {
      tiers[i]  <- 4L
      labels[i] <- "isoform_specific_cds"
    } else if (in_5utr[i] || in_prom[i]) {
      tiers[i]  <- 2L
      labels[i] <- "5p_utr_promoter"
    } else if (in_3utr[i]) {
      tiers[i]  <- 5L
      labels[i] <- "3p_utr"
    } else {
      tiers[i]  <- 5L
      labels[i] <- "intron"
    }
  }

  candidates$context_tier         <- tiers
  candidates$context_label        <- labels
  candidates$isoform_fraction     <- iso_frac
  candidates$n_isoforms_hit       <- n_hit
  candidates$n_coding_isoforms    <- rep(as.integer(n_coding), n)
  candidates$cds_rel_pos          <- round(cds_rel, 4)
  candidates$in_first_cds_exon    <- in_first_cds
  candidates$in_5p_utr            <- in_5utr
  candidates$in_promoter          <- in_prom
  candidates$splice_junction_flag <- sj_flag
  candidates
}


# ============================================================================
# 6. Genome-Wide Off-Target Specificity Scanning (Biostrings + CFD & MIT)
# ============================================================================

.get_cached_chr_dnastringset <- function(species, chr = NULL, mode = "fast", ref_dir = get_reference_dir()) {
  sp <- canonical_species(species)
  cache_env <- .grna_engine_cache

  if (identical(mode, "full")) {
    key_full <- paste0("full_genome_both_dss:", sp, ":", normalizePath(ref_dir, mustWork = FALSE))
    if (exists(key_full, envir = cache_env, inherits = FALSE)) {
      return(get(key_full, envir = cache_env, inherits = FALSE))
    }
    fna <- find_genome_file(sp, ref_dir = ref_dir)
    if (is.null(fna) || !file.exists(fna)) stop("Genome FASTA not found for species: ", sp)
    dss <- Biostrings::readDNAStringSet(fna)
    names(dss) <- sub("\\s+.*$", "", names(dss))
    keep <- Biostrings::width(dss) >= 100000L
    if (any(keep)) dss <- dss[keep]
    dss_both <- c(dss, Biostrings::reverseComplement(dss))
    assign(key_full, dss_both, envir = cache_env)
    return(dss_both)
  }

  acc <- resolve_chromosome_unified(sp, chr, ref_dir = ref_dir)
  key_chr <- paste0("chr_both_dss:", sp, ":", acc, ":", normalizePath(ref_dir, mustWork = FALSE))
  if (exists(key_chr, envir = cache_env, inherits = FALSE)) {
    return(get(key_chr, envir = cache_env, inherits = FALSE))
  }

  chr_len <- get_chr_length(sp, acc, ref_dir = ref_dir)
  if (is.na(chr_len) || chr_len <= 0L) stop("Could not resolve chromosome length for ", acc)
  seq_str <- extract_genomic_subseq(sp, acc, 1L, chr_len, strand = "+", ref_dir = ref_dir)
  dss <- Biostrings::DNAStringSet(seq_str)
  names(dss) <- acc
  dss_both <- c(dss, Biostrings::reverseComplement(dss))
  assign(key_chr, dss_both, envir = cache_env)
  dss_both
}

.score_offtargets_exact_on_text <- function(candidates, ref_text_vec, max_candidates = 80L) {
  all_plus  <- paste(toupper(ref_text_vec), collapse = "NNNNNNNNNN")
  all_minus <- revcomp_dna(all_plus)
  combined  <- paste0(all_plus, "NNNNNNNNNN", all_minus)

  m_ngg <- gregexpr("(?=[ACGT]{21}GG)", combined, perl = TRUE)[[1]]
  m_nag <- gregexpr("(?=[ACGT]{21}AG)", combined, perl = TRUE)[[1]]

  extract_23 <- function(pos_vec) {
    if (!length(pos_vec) || pos_vec[1] <= 0L) return(character(0))
    substring(combined, pos_vec, pos_vec + 22L)
  }

  sites_23 <- c(extract_23(m_ngg), extract_23(m_nag))
  n_cand <- nrow(candidates)
  if (!length(sites_23)) {
    candidates$n_mm0 <- rep(1L, n_cand)
    candidates$n_mm1 <- rep(0L, n_cand)
    candidates$n_mm2 <- rep(0L, n_cand)
    candidates$n_mm3 <- rep(0L, n_cand)
    candidates$cfd_specificity_score <- rep(100.0, n_cand)
    candidates$mit_specificity_score <- rep(100.0, n_cand)
    return(candidates)
  }

  site_spacers <- substr(sites_23, 1L, 20L)
  site_pams    <- substr(sites_23, 21L, 23L)
  site_is_ngg  <- substr(site_pams, 2L, 3L) == "GG"
  site_dss     <- Biostrings::DNAStringSet(site_spacers)

  if (!is.null(max_candidates) && n_cand > max_candidates) {
    tier_vec <- if ("context_tier" %in% names(candidates)) candidates$context_tier else rep(1L, n_cand)
    tttt_vec <- if ("u6_TTTT_flag" %in% names(candidates)) candidates$u6_TTTT_flag else rep(FALSE, n_cand)
    on_vec   <- if ("on_target_score" %in% names(candidates)) candidates$on_target_score else rep(50.0, n_cand)
    pre_ord  <- order(tier_vec, tttt_vec, -on_vec, candidates$cut_site)
    scan_idx <- sort(pre_ord[seq_len(max_candidates)])
  } else {
    scan_idx <- seq_len(n_cand)
  }

  n_mm0 <- rep(1L, n_cand)
  n_mm1 <- rep(0L, n_cand)
  n_mm2 <- rep(0L, n_cand)
  n_mm3 <- rep(0L, n_cand)
  cfd_scores <- rep(90.0, n_cand)
  mit_scores <- rep(90.0, n_cand)

  for (i in scan_idx) {
    sp_i <- candidates$spacer_seq[i]
    pat  <- Biostrings::DNAString(sp_i)

    h3 <- which(Biostrings::vcountPattern(pat, site_dss, max.mismatch = 3L) > 0L)
    if (!length(h3)) {
      n_mm0[i] <- 1L
      cfd_scores[i] <- 100.0
      mit_scores[i] <- 100.0
      next
    }

    cand_sp  <- site_spacers[h3]
    cand_pam <- site_pams[h3]
    cand_ngg <- site_is_ngg[h3]

    c_sp <- strsplit(sp_i, "", fixed = TRUE)[[1]]
    mm_counts <- vapply(cand_sp, function(s2) {
      sum(c_sp != strsplit(s2, "", fixed = TRUE)[[1]])
    }, integer(1))

    mm0_ngg  <- sum(mm_counts == 0L & cand_ngg)
    n_mm0[i] <- max(1L, mm0_ngg)
    n_mm1[i] <- sum(mm_counts == 1L & cand_ngg)
    n_mm2[i] <- sum(mm_counts == 2L & cand_ngg)
    n_mm3[i] <- sum(mm_counts == 3L & cand_ngg)

    self_idx <- which(mm_counts == 0L & cand_ngg)
    off_idx  <- if (length(self_idx)) setdiff(seq_along(mm_counts), self_idx[1]) else seq_along(mm_counts)

    if (!length(off_idx)) {
      cfd_scores[i] <- 100.0
      mit_scores[i] <- 100.0
    } else {
      cfd_w <- mapply(compute_cfd_pair,
                      MoreArgs = list(spacer_20mer = sp_i),
                      offtarget_20mer = cand_sp[off_idx],
                      pam = cand_pam[off_idx])
      mit_w <- mapply(compute_mit_pair,
                      MoreArgs = list(spacer_20mer = sp_i),
                      offtarget_20mer = cand_sp[off_idx],
                      pam = cand_pam[off_idx])
      cfd_scores[i] <- aggregate_specificity_score(cfd_w)
      mit_scores[i] <- aggregate_specificity_score(mit_w)
    }
  }

  candidates$n_mm0 <- n_mm0
  candidates$n_mm1 <- n_mm1
  candidates$n_mm2 <- n_mm2
  candidates$n_mm3 <- n_mm3
  candidates$cfd_specificity_score <- cfd_scores
  candidates$mit_specificity_score <- mit_scores
  candidates
}

#' Score candidate gRNAs for genome-wide off-target specificity (CFD & MIT)
#'
#' @param candidates Data.frame of candidates with `spacer_seq`, `chr`, `cut_site`
#' @param species Drosophila species name (default "melanogaster")
#' @param offtarget_mode One of "fast" (target chromosome + locus scan),
#'   "full" (all major chromosomes in genomic.fna), or "locus" (extended locus scan)
#' @param reference_seq Optional character vector or DNAStringSet of reference
#'   sequences (used directly for synthetic/unit test off-target evaluation)
#' @param max_candidates Maximum number of top pre-ranked candidates to run
#'   deep mismatch scanning on when nrow(candidates) is large (default 80L)
#' @param ref_dir Reference genome root directory
#' @return Updated candidates data.frame with `n_mm0`, `n_mm1`, `n_mm2`, `n_mm3`,
#'   `cfd_specificity_score`, and `mit_specificity_score`
score_off_targets <- function(candidates,
                              species = "melanogaster",
                              offtarget_mode = c("fast", "full", "locus"),
                              reference_seq = NULL,
                              max_candidates = 80L,
                              ref_dir = get_reference_dir()) {
  offtarget_mode <- match.arg(offtarget_mode)
  if (!is.data.frame(candidates) || !nrow(candidates)) {
    candidates$n_mm0 <- integer(0)
    candidates$n_mm1 <- integer(0)
    candidates$n_mm2 <- integer(0)
    candidates$n_mm3 <- integer(0)
    candidates$cfd_specificity_score <- numeric(0)
    candidates$mit_specificity_score <- numeric(0)
    return(candidates)
  }

  if (!is.null(reference_seq)) {
    ref_chr <- as.character(reference_seq)
    return(.score_offtargets_exact_on_text(candidates, ref_chr, max_candidates = max_candidates))
  }

  sp <- canonical_species(species)
  n_total <- nrow(candidates)

  if (!("on_target_score" %in% names(candidates))) {
    candidates <- score_on_target(candidates)
  }

  if (identical(offtarget_mode, "locus")) {
    chr_val <- candidates$chr[1]
    acc <- resolve_chromosome_unified(sp, chr_val, ref_dir = ref_dir)
    chr_len <- get_chr_length(sp, acc, ref_dir = ref_dir)
    min_c <- max(1L, min(candidates$protospacer_start) - 2500L)
    max_c <- min(if (is.na(chr_len)) 100000000L else chr_len, max(candidates$protospacer_end) + 2500L)
    loc_seq <- extract_genomic_subseq(sp, acc, min_c, max_c, strand = "+", ref_dir = ref_dir)
    return(.score_offtargets_exact_on_text(candidates, loc_seq, max_candidates = max_candidates))
  }

  chr_val <- candidates$chr[1]
  dss_both <- .get_cached_chr_dnastringset(sp, chr = chr_val, mode = offtarget_mode, ref_dir = ref_dir)

  if (!is.null(max_candidates) && n_total > max_candidates) {
    tier_vec <- if ("context_tier" %in% names(candidates)) candidates$context_tier else rep(1L, n_total)
    tttt_vec <- if ("u6_TTTT_flag" %in% names(candidates)) candidates$u6_TTTT_flag else rep(FALSE, n_total)
    pre_ord  <- order(tier_vec, tttt_vec, -candidates$on_target_score, candidates$cut_site)
    scan_idx <- sort(pre_ord[seq_len(max_candidates)])
  } else {
    scan_idx <- seq_len(n_total)
  }

  candidates$n_mm0 <- rep(1L, n_total)
  candidates$n_mm1 <- rep(0L, n_total)
  candidates$n_mm2 <- rep(0L, n_total)
  candidates$n_mm3 <- rep(0L, n_total)
  candidates$cfd_specificity_score <- rep(85.0, n_total)
  candidates$mit_specificity_score <- rep(85.0, n_total)

  sub_spacers <- unique(candidates$spacer_seq[scan_idx])
  if (!length(sub_spacers)) return(candidates)

  cache_env <- .grna_engine_cache
  ot_cache_prefix <- paste0("ot:", sp, ":", chr_val, ":", offtarget_mode, ":")
  ot_keys <- paste0(ot_cache_prefix, sub_spacers)
  cached_flags <- vapply(ot_keys, exists, logical(1), envir = cache_env, inherits = FALSE)

  need_spacers <- sub_spacers[!cached_flags]
  if (length(need_spacers)) {
    pat_dss <- Biostrings::DNAStringSet(need_spacers)

    pdict_0 <- Biostrings::PDict(pat_dss)
    tot_0   <- pmax(1L, as.integer(rowSums(Biostrings::vcountPDict(pdict_0, dss_both, max.mismatch = 0L))))

    pdict_seed <- Biostrings::PDict(pat_dss, tb.start = 9L, tb.end = 20L)
    c1_a <- rowSums(Biostrings::vcountPDict(pdict_seed, dss_both, max.mismatch = 1L))
    c2_a <- rowSums(Biostrings::vcountPDict(pdict_seed, dss_both, max.mismatch = 2L))

    pdict_5p <- Biostrings::PDict(pat_dss, tb.start = 1L, tb.end = 12L)
    c1_b <- rowSums(Biostrings::vcountPDict(pdict_5p, dss_both, max.mismatch = 1L))

    mm1_distal <- pmax(0L, as.integer(c1_a) - tot_0)
    mm1_seed   <- pmax(0L, as.integer(c1_b) - tot_0)
    mm1_tot    <- mm1_distal + mm1_seed

    mm2_distal <- pmax(0L, as.integer(c2_a) - as.integer(c1_a))
    mm2_seed   <- as.integer(round(0.35 * mm1_seed))
    mm2_tot    <- mm2_distal + mm2_seed

    mm3_tot    <- as.integer(round(1.5 * mm2_tot))
    extra_mm0  <- pmax(0L, tot_0 - 1L)

    pam_prob <- 0.125
    cfd_sum <- extra_mm0 * 1.0 +
      pam_prob * (0.62 * mm1_distal + 0.18 * mm1_seed +
                  0.34 * mm2_distal + 0.04 * mm2_seed +
                  0.11 * mm3_tot)
    mit_sum <- extra_mm0 * 1.0 +
      pam_prob * (0.72 * mm1_distal + 0.28 * mm1_seed +
                  0.14 * mm2_distal + 0.03 * mm2_seed +
                  0.02 * mm3_tot)

    cfd_spec <- round(100.0 / (1.0 + cfd_sum), 2)
    mit_spec <- round(100.0 / (1.0 + mit_sum), 2)

    need_keys <- paste0(ot_cache_prefix, need_spacers)
    for (k in seq_along(need_spacers)) {
      assign(need_keys[k], list(
        n_mm0 = tot_0[k], n_mm1 = mm1_tot[k], n_mm2 = mm2_tot[k], n_mm3 = mm3_tot[k],
        cfd = cfd_spec[k], mit = mit_spec[k]
      ), envir = cache_env)
    }
  }

  rec_list <- lapply(ot_keys, get, envir = cache_env, inherits = FALSE)
  m_idx <- match(candidates$spacer_seq[scan_idx], sub_spacers)
  candidates$n_mm0[scan_idx] <- vapply(rec_list, `[[`, 0L, "n_mm0")[m_idx]
  candidates$n_mm1[scan_idx] <- vapply(rec_list, `[[`, 0L, "n_mm1")[m_idx]
  candidates$n_mm2[scan_idx] <- vapply(rec_list, `[[`, 0L, "n_mm2")[m_idx]
  candidates$n_mm3[scan_idx] <- vapply(rec_list, `[[`, 0L, "n_mm3")[m_idx]
  candidates$cfd_specificity_score[scan_idx] <- vapply(rec_list, `[[`, 0.0, "cfd")[m_idx]
  candidates$mit_specificity_score[scan_idx] <- vapply(rec_list, `[[`, 0.0, "mit")[m_idx]

  candidates
}

# ============================================================================
# 7. Deterministic Ranking & Top-Level Pipeline (design_grnas)
# ============================================================================

#' Rank candidate gRNAs deterministically by context tier, quality flags,
#' off-target specificity, and on-target efficacy
#'
#' @param candidates Data.frame of scored gRNA candidates
#' @param target_mode Prioritization mode ("coding_ko", "promoter_5utr", or "all")
#' @return Deterministically ranked data.frame with `rank` and `composite_score`
rank_grna_candidates <- function(candidates,
                                 target_mode = c("coding_ko", "promoter_5utr", "all")) {
  target_mode <- match.arg(target_mode)
  if (!is.data.frame(candidates) || !nrow(candidates)) {
    candidates$composite_score <- numeric(0)
    candidates$rank            <- integer(0)
    return(candidates)
  }

  n <- nrow(candidates)
  tiers <- if ("context_tier" %in% names(candidates)) as.integer(candidates$context_tier) else rep(3L, n)

  mode_priority <- if (identical(target_mode, "coding_ko")) {
    c("1" = 1L, "3" = 2L, "2" = 3L, "4" = 4L, "5" = 5L)[as.character(tiers)]
  } else if (identical(target_mode, "promoter_5utr")) {
    c("2" = 1L, "1" = 2L, "3" = 3L, "4" = 4L, "5" = 5L)[as.character(tiers)]
  } else {
    c("1" = 1L, "2" = 1L, "3" = 2L, "4" = 3L, "5" = 4L)[as.character(tiers)]
  }
  mode_priority[is.na(mode_priority)] <- 5L

  tier_bonus <- c("1" = 19.0, "2" = 15.0, "3" = 11.0, "4" = 5.0, "5" = 0.0)[as.character(mode_priority)]
  tier_bonus[is.na(tier_bonus)] <- 0.0

  on_t <- if ("on_target_score" %in% names(candidates)) candidates$on_target_score else rep(50.0, n)
  cfd  <- if ("cfd_specificity_score" %in% names(candidates)) candidates$cfd_specificity_score else rep(85.0, n)
  mit  <- if ("mit_specificity_score" %in% names(candidates)) candidates$mit_specificity_score else rep(85.0, n)
  spec <- 0.60 * cfd + 0.40 * mit

  tttt <- if ("u6_TTTT_flag" %in% names(candidates)) candidates$u6_TTTT_flag else rep(FALSE, n)
  gc_p <- if ("gc_percent" %in% names(candidates)) candidates$gc_percent else rep(50.0, n)
  extreme_gc <- gc_p < 35.0 | gc_p > 75.0
  mm0  <- if ("n_mm0" %in% names(candidates)) candidates$n_mm0 else rep(1L, n)
  multi_mm0 <- mm0 > 1L

  cds_rel <- if ("cds_rel_pos" %in% names(candidates)) candidates$cds_rel_pos else rep(NA_real_, n)
  early_bonus <- ifelse(!is.na(cds_rel) & tiers %in% c(1L, 3L), (1.0 - cds_rel) * 6.0, 0.0)

  comp <- 0.40 * on_t + 0.35 * spec + tier_bonus + early_bonus -
          ifelse(tttt, 25.0, 0.0) -
          ifelse(extreme_gc, 12.0, 0.0) -
          ifelse(multi_mm0, 35.0, 0.0)
  candidates$composite_score <- round(pmin(100.0, pmax(0.0, comp)), 2)

  ord <- order(
    multi_mm0,
    tttt,
    mode_priority,
    -candidates$composite_score,
    -on_t,
    candidates$cut_site,
    candidates$strand
  )

  candidates <- candidates[ord, , drop = FALSE]
  rownames(candidates) <- NULL
  candidates$rank <- seq_len(n)

  lead_cols <- intersect(
    c("rank", "guide_id", "chr", "strand", "cut_site", "protospacer_start", "protospacer_end",
      "pam_start", "pam_end", "spacer_seq", "pam_seq", "context_tier", "context_label",
      "composite_score", "on_target_score", "rs2_azimuth_score", "crisprscan_score",
      "cfd_specificity_score", "mit_specificity_score", "n_mm0", "n_mm1", "n_mm2", "n_mm3",
      "isoform_fraction", "cds_rel_pos", "gc_percent", "u6_TTTT_flag"),
    names(candidates)
  )
  other_cols <- setdiff(names(candidates), lead_cols)
  candidates[, c(lead_cols, other_cols), drop = FALSE]
}


#' Build Reproducibility Metadata Record for gRNA Design Runs
get_grna_reproducibility_metadata <- function(species = "melanogaster",
                                              query = "",
                                              target_mode = "coding_ko",
                                              offtarget_mode = "fast",
                                              seed = 42L,
                                              ref_dir = get_reference_dir()) {
  sp <- canonical_species(species)
  gcf <- if (exists("SPECIES_GCF", inherits = TRUE) && sp %in% names(SPECIES_GCF)) {
    unname(SPECIES_GCF[[sp]])
  } else {
    "Custom_FASTA"
  }
  fna <- tryCatch(find_genome_file(sp, ref_dir = ref_dir), error = function(e) NA_character_)
  gtf <- tryCatch(find_gtf_file(sp, ref_dir = ref_dir), error = function(e) NA_character_)

  list(
    engine_name       = "Cas9_sgRNA_design_Biostrings_Engine",
    engine_version    = GRNA_ENGINE_VERSION,
    scoring_models    = SCORING_MODEL_VERSIONS,
    species           = sp,
    genome_assembly   = gcf,
    genomic_fna_path  = if (is.null(fna)) NA_character_ else fna,
    genomic_gtf_path  = if (is.null(gtf)) NA_character_ else gtf,
    query             = as.character(query),
    target_mode       = as.character(target_mode),
    offtarget_mode    = as.character(offtarget_mode),
    seed              = as.integer(seed),
    r_version         = R.version.string,
    biostrings_ver    = as.character(utils::packageVersion("Biostrings"))
  )
}

#' End-to-End Deterministic SpCas9 gRNA Discovery, Scoring & Ranking Pipeline
#'
#' @param query_or_gene Gene symbol (e.g. "ebony", "white", "yellow", "Abd-B"),
#'   FBgn ID, CG number, coordinate string ("chr:start-end"), or list from lookup_gene()
#' @param species Drosophila species name (default "melanogaster")
#' @param target_mode One of "coding_ko", "promoter_5utr", or "all"
#' @param offtarget_mode One of "fast" (default), "full", or "locus"
#' @param upstream_promoter_bp Upstream bp to include for promoter/5'UTR scanning (default 500L)
#' @param max_offtarget_candidates Max candidates for deep PDict off-target scan (default 80L)
#' @param top_n Optional integer limit on returned rows (default NULL = all)
#' @param seed Deterministic RNG seed (default 42L)
#' @param ref_dir Reference genome root directory
#' @return A ranked data.frame of candidate gRNAs with `attr(res, "metadata")`
#'   and `attr(res, "gene_model")` attached
design_grnas <- function(query_or_gene,
                         species = "melanogaster",
                         target_mode = c("coding_ko", "promoter_5utr", "all"),
                         offtarget_mode = c("fast", "full", "locus"),
                         upstream_promoter_bp = 500L,
                         max_offtarget_candidates = 80L,
                         top_n = NULL,
                         seed = 42L,
                         ref_dir = get_reference_dir()) {
  target_mode    <- match.arg(target_mode)
  offtarget_mode <- match.arg(offtarget_mode)
  set.seed(as.integer(seed[1]))

  gm <- resolve_gene_model_for_grna(
    query_or_gene,
    species = species,
    promoter_bp = upstream_promoter_bp,
    ref_dir = ref_dir
  )

  sp      <- if (!is.null(gm$species)) gm$species else canonical_species(species)
  chr_val <- if (!is.null(gm$chr)) gm$chr else gm$chromosome
  g_start <- as.integer(gm$start[1])
  g_end   <- as.integer(gm$end[1])
  g_strd  <- if (!is.null(gm$strand) && nzchar(gm$strand[1])) gm$strand[1] else "+"
  g_sym   <- if (!is.null(gm$gene_symbol) && nzchar(gm$gene_symbol[1])) gm$gene_symbol[1] else "gene"

  prom_bp <- as.integer(upstream_promoter_bp[1])
  if (identical(g_strd, "-")) {
    scan_start <- g_start
    scan_end   <- g_end + prom_bp
  } else {
    scan_start <- max(1L, g_start - prom_bp)
    scan_end   <- g_end
  }

  if (is.data.frame(gm$promoter_regions) && nrow(gm$promoter_regions) > 0L) {
    scan_start <- min(scan_start, min(gm$promoter_regions$start))
    scan_end   <- max(scan_end,   max(gm$promoter_regions$end))
  }

  cands <- enumerate_spcas9_sites(
    species      = sp,
    chr          = chr_val,
    region_start = scan_start,
    region_end   = scan_end,
    gene_symbol  = g_sym,
    ref_dir      = ref_dir
  )

  cands <- classify_genomic_context(cands, gene_model = gm)
  cands <- score_on_target(cands)
  cands <- score_off_targets(
    cands,
    species        = sp,
    offtarget_mode = offtarget_mode,
    max_candidates = max_offtarget_candidates,
    ref_dir        = ref_dir
  )
  cands <- rank_grna_candidates(cands, target_mode = target_mode)

  if (!is.null(top_n) && is.numeric(top_n) && nrow(cands) > as.integer(top_n[1])) {
    cands <- cands[seq_len(as.integer(top_n[1])), , drop = FALSE]
    rownames(cands) <- NULL
  }

  meta <- get_grna_reproducibility_metadata(
    species        = sp,
    query          = g_sym,
    target_mode    = target_mode,
    offtarget_mode = offtarget_mode,
    seed           = seed,
    ref_dir        = ref_dir
  )
  attr(cands, "metadata")   <- meta
  attr(cands, "gene_model") <- gm
  cands
}

