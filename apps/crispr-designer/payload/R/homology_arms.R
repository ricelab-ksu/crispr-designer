# ============================================================================
# R/homology_arms.R — HDR Homology Arm Extraction, Cut-Site Geometry,
# Sequence QC, PAM/Seed Shielding & Multi-Format Export Engine (Ticket T5)
# ============================================================================
# Implements exact SpCas9 blunt cut-site geometry (3 bp 5' of PAM on the target
# strand) for both '+' and '-' strands, extracts 5' and 3' homology arms
# (default 1000 bp, configurable 1000-2000 bp or any positive integer),
# performs sequence QC (N tracts, homopolymers, low-complexity entropy windows,
# extreme GC, chromosome boundary truncation), detects donor re-cleavage risk
# and designs silent/PAM-blocking mutations, builds nucleotide-level cut-site
# zoom views, and exports FASTA, CSV, JSON, reproducibility logs, and ZIP bundles.
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
HOMOLOGY_ARMS_VERSION <- "1.0.0"

# Source R/genome_registry.R, R/gene_lookup.R, and R/grna_engine.R safely
local({
  roots <- unique(c(
    getwd(),
    app_root()
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
    ge_p <- file.path(r, "R", "grna_engine.R")
    if (file.exists(ge_p) && !exists("design_grnas", mode = "function")) {
      tryCatch(source(ge_p, local = FALSE), error = function(e) NULL)
    }
    pdb_p <- file.path(r, "R", "plasmid_designer_bundle.R")
    if (file.exists(pdb_p) && !exists("build_plasmid_designer_bundle", mode = "function")) {
      tryCatch(source(pdb_p, local = FALSE), error = function(e) NULL)
    }
  }
  if (exists("resolve_chromosome_unified", mode = "function")) {
    assign("resolve_chromosome", resolve_chromosome_unified, envir = .GlobalEnv)
  }
})

# ============================================================================
# 1. Standard Genetic Code & Synonymous Codon Table (for CDS PAM Shielding)
# ============================================================================

.STD_GENETIC_CODE <- c(
  "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
  "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
  "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
  "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
  "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
  "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
  "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
  "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
  "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
  "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
  "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
  "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
  "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
  "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
  "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
  "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
)

.DMEL_CODON_PREF <- c(
  "TTC" = 0.62, "TTT" = 0.38, "CTG" = 0.43, "CTC" = 0.15, "CTT" = 0.10, "TTG" = 0.18, "CTA" = 0.09, "TTA" = 0.05,
  "ATC" = 0.64, "ATT" = 0.31, "ATA" = 0.05, "ATG" = 1.00,
  "GTG" = 0.47, "GTC" = 0.24, "GTT" = 0.18, "GTA" = 0.11,
  "TCC" = 0.24, "TCG" = 0.20, "AGC" = 0.25, "TCT" = 0.09, "AGT" = 0.13, "TCA" = 0.09,
  "CCC" = 0.33, "CCA" = 0.25, "CCG" = 0.29, "CCT" = 0.13,
  "ACC" = 0.38, "ACG" = 0.26, "ACA" = 0.19, "ACT" = 0.17,
  "GCC" = 0.45, "GCT" = 0.19, "GCA" = 0.17, "GCG" = 0.19,
  "TAC" = 0.63, "TAT" = 0.37, "TAA" = 0.42, "TAG" = 0.32, "TGA" = 0.26,
  "CAC" = 0.60, "CAT" = 0.40, "CAG" = 0.70, "CAA" = 0.30,
  "AAC" = 0.56, "AAT" = 0.44, "AAG" = 0.71, "AAA" = 0.29,
  "GAC" = 0.47, "GAT" = 0.53, "GAG" = 0.67, "GAA" = 0.33,
  "TGC" = 0.71, "TGT" = 0.29, "TGG" = 1.00,
  "CGC" = 0.33, "CGG" = 0.15, "CGT" = 0.16, "AGA" = 0.09, "AGG" = 0.11, "CGA" = 0.16,
  "GGC" = 0.43, "GGA" = 0.29, "GGT" = 0.21, "GGG" = 0.07
)

# ============================================================================
# 2. Sequence QC Helpers (N-tracts, Homopolymers, Shannon Entropy, GC%)
# ============================================================================

.compute_gc_pct <- function(seq_str) {
  if (is.null(seq_str) || !length(seq_str) || is.na(seq_str[1]) || !nzchar(seq_str[1])) {
    return(0.0)
  }
  s <- toupper(seq_str[1])
  acgt_only <- gsub("[^ACGT]", "", s)
  denom <- nchar(acgt_only)
  if (denom == 0L) return(0.0)
  gc_only <- gsub("[^GC]", "", acgt_only)
  round(100.0 * nchar(gc_only) / as.numeric(denom), 2)
}

.max_char_run <- function(seq_str, chars = c("A", "C", "G", "T")) {
  if (is.null(seq_str) || !length(seq_str) || is.na(seq_str[1]) || !nzchar(seq_str[1])) {
    return(0L)
  }
  s <- toupper(seq_str[1])
  chars_vec <- strsplit(s, "", fixed = TRUE)[[1]]
  r <- rle(chars_vec)
  keep <- r$values %in% chars
  if (!any(keep)) return(0L)
  as.integer(max(r$lengths[keep]))
}

.max_dinuc_repeat_bp <- function(seq_str) {
  if (is.null(seq_str) || !length(seq_str) || is.na(seq_str[1]) || nchar(seq_str[1]) < 4L) {
    return(0L)
  }
  s <- toupper(seq_str[1])
  dinucs <- c("AT", "TA", "AC", "CA", "AG", "GA", "CT", "TC", "CG", "GC", "GT", "TG")
  max_bp <- 0L
  for (dn in dinucs) {
    m <- gregexpr(paste0("(", dn, "){3,}"), s, perl = TRUE)[[1]]
    if (length(m) && m[1] > 0L) {
      ml <- attr(m, "match.length")
      if (max(ml) > max_bp) max_bp <- as.integer(max(ml))
    }
  }
  max_bp
}

.window_shannon_entropy <- function(seq_str, window_bp = 50L, step_bp = 10L) {
  if (is.null(seq_str) || !length(seq_str) || is.na(seq_str[1]) || !nzchar(seq_str[1])) {
    return(list(min_entropy = 0.0, low_complexity_pct = 0.0))
  }
  s <- toupper(seq_str[1])
  len <- nchar(s)
  w <- min(as.integer(window_bp), len)
  if (w <= 0L) return(list(min_entropy = 0.0, low_complexity_pct = 0.0))

  starts <- if (len <= w) 1L else seq.int(1L, len - w + 1L, by = max(1L, as.integer(step_bp)))
  entropies <- numeric(length(starts))

  for (i in seq_along(starts)) {
    sub_s <- substr(s, starts[i], starts[i] + w - 1L)
    ch <- strsplit(sub_s, "", fixed = TRUE)[[1]]
    ch <- ch[ch %in% c("A", "C", "G", "T")]
    if (!length(ch)) {
      entropies[i] <- 0.0
    } else {
      counts <- tabulate(match(ch, c("A", "C", "G", "T")), nbins = 4L)
      probs <- counts[counts > 0L] / length(ch)
      entropies[i] <- -sum(probs * log2(probs))
    }
  }

  min_ent <- round(min(entropies), 3)
  low_pct <- round(100.0 * mean(entropies < 1.5), 2)
  list(min_entropy = min_ent, low_complexity_pct = low_pct)
}

.analyze_single_arm_qc <- function(seq_str, requested_len) {
  s <- if (is.null(seq_str) || is.na(seq_str[1])) "" else toupper(as.character(seq_str[1]))
  obs_len <- nchar(s)
  req_len <- as.integer(requested_len[1])

  n_only <- gsub("[^N]", "", s)
  n_count <- as.integer(nchar(n_only))
  max_n_run <- .max_char_run(s, chars = "N")
  has_n_tract <- (max_n_run >= 5L) || (n_count > 0L)

  max_hp <- .max_char_run(s, chars = c("A", "C", "G", "T"))
  max_dinuc <- .max_dinuc_repeat_bp(s)
  ent <- .window_shannon_entropy(s, window_bp = 50L, step_bp = 10L)

  repeat_warn <- (max_hp >= 8L) || (ent$min_entropy < 1.5) || (max_dinuc >= 12L)
  gc_pct <- .compute_gc_pct(s)
  gc_warn <- (obs_len > 0L) && (gc_pct < 25.0 || gc_pct > 75.0)
  truncated <- obs_len < req_len

  list(
    observed_len       = obs_len,
    requested_len      = req_len,
    gc_pct             = gc_pct,
    gc_warning         = gc_warn,
    n_count            = n_count,
    max_n_run          = max_n_run,
    has_n_tract        = has_n_tract,
    max_homopolymer    = max_hp,
    max_dinuc_repeat   = max_dinuc,
    min_window_entropy = ent$min_entropy,
    low_complexity_pct = ent$low_complexity_pct,
    repeat_warning     = repeat_warn,
    boundary_truncated = truncated
  )
}

# ============================================================================
# 3. Re-Cleavage Detection & Reading-Frame-Aware PAM/Seed Shielding Designer
# ============================================================================

.find_cds_codon_at_pos <- function(genomic_pos, gene_model, species, acc, cache = NULL, ref_dir = get_reference_dir()) {
  if (is.null(gene_model) || !is.list(gene_model) || is.null(gene_model$cds_by_tx)) {
    return(NULL)
  }
  cds_list <- Filter(function(d) is.data.frame(d) && nrow(d) > 0L, gene_model$cds_by_tx)
  if (!length(cds_list)) return(NULL)

  g_strand <- if (!is.null(gene_model$strand) && nzchar(gene_model$strand[1])) gene_model$strand[1] else "+"
  pos <- as.integer(genomic_pos[1])

  for (tx_name in names(cds_list)) {
    df <- cds_list[[tx_name]]
    hit_rows <- which(df$start <= pos & df$end >= pos)
    if (!length(hit_rows)) next

    if (identical(g_strand, "+")) {
      df <- df[order(df$start, df$end), , drop = FALSE]
      cds_coords <- unlist(Map(seq.int, as.integer(df$start), as.integer(df$end)), use.names = FALSE)
    } else {
      df <- df[order(-df$end, -df$start), , drop = FALSE]
      cds_coords <- unlist(Map(function(s, e) seq.int(e, s, by = -1L),
                               as.integer(df$start), as.integer(df$end)), use.names = FALSE)
    }

    phase0 <- if ("phase" %in% names(df)) suppressWarnings(as.integer(df$phase[1])) else 0L
    if (is.na(phase0) || phase0 < 0L || phase0 > 2L) phase0 <- 0L
    if (phase0 > 0L && length(cds_coords) > phase0) {
      cds_coords <- cds_coords[(phase0 + 1L):length(cds_coords)]
    }

    idx <- match(pos, cds_coords)
    if (is.na(idx)) next
    codon_idx <- ((idx - 1L) %/% 3L) + 1L
    codon_start_i <- (codon_idx - 1L) * 3L + 1L
    if ((codon_start_i + 2L) > length(cds_coords)) next

    codon_genomic <- cds_coords[codon_start_i:(codon_start_i + 2L)]
    frame_pos <- ((idx - 1L) %% 3L) + 1L

    bases <- vapply(codon_genomic, function(gp) {
      extract_genomic_subseq(species, acc, gp, gp, strand = g_strand, cache = cache, ref_dir = ref_dir)
    }, character(1))
    wt_codon <- toupper(paste(bases, collapse = ""))
    if (nchar(wt_codon) != 3L || grepl("[^ACGT]", wt_codon)) next

    return(list(
      transcript_id = tx_name,
      gene_strand   = g_strand,
      codon_index   = as.integer(codon_idx),
      frame_pos     = as.integer(frame_pos),
      codon_genomic = as.integer(codon_genomic),
      wt_codon      = wt_codon,
      wt_aa         = unname(.STD_GENETIC_CODE[wt_codon])
    ))
  }
  NULL
}

.comp_single_base <- function(b) {
  chartr("ACGTacgt", "TGCAtgca", b)
}

.eval_shield_codon <- function(codon_info, alt_target_base, target_strand) {
  if (is.null(codon_info)) {
    return(list(
      is_synonymous = FALSE,
      in_cds        = FALSE,
      wt_codon      = NA_character_,
      mut_codon     = NA_character_,
      aa_change     = "noncoding",
      codon_pref    = 0.50
    ))
  }
  wt_cod <- codon_info$wt_codon
  cod_base <- if (identical(codon_info$gene_strand, target_strand)) alt_target_base else .comp_single_base(alt_target_base)
  mut_cod <- wt_cod
  substr(mut_cod, codon_info$frame_pos, codon_info$frame_pos) <- cod_base
  wt_aa  <- codon_info$wt_aa
  mut_aa <- unname(.STD_GENETIC_CODE[mut_cod])
  is_syn <- !is.na(wt_aa) && !is.na(mut_aa) && identical(wt_aa, mut_aa)
  aa_str <- if (is_syn) {
    paste0("p.", wt_aa, codon_info$codon_index, "=")
  } else {
    paste0("p.", wt_aa, codon_info$codon_index, mut_aa)
  }
  cpref <- if (!is.na(.DMEL_CODON_PREF[mut_cod])) unname(.DMEL_CODON_PREF[mut_cod]) else 0.50
  list(
    is_synonymous = is_syn,
    in_cds        = TRUE,
    wt_codon      = wt_cod,
    mut_codon     = mut_cod,
    aa_change     = aa_str,
    codon_pref    = cpref
  )
}

.design_pam_shielding <- function(spacer_seq,
                                  pam_seq,
                                  strand,
                                  protospacer_start,
                                  protospacer_end,
                                  pam_start,
                                  pam_end,
                                  arm_3p_seq,
                                  arm_3p_anchor = "cut_site",
                                  gene_model = NULL,
                                  species = "melanogaster",
                                  acc = NULL,
                                  cache = NULL,
                                  ref_dir = get_reference_dir()) {
  pam_u <- toupper(as.character(pam_seq[1]))
  sp_u  <- toupper(as.character(spacer_seq[1]))
  strd  <- as.character(strand[1])

  if (identical(strd, "+")) {
    pam_genomic  <- c(as.integer(pam_start), as.integer(pam_start) + 1L, as.integer(pam_start) + 2L)
    seed_genomic <- c(as.integer(pam_start) - 3L, as.integer(pam_start) - 2L, as.integer(pam_start) - 1L)
  } else {
    pam_genomic  <- c(as.integer(pam_end), as.integer(pam_end) - 1L, as.integer(pam_end) - 2L)
    seed_genomic <- c(as.integer(pam_end) + 3L, as.integer(pam_end) + 2L, as.integer(pam_end) + 1L)
  }

  candidates <- list()

  for (pam_pos_idx in c(3L, 2L)) {
    gpos <- pam_genomic[pam_pos_idx]
    wt_target_base <- substr(pam_u, pam_pos_idx, pam_pos_idx)
    alt_target_bases <- setdiff(c("T", "C", "A"), wt_target_base)
    codon_info <- if (!is.null(acc)) {
      tryCatch(.find_cds_codon_at_pos(gpos, gene_model, species, acc, cache = cache, ref_dir = ref_dir),
               error = function(e) NULL)
    } else NULL

    for (alt_b in alt_target_bases) {
      if (pam_pos_idx == 2L && alt_b == "A") next
      mut_pam <- pam_u
      substr(mut_pam, pam_pos_idx, pam_pos_idx) <- alt_b
      ev <- .eval_shield_codon(codon_info, alt_b, strd)
      candidates[[length(candidates) + 1L]] <- list(
        mutation_class   = "PAM_disruption",
        target_element   = paste0("PAM_pos", pam_pos_idx),
        genomic_pos      = as.integer(gpos),
        wt_target_base   = wt_target_base,
        mut_target_base  = alt_b,
        wt_pam           = pam_u,
        mut_pam          = mut_pam,
        wt_spacer        = sp_u,
        mut_spacer       = sp_u,
        is_synonymous    = ev$is_synonymous,
        in_cds           = ev$in_cds,
        wt_codon         = ev$wt_codon,
        mut_codon        = ev$mut_codon,
        aa_change        = ev$aa_change,
        codon_preference = ev$codon_pref,
        arm_3p_offset    = if (identical(arm_3p_anchor, "cut_site")) 3L + pam_pos_idx else pam_pos_idx
      )
    }
  }

  for (seed_i in c(3L, 2L, 1L)) {
    sp_pos <- 17L + seed_i
    gpos   <- seed_genomic[seed_i]
    wt_b   <- substr(sp_u, sp_pos, sp_pos)
    alt_bs <- setdiff(c("A", "C", "G", "T"), wt_b)
    codon_info <- if (!is.null(acc)) {
      tryCatch(.find_cds_codon_at_pos(gpos, gene_model, species, acc, cache = cache, ref_dir = ref_dir),
               error = function(e) NULL)
    } else NULL

    for (alt_b in alt_bs) {
      ev <- .eval_shield_codon(codon_info, alt_b, strd)
      mut_sp <- sp_u
      substr(mut_sp, sp_pos, sp_pos) <- alt_b
      candidates[[length(candidates) + 1L]] <- list(
        mutation_class   = "seed_disruption",
        target_element   = paste0("spacer_pos", sp_pos),
        genomic_pos      = as.integer(gpos),
        wt_target_base   = wt_b,
        mut_target_base  = alt_b,
        wt_pam           = pam_u,
        mut_pam          = pam_u,
        wt_spacer        = sp_u,
        mut_spacer       = mut_sp,
        is_synonymous    = ev$is_synonymous,
        in_cds           = ev$in_cds,
        wt_codon         = ev$wt_codon,
        mut_codon        = ev$mut_codon,
        aa_change        = ev$aa_change,
        codon_preference = ev$codon_pref,
        arm_3p_offset    = if (identical(arm_3p_anchor, "cut_site")) seed_i else NA_integer_
      )
    }
  }

  scores <- vapply(candidates, function(x) {
    syn_bonus <- if (x$in_cds) (if (x$is_synonymous) 100.0 else 0.0) else 80.0
    pam_bonus <- if (identical(x$mutation_class, "PAM_disruption")) 30.0 else 10.0
    pos_bonus <- if (identical(x$target_element, "PAM_pos3")) 5.0 else 0.0
    syn_bonus + pam_bonus + pos_bonus + x$codon_preference
  }, numeric(1))

  ord <- order(-scores)
  best <- candidates[[ord[1L]]]

  shielded_arm_3p <- arm_3p_seq
  if (!is.na(best$arm_3p_offset) && best$arm_3p_offset >= 1L && best$arm_3p_offset <= nchar(shielded_arm_3p)) {
    substr(shielded_arm_3p, best$arm_3p_offset, best$arm_3p_offset) <- best$mut_target_base
  }

  best$shielded_arm_3p_seq <- shielded_arm_3p
  best$summary <- if (identical(best$mutation_class, "PAM_disruption")) {
    sprintf("%s (%s -> %s at genomic %d, %s%s)",
            best$target_element, best$wt_pam, best$mut_pam, best$genomic_pos,
            best$aa_change,
            if (!is.na(best$wt_codon)) sprintf(", %s->%s", best$wt_codon, best$mut_codon) else "")
  } else {
    sprintf("%s (%s>%s at genomic %d, %s%s)",
            best$target_element, best$wt_target_base, best$mut_target_base, best$genomic_pos,
            best$aa_change,
            if (!is.na(best$wt_codon)) sprintf(", %s->%s", best$wt_codon, best$mut_codon) else "")
  }
  best$alternatives <- lapply(candidates[ord[seq_len(min(5L, length(ord)))]], function(x) {
    x[c("mutation_class", "target_element", "genomic_pos", "wt_target_base",
        "mut_target_base", "wt_pam", "mut_pam", "is_synonymous", "wt_codon",
        "mut_codon", "aa_change")]
  })
  best
}

.has_spacer_pam_match <- function(seq_str, spacer_seq) {
  if (is.null(seq_str) || !nzchar(seq_str[1]) || is.null(spacer_seq) || nchar(spacer_seq[1]) != 20L) {
    return(FALSE)
  }
  s  <- toupper(seq_str[1])
  sp <- toupper(spacer_seq[1])
  pat_fwd <- paste0(sp, "[ACGT][AG]G")
  rc_s <- revcomp_dna(s)
  grepl(pat_fwd, s, perl = TRUE) || grepl(pat_fwd, rc_s, perl = TRUE)
}

# ============================================================================
# 4. Guide Parameter Normalization & Core design_homology_arms() Engine
# ============================================================================

.normalize_guide_params <- function(guide_row_or_params,
                                    species = "melanogaster",
                                    guide_index = 1L,
                                    cache = NULL,
                                    ref_dir = get_reference_dir()) {
  if (is.data.frame(guide_row_or_params)) {
    if (nrow(guide_row_or_params) == 0L) {
      stop("guide_row_or_params data.frame is empty (0 rows).")
    }
    r_idx <- max(1L, min(nrow(guide_row_or_params), as.integer(guide_index[1])))
    row_l <- as.list(guide_row_or_params[r_idx, , drop = FALSE])
    attr_gm   <- attr(guide_row_or_params, "gene_model")
    attr_meta <- attr(guide_row_or_params, "metadata")
  } else if (is.list(guide_row_or_params)) {
    row_l <- guide_row_or_params
    attr_gm   <- row_l$gene_model
    attr_meta <- row_l$metadata
  } else {
    stop("guide_row_or_params must be a data.frame row or named list.")
  }

  sp_raw <- if (!is.null(row_l$species) && !is.na(row_l$species[1]) && nzchar(row_l$species[1])) {
    row_l$species[1]
  } else if (!is.null(attr_meta$species) && nzchar(attr_meta$species[1])) {
    attr_meta$species[1]
  } else {
    species
  }
  sp <- canonical_species(sp_raw)
  if (is.na(sp) || !nzchar(sp)) sp <- "melanogaster"

  chr_in <- if (!is.null(row_l$chr) && !is.na(row_l$chr[1]) && nzchar(row_l$chr[1])) {
    as.character(row_l$chr[1])
  } else if (!is.null(row_l$accession) && !is.na(row_l$accession[1]) && nzchar(row_l$accession[1])) {
    as.character(row_l$accession[1])
  } else {
    stop("Missing 'chr' or 'accession' in guide_row_or_params.")
  }

  acc <- if (!is.null(row_l$accession) && !is.na(row_l$accession[1]) && nzchar(row_l$accession[1])) {
    resolve_chromosome_unified(sp, as.character(row_l$accession[1]), cache = cache, ref_dir = ref_dir)
  } else {
    resolve_chromosome_unified(sp, chr_in, cache = cache, ref_dir = ref_dir)
  }

  strd <- if (!is.null(row_l$strand) && !is.na(row_l$strand[1]) && nzchar(row_l$strand[1])) {
    as.character(row_l$strand[1])
  } else {
    "+"
  }
  if (!strd %in% c("+", "-")) stop("strand must be '+' or '-'.")

  p_start <- if (!is.null(row_l$protospacer_start)) as.integer(row_l$protospacer_start[1]) else NA_integer_
  p_end   <- if (!is.null(row_l$protospacer_end))   as.integer(row_l$protospacer_end[1])   else NA_integer_
  pam_s   <- if (!is.null(row_l$pam_start))         as.integer(row_l$pam_start[1])         else NA_integer_
  pam_e   <- if (!is.null(row_l$pam_end))           as.integer(row_l$pam_end[1])           else NA_integer_

  if (identical(strd, "+")) {
    if (is.na(p_start) && !is.na(pam_s)) p_start <- pam_s - 20L
    if (is.na(p_start) && !is.null(row_l$cut_site)) p_start <- as.integer(row_l$cut_site[1]) - 17L
    if (is.na(p_start)) stop("Must provide protospacer_start, pam_start, or cut_site.")
    if (is.na(p_end))   p_end <- p_start + 19L
    if (is.na(pam_s))   pam_s <- p_end + 1L
    if (is.na(pam_e))   pam_e <- pam_s + 2L
  } else {
    if (is.na(p_start) && !is.na(pam_e)) p_start <- pam_e + 1L
    if (is.na(p_start) && !is.na(pam_s)) p_start <- pam_s + 3L
    if (is.na(p_start) && !is.null(row_l$cut_site)) p_start <- as.integer(row_l$cut_site[1]) - 2L
    if (is.na(p_start)) stop("Must provide protospacer_start, pam_start, pam_end, or cut_site.")
    if (is.na(p_end))   p_end <- p_start + 19L
    if (is.na(pam_e))   pam_e <- p_start - 1L
    if (is.na(pam_s))   pam_s <- pam_e - 2L
  }

  sp_seq <- if (!is.null(row_l$spacer_seq) && !is.na(row_l$spacer_seq[1]) && nzchar(row_l$spacer_seq[1])) {
    toupper(as.character(row_l$spacer_seq[1]))
  } else {
    extract_genomic_subseq(sp, acc, p_start, p_end, strand = strd, cache = cache, ref_dir = ref_dir)
  }

  pm_seq <- if (!is.null(row_l$pam_seq) && !is.na(row_l$pam_seq[1]) && nzchar(row_l$pam_seq[1])) {
    toupper(as.character(row_l$pam_seq[1]))
  } else {
    extract_genomic_subseq(sp, acc, pam_s, pam_e, strand = strd, cache = cache, ref_dir = ref_dir)
  }

  gid <- if (!is.null(row_l$guide_id) && !is.na(row_l$guide_id[1]) && nzchar(row_l$guide_id[1])) {
    as.character(row_l$guide_id[1])
  } else {
    paste0(chr_in, ":", p_start, "-", p_end, "(", strd, ")")
  }

  gsym <- if (!is.null(row_l$gene_symbol) && !is.na(row_l$gene_symbol[1]) && nzchar(row_l$gene_symbol[1])) {
    as.character(row_l$gene_symbol[1])
  } else if (!is.null(attr_gm$gene_symbol) && nzchar(attr_gm$gene_symbol[1])) {
    as.character(attr_gm$gene_symbol[1])
  } else {
    "target_locus"
  }

  list(
    species           = sp,
    chr               = chr_in,
    accession         = acc,
    strand            = strd,
    protospacer_start = as.integer(p_start),
    protospacer_end   = as.integer(p_end),
    pam_start         = as.integer(pam_s),
    pam_end           = as.integer(pam_e),
    spacer_seq        = sp_seq,
    pam_seq           = pm_seq,
    guide_id          = gid,
    gene_symbol       = gsym,
    raw_row           = row_l,
    attr_gene_model   = attr_gm,
    attr_metadata     = attr_meta
  )
}

design_homology_arms <- function(guide_row_or_params,
                                 species = "melanogaster",
                                 arm_5p_len = 1000L,
                                 arm_3p_len = 1000L,
                                 arm_3p_anchor = c("cut_site", "pam_start"),
                                 donor_type = c("cargo_insertion", "scarless_point_edit"),
                                 gene_model = NULL,
                                 cache = NULL,
                                 ref_dir = get_reference_dir()) {
  arm_3p_anchor <- match.arg(arm_3p_anchor)
  donor_type    <- match.arg(donor_type)

  arm_5p_len <- as.integer(round(as.numeric(arm_5p_len[1])))
  arm_3p_len <- as.integer(round(as.numeric(arm_3p_len[1])))
  if (is.na(arm_5p_len) || arm_5p_len < 1L) stop("arm_5p_len must be a positive integer.")
  if (is.na(arm_3p_len) || arm_3p_len < 1L) stop("arm_3p_len must be a positive integer.")

  norm <- .normalize_guide_params(
    guide_row_or_params,
    species = species,
    cache   = cache,
    ref_dir = ref_dir
  )
  if (is.null(gene_model) && !is.null(norm$attr_gene_model)) {
    gene_model <- norm$attr_gene_model
  }

  sp      <- norm$species
  chr     <- norm$chr
  acc     <- norm$accession
  strd    <- norm$strand
  p_start <- norm$protospacer_start
  p_end   <- norm$protospacer_end
  pam_s   <- norm$pam_start
  pam_e   <- norm$pam_end
  sp_seq  <- norm$spacer_seq
  pm_seq  <- norm$pam_seq

  chr_len <- get_chr_length(sp, acc, cache = cache, ref_dir = ref_dir)
  if (is.na(chr_len) || chr_len <= 0L) {
    stop("Could not resolve chromosome length for ", chr, " (", acc, ") in species ", sp)
  }

  if (identical(strd, "+")) {
    cut_5p_genomic <- as.integer(pam_s - 4L)
    cut_3p_genomic <- as.integer(pam_s - 3L)
    cut_site       <- cut_3p_genomic
    cut_left       <- cut_5p_genomic
    cut_right      <- cut_3p_genomic

    arm_5p_start     <- as.integer(max(1L, cut_5p_genomic - arm_5p_len + 1L))
    arm_5p_end       <- as.integer(cut_5p_genomic)
    arm_3p_cut_start <- as.integer(cut_3p_genomic)
    arm_3p_cut_end   <- as.integer(min(chr_len, cut_3p_genomic + arm_3p_len - 1L))
    arm_3p_pam_start <- as.integer(pam_s)
    arm_3p_pam_end   <- as.integer(min(chr_len, pam_s + arm_3p_len - 1L))

    if (identical(arm_3p_anchor, "cut_site")) {
      arm_3p_start <- arm_3p_cut_start
      arm_3p_end   <- arm_3p_cut_end
      pam_offset_in_arm_3p <- 4L
    } else {
      arm_3p_start <- arm_3p_pam_start
      arm_3p_end   <- arm_3p_pam_end
      pam_offset_in_arm_3p <- 1L
    }

    arm_5p_seq     <- extract_genomic_subseq(sp, acc, arm_5p_start, arm_5p_end, strand = "+", cache = cache, ref_dir = ref_dir)
    arm_3p_cut_seq <- extract_genomic_subseq(sp, acc, arm_3p_cut_start, arm_3p_cut_end, strand = "+", cache = cache, ref_dir = ref_dir)
    arm_3p_pam_seq <- extract_genomic_subseq(sp, acc, arm_3p_pam_start, arm_3p_pam_end, strand = "+", cache = cache, ref_dir = ref_dir)
    arm_3p_seq     <- if (identical(arm_3p_anchor, "cut_site")) arm_3p_cut_seq else arm_3p_pam_seq

    arm_5p_seq_plus <- arm_5p_seq
    arm_3p_seq_plus <- arm_3p_seq

    genomic_left_arm_start     <- arm_5p_start
    genomic_left_arm_end       <- arm_5p_end
    genomic_left_arm_seq_plus  <- arm_5p_seq
    genomic_right_arm_start    <- arm_3p_start
    genomic_right_arm_end      <- arm_3p_end
    genomic_right_arm_seq_plus <- arm_3p_seq
  } else {
    cut_5p_genomic <- as.integer(pam_e + 4L)
    cut_3p_genomic <- as.integer(pam_e + 3L)
    cut_site       <- cut_3p_genomic
    cut_left       <- cut_3p_genomic
    cut_right      <- cut_5p_genomic

    arm_5p_start     <- as.integer(cut_5p_genomic)
    arm_5p_end       <- as.integer(min(chr_len, cut_5p_genomic + arm_5p_len - 1L))
    arm_3p_cut_start <- as.integer(max(1L, cut_3p_genomic - arm_3p_len + 1L))
    arm_3p_cut_end   <- as.integer(cut_3p_genomic)
    arm_3p_pam_start <- as.integer(max(1L, pam_e - arm_3p_len + 1L))
    arm_3p_pam_end   <- as.integer(pam_e)

    if (identical(arm_3p_anchor, "cut_site")) {
      arm_3p_start <- arm_3p_cut_start
      arm_3p_end   <- arm_3p_cut_end
      pam_offset_in_arm_3p <- 4L
    } else {
      arm_3p_start <- arm_3p_pam_start
      arm_3p_end   <- arm_3p_pam_end
      pam_offset_in_arm_3p <- 1L
    }

    arm_5p_seq     <- extract_genomic_subseq(sp, acc, arm_5p_start, arm_5p_end, strand = "-", cache = cache, ref_dir = ref_dir)
    arm_3p_cut_seq <- extract_genomic_subseq(sp, acc, arm_3p_cut_start, arm_3p_cut_end, strand = "-", cache = cache, ref_dir = ref_dir)
    arm_3p_pam_seq <- extract_genomic_subseq(sp, acc, arm_3p_pam_start, arm_3p_pam_end, strand = "-", cache = cache, ref_dir = ref_dir)
    arm_3p_seq     <- if (identical(arm_3p_anchor, "cut_site")) arm_3p_cut_seq else arm_3p_pam_seq

    arm_5p_seq_plus <- revcomp_dna(arm_5p_seq)
    arm_3p_seq_plus <- revcomp_dna(arm_3p_seq)

    genomic_left_arm_start     <- arm_3p_start
    genomic_left_arm_end       <- arm_3p_end
    genomic_left_arm_seq_plus  <- arm_3p_seq_plus
    genomic_right_arm_start    <- arm_5p_start
    genomic_right_arm_end      <- arm_5p_end
    genomic_right_arm_seq_plus <- arm_5p_seq_plus
  }

  qc_5p <- .analyze_single_arm_qc(arm_5p_seq, arm_5p_len)
  qc_3p <- .analyze_single_arm_qc(arm_3p_seq, arm_3p_len)

  min_coord <- min(arm_5p_start, arm_3p_start)
  max_coord <- max(arm_5p_end, arm_3p_end)
  boundary_dist <- as.integer(min(min_coord - 1L, chr_len - max_coord))

  internal_5p_match <- .has_spacer_pam_match(arm_5p_seq, sp_seq)
  internal_3p_match <- .has_spacer_pam_match(arm_3p_seq, sp_seq)
  internal_recleave <- internal_5p_match || internal_3p_match
  junction_recleave <- .has_spacer_pam_match(paste0(arm_5p_seq, arm_3p_seq), sp_seq)
  cargo_disrupts    <- identical(arm_3p_anchor, "cut_site") && !internal_recleave
  recleavage_risk   <- internal_recleave || (identical(donor_type, "scarless_point_edit") && junction_recleave)

  shield <- .design_pam_shielding(
    spacer_seq        = sp_seq,
    pam_seq           = pm_seq,
    strand            = strd,
    protospacer_start = p_start,
    protospacer_end   = p_end,
    pam_start         = pam_s,
    pam_end           = pam_e,
    arm_3p_seq        = arm_3p_seq,
    arm_3p_anchor     = arm_3p_anchor,
    gene_model        = gene_model,
    species           = sp,
    acc               = acc,
    cache             = cache,
    ref_dir           = ref_dir
  )

  n_count_total   <- as.integer(qc_5p$n_count + qc_3p$n_count)
  max_n_run_total <- as.integer(max(qc_5p$max_n_run, qc_3p$max_n_run))
  has_n_tract     <- qc_5p$has_n_tract || qc_3p$has_n_tract
  max_hp_total    <- as.integer(max(qc_5p$max_homopolymer, qc_3p$max_homopolymer))
  low_comp_pct    <- round(max(qc_5p$low_complexity_pct, qc_3p$low_complexity_pct), 2)
  min_ent_total   <- round(min(qc_5p$min_window_entropy, qc_3p$min_window_entropy), 3)
  repeat_warn     <- qc_5p$repeat_warning || qc_3p$repeat_warning
  gc_warn         <- qc_5p$gc_warning || qc_3p$gc_warning
  trunc_any       <- qc_5p$boundary_truncated || qc_3p$boundary_truncated

  warnings_vec <- character(0)
  if (has_n_tract) {
    warnings_vec <- c(warnings_vec, sprintf("N-tract detected (total N=%d, max run=%d bp).", n_count_total, max_n_run_total))
  }
  if (repeat_warn) {
    warnings_vec <- c(warnings_vec, sprintf("Low-complexity/repeat region detected (max homopolymer=%d bp, min 50bp entropy=%.2f bits).", max_hp_total, min_ent_total))
  }
  if (gc_warn) {
    warnings_vec <- c(warnings_vec, sprintf("Extreme GC content (5' arm=%.1f%%, 3' arm=%.1f%%).", qc_5p$gc_pct, qc_3p$gc_pct))
  }
  if (trunc_any) {
    warnings_vec <- c(warnings_vec, sprintf("Chromosome scaffold boundary truncation (5'=%d/%d bp, 3'=%d/%d bp).",
                                            nchar(arm_5p_seq), arm_5p_len, nchar(arm_3p_seq), arm_3p_len))
  }
  if (recleavage_risk) {
    warnings_vec <- c(warnings_vec, paste0("Donor re-cleavage risk detected; recommended shield: ", shield$summary))
  }

  out <- list(
    engine_version            = HOMOLOGY_ARMS_VERSION,
    species                   = sp,
    gene_symbol               = norm$gene_symbol,
    guide_id                  = norm$guide_id,
    chr                       = chr,
    accession                 = acc,
    chr_len                   = as.integer(chr_len),
    strand                    = strd,
    protospacer_start         = p_start,
    protospacer_end           = p_end,
    pam_start                 = pam_s,
    pam_end                   = pam_e,
    spacer_seq                = sp_seq,
    pam_seq                   = pm_seq,
    cut_site                  = cut_site,
    cut_5p_genomic            = cut_5p_genomic,
    cut_3p_genomic            = cut_3p_genomic,
    cut_left                  = cut_left,
    cut_right                 = cut_right,
    arm_3p_anchor             = arm_3p_anchor,
    donor_type                = donor_type,
    arm_5p_len                = as.integer(nchar(arm_5p_seq)),
    arm_3p_len                = as.integer(nchar(arm_3p_seq)),
    arm_5p_requested_len      = arm_5p_len,
    arm_3p_requested_len      = arm_3p_len,
    arm_5p_start              = arm_5p_start,
    arm_5p_end                = arm_5p_end,
    arm_5p_seq                = arm_5p_seq,
    arm_5p_seq_plus           = arm_5p_seq_plus,
    arm_3p_start              = arm_3p_start,
    arm_3p_end                = arm_3p_end,
    arm_3p_seq                = arm_3p_seq,
    arm_3p_seq_plus           = arm_3p_seq_plus,
    arm_3p_cut_start          = arm_3p_cut_start,
    arm_3p_cut_end            = arm_3p_cut_end,
    arm_3p_cut_seq            = arm_3p_cut_seq,
    arm_3p_pam_start          = arm_3p_pam_start,
    arm_3p_pam_end            = arm_3p_pam_end,
    arm_3p_pam_seq            = arm_3p_pam_seq,
    pam_offset_in_arm_3p      = as.integer(pam_offset_in_arm_3p),
    genomic_left_arm_start    = genomic_left_arm_start,
    genomic_left_arm_end      = genomic_left_arm_end,
    genomic_left_arm_seq_plus = genomic_left_arm_seq_plus,
    genomic_right_arm_start   = genomic_right_arm_start,
    genomic_right_arm_end     = genomic_right_arm_end,
    genomic_right_arm_seq_plus = genomic_right_arm_seq_plus,
    arm_5p_gc                 = qc_5p$gc_pct,
    arm_3p_gc                 = qc_3p$gc_pct,
    gc_warning_5p             = qc_5p$gc_warning,
    gc_warning_3p             = qc_3p$gc_warning,
    gc_warning                = gc_warn,
    n_count_5p                = qc_5p$n_count,
    n_count_3p                = qc_3p$n_count,
    n_count                   = n_count_total,
    max_n_run_5p              = qc_5p$max_n_run,
    max_n_run_3p              = qc_3p$max_n_run,
    max_n_run                 = max_n_run_total,
    has_n_tract_5p            = qc_5p$has_n_tract,
    has_n_tract_3p            = qc_3p$has_n_tract,
    has_n_tract               = has_n_tract,
    max_homopolymer_5p        = qc_5p$max_homopolymer,
    max_homopolymer_3p        = qc_3p$max_homopolymer,
    max_homopolymer           = max_hp_total,
    low_complexity_pct_5p     = qc_5p$low_complexity_pct,
    low_complexity_pct_3p     = qc_3p$low_complexity_pct,
    low_complexity_pct        = low_comp_pct,
    min_window_entropy_5p     = qc_5p$min_window_entropy,
    min_window_entropy_3p     = qc_3p$min_window_entropy,
    min_window_entropy        = min_ent_total,
    repeat_warning_5p         = qc_5p$repeat_warning,
    repeat_warning_3p         = qc_3p$repeat_warning,
    repeat_warning            = repeat_warn,
    boundary_truncated_5p     = qc_5p$boundary_truncated,
    boundary_truncated_3p     = qc_3p$boundary_truncated,
    boundary_truncated        = trunc_any,
    boundary_distance_bp      = boundary_dist,
    recleavage_risk           = recleavage_risk,
    scarless_recleavage_risk  = junction_recleave,
    internal_arm_recleavage   = internal_recleave,
    cargo_disrupts_target     = cargo_disrupts,
    pam_shielding_mutation    = shield,
    pam_shielding_summary     = shield$summary,
    shielded_arm_3p_seq       = shield$shielded_arm_3p_seq,
    qc_pass                   = !has_n_tract && !repeat_warn && !gc_warn && !trunc_any,
    qc_warnings               = warnings_vec
  )
  class(out) <- c("cas9_homology_arms", "list")
  out
}

# ============================================================================
# 5. Nucleotide-Level Cut-Site Zoom Helper (build_cut_site_zoom)
# ============================================================================

#' Build a Nucleotide-Level Cut-Site Zoom Window Around the SpCas9 Blunt Cleavage Site
#'
#' @param arm_design Output list from `design_homology_arms()`
#' @param flank_bp Number of nucleotides on each side of the blunt cut site (default `30L`)
#' @return A structured list with `top_strand_5to3`, `bottom_strand_3to5`,
#'   `genomic_positions`, `spacer_range`, `pam_range`, `cut_index`, and `formatted_ascii`
build_cut_site_zoom <- function(arm_design, flank_bp = 30L) {
  if (!is.list(arm_design) || is.null(arm_design$arm_5p_seq)) {
    stop("arm_design must be a list returned by design_homology_arms().")
  }
  flank_bp <- max(10L, as.integer(round(as.numeric(flank_bp[1]))))
  strd     <- arm_design$strand

  # Extract +/- flank_bp on the '+' genomic reference strand around [cut_left, cut_right]
  if (identical(strd, "+")) {
    left_plus_seq  <- arm_design$arm_5p_seq
    right_plus_seq <- arm_design$arm_3p_cut_seq
  } else {
    left_plus_seq  <- revcomp_dna(arm_design$arm_3p_cut_seq)
    right_plus_seq <- revcomp_dna(arm_design$arm_5p_seq)
  }

  l_len <- min(flank_bp, nchar(left_plus_seq))
  r_len <- min(flank_bp, nchar(right_plus_seq))

  top_left  <- if (l_len > 0L) substr(left_plus_seq, nchar(left_plus_seq) - l_len + 1L, nchar(left_plus_seq)) else ""
  top_right <- if (r_len > 0L) substr(right_plus_seq, 1L, r_len) else ""
  top_strand_5to3    <- paste0(top_left, top_right)
  bottom_strand_3to5 <- chartr("ACGTNacgtn", "TGCANtgcan", top_strand_5to3)

  win_start_genomic <- as.integer(arm_design$cut_left - l_len + 1L)
  win_end_genomic   <- as.integer(arm_design$cut_right + r_len - 1L)
  genomic_positions <- seq.int(win_start_genomic, win_end_genomic)
  cut_index         <- as.integer(l_len)

  sp_idx_start  <- match(arm_design$protospacer_start, genomic_positions)
  sp_idx_end    <- match(arm_design$protospacer_end,   genomic_positions)
  pam_idx_start <- match(arm_design$pam_start,         genomic_positions)
  pam_idx_end   <- match(arm_design$pam_end,           genomic_positions)

  spacer_range <- c(start_idx = as.integer(sp_idx_start), end_idx = as.integer(sp_idx_end))
  pam_range    <- c(start_idx = as.integer(pam_idx_start), end_idx = as.integer(pam_idx_end))

  # Target-strand 5'->3' junction sequences
  t5_len <- min(flank_bp, nchar(arm_design$arm_5p_seq))
  t3_len <- min(flank_bp, nchar(arm_design$arm_3p_cut_seq))
  t_left  <- if (t5_len > 0L) substr(arm_design$arm_5p_seq, nchar(arm_design$arm_5p_seq) - t5_len + 1L, nchar(arm_design$arm_5p_seq)) else ""
  t_right <- if (t3_len > 0L) substr(arm_design$arm_3p_cut_seq, 1L, t3_len) else ""
  target_strand_5to3    <- paste0(t_left, t_right)
  non_target_strand_3to5 <- chartr("ACGTNacgtn", "TGCANtgcan", target_strand_5to3)

  # Build annotation track above the split sequence (with '|' inserted after cut_index)
  n_tot <- nchar(top_strand_5to3)
  annot_chars <- rep(" ", n_tot)
  if (!is.na(sp_idx_start) && !is.na(sp_idx_end)) {
    annot_chars[sp_idx_start:sp_idx_end] <- "-"
    annot_chars[sp_idx_start] <- if (identical(strd, "+")) "5" else "3"
    annot_chars[sp_idx_end]   <- if (identical(strd, "+")) ">" else "<"
  }
  if (!is.na(pam_idx_start) && !is.na(pam_idx_end)) {
    annot_chars[pam_idx_start:pam_idx_end] <- c("P", "A", "M")
  }
  annot_left  <- if (cut_index > 0L) paste(annot_chars[seq_len(cut_index)], collapse = "") else ""
  annot_right <- if (cut_index < n_tot) paste(annot_chars[(cut_index + 1L):n_tot], collapse = "") else ""
  annot_line  <- paste0("  Track   : ", annot_left, "|", annot_right)

  bot_left  <- if (l_len > 0L) substr(bottom_strand_3to5, 1L, l_len) else ""
  bot_right <- if (r_len > 0L) substr(bottom_strand_3to5, l_len + 1L, n_tot) else ""

  arm_label_left  <- if (identical(strd, "+")) "5' Homology Arm (arm_5p)" else "3' Homology Arm (arm_3p, '-' strand)"
  arm_label_right <- if (identical(strd, "+")) "3' Homology Arm (arm_3p)" else "5' Homology Arm (arm_5p, '-' strand)"

  ascii_lines <- c(
    sprintf("SpCas9 Blunt Cut-Site Zoom (%s %s:%d|%d, target strand '%s', cut 3 bp 5' of PAM):",
            arm_design$gene_symbol, arm_design$chr, arm_design$cut_left, arm_design$cut_right, strd),
    sprintf("  Arms    : [<- %s ... %d] | [%d ... %s ->]",
            arm_label_left, arm_design$cut_left, arm_design$cut_right, arm_label_right),
    annot_line,
    sprintf("  5'->3' +: %s|%s   (%d..%d)", top_left, top_right, win_start_genomic, win_end_genomic),
    sprintf("  3'->5' -: %s|%s", bot_left, bot_right),
    sprintf("  Target  : 5'-%s|%s-3' (spacer=%s, PAM=%s)",
            t_left, t_right, arm_design$spacer_seq, arm_design$pam_seq)
  )

  list(
    flank_bp               = flank_bp,
    strand                 = strd,
    top_strand_5to3        = top_strand_5to3,
    bottom_strand_3to5     = bottom_strand_3to5,
    target_strand_5to3     = target_strand_5to3,
    non_target_strand_3to5 = non_target_strand_3to5,
    genomic_positions      = genomic_positions,
    spacer_range           = spacer_range,
    pam_range              = pam_range,
    cut_index              = cut_index,
    cut_left_genomic       = arm_design$cut_left,
    cut_right_genomic      = arm_design$cut_right,
    arm_5p_junction_seq    = t_left,
    arm_3p_junction_seq    = t_right,
    formatted_ascii        = paste(ascii_lines, collapse = "\n")
  )
}

# ============================================================================
# 6. Export Functions (FASTA, CSV, JSON, Design Log, and ZIP Bundle)
# ============================================================================

.wrap_fasta_seq <- function(seq_str, width = 80L) {
  s <- as.character(seq_str[1])
  n <- nchar(s)
  if (n == 0L) return("")
  starts <- seq.int(1L, n, by = as.integer(width))
  ends   <- pmin(n, starts + as.integer(width) - 1L)
  substring(s, starts, ends)
}

export_homology_arms_fasta <- function(arm_design, file = NULL, include_plus_strand = TRUE) {
  if (!is.list(arm_design) || is.null(arm_design$arm_5p_seq)) {
    stop("arm_design must be a list returned by design_homology_arms().")
  }

  hdr_5p <- sprintf(
    ">%s|%s|arm_5p|species=%s|chr=%s|accession=%s|coords=%d-%d|strand=%s|len=%dbp|gc=%.2f%%|cut_site=%d",
    arm_design$gene_symbol, arm_design$guide_id, arm_design$species,
    arm_design$chr, arm_design$accession, arm_design$arm_5p_start, arm_design$arm_5p_end,
    arm_design$strand, arm_design$arm_5p_len, arm_design$arm_5p_gc, arm_design$cut_site
  )

  hdr_3p <- sprintf(
    ">%s|%s|arm_3p|species=%s|chr=%s|accession=%s|coords=%d-%d|strand=%s|anchor=%s|len=%dbp|gc=%.2f%%|cut_site=%d",
    arm_design$gene_symbol, arm_design$guide_id, arm_design$species,
    arm_design$chr, arm_design$accession, arm_design$arm_3p_start, arm_design$arm_3p_end,
    arm_design$strand, arm_design$arm_3p_anchor, arm_design$arm_3p_len,
    arm_design$arm_3p_gc, arm_design$cut_site
  )

  lines_out <- c(
    hdr_5p, .wrap_fasta_seq(arm_design$arm_5p_seq),
    hdr_3p, .wrap_fasta_seq(arm_design$arm_3p_seq)
  )

  if (!is.null(arm_design$shielded_arm_3p_seq) && nzchar(arm_design$shielded_arm_3p_seq)) {
    hdr_shield <- sprintf(
      ">%s|%s|arm_3p_pam_shielded|species=%s|chr=%s|coords=%d-%d|strand=%s|shield=%s",
      arm_design$gene_symbol, arm_design$guide_id, arm_design$species,
      arm_design$chr, arm_design$arm_3p_start, arm_design$arm_3p_end,
      arm_design$strand, gsub("\\s+", "_", arm_design$pam_shielding_summary)
    )
    lines_out <- c(lines_out, hdr_shield, .wrap_fasta_seq(arm_design$shielded_arm_3p_seq))
  }

  if (isTRUE(include_plus_strand) && identical(arm_design$strand, "-")) {
    hdr_left_plus <- sprintf(
      ">%s|%s|genomic_left_arm_plus_strand|chr=%s|coords=%d-%d|strand=+|len=%dbp",
      arm_design$gene_symbol, arm_design$guide_id, arm_design$chr,
      arm_design$genomic_left_arm_start, arm_design$genomic_left_arm_end,
      nchar(arm_design$genomic_left_arm_seq_plus)
    )
    hdr_right_plus <- sprintf(
      ">%s|%s|genomic_right_arm_plus_strand|chr=%s|coords=%d-%d|strand=+|len=%dbp",
      arm_design$gene_symbol, arm_design$guide_id, arm_design$chr,
      arm_design$genomic_right_arm_start, arm_design$genomic_right_arm_end,
      nchar(arm_design$genomic_right_arm_seq_plus)
    )
    lines_out <- c(
      lines_out,
      hdr_left_plus,  .wrap_fasta_seq(arm_design$genomic_left_arm_seq_plus),
      hdr_right_plus, .wrap_fasta_seq(arm_design$genomic_right_arm_seq_plus)
    )
  }

  fasta_txt <- paste(lines_out, collapse = "\n")
  if (!is.null(file) && nzchar(file[1])) {
    writeLines(lines_out, con = file[1])
  }
  fasta_txt
}

export_design_csv <- function(grna_df, arm_design = NULL, file = NULL) {
  if (!is.data.frame(grna_df)) {
    stop("grna_df must be a data.frame.")
  }
  df <- grna_df

  if (!is.null(arm_design) && is.list(arm_design) && nrow(df) > 0L) {
    match_idx <- rep(FALSE, nrow(df))
    if ("guide_id" %in% names(df) && !is.null(arm_design$guide_id)) {
      match_idx <- df$guide_id == arm_design$guide_id
    }
    if (!any(match_idx) && all(c("protospacer_start", "strand") %in% names(df))) {
      match_idx <- (df$protospacer_start == arm_design$protospacer_start) &
                   (df$strand == arm_design$strand)
    }
    if (!any(match_idx)) match_idx[1L] <- TRUE

    df$selected_for_hdr      <- match_idx
    df$arm_5p_start          <- ifelse(match_idx, arm_design$arm_5p_start, NA_integer_)
    df$arm_5p_end            <- ifelse(match_idx, arm_design$arm_5p_end,   NA_integer_)
    df$arm_5p_len            <- ifelse(match_idx, arm_design$arm_5p_len,   NA_integer_)
    df$arm_5p_gc             <- ifelse(match_idx, arm_design$arm_5p_gc,    NA_real_)
    df$arm_3p_start          <- ifelse(match_idx, arm_design$arm_3p_start, NA_integer_)
    df$arm_3p_end            <- ifelse(match_idx, arm_design$arm_3p_end,   NA_integer_)
    df$arm_3p_len            <- ifelse(match_idx, arm_design$arm_3p_len,   NA_integer_)
    df$arm_3p_gc             <- ifelse(match_idx, arm_design$arm_3p_gc,    NA_real_)
    df$pam_shielding_summary <- ifelse(match_idx, arm_design$pam_shielding_summary, NA_character_)
    df$arm_5p_seq            <- ifelse(match_idx, arm_design$arm_5p_seq,   NA_character_)
    df$arm_3p_seq            <- ifelse(match_idx, arm_design$arm_3p_seq,   NA_character_)
  }

  tc <- textConnection("csv_buf", open = "w", local = TRUE)
  utils::write.csv(df, file = tc, row.names = FALSE)
  close(tc)
  csv_txt <- paste(csv_buf, collapse = "\n")

  if (!is.null(file) && nzchar(file[1])) {
    writeLines(csv_buf, con = file[1])
  }
  attr(csv_txt, "data") <- df
  csv_txt
}

.escape_json_str <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n",  s, fixed = TRUE)
  s <- gsub("\r", "\\r",  s, fixed = TRUE)
  s <- gsub("\t", "\\t",  s, fixed = TRUE)
  paste0("\"", s, "\"")
}

.to_json_str <- function(x, indent = 0L) {
  pad     <- paste(rep("  ", indent), collapse = "")
  pad_in  <- paste(rep("  ", indent + 1L), collapse = "")

  if (is.null(x)) return("null")

  if (is.data.frame(x)) {
    if (nrow(x) == 0L) return("[]")
    rows <- vapply(seq_len(nrow(x)), function(i) {
      .to_json_str(as.list(x[i, , drop = FALSE]), indent = indent + 1L)
    }, character(1))
    return(paste0("[\n", paste0(pad_in, rows, collapse = ",\n"), "\n", pad, "]"))
  }

  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      if (length(x) == 0L) return("[]")
      elems <- vapply(x, .to_json_str, character(1), indent = indent + 1L)
      return(paste0("[\n", paste0(pad_in, elems, collapse = ",\n"), "\n", pad, "]"))
    }
    if (length(x) == 0L) return("{}")
    pairs <- vapply(seq_along(x), function(i) {
      paste0(pad_in, .escape_json_str(nms[i]), ": ", .to_json_str(x[[i]], indent = indent + 1L))
    }, character(1))
    return(paste0("{\n", paste(pairs, collapse = ",\n"), "\n", pad, "}"))
  }

  if (is.logical(x)) {
    if (length(x) == 0L) return("[]")
    vals <- ifelse(is.na(x), "null", ifelse(x, "true", "false"))
    return(if (length(vals) == 1L) vals else paste0("[", paste(vals, collapse = ", "), "]"))
  }

  if (is.numeric(x)) {
    if (length(x) == 0L) return("[]")
    vals <- ifelse(is.na(x) | !is.finite(x), "null", as.character(x))
    return(if (length(vals) == 1L) vals else paste0("[", paste(vals, collapse = ", "), "]"))
  }

  ch <- as.character(x)
  if (length(ch) == 0L) return("[]")
  vals <- ifelse(is.na(ch), "null", vapply(ch, .escape_json_str, character(1), USE.NAMES = FALSE))
  if (length(vals) == 1L) vals else paste0("[", paste(vals, collapse = ", "), "]")
}

export_design_json <- function(grna_df, arm_design = NULL, metadata = NULL, file = NULL) {
  if (is.null(metadata) && is.data.frame(grna_df)) {
    metadata <- attr(grna_df, "metadata")
  }
  if (is.null(metadata)) {
    metadata <- list(
      engine_version        = if (exists("GRNA_ENGINE_VERSION")) GRNA_ENGINE_VERSION else "1.0.0",
      homology_arms_version = HOMOLOGY_ARMS_VERSION,
      timestamp_utc         = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
  } else {
    metadata$homology_arms_version <- HOMOLOGY_ARMS_VERSION
  }

  zoom <- if (!is.null(arm_design) && is.list(arm_design)) {
    tryCatch(build_cut_site_zoom(arm_design, flank_bp = 30L), error = function(e) NULL)
  } else NULL

  payload <- list(
    metadata       = metadata,
    homology_arms  = if (!is.null(arm_design)) unclass(arm_design) else NULL,
    cut_site_zoom  = zoom,
    candidates     = if (is.data.frame(grna_df)) grna_df else NULL
  )

  json_txt <- .to_json_str(payload, indent = 0L)
  if (!is.null(file) && nzchar(file[1])) {
    writeLines(json_txt, con = file[1])
  }
  json_txt
}

export_design_log <- function(grna_df, arm_design = NULL, metadata = NULL, file = NULL) {
  if (is.null(metadata) && is.data.frame(grna_df)) {
    metadata <- attr(grna_df, "metadata")
  }
  ts <- if (!is.null(metadata$timestamp_utc)) metadata$timestamp_utc else format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  sp <- if (!is.null(arm_design$species)) arm_design$species else if (!is.null(metadata$species)) metadata$species else "melanogaster"
  asm <- if (!is.null(metadata$genome_assembly)) metadata$genome_assembly else unname(SPECIES_GCF[canonical_species(sp)])
  n_cands <- if (is.data.frame(grna_df)) nrow(grna_df) else 0L

  lines <- c(
    "============================================================================",
    "Cas9 sgRNA & HDR Homology Arm Reproducibility Design Log",
    "============================================================================",
    sprintf("Timestamp (UTC)        : %s", ts),
    sprintf("Homology Arms Version  : %s", HOMOLOGY_ARMS_VERSION),
    sprintf("gRNA Engine Version    : %s", if (exists("GRNA_ENGINE_VERSION")) GRNA_ENGINE_VERSION else "1.0.0"),
    sprintf("Species                : %s (%s)", sp, ifelse(is.na(asm), "custom", asm)),
    sprintf("Total Candidate gRNAs  : %d", n_cands)
  )

  if (!is.null(metadata$scoring_models)) {
    lines <- c(
      lines,
      "----------------------------------------------------------------------------",
      "Scoring Models:",
      sprintf("  On-Target Primary    : %s", metadata$scoring_models$on_target_primary),
      sprintf("  On-Target Secondary  : %s", metadata$scoring_models$on_target_secondary),
      sprintf("  Off-Target CFD       : %s", metadata$scoring_models$off_target_cfd),
      sprintf("  Off-Target MIT       : %s", metadata$scoring_models$off_target_mit)
    )
  }

  if (!is.null(arm_design) && is.list(arm_design)) {
    zoom <- tryCatch(build_cut_site_zoom(arm_design, flank_bp = 30L), error = function(e) NULL)
    lines <- c(
      lines,
      "----------------------------------------------------------------------------",
      "Selected Guide & SpCas9 Blunt Cut-Site Geometry (3 bp 5' of PAM):",
      sprintf("  Gene / Locus         : %s", arm_design$gene_symbol),
      sprintf("  Guide ID             : %s", arm_design$guide_id),
      sprintf("  Chromosome           : %s (%s, len=%d bp)", arm_design$chr, arm_design$accession, arm_design$chr_len),
      sprintf("  Target Strand        : %s", arm_design$strand),
      sprintf("  Protospacer (20 bp)  : %d..%d (%s)", arm_design$protospacer_start, arm_design$protospacer_end, arm_design$spacer_seq),
      sprintf("  PAM (3 bp)           : %d..%d (%s)", arm_design$pam_start, arm_design$pam_end, arm_design$pam_seq),
      sprintf("  Blunt Cut Site       : between %d (5' side) and %d (3' side) [cut_site=%d]",
              arm_design$cut_5p_genomic, arm_design$cut_3p_genomic, arm_design$cut_site),
      "----------------------------------------------------------------------------",
      "Homology Arms & Sequence Quality Control:",
      sprintf("  5' Arm (arm_5p)      : %s:%d..%d (strand=%s, len=%d bp, GC=%.2f%%)",
              arm_design$chr, arm_design$arm_5p_start, arm_design$arm_5p_end,
              arm_design$strand, arm_design$arm_5p_len, arm_design$arm_5p_gc),
      sprintf("  3' Arm (arm_3p)      : %s:%d..%d (strand=%s, anchor=%s, len=%d bp, GC=%.2f%%)",
              arm_design$chr, arm_design$arm_3p_start, arm_design$arm_3p_end,
              arm_design$strand, arm_design$arm_3p_anchor, arm_design$arm_3p_len, arm_design$arm_3p_gc),
      sprintf("  PAM-Anchored 3' Arm  : %s:%d..%d (pam_offset_in_arm_3p=%d)",
              arm_design$chr, arm_design$arm_3p_pam_start, arm_design$arm_3p_pam_end, arm_design$pam_offset_in_arm_3p),
      sprintf("  QC Pass              : %s", ifelse(arm_design$qc_pass, "PASS", "WARN")),
      sprintf("  N-Tract Check        : total N=%d, max N run=%d bp (has_n_tract=%s)",
              arm_design$n_count, arm_design$max_n_run, arm_design$has_n_tract),
      sprintf("  Repeat/Entropy Check : max homopolymer=%d bp, min 50bp entropy=%.2f bits, low_comp=%.1f%% (repeat_warning=%s)",
              arm_design$max_homopolymer, arm_design$min_window_entropy, arm_design$low_complexity_pct, arm_design$repeat_warning),
      sprintf("  Boundary Check       : truncated_5p=%s, truncated_3p=%s, min_dist=%d bp",
              arm_design$boundary_truncated_5p, arm_design$boundary_truncated_3p, arm_design$boundary_distance_bp),
      sprintf("  Re-Cleavage & Shield : recleavage_risk=%s | %s",
              arm_design$recleavage_risk, arm_design$pam_shielding_summary)
    )
    if (!is.null(zoom)) {
      lines <- c(
        lines,
        "----------------------------------------------------------------------------",
        zoom$formatted_ascii
      )
    }
  }

  lines <- c(lines, "============================================================================")
  log_txt <- paste(lines, collapse = "\n")
  if (!is.null(file) && nzchar(file[1])) {
    writeLines(lines, con = file[1])
  }
  log_txt
}

export_design_bundle_zip <- function(zip_path,
                                     grna_df,
                                     arm_design = NULL,
                                     gene_model = NULL,
                                     metadata = NULL,
                                     plot_fn = NULL) {
  if (is.null(zip_path) || !nzchar(zip_path[1])) {
    stop("zip_path must be a non-empty file path.")
  }
  if (is.null(arm_design) && is.data.frame(grna_df) && nrow(grna_df) > 0L) {
    arm_design <- tryCatch(
      design_homology_arms(grna_df[1L, , drop = FALSE], gene_model = gene_model),
      error = function(e) NULL
    )
  }

  stage_dir <- tempfile("cas9_hdr_bundle_")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)

  csv_file   <- file.path(stage_dir, "grna_candidates.csv")
  fasta_file <- file.path(stage_dir, "homology_arms.fasta")
  json_file  <- file.path(stage_dir, "design_bundle.json")
  log_file   <- file.path(stage_dir, "design_reproducibility.log")

  export_design_csv(grna_df, arm_design = arm_design, file = csv_file)
  if (!is.null(arm_design)) {
    export_homology_arms_fasta(arm_design, file = fasta_file, include_plus_strand = TRUE)
  }
  export_design_json(grna_df, arm_design = arm_design, metadata = metadata, file = json_file)
  export_design_log(grna_df, arm_design = arm_design, metadata = metadata, file = log_file)

  if (is.function(plot_fn)) {
    plot_file <- file.path(stage_dir, "locus_design_plot.pdf")
    tryCatch(plot_fn(plot_file), error = function(e) NULL)
  }

  staged_files <- list.files(stage_dir, full.names = FALSE)
  out_dir <- dirname(normalizePath(zip_path[1], mustWork = FALSE))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  target_zip <- normalizePath(zip_path[1], mustWork = FALSE)
  if (file.exists(target_zip)) unlink(target_zip)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(stage_dir)
  suppressWarnings(utils::zip(zipfile = target_zip, files = staged_files, flags = "-q"))
  setwd(old_wd)

  if (!file.exists(target_zip) || is.na(file.info(target_zip)$size) || file.info(target_zip)$size == 0) {
    stop("Failed to create ZIP archive at: ", target_zip)
  }

  attr(target_zip, "files") <- staged_files
  target_zip
}







