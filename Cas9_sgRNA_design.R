#!/usr/bin/env Rscript
# ============================================================================
# CRISPR Cas9 gRNA & Homology Arm Designer for Drosophila (v1.0.0)
# ============================================================================
# Standalone Shiny GUI + Headless CLI application modeled after the lab's
# HCR Probe Designer (HCR_probe_design_v41.R).
#
# Core Capabilities:
#   1. Multi-species Drosophila reference genome registry (7 species:
#      d_melanogaster, d_simulans, d_eugracilis, d_erecta, d_teissieri,
#      d_ananassae, d_virilis) with automatic detection and one-click NCBI
#      RefSeq lazy-download / decompression / .fai indexing.
#   2. Multi-modal gene lookup (Symbol, FlyBase FBgn ID, CG number, NCBI LOC ID,
#      synonym/alias, protein product name, or genomic coordinates chr:start-end)
#      with <1s GTF indexing and isoform exon/CDS/5'UTR/promoter extraction.
#   3. Reproducible SpCas9 (20-nt protospacer + NGG PAM) discovery and ranking:
#      - Genomic context tiering (Tier 1: 5'-most shared coding exons;
#        Tier 2: 5' UTR / upstream promoter; Tier 3: late shared CDS;
#        Tier 4: isoform-specific CDS; Tier 5: 3' UTR / intron).
#      - On-target efficacy scoring: Doench 2016 Rule Set 2 / Azimuth +
#        Moreno-Mateos 2015 CRISPRscan.
#      - Off-target specificity scoring: Doench 2016 CFD + Hsu 2013 MIT.
#   4. Exact SpCas9 blunt cut-site geometry (3 bp upstream / 5' of the first
#      nucleotide of the PAM on the target strand) and adjustable 1-2 kb
#      flanking 5' and 3' homology arms on both '+' and '-' strands, with
#      donor re-cleavage / PAM-shielding check, N-tract / repeat / GC QC,
#      and nucleotide-level cut-site duplex zoom view.
#   5. Isoform + Cut-Site & Homology-Arm track schematic plot modeled on
#      plot_isoform_tracks() from HCR_probe_design_v41.R.
#   6. Reproducible one-click export (ZIP bundle, Homology Arms FASTA,
#      Ranked gRNAs CSV, Design Summary JSON, and Reproducibility Log TXT).
# ============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("shiny", quietly = TRUE)) library(shiny)
})

CAS9_APP_VERSION <- "1.0.0"
PORT <- as.integer(Sys.getenv("CAS9_PORT", unset = "7789"))
if (is.na(PORT) || PORT <= 0L) PORT <- 7789L

# ---- Resolve project root and source core R modules ------------------------

.resolve_cas9_app_root <- function() {
  candidates <- unique(c(
    getwd(),
    "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design"
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
    if (dir.exists(file.path(cand, "R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  candidates[1]
}

CAS9_APP_ROOT <- .resolve_cas9_app_root()

.source_cas9_modules <- function(root = CAS9_APP_ROOT, force_homology = FALSE) {
  mod_files <- c(
    "genome_registry.R",
    "gene_lookup.R",
    "grna_engine.R",
    "homology_arms.R",
    "pjat_plasmids.R",
    "donor_payload.R",
    "plasmid_designer_bundle.R"
  )
  for (mf in mod_files) {
    p <- file.path(root, "R", mf)
    if (file.exists(p)) {
      if (identical(mf, "homology_arms.R") && !isTRUE(force_homology) &&
          exists("design_homology_arms", mode = "function", envir = .GlobalEnv) &&
          exists("export_design_bundle_zip", mode = "function", envir = .GlobalEnv)) {
        next
      }
      tryCatch(source(p, local = FALSE), error = function(e) {
        warning("Could not source ", p, ": ", conditionMessage(e))
      })
    }
  }
  if (exists("resolve_chromosome_unified", mode = "function")) {
    assign("resolve_chromosome", resolve_chromosome_unified, envir = .GlobalEnv)
  }
}

.source_cas9_modules(CAS9_APP_ROOT)

# ============================================================================
# 1. Defensive Homology-Arm Normalization & Fallback Wrappers
# ============================================================================

.ensure_homology_module <- function() {
  if (!exists("export_design_bundle_zip", mode = "function", envir = .GlobalEnv)) {
    .source_cas9_modules(CAS9_APP_ROOT, force_homology = TRUE)
  }
}

.fallback_gc_pct <- function(seq_str) {
  if (is.null(seq_str) || !length(seq_str) || is.na(seq_str[1]) || !nzchar(seq_str[1])) return(0.0)
  s <- gsub("[^ACGT]", "", toupper(seq_str[1]))
  if (!nzchar(s)) return(0.0)
  round(100.0 * nchar(gsub("[^GC]", "", s)) / nchar(s), 2)
}

.fallback_design_homology_arms <- function(guide_row_or_params,
                                           species = "melanogaster",
                                           arm_5p_len = 1000L,
                                           arm_3p_len = 1000L,
                                           arm_3p_anchor = c("cut_site", "pam_start"),
                                           gene_model = NULL,
                                           cache = NULL,
                                           ref_dir = get_reference_dir()) {
  arm_3p_anchor <- match.arg(arm_3p_anchor)
  sp <- canonical_species(species)
  r <- if (is.data.frame(guide_row_or_params)) guide_row_or_params[1, , drop = FALSE] else as.data.frame(guide_row_or_params, stringsAsFactors = FALSE)

  chr_val <- as.character(r$chr[1])
  acc_val <- resolve_chromosome_unified(sp, chr_val, ref_dir = ref_dir)
  if (is.na(acc_val) || !nzchar(acc_val)) acc_val <- chr_val
  strand_val <- as.character(r$strand[1])
  if (!strand_val %in% c("+", "-")) strand_val <- "+"

  p_start <- as.integer(r$protospacer_start[1])
  p_end   <- as.integer(r$protospacer_end[1])
  pam_s   <- as.integer(r$pam_start[1])
  pam_e   <- as.integer(r$pam_end[1])
  cut_pos <- as.integer(r$cut_site[1])
  gid     <- if ("guide_id" %in% names(r)) as.character(r$guide_id[1]) else "guide_1"
  spacer  <- as.character(r$spacer_seq[1])
  pam_sq  <- as.character(r$pam_seq[1])

  L5 <- as.integer(round(arm_5p_len[1]))
  L3 <- as.integer(round(arm_3p_len[1]))
  chr_len <- get_chr_length(sp, acc_val, ref_dir = ref_dir)
  if (is.na(chr_len) || chr_len <= 0L) chr_len <- 100000000L

  if (identical(strand_val, "+")) {
    cut_5p_gen <- cut_pos
    cut_3p_gen <- cut_pos + 1L
    a5_start <- max(1L, cut_pos - L5 + 1L)
    a5_end   <- cut_pos
    a3_start <- if (identical(arm_3p_anchor, "pam_start")) pam_s else (cut_pos + 1L)
    a3_end   <- min(chr_len, a3_start + L3 - 1L)
  } else {
    cut_5p_gen <- cut_pos + 1L
    cut_3p_gen <- cut_pos
    a5_start <- cut_pos + 1L
    a5_end   <- min(chr_len, a5_start + L5 - 1L)
    a3_end   <- if (identical(arm_3p_anchor, "pam_start")) pam_e else cut_pos
    a3_start <- max(1L, a3_end - L3 + 1L)
  }

  seq5_plus <- extract_genomic_subseq(sp, acc_val, a5_start, a5_end, strand = "+", ref_dir = ref_dir)
  seq3_plus <- extract_genomic_subseq(sp, acc_val, a3_start, a3_end, strand = "+", ref_dir = ref_dir)
  seq5_tgt  <- if (identical(strand_val, "-")) revcomp_dna(seq5_plus) else seq5_plus
  seq3_tgt  <- if (identical(strand_val, "-")) revcomp_dna(seq3_plus) else seq3_plus

  recleave_flag <- identical(arm_3p_anchor, "cut_site")
  pam_mut_note  <- if (recleave_flag) {
    paste0("3' arm retains the 3-nt spacer tail + ", pam_sq,
           " PAM; introduce a silent mutation in the NGG PAM or seed region (",
           substr(spacer, 18, 20), ") in the donor template to prevent Cas9 re-cleavage.")
  } else {
    "3' arm anchored at PAM start; protospacer is split at the cut site."
  }

  list(
    species             = sp,
    guide_id            = gid,
    chr                 = chr_val,
    accession           = acc_val,
    strand              = strand_val,
    spacer_seq          = spacer,
    pam_seq             = pam_sq,
    protospacer_start   = p_start,
    protospacer_end     = p_end,
    pam_start           = pam_s,
    pam_end             = pam_e,
    cut_site            = cut_pos,
    cut_left            = cut_pos,
    cut_right           = cut_pos + 1L,
    cut_5p_genomic      = cut_5p_gen,
    cut_3p_genomic      = cut_3p_gen,
    arm_3p_anchor       = arm_3p_anchor,
    arm_5p_start        = as.integer(a5_start),
    arm_5p_end          = as.integer(a5_end),
    arm_5p_length       = as.integer(nchar(seq5_tgt)),
    arm_5p_requested    = L5,
    arm_5p_seq          = seq5_tgt,
    arm_5p_seq_plus     = seq5_plus,
    arm_5p_gc           = .fallback_gc_pct(seq5_tgt),
    arm_3p_start        = as.integer(a3_start),
    arm_3p_end          = as.integer(a3_end),
    arm_3p_length       = as.integer(nchar(seq3_tgt)),
    arm_3p_requested    = L3,
    arm_3p_seq          = seq3_tgt,
    arm_3p_seq_plus     = seq3_plus,
    arm_3p_gc           = .fallback_gc_pct(seq3_tgt),
    recleavage_flag     = recleave_flag,
    pam_shield_mutation = pam_mut_note,
    n_tract_flag        = grepl("NNNN", paste0(seq5_tgt, seq3_tgt), fixed = TRUE),
    boundary_truncated  = (nchar(seq5_tgt) < L5) || (nchar(seq3_tgt) < L3),
    qc_pass             = !grepl("NNNN", paste0(seq5_tgt, seq3_tgt), fixed = TRUE) &&
                          (nchar(seq5_tgt) == L5) && (nchar(seq3_tgt) == L3),
    warnings            = character(0)
  )
}

normalize_arm_design <- function(arm_obj,
                                 guide_row,
                                 species = "melanogaster",
                                 arm_5p_len = 1000L,
                                 arm_3p_len = 1000L,
                                 arm_3p_anchor = "cut_site",
                                 gene_model = NULL,
                                 ref_dir = get_reference_dir()) {
  fb <- .fallback_design_homology_arms(
    guide_row_or_params = guide_row,
    species             = species,
    arm_5p_len          = arm_5p_len,
    arm_3p_len          = arm_3p_len,
    arm_3p_anchor       = arm_3p_anchor,
    gene_model          = gene_model,
    ref_dir             = ref_dir
  )
  if (is.null(arm_obj) || !is.list(arm_obj)) return(fb)

  pick_val <- function(primary, secondary = NULL, default = NULL) {
    if (!is.null(primary) && length(primary) > 0 && !all(is.na(primary))) return(primary[1])
    if (!is.null(secondary) && length(secondary) > 0 && !all(is.na(secondary))) return(secondary[1])
    default
  }

  a5 <- arm_obj$arm_5p
  a3 <- arm_obj$arm_3p

  arm_obj$species           <- pick_val(arm_obj$species, NULL, fb$species)
  arm_obj$guide_id          <- pick_val(arm_obj$guide_id, NULL, fb$guide_id)
  arm_obj$chr               <- pick_val(arm_obj$chr, NULL, fb$chr)
  arm_obj$accession         <- pick_val(arm_obj$accession, arm_obj$chr_accession, fb$accession)
  arm_obj$strand            <- pick_val(arm_obj$strand, NULL, fb$strand)
  arm_obj$spacer_seq        <- pick_val(arm_obj$spacer_seq, NULL, fb$spacer_seq)
  arm_obj$pam_seq           <- pick_val(arm_obj$pam_seq, NULL, fb$pam_seq)
  arm_obj$protospacer_start <- as.integer(pick_val(arm_obj$protospacer_start, NULL, fb$protospacer_start))
  arm_obj$protospacer_end   <- as.integer(pick_val(arm_obj$protospacer_end, NULL, fb$protospacer_end))
  arm_obj$pam_start         <- as.integer(pick_val(arm_obj$pam_start, NULL, fb$pam_start))
  arm_obj$pam_end           <- as.integer(pick_val(arm_obj$pam_end, NULL, fb$pam_end))
  arm_obj$cut_site          <- as.integer(pick_val(arm_obj$cut_site, NULL, fb$cut_site))
  arm_obj$cut_left          <- as.integer(pick_val(arm_obj$cut_left, NULL, fb$cut_left))
  arm_obj$cut_right         <- as.integer(pick_val(arm_obj$cut_right, NULL, fb$cut_right))
  arm_obj$cut_5p_genomic    <- as.integer(pick_val(arm_obj$cut_5p_genomic, NULL, fb$cut_5p_genomic))
  arm_obj$cut_3p_genomic    <- as.integer(pick_val(arm_obj$cut_3p_genomic, NULL, fb$cut_3p_genomic))
  arm_obj$arm_3p_anchor     <- pick_val(arm_obj$arm_3p_anchor, NULL, fb$arm_3p_anchor)

  arm_obj$arm_5p_start      <- as.integer(pick_val(arm_obj$arm_5p_start, if (is.list(a5)) a5$start else NULL, fb$arm_5p_start))
  arm_obj$arm_5p_end        <- as.integer(pick_val(arm_obj$arm_5p_end,   if (is.list(a5)) a5$end   else NULL, fb$arm_5p_end))
  arm_obj$arm_5p_seq        <- as.character(pick_val(arm_obj$arm_5p_seq, if (is.list(a5)) pick_val(a5$seq_target_strand, a5$seq) else NULL, fb$arm_5p_seq))
  arm_obj$arm_5p_seq_plus   <- as.character(pick_val(arm_obj$arm_5p_seq_plus, if (is.list(a5)) pick_val(a5$seq_plus_strand, a5$seq_genomic) else NULL, fb$arm_5p_seq_plus))
  arm_obj$arm_5p_length     <- as.integer(nchar(arm_obj$arm_5p_seq))
  arm_obj$arm_5p_gc         <- as.numeric(pick_val(arm_obj$arm_5p_gc, if (is.list(a5)) pick_val(a5$gc_percent, a5$gc_pct) else NULL, fb$arm_5p_gc))

  arm_obj$arm_3p_start      <- as.integer(pick_val(arm_obj$arm_3p_start, if (is.list(a3)) a3$start else NULL, fb$arm_3p_start))
  arm_obj$arm_3p_end        <- as.integer(pick_val(arm_obj$arm_3p_end,   if (is.list(a3)) a3$end   else NULL, fb$arm_3p_end))
  arm_obj$arm_3p_seq        <- as.character(pick_val(arm_obj$arm_3p_seq, if (is.list(a3)) pick_val(a3$seq_target_strand, a3$seq) else NULL, fb$arm_3p_seq))
  arm_obj$arm_3p_seq_plus   <- as.character(pick_val(arm_obj$arm_3p_seq_plus, if (is.list(a3)) pick_val(a3$seq_plus_strand, a3$seq_genomic) else NULL, fb$arm_3p_seq_plus))
  arm_obj$arm_3p_length     <- as.integer(nchar(arm_obj$arm_3p_seq))
  arm_obj$arm_3p_gc         <- as.numeric(pick_val(arm_obj$arm_3p_gc, if (is.list(a3)) pick_val(a3$gc_percent, a3$gc_pct) else NULL, fb$arm_3p_gc))

  arm_obj$recleavage_flag   <- isTRUE(pick_val(arm_obj$recleavage_flag, arm_obj$recleavage_risk, fb$recleavage_flag))
  arm_obj$pam_shield_mutation <- as.character(pick_val(
    arm_obj$pam_shield_mutation,
    if (is.list(arm_obj$shielding) && !is.null(arm_obj$shielding$summary)) arm_obj$shielding$summary else arm_obj$pam_blocking_suggestion,
    fb$pam_shield_mutation
  ))
  arm_obj$n_tract_flag      <- isTRUE(pick_val(arm_obj$n_tract_flag, NULL, fb$n_tract_flag))
  arm_obj$boundary_truncated <- isTRUE(pick_val(arm_obj$boundary_truncated, NULL, fb$boundary_truncated))
  arm_obj$qc_pass           <- isTRUE(pick_val(arm_obj$qc_pass, NULL, fb$qc_pass))
  if (is.null(arm_obj$warnings)) arm_obj$warnings <- fb$warnings
  arm_obj
}

compute_cut_site_zoom_view <- function(arm_design, flank_bp = 30L) {
  .ensure_homology_module()
  zoom_raw <- NULL
  if (exists("build_cut_site_zoom", mode = "function")) {
    zoom_raw <- tryCatch(build_cut_site_zoom(arm_design, flank_bp = flank_bp), error = function(e) NULL)
  }

  spacer <- as.character(arm_design$spacer_seq[1])
  pam    <- as.character(arm_design$pam_seq[1])
  strand <- as.character(arm_design$strand[1])

  arm5_tail <- substr(arm_design$arm_5p_seq, max(1L, nchar(arm_design$arm_5p_seq) - flank_bp + 1L), nchar(arm_design$arm_5p_seq))
  spacer_5p_17 <- substr(spacer, 1L, 17L)
  spacer_3p_3  <- substr(spacer, 18L, 20L)
  arm3_head <- substr(arm_design$arm_3p_seq, 1L, min(nchar(arm_design$arm_3p_seq), flank_bp))

  pre_proto <- if (nchar(arm_design$arm_5p_seq) > 17L) {
    s_idx <- max(1L, nchar(arm_design$arm_5p_seq) - 17L - flank_bp + 1L)
    substr(arm_design$arm_5p_seq, s_idx, nchar(arm_design$arm_5p_seq) - 17L)
  } else ""

  post_pam <- if (identical(arm_design$arm_3p_anchor, "cut_site") && nchar(arm_design$arm_3p_seq) > 6L) {
    substr(arm_design$arm_3p_seq, 7L, min(nchar(arm_design$arm_3p_seq), 6L + flank_bp))
  } else if (nchar(arm_design$arm_3p_seq) > 3L) {
    substr(arm_design$arm_3p_seq, 4L, min(nchar(arm_design$arm_3p_seq), 3L + flank_bp))
  } else ""

  comp_base <- function(s) {
    if (is.null(s) || !nzchar(s)) return("")
    chartr("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn", s)
  }

  target_line <- paste0(pre_proto, spacer_5p_17, "|", spacer_3p_3, pam, post_pam)
  comp_line   <- paste0(comp_base(pre_proto), comp_base(spacer_5p_17), "|",
                        comp_base(spacer_3p_3), comp_base(pam), comp_base(post_pam))

  list(
    raw              = zoom_raw,
    strand           = strand,
    cut_site         = arm_design$cut_site,
    cut_5p_genomic   = arm_design$cut_5p_genomic,
    cut_3p_genomic   = arm_design$cut_3p_genomic,
    pre_protospacer  = pre_proto,
    spacer_5p_17     = spacer_5p_17,
    spacer_3p_3      = spacer_3p_3,
    pam_seq          = pam,
    post_pam         = post_pam,
    comp_pre         = comp_base(pre_proto),
    comp_spacer_17   = comp_base(spacer_5p_17),
    comp_spacer_3    = comp_base(spacer_3p_3),
    comp_pam         = comp_base(pam),
    comp_post        = comp_base(post_pam),
    arm5_tail        = arm5_tail,
    arm3_head        = arm3_head,
    target_5to3_text = target_line,
    comp_3to5_text   = comp_line
  )
}

safe_design_homology_arms <- function(guide_row,
                                      species = "melanogaster",
                                      arm_5p_len = 1000L,
                                      arm_3p_len = 1000L,
                                      arm_3p_anchor = "cut_site",
                                      gene_model = NULL,
                                      ref_dir = get_reference_dir()) {
  .ensure_homology_module()
  raw_arm <- NULL
  if (exists("design_homology_arms", mode = "function")) {
    raw_arm <- tryCatch(
      design_homology_arms(
        guide_row_or_params = guide_row,
        species             = species,
        arm_5p_len          = arm_5p_len,
        arm_3p_len          = arm_3p_len,
        arm_3p_anchor       = arm_3p_anchor,
        gene_model          = gene_model,
        ref_dir             = ref_dir
      ),
      error = function(e) NULL
    )
  }
  normalize_arm_design(
    arm_obj       = raw_arm,
    guide_row     = guide_row,
    species       = species,
    arm_5p_len    = arm_5p_len,
    arm_3p_len    = arm_3p_len,
    arm_3p_anchor = arm_3p_anchor,
    gene_model    = gene_model,
    ref_dir       = ref_dir
  )
}

safe_export_fasta <- function(arm_design, file = NULL, include_plus_strand = TRUE) {
  .ensure_homology_module()
  if (exists("export_homology_arms_fasta", mode = "function")) {
    res <- tryCatch(
      export_homology_arms_fasta(arm_design, file = file, include_plus_strand = include_plus_strand),
      error = function(e) NULL
    )
    if (!is.null(res) && (!is.null(file) && file.exists(file) || is.character(res))) {
      return(res)
    }
  }
  lines <- c(
    sprintf(">%s|5p_homology_arm|species=%s|chr=%s:%d-%d|strand=%s|len=%d|cut_site=%d",
            arm_design$guide_id, arm_design$species, arm_design$chr,
            arm_design$arm_5p_start, arm_design$arm_5p_end, arm_design$strand,
            arm_design$arm_5p_length, arm_design$cut_site),
    arm_design$arm_5p_seq,
    sprintf(">%s|3p_homology_arm|species=%s|chr=%s:%d-%d|strand=%s|len=%d|anchor=%s",
            arm_design$guide_id, arm_design$species, arm_design$chr,
            arm_design$arm_3p_start, arm_design$arm_3p_end, arm_design$strand,
            arm_design$arm_3p_length, arm_design$arm_3p_anchor),
    arm_design$arm_3p_seq
  )
  if (isTRUE(include_plus_strand)) {
    lines <- c(
      lines,
      sprintf(">%s|5p_homology_arm_plus_strand|chr=%s:%d-%d(+)",
              arm_design$guide_id, arm_design$chr, arm_design$arm_5p_start, arm_design$arm_5p_end),
      arm_design$arm_5p_seq_plus,
      sprintf(">%s|3p_homology_arm_plus_strand|chr=%s:%d-%d(+)",
              arm_design$guide_id, arm_design$chr, arm_design$arm_3p_start, arm_design$arm_3p_end),
      arm_design$arm_3p_seq_plus
    )
  }
  txt <- paste(lines, collapse = "\n")
  if (!is.null(file)) writeLines(txt, file)
  txt
}

safe_export_csv <- function(grna_df, arm_design = NULL, file = NULL) {
  .ensure_homology_module()
  if (exists("export_design_csv", mode = "function")) {
    res <- tryCatch(export_design_csv(grna_df, arm_design = arm_design, file = file), error = function(e) NULL)
    if (!is.null(res) && (!is.null(file) && file.exists(file) || is.character(res) || is.data.frame(res))) {
      return(res)
    }
  }
  df_out <- grna_df
  if (!is.null(arm_design) && nrow(df_out) > 0L) {
    sel_mask <- df_out$guide_id == arm_design$guide_id
    df_out$selected_for_arms <- sel_mask
    df_out$arm_5p_coords <- ifelse(sel_mask, sprintf("%s:%d-%d", arm_design$chr, arm_design$arm_5p_start, arm_design$arm_5p_end), NA_character_)
    df_out$arm_3p_coords <- ifelse(sel_mask, sprintf("%s:%d-%d", arm_design$chr, arm_design$arm_3p_start, arm_design$arm_3p_end), NA_character_)
  }
  if (!is.null(file)) utils::write.csv(df_out, file = file, row.names = FALSE)
  invisible(df_out)
}

safe_export_json <- function(grna_df, arm_design = NULL, metadata = NULL, file = NULL) {
  .ensure_homology_module()
  if (exists("export_design_json", mode = "function")) {
    res <- tryCatch(export_design_json(grna_df, arm_design = arm_design, metadata = metadata, file = file), error = function(e) NULL)
    if (!is.null(res) && (!is.null(file) && file.exists(file) || is.character(res))) {
      return(res)
    }
  }
  payload <- list(
    app_version   = CAS9_APP_VERSION,
    metadata      = if (!is.null(metadata)) metadata else attr(grna_df, "metadata"),
    selected_arms = if (!is.null(arm_design)) arm_design[setdiff(names(arm_design), c("raw_row", "attr_gene_model"))] else NULL,
    top_guides    = if (is.data.frame(grna_df) && nrow(grna_df) > 0L) utils::head(grna_df, 25L) else list()
  )
  json_txt <- jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (!is.null(file)) writeLines(as.character(json_txt), file)
  as.character(json_txt)
}

safe_export_log <- function(grna_df, arm_design = NULL, metadata = NULL, file = NULL) {
  .ensure_homology_module()
  if (exists("export_design_log", mode = "function")) {
    res <- tryCatch(export_design_log(grna_df, arm_design = arm_design, metadata = metadata, file = file), error = function(e) NULL)
    if (!is.null(res) && (!is.null(file) && file.exists(file) || is.character(res))) {
      return(res)
    }
  }
  meta <- if (!is.null(metadata)) metadata else attr(grna_df, "metadata")
  lines <- c(
    "============================================================================",
    "CRISPR Cas9 gRNA & Homology Arm Design — Reproducibility Log",
    "============================================================================",
    paste0("App Version          : ", CAS9_APP_VERSION),
    paste0("gRNA Engine Version  : ", if (!is.null(meta$engine_version)) meta$engine_version else "1.0.0"),
    paste0("R Version            : ", R.version.string),
    paste0("Species              : ", if (!is.null(meta$species)) meta$species else if (!is.null(arm_design)) arm_design$species else "melanogaster"),
    paste0("Assembly Accession   : ", if (!is.null(meta$assembly_accession)) meta$assembly_accession else "RefSeq"),
    paste0("Query                : ", if (!is.null(meta$query)) meta$query else ""),
    paste0("Target Priority Mode : ", if (!is.null(meta$target_mode)) meta$target_mode else "coding_ko"),
    paste0("Off-Target Mode      : ", if (!is.null(meta$offtarget_mode)) meta$offtarget_mode else "fast"),
    paste0("RNG Seed             : ", if (!is.null(meta$seed)) meta$seed else 42L),
    paste0("Total Candidates     : ", if (is.data.frame(grna_df)) nrow(grna_df) else 0L),
    "----------------------------------------------------------------------------"
  )
  if (!is.null(arm_design)) {
    lines <- c(
      lines,
      paste0("Selected Guide ID    : ", arm_design$guide_id),
      paste0("Spacer (20-nt)       : ", arm_design$spacer_seq),
      paste0("PAM (NGG)            : ", arm_design$pam_seq),
      paste0("Target Strand        : ", arm_design$strand),
      paste0("Blunt Cut Site (1b)  : ", arm_design$chr, ":", arm_design$cut_site, " (between ",
             arm_design$cut_5p_genomic, " and ", arm_design$cut_3p_genomic, ")"),
      paste0("5' Homology Arm      : ", arm_design$chr, ":", arm_design$arm_5p_start, "-", arm_design$arm_5p_end,
             " (", arm_design$arm_5p_length, " bp, GC=", arm_design$arm_5p_gc, "%)"),
      paste0("3' Homology Arm      : ", arm_design$chr, ":", arm_design$arm_3p_start, "-", arm_design$arm_3p_end,
             " (", arm_design$arm_3p_length, " bp, GC=", arm_design$arm_3p_gc, "%, anchor=", arm_design$arm_3p_anchor, ")"),
      paste0("PAM Shielding Note   : ", arm_design$pam_shield_mutation),
      "============================================================================"
    )
  }
  txt <- paste(lines, collapse = "\n")
  if (!is.null(file)) writeLines(txt, file)
  txt
}

safe_export_zip <- function(zip_path, grna_df, arm_design = NULL, gene_model = NULL, metadata = NULL, plot_fn = NULL) {
  .ensure_homology_module()
  tmp_dir <- tempfile("crispr_bundle_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  csv_file   <- file.path(tmp_dir, "gRNA_candidates.csv")
  fasta_file <- file.path(tmp_dir, "homology_arms.fasta")
  json_file  <- file.path(tmp_dir, "design_summary.json")
  log_file   <- file.path(tmp_dir, "CRISPR_design_log.txt")
  pdf_file   <- file.path(tmp_dir, "isoform_track_plot.pdf")

  safe_export_csv(grna_df, arm_design = arm_design, file = csv_file)
  if (!is.null(arm_design)) {
    safe_export_fasta(arm_design, file = fasta_file, include_plus_strand = TRUE)
  }
  safe_export_json(grna_df, arm_design = arm_design, metadata = metadata, file = json_file)
  safe_export_log(grna_df, arm_design = arm_design, metadata = metadata, file = log_file)

  if (is.function(plot_fn)) {
    tryCatch({
      grDevices::pdf(pdf_file, width = 11, height = 6)
      plot_fn()
      grDevices::dev.off()
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
    })
  }

  files_to_zip <- list.files(tmp_dir, full.names = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  utils::zip(zipfile = zip_path, files = files_to_zip, flags = "-q")
  invisible(zip_path)
}

# ============================================================================
# 2. Isoform + Cut-Site & Homology-Arm Track Schematic Plot
#    (Modeled directly on plot_isoform_tracks() in HCR_probe_design_v41.R)
# ============================================================================

TIER_COLORS <- c(
  "1" = "#059669",  # Early Shared CDS (Emerald)
  "2" = "#7c3aed",  # 5' UTR / Promoter (Violet)
  "3" = "#0284c7",  # Late Shared CDS (Sky)
  "4" = "#d97706",  # Isoform-Specific CDS (Amber)
  "5" = "#64748b"   # 3' UTR / Intron (Slate)
)

TIER_LABELS <- c(
  "1" = "Tier 1: Early Shared CDS",
  "2" = "Tier 2: 5' UTR / Promoter",
  "3" = "Tier 3: Late Shared CDS",
  "4" = "Tier 4: Isoform-Specific CDS",
  "5" = "Tier 5: 3' UTR / Intron"
)

flip_coord_se <- function(s, e, lo, hi) {
  list(s = lo + hi - e, e = lo + hi - s)
}

plot_crispr_isoform_tracks <- function(gene_model,
                                       candidates_df = NULL,
                                       selected_guide = NULL,
                                       arm_design = NULL,
                                       orient_5p_to_3p = FALSE,
                                       show_all_candidates = TRUE) {
  if (is.null(gene_model)) return(invisible(NULL))

  exons_by_tx <- gene_model$exons_by_tx
  cds_by_tx   <- gene_model$cds_by_tx
  shared_cds  <- gene_model$shared_cds
  five_utr    <- gene_model$five_prime_utr
  prom_df     <- gene_model$promoter_regions

  if (is.null(exons_by_tx) || !length(exons_by_tx)) {
    exons_by_tx <- list(
      locus = data.frame(
        chr    = as.character(gene_model$chr[1]),
        start  = as.integer(gene_model$start[1]),
        end    = as.integer(gene_model$end[1]),
        strand = as.character(gene_model$strand[1]),
        stringsAsFactors = FALSE
      )
    )
  }

  tx_ids <- names(exons_by_tx)
  if (is.null(tx_ids) || any(!nzchar(tx_ids))) {
    tx_ids <- paste0("tx_", seq_along(exons_by_tx))
    names(exons_by_tx) <- tx_ids
  }

  tx_starts <- vapply(exons_by_tx, function(df) if (nrow(df)) min(df$start) else Inf, numeric(1))
  tx_ends   <- vapply(exons_by_tx, function(df) if (nrow(df)) max(df$end)   else -Inf, numeric(1))
  ord <- order(tx_starts, tx_ends)
  exons_by_tx <- exons_by_tx[ord]
  tx_ids      <- tx_ids[ord]

  tx_meta <- gene_model$transcripts
  row_labels <- tx_ids
  if (is.data.frame(tx_meta) && nrow(tx_meta) > 0L) {
    m <- match(tx_ids, tx_meta$transcript_id)
    v_lbl <- tx_meta$variant[m]
    row_labels <- ifelse(!is.na(v_lbl) & nzchar(v_lbl) & v_lbl != tx_ids,
                         paste0(tx_ids, " (", v_lbl, ")"),
                         tx_ids)
  }
  names(row_labels) <- tx_ids

  has_shared <- is.data.frame(shared_cds) && nrow(shared_cds) > 0L
  utr_prom_df <- rbind(
    if (is.data.frame(five_utr) && nrow(five_utr)) five_utr[, c("chr", "start", "end", "strand"), drop = FALSE] else NULL,
    if (is.data.frame(prom_df)  && nrow(prom_df))  prom_df[,  c("chr", "start", "end", "strand"), drop = FALSE] else NULL
  )
  has_utr_prom <- !is.null(utr_prom_df) && nrow(utr_prom_df) > 0L
  has_cands    <- isTRUE(show_all_candidates) && is.data.frame(candidates_df) && nrow(candidates_df) > 0L
  has_arms     <- !is.null(arm_design)

  x_min <- min(tx_starts[is.finite(tx_starts)])
  x_max <- max(tx_ends[is.finite(tx_ends)])
  if (!is.finite(x_min) || !is.finite(x_max)) {
    x_min <- as.numeric(gene_model$start[1])
    x_max <- as.numeric(gene_model$end[1])
  }
  if (has_shared) {
    x_min <- min(x_min, min(shared_cds$start))
    x_max <- max(x_max, max(shared_cds$end))
  }
  if (has_utr_prom) {
    x_min <- min(x_min, min(utr_prom_df$start))
    x_max <- max(x_max, max(utr_prom_df$end))
  }
  if (has_arms) {
    x_min <- min(x_min, arm_design$arm_5p_start, arm_design$arm_3p_start)
    x_max <- max(x_max, arm_design$arm_5p_end,   arm_design$arm_3p_end)
  }

  locus_strand <- if (!is.null(gene_model$strand) && gene_model$strand[1] %in% c("+", "-")) {
    gene_model$strand[1]
  } else "+"

  flip_view <- isTRUE(orient_5p_to_3p) && identical(locus_strand, "-")
  lo_ref <- x_min
  hi_ref <- x_max

  map_x <- function(pos) {
    if (flip_view) (lo_ref + hi_ref - as.numeric(pos)) else as.numeric(pos)
  }
  map_interval <- function(s, e) {
    if (flip_view) flip_coord_se(s, e, lo_ref, hi_ref) else list(s = as.numeric(s), e = as.numeric(e))
  }

  n_tx <- length(exons_by_tx)
  extra_rows <- (if (has_shared) 1L else 0L) +
                (if (has_utr_prom) 1L else 0L) +
                (if (has_cands) 1L else 0L) +
                (if (has_arms) 1L else 0L)
  n_rows <- n_tx + extra_rows

  pad <- max(50, 0.03 * (x_max - x_min))
  max_label_chars <- max(c(nchar(row_labels), 18L), na.rm = TRUE)
  left_in <- max(1.8, min(3.8, 0.45 + 0.068 * max_label_chars))

  old_par <- graphics::par(mai = c(0.65, left_in, 0.48, 0.25))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot(
    NA,
    xlim = c(x_min - pad, x_max + pad),
    ylim = c(0.35, n_rows + 0.65),
    yaxt = "n",
    ylab = "",
    xlab = "",
    bty  = "n",
    xaxt = "n"
  )

  # Top header banner
  gene_sym <- if (!is.null(gene_model$gene_symbol)) gene_model$gene_symbol[1] else "Locus"
  chr_lbl  <- if (!is.null(gene_model$chr)) gene_model$chr[1] else ""
  hdr_txt  <- sprintf("%s (%s, %s:%s-%s, strand %s)%s",
                      gene_sym,
                      gene_model$species,
                      chr_lbl,
                      format(as.integer(gene_model$start[1]), big.mark = ","),
                      format(as.integer(gene_model$end[1]),   big.mark = ","),
                      locus_strand,
                      if (flip_view) "  [5' -> 3' oriented view]" else "  [genomic coordinate view]")
  graphics::mtext(hdr_txt, side = 3, line = 1.1, adj = 0, cex = 0.82, font = 2, col = "#0f172a")

  # X-axis ticks
  at_ticks <- graphics::axTicks(1)
  lbl_ticks <- if (flip_view) {
    format(round(lo_ref + hi_ref - at_ticks), big.mark = ",", scientific = FALSE, trim = TRUE)
  } else {
    format(round(at_ticks), big.mark = ",", scientific = FALSE, trim = TRUE)
  }
  graphics::axis(1, at = at_ticks, labels = lbl_ticks, cex.axis = 0.75, col = "#cbd5e1", col.axis = "#475569")
  graphics::mtext(
    if (flip_view) paste0("Genomic coordinate on ", chr_lbl, " (bp, minus-strand flipped 5' -> 3')")
    else paste0("Genomic coordinate on ", chr_lbl, " (bp)"),
    side = 1, line = 2.0, cex = 0.74, col = "#475569"
  )

  draw_interval_track <- function(y, ex_df, col_fill, strand_char = "+", cds_df = NULL, border_col = NA) {
    if (is.null(ex_df) || !nrow(ex_df)) return(invisible(NULL))
    iv <- map_interval(ex_df$start, ex_df$end)
    s <- iv$s
    e <- iv$e
    graphics::segments(min(s), y, max(e), y, col = "#94a3b8", lwd = 2, lend = 2)

    disp_strand <- if (flip_view && strand_char == "-") "+" else if (flip_view && strand_char == "+") "-" else strand_char
    if (!is.null(disp_strand) && disp_strand %in% c("+", "-") && (max(e) - min(s)) > 40) {
      xs <- seq(min(s), max(e), length.out = 7)[2:6]
      graphics::text(xs, rep(y, length(xs)), labels = if (disp_strand == "+") ">" else "<",
                     col = "#64748b", cex = 0.62)
    }

    if (!is.null(cds_df) && nrow(cds_df) > 0L) {
      graphics::rect(s, y - 0.10, e, y + 0.10, col = "#cbd5e1", border = NA)
      iv_cds <- map_interval(cds_df$start, cds_df$end)
      graphics::rect(iv_cds$s, y - 0.21, iv_cds$e, y + 0.21, col = col_fill, border = border_col)
    } else {
      graphics::rect(s, y - 0.19, e, y + 0.19, col = col_fill, border = border_col)
    }
  }

  # 1. Transcript isoform rows (top to bottom)
  for (i in seq_along(tx_ids)) {
    y <- n_rows - i + 1L
    tid <- tx_ids[i]
    ex_df <- exons_by_tx[[tid]]
    cds_df <- if (!is.null(cds_by_tx) && tid %in% names(cds_by_tx)) cds_by_tx[[tid]] else NULL
    st_ch <- if (nrow(ex_df) && "strand" %in% names(ex_df)) ex_df$strand[1] else locus_strand
    draw_interval_track(y, ex_df, col_fill = "#0ea5e9", strand_char = st_ch, cds_df = cds_df)
    graphics::text(graphics::par("usr")[1], y, labels = row_labels[[tid]],
                   adj = c(1, 0.5), xpd = NA, cex = 0.70, font = 1, col = "#1e293b")
  }

  cur_y <- n_rows - n_tx

  # 2. Shared CDS track (emerald green)
  if (has_shared) {
    draw_interval_track(cur_y, shared_cds, col_fill = "#10b981", strand_char = locus_strand)
    graphics::text(graphics::par("usr")[1], cur_y, labels = "shared CDS",
                   adj = c(1, 0.5), xpd = NA, cex = 0.72, font = 2, col = "#059669")
    cur_y <- cur_y - 1L
  }

  # 3. 5' UTR / Promoter track (violet)
  if (has_utr_prom) {
    draw_interval_track(cur_y, utr_prom_df, col_fill = "#8b5cf6", strand_char = locus_strand)
    graphics::text(graphics::par("usr")[1], cur_y, labels = "5' UTR / promoter",
                   adj = c(1, 0.5), xpd = NA, cex = 0.72, font = 2, col = "#6d28d9")
    cur_y <- cur_y - 1L
  }

  # 4. Candidate gRNA cut-site ticks track
  if (has_cands) {
    graphics::segments(x_min, cur_y, x_max, cur_y, col = "#e2e8f0", lwd = 1.5)
    c_sub <- if (nrow(candidates_df) > 150L) candidates_df[seq_len(150L), , drop = FALSE] else candidates_df
    cx <- map_x(c_sub$cut_site)
    t_cols <- TIER_COLORS[as.character(c_sub$context_tier)]
    t_cols[is.na(t_cols)] <- "#64748b"
    graphics::segments(cx, cur_y - 0.17, cx, cur_y + 0.17, col = grDevices::adjustcolor(t_cols, alpha.f = 0.78), lwd = 1.6)
    graphics::text(graphics::par("usr")[1], cur_y, labels = paste0("gRNA cuts (n=", nrow(candidates_df), ")"),
                   adj = c(1, 0.5), xpd = NA, cex = 0.70, font = 2, col = "#334155")
    cur_y <- cur_y - 1L
  }

  # 5. Selected gRNA Cut Site & Flanking 5'/3' Homology Arms track
  if (has_arms) {
    iv5 <- map_interval(arm_design$arm_5p_start, arm_design$arm_5p_end)
    iv3 <- map_interval(arm_design$arm_3p_start, arm_design$arm_3p_end)
    cut_x <- map_x(arm_design$cut_site)

    # Shaded background band across all tracks for 5' and 3' homology arms
    graphics::rect(iv5$s, 0.45, iv5$e, n_rows + 0.45,
                   col = grDevices::adjustcolor("#2563eb", alpha.f = 0.08), border = NA)
    graphics::rect(iv3$s, 0.45, iv3$e, n_rows + 0.45,
                   col = grDevices::adjustcolor("#d97706", alpha.f = 0.08), border = NA)

    # Dedicated homology arm blocks on bottom track
    graphics::rect(iv5$s, cur_y - 0.22, iv5$e, cur_y + 0.22,
                   col = grDevices::adjustcolor("#2563eb", alpha.f = 0.85), border = "#1d4ed8", lwd = 1.2)
    graphics::rect(iv3$s, cur_y - 0.22, iv3$e, cur_y + 0.22,
                   col = grDevices::adjustcolor("#f59e0b", alpha.f = 0.85), border = "#b45309", lwd = 1.2)

    graphics::text(mean(c(iv5$s, iv5$e)), cur_y,
                   labels = paste0("5' Arm (", arm_design$arm_5p_length, " bp)"),
                   col = "#ffffff", cex = 0.65, font = 2)
    graphics::text(mean(c(iv3$s, iv3$e)), cur_y,
                   labels = paste0("3' Arm (", arm_design$arm_3p_length, " bp)"),
                   col = "#ffffff", cex = 0.65, font = 2)

    # Vertical red dashed cut line across entire plot + red triangle marker
    graphics::abline(v = cut_x, col = "#dc2626", lwd = 1.8, lty = 2)
    graphics::points(cut_x, cur_y + 0.28, pch = 25, bg = "#dc2626", col = "#991b1b", cex = 1.15)

    graphics::text(graphics::par("usr")[1], cur_y,
                   labels = paste0(arm_design$guide_id, " + arms"),
                   adj = c(1, 0.5), xpd = NA, cex = 0.72, font = 2, col = "#dc2626")

    sub_txt <- sprintf(
      "Selected: %s (%s%s, strand %s) | Blunt cut @ %s:%s (3 bp 5' of PAM) | 5' Arm: %s-%s | 3' Arm: %s-%s",
      arm_design$guide_id,
      arm_design$spacer_seq,
      arm_design$pam_seq,
      arm_design$strand,
      arm_design$chr,
      format(arm_design$cut_site, big.mark = ","),
      format(arm_design$arm_5p_start, big.mark = ","),
      format(arm_design$arm_5p_end,   big.mark = ","),
      format(arm_design$arm_3p_start, big.mark = ","),
      format(arm_design$arm_3p_end,   big.mark = ",")
    )
    graphics::mtext(sub_txt, side = 3, line = 0.05, adj = 0, cex = 0.68, col = "#dc2626", font = 2)
  }

  invisible(TRUE)
}

# ============================================================================
# 3. End-to-End Headless & Interactive Design Workflow Engine
# ============================================================================

filter_grna_candidates <- function(candidates_df,
                                   allowed_tiers = 1:5,
                                   min_on_target = 0,
                                   min_specificity = 0,
                                   exclude_tttt = FALSE) {
  if (!is.data.frame(candidates_df) || nrow(candidates_df) == 0L) return(candidates_df)
  keep <- rep(TRUE, nrow(candidates_df))
  if (!is.null(allowed_tiers) && length(allowed_tiers) > 0L && "context_tier" %in% names(candidates_df)) {
    keep <- keep & (candidates_df$context_tier %in% as.integer(allowed_tiers))
  }
  if (!is.null(min_on_target) && is.finite(min_on_target) && min_on_target > 0) {
    keep <- keep & (candidates_df$on_target_score >= as.numeric(min_on_target))
  }
  if (!is.null(min_specificity) && is.finite(min_specificity) && min_specificity > 0) {
    keep <- keep & (candidates_df$cfd_specificity_score >= as.numeric(min_specificity))
  }
  if (isTRUE(exclude_tttt) && "u6_TTTT_flag" %in% names(candidates_df)) {
    keep <- keep & (!candidates_df$u6_TTTT_flag)
  }
  out <- candidates_df[keep, , drop = FALSE]
  if (nrow(out) > 0L) {
    out$filtered_rank <- seq_len(nrow(out))
    rownames(out) <- NULL
  }
  attr(out, "metadata")   <- attr(candidates_df, "metadata")
  attr(out, "gene_model") <- attr(candidates_df, "gene_model")
  out
}

run_crispr_design_workflow <- function(query,
                                       species = "melanogaster",
                                       target_mode = c("coding_ko", "promoter_5utr", "all"),
                                       offtarget_mode = c("fast", "full", "locus"),
                                       upstream_promoter_bp = 500L,
                                       max_offtarget_candidates = 80L,
                                       allowed_tiers = 1:5,
                                       min_on_target = 0,
                                       min_specificity = 0,
                                       exclude_tttt = FALSE,
                                       arm_5p_len = 1000L,
                                       arm_3p_len = 1000L,
                                       arm_3p_anchor = c("cut_site", "pam_start"),
                                       selected_rank = 1L,
                                       selected_guide_id = NULL,
                                       seed = 42L,
                                       ref_dir = get_reference_dir()) {
  target_mode    <- match.arg(target_mode)
  offtarget_mode <- match.arg(offtarget_mode)
  arm_3p_anchor  <- match.arg(arm_3p_anchor)
  sp             <- canonical_species(species)
  set.seed(as.integer(seed[1]))

  warnings_vec <- character(0)

  gm <- lookup_gene(
    query       = query,
    species     = sp,
    promoter_bp = as.integer(upstream_promoter_bp[1]),
    ref_dir     = ref_dir
  )

  all_cands <- design_grnas(
    query_or_gene            = gm,
    species                  = sp,
    target_mode              = target_mode,
    offtarget_mode           = offtarget_mode,
    upstream_promoter_bp     = as.integer(upstream_promoter_bp[1]),
    max_offtarget_candidates = as.integer(max_offtarget_candidates[1]),
    seed                     = as.integer(seed[1]),
    ref_dir                  = ref_dir
  )

  filt_cands <- filter_grna_candidates(
    candidates_df   = all_cands,
    allowed_tiers   = allowed_tiers,
    min_on_target   = min_on_target,
    min_specificity = min_specificity,
    exclude_tttt    = exclude_tttt
  )

  if (nrow(filt_cands) == 0L && nrow(all_cands) > 0L) {
    warnings_vec <- c(warnings_vec, "No candidates matched all active score/tier filters; showing unfiltered candidates.")
    filt_cands <- all_cands
    filt_cands$filtered_rank <- seq_len(nrow(filt_cands))
  }

  sel_row <- NULL
  if (nrow(filt_cands) > 0L) {
    idx <- 1L
    if (!is.null(selected_guide_id) && nzchar(selected_guide_id[1])) {
      m <- match(selected_guide_id[1], filt_cands$guide_id)
      if (!is.na(m)) idx <- m
    } else if (!is.null(selected_rank) && is.finite(selected_rank[1])) {
      idx <- max(1L, min(nrow(filt_cands), as.integer(selected_rank[1])))
    }
    sel_row <- filt_cands[idx, , drop = FALSE]
  }

  arm_obj  <- NULL
  zoom_obj <- NULL
  if (!is.null(sel_row) && nrow(sel_row) > 0L) {
    arm_obj <- safe_design_homology_arms(
      guide_row     = sel_row,
      species       = sp,
      arm_5p_len    = as.integer(arm_5p_len[1]),
      arm_3p_len    = as.integer(arm_3p_len[1]),
      arm_3p_anchor = arm_3p_anchor,
      gene_model    = gm,
      ref_dir       = ref_dir
    )
    zoom_obj <- compute_cut_site_zoom_view(arm_obj, flank_bp = 30L)
  }

  meta <- attr(all_cands, "metadata")
  if (is.null(meta)) meta <- list(species = sp, query = query, target_mode = target_mode, offtarget_mode = offtarget_mode, seed = seed)

  list(
    status              = "ok",
    species             = sp,
    query               = query,
    gene_model          = gm,
    all_candidates      = all_cands,
    filtered_candidates = filt_cands,
    selected_guide      = sel_row,
    arm_design          = arm_obj,
    cut_site_zoom       = zoom_obj,
    metadata            = meta,
    warnings            = warnings_vec
  )
}

write_crispr_output_bundle <- function(out_dir, workflow_result, orient_5p_to_3p = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cands <- workflow_result$filtered_candidates
  arm   <- workflow_result$arm_design
  gm    <- workflow_result$gene_model
  meta  <- workflow_result$metadata

  csv_path   <- file.path(out_dir, "gRNA_candidates.csv")
  fasta_path <- file.path(out_dir, "homology_arms.fasta")
  json_path  <- file.path(out_dir, "design_summary.json")
  log_path   <- file.path(out_dir, "CRISPR_design_log.txt")
  pdf_path   <- file.path(out_dir, "isoform_track_plot.pdf")
  zip_path   <- file.path(out_dir, "CRISPR_design_bundle.zip")

  safe_export_csv(cands, arm_design = arm, file = csv_path)
  if (!is.null(arm)) {
    safe_export_fasta(arm, file = fasta_path, include_plus_strand = TRUE)
  }
  safe_export_json(cands, arm_design = arm, metadata = meta, file = json_path)
  safe_export_log(cands, arm_design = arm, metadata = meta, file = log_path)

  plot_fn <- function() {
    plot_crispr_isoform_tracks(
      gene_model          = gm,
      candidates_df       = cands,
      selected_guide      = workflow_result$selected_guide,
      arm_design          = arm,
      orient_5p_to_3p     = orient_5p_to_3p,
      show_all_candidates = TRUE
    )
  }

  tryCatch({
    grDevices::pdf(pdf_path, width = 11, height = 6)
    plot_fn()
    grDevices::dev.off()
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
  })

  safe_export_zip(
    zip_path   = zip_path,
    grna_df    = cands,
    arm_design = arm,
    gene_model = gm,
    metadata   = meta,
    plot_fn    = plot_fn
  )

  invisible(list(
    out_dir    = normalizePath(out_dir, mustWork = FALSE),
    csv_path   = csv_path,
    fasta_path = fasta_path,
    json_path  = json_path,
    log_path   = log_path,
    pdf_path   = pdf_path,
    zip_path   = zip_path
  ))
}

# ============================================================================
# 4. UI Helpers, Modern Slate/Sky/Emerald CSS Theme & Shiny UI Definition
# ============================================================================

species_choices_labeled <- function(ref_dir = get_reference_dir()) {
  reg <- tryCatch(get_species_registry(ref_dir = ref_dir), error = function(e) NULL)
  if (is.null(reg) || !nrow(reg)) {
    return(c(
      "D. melanogaster (Ready)" = "melanogaster",
      "D. simulans (Ready)"     = "simulans",
      "D. eugracilis (Ready)"   = "eugracilis",
      "D. erecta (.gz)"         = "erecta",
      "D. teissieri (.gz)"      = "teissieri",
      "D. ananassae (NCBI)"     = "ananassae",
      "D. virilis (NCBI)"       = "virilis"
    ))
  }
  status_tag <- ifelse(
    reg$status == "ready", "Ready",
    ifelse(reg$status == "compressed", "Compressed (.gz)", "Needs NCBI Download")
  )
  lbls <- sprintf("D. %s [%s — %s]", reg$species_key, reg$folder_name, status_tag)
  stats::setNames(reg$species_key, lbls)
}

build_gene_summary_banner <- function(gm, all_cands, filt_cands, selected_guide) {
  if (is.null(gm)) return(NULL)
  fb_url <- if (!is.na(gm$flybase_id) && grepl("^FBgn", gm$flybase_id)) {
    paste0("https://flybase.org/reports/", gm$flybase_id)
  } else NULL
  ncbi_url <- if (!is.na(gm$ncbi_gene_id) && nzchar(gm$ncbi_gene_id)) {
    paste0("https://www.ncbi.nlm.nih.gov/gene/", sub("^LOC", "", gm$ncbi_gene_id))
  } else NULL
  n_tx <- if (is.data.frame(gm$transcripts)) nrow(gm$transcripts) else length(gm$exons_by_tx)
  n_sh <- if (is.data.frame(gm$shared_cds)) nrow(gm$shared_cds) else 0L
  sh_bp <- if (n_sh > 0L) sum(gm$shared_cds$end - gm$shared_cds$start + 1L) else 0L

  div(class = "cas9-card cas9-banner",
    div(style = "display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px;",
      div(
        span(style = "font-size:20px; font-weight:700; color:#0f172a;", gm$gene_symbol),
        span(style = "margin-left:8px; font-style:italic; color:#475569; font-size:14px;",
             paste0("Drosophila ", gm$species)),
        if (!is.na(gm$flybase_id) && nzchar(gm$flybase_id))
          span(class = "cas9-badge cas9-badge-sky", style = "margin-left:8px;", gm$flybase_id),
        if (!is.na(gm$cg_number) && nzchar(gm$cg_number))
          span(class = "cas9-badge cas9-badge-slate", style = "margin-left:6px;", gm$cg_number),
        if (!is.na(gm$ncbi_gene_id) && nzchar(gm$ncbi_gene_id))
          span(class = "cas9-badge cas9-badge-slate", style = "margin-left:6px;", paste0("NCBI: ", gm$ncbi_gene_id)),
        span(class = "cas9-badge cas9-badge-emerald", style = "margin-left:6px;",
             paste0(nrow(filt_cands), " / ", nrow(all_cands), " gRNAs passing filters"))
      ),
      div(style = "display:flex; gap:8px; align-items:center;",
        if (!is.null(fb_url))
          tags$a(href = fb_url, target = "_blank", class = "btn btn-xs btn-default", "FlyBase Report ↗"),
        if (!is.null(ncbi_url))
          tags$a(href = ncbi_url, target = "_blank", class = "btn btn-xs btn-default", "NCBI Gene ↗")
      )
    ),
    div(style = "margin-top:8px; font-size:12.5px; color:#334155; display:flex; flex-wrap:wrap; gap:16px;",
      div(strong("Locus: "), code(sprintf("%s:%s-%s (%s)",
                                          gm$chr,
                                          format(gm$start, big.mark = ","),
                                          format(gm$end, big.mark = ","),
                                          gm$strand))),
      div(strong("Span: "), paste0(format(gm$end - gm$start + 1L, big.mark = ","), " bp")),
      div(strong("Isoforms: "), paste0(n_tx, " transcript(s)")),
      div(strong("Shared CDS: "), paste0(n_sh, " exon interval(s) (", format(sh_bp, big.mark = ","), " bp)")),
      if (!is.na(gm$product) && nzchar(gm$product))
        div(strong("Product: "), gm$product),
      if (!is.na(gm$synonyms) && nzchar(gm$synonyms))
        div(strong("Synonyms: "), gm$synonyms)
    )
  )
}
build_cut_site_zoom_ui <- function(arm, zoom, sel_guide) {
  if (is.null(arm) || is.null(zoom)) return(NULL)

  qc_badges <- tagList(
    if (isTRUE(arm$qc_pass))
      span(class = "cas9-badge cas9-badge-emerald", "✓ Arm Sequence QC Pass")
    else
      span(class = "cas9-badge cas9-badge-amber", "⚠ Check Sequence QC Warnings"),
    if (isTRUE(arm$recleavage_flag))
      span(class = "cas9-badge cas9-badge-amber", "⚠ Donor Re-Cleavage Risk (Mutate PAM/Seed in Donor)")
    else
      span(class = "cas9-badge cas9-badge-emerald", "✓ Split Protospacer (Low Re-Cleavage Risk)"),
    if (isTRUE(arm$n_tract_flag))
      span(class = "cas9-badge cas9-badge-red", "⚠ Contains N Tract")
    else
      span(class = "cas9-badge cas9-badge-slate", "0 N-tracts"),
    if (isTRUE(arm$boundary_truncated))
      span(class = "cas9-badge cas9-badge-red", "⚠ Truncated at Scaffold Boundary")
    else
      span(class = "cas9-badge cas9-badge-slate", "Full Requested Arm Length")
  )

  duplex_box <- div(
    style = "background:#0f172a; color:#f8fafc; padding:12px 14px; border-radius:6px; font-family:'Courier New',monospace; font-size:12.5px; overflow-x:auto; margin-bottom:12px; line-height:1.55;",
    div(style = "color:#94a3b8; font-size:11px; margin-bottom:4px;",
      sprintf("Target Strand (5' -> 3', %s genomic strand) | Spacer = %s (17 nt 5' of cut + 3 nt 3' of cut) | PAM = %s:",
              arm$strand, arm$spacer_seq, arm$pam_seq)),
    div(
      span(style = "color:#94a3b8;", "5'- "),
      span(style = "color:#60a5fa;", zoom$pre_protospacer),
      span(style = "color:#38bdf8; font-weight:bold; text-decoration:underline;", zoom$spacer_5p_17),
      span(style = "color:#ef4444; font-weight:bold; font-size:15px;", " | "),
      span(style = "color:#fbbf24; font-weight:bold; text-decoration:underline;", zoom$spacer_3p_3),
      span(style = "color:#f43f5e; font-weight:bold; background:#3f1d2e; padding:0 2px; border-radius:2px;", zoom$pam_seq),
      span(style = "color:#34d399;", zoom$post_pam),
      span(style = "color:#94a3b8;", " -3'")
    ),
    div(
      span(style = "color:#94a3b8;", "3'- "),
      span(style = "color:#64748b;", zoom$comp_pre),
      span(style = "color:#cbd5e1;", zoom$comp_spacer_17),
      span(style = "color:#ef4444; font-weight:bold; font-size:15px;", " | "),
      span(style = "color:#cbd5e1;", zoom$comp_spacer_3),
      span(style = "color:#fda4af;", zoom$comp_pam),
      span(style = "color:#64748b;", zoom$comp_post),
      span(style = "color:#94a3b8;", " -5'")
    ),
    div(style = "margin-top:6px; font-size:11px; color:#cbd5e1; display:flex; flex-wrap:wrap; gap:14px;",
      span(style = "color:#60a5fa;", "■ 5' Flanking Arm"),
      span(style = "color:#38bdf8; font-weight:bold;", "■ Protospacer 5' (nt 1–17)"),
      span(style = "color:#ef4444; font-weight:bold;", "| Blunt Cas9 Cut (3 bp 5' of PAM)"),
      span(style = "color:#fbbf24; font-weight:bold;", "■ Protospacer 3' Seed Tail (nt 18–20)"),
      span(style = "color:#f43f5e; font-weight:bold;", "■ NGG PAM"),
      span(style = "color:#34d399;", "■ 3' Flanking Arm")
    )
  )

  div(class = "cas9-card",
    div(style = "display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; margin-bottom:8px;",
      h4(style = "margin:0; color:#0f172a; font-weight:700;",
         paste0("Nucleotide Cut-Site Zoom & Homology Arms — ", arm$guide_id)),
      div(qc_badges)
    ),
    p(style = "font-size:12px; color:#475569; margin-bottom:8px;",
      "SpCas9 induces a blunt double-strand break exactly ",
      strong("3 bp upstream (5') of the first nucleotide of the NGG PAM"),
      " on the target strand (between spacer positions 17 and 18, marked by ",
      span(style = "color:#dc2626; font-weight:700;", "|"), "). ",
      sprintf("Genomic cut coordinates on %s (%s strand): between %s and %s (1-based cut_site = %s).",
              arm$chr, arm$strand,
              format(arm$cut_5p_genomic, big.mark = ","),
              format(arm$cut_3p_genomic, big.mark = ","),
              format(arm$cut_site, big.mark = ","))),
    duplex_box,
    div(style = "background:#fffbeb; border:1px solid #fde68a; border-radius:6px; padding:8px 12px; margin-bottom:12px; font-size:12px; color:#92400e;",
      strong("Donor Re-Cleavage & PAM Shielding Recommendation: "),
      arm$pam_shield_mutation
    ),
    fluidRow(
      column(6,
        div(style = "background:#eff6ff; border:1px solid #bfdbfe; border-radius:6px; padding:10px 12px;",
          div(style = "font-weight:700; color:#1e40af; font-size:13px; margin-bottom:4px;",
              sprintf("5' Homology Arm (%s bp, GC: %s%%)",
                      format(arm$arm_5p_length, big.mark = ","), arm$arm_5p_gc)),
          div(style = "font-size:12px; color:#1e3a8a; margin-bottom:6px;",
              code(sprintf("%s:%s-%s (strand %s)",
                           arm$chr,
                           format(arm$arm_5p_start, big.mark = ","),
                           format(arm$arm_5p_end, big.mark = ","),
                           arm$strand))),
          tags$label(style = "font-size:11px; color:#475569;", "5' Arm Sequence (5' -> 3' on target strand):"),
          tags$textarea(readonly = "readonly", rows = 4,
                        style = "width:100%; font-family:'Courier New',monospace; font-size:11px; background:#ffffff; border:1px solid #cbd5e1; border-radius:4px; padding:6px;",
                        arm$arm_5p_seq)
        )
      ),
      column(6,
        div(style = "background:#fffbeb; border:1px solid #fde68a; border-radius:6px; padding:10px 12px;",
          div(style = "font-weight:700; color:#b45309; font-size:13px; margin-bottom:4px;",
              sprintf("3' Homology Arm (%s bp, GC: %s%%, anchor: %s)",
                      format(arm$arm_3p_length, big.mark = ","), arm$arm_3p_gc, arm$arm_3p_anchor)),
          div(style = "font-size:12px; color:#78350f; margin-bottom:6px;",
              code(sprintf("%s:%s-%s (strand %s)",
                           arm$chr,
                           format(arm$arm_3p_start, big.mark = ","),
                           format(arm$arm_3p_end, big.mark = ","),
                           arm$strand))),
          tags$label(style = "font-size:11px; color:#475569;", "3' Arm Sequence (5' -> 3' on target strand):"),
          tags$textarea(readonly = "readonly", rows = 4,
                        style = "width:100%; font-family:'Courier New',monospace; font-size:11px; background:#ffffff; border:1px solid #cbd5e1; border-radius:4px; padding:6px;",
                        arm$arm_3p_seq)
        )
      )
    )
  )
}

CAS9_CUSTOM_CSS <- "
  body { background-color: #f8fafc; color: #0f172a; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
  .well { background-color: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; box-shadow: 0 1px 2px rgba(15,23,42,0.04); }
  .cas9-card { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 16px; margin-bottom: 14px; box-shadow: 0 1px 2px rgba(15,23,42,0.04); }
  .cas9-banner { border-left: 4px solid #0ea5e9; }
  .cas9-section-hdr { font-size: 13px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em; color: #0284c7; margin-top: 8px; margin-bottom: 8px; }
  .cas9-badge { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px; font-weight: 600; }
  .cas9-badge-sky { background: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd; }
  .cas9-badge-emerald { background: #d1fae5; color: #047857; border: 1px solid #a7f3d0; }
  .cas9-badge-amber { background: #fef3c7; color: #b45309; border: 1px solid #fde68a; }
  .cas9-badge-red { background: #fee2e2; color: #b91c1c; border: 1px solid #fecaca; }
  .cas9-badge-slate { background: #f1f5f9; color: #334155; border: 1px solid #cbd5e1; }
  .cas9-dl-bar { background: #f0f9ff; border: 1px solid #bae6fd; border-radius: 8px; padding: 10px 14px; margin-bottom: 14px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 8px; }
"

ui <- fluidPage(
  tags$head(tags$style(HTML(CAS9_CUSTOM_CSS))),
  div(style = "padding: 12px 4px 6px 4px; border-bottom: 1px solid #e2e8f0; margin-bottom: 14px; display:flex; justify-content:space-between; align-items:baseline; flex-wrap:wrap;",
    div(
      span(style = "font-size:24px; font-weight:800; color:#0f172a;", "CRISPR Cas9 gRNA & Homology Arm Designer"),
      span(class = "cas9-badge cas9-badge-sky", style = "margin-left:10px;", paste0("v", CAS9_APP_VERSION, " · Drosophila Edition"))
    ),
    div(style = "font-size:12px; color:#64748b;",
      "SpCas9 (20-nt + NGG) · Rule Set 2 / Azimuth + CRISPRscan · CFD + MIT Off-Target · 1–2 kb HDR Arms"
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      # 1. Species & Reference Genome Management
      div(class = "cas9-section-hdr", "1. Drosophila Species & Reference"),
      selectInput(
        "species",
        "Target Drosophila species",
        choices  = species_choices_labeled(),
        selected = "melanogaster"
      ),
      uiOutput("species_status_box"),
      div(style = "margin-bottom:8px;",
        actionButton(
          "dl_prepare_species_btn",
          "Download / Prepare Selected Species Genome (NCBI RefSeq)",
          class = "btn-default btn-xs",
          style = "width:100%;"
        )
      ),

      hr(style = "margin:10px 0;"),

      # 2. Gene Lookup Search Bar
      div(class = "cas9-section-hdr", "2. Gene Lookup or Genomic Locus"),
      textInput(
        "gene_query",
        "Gene Symbol, FBgn ID, CG #, LOC ID, Product, or chr:start-end",
        value       = "ebony",
        placeholder = "e.g. ebony, white, yellow, Abd-B, FBgn0000527, 3R:21231000-21236000"
      ),
      div(style = "margin-top:-6px; margin-bottom:8px; display:flex; flex-wrap:wrap; gap:4px; align-items:center;",
        span(style = "font-size:11px; color:#64748b; margin-right:4px;", "Presets:"),
        actionButton("preset_ebony",  "ebony",  class = "btn-default btn-xs"),
        actionButton("preset_white",  "white",  class = "btn-default btn-xs"),
        actionButton("preset_yellow", "yellow", class = "btn-default btn-xs"),
        actionButton("preset_abdb",   "Abd-B",  class = "btn-default btn-xs"),
        actionButton("preset_plexb",  "PlexB",  class = "btn-default btn-xs"),
        actionButton("search_disambig_btn", "Search Matches", class = "btn-info btn-xs")
      ),
      uiOutput("gene_disambig_ui"),

      hr(style = "margin:10px 0;"),

      # 3. Target Context & Scoring Filters
      div(class = "cas9-section-hdr", "3. Target Context & Scoring Filters"),
      radioButtons(
        "target_mode",
        "Target Region Priority Mode",
        choices = c(
          "Coding Disruption (5'-most shared CDS exons prioritized)" = "coding_ko",
          "5' UTR / Promoter Prioritized (regulatory / knock-in)"    = "promoter_5utr",
          "All Genomic Regions (unweighted context)"                 = "all"
        ),
        selected = "coding_ko"
      ),
      checkboxGroupInput(
        "allowed_tiers",
        "Include Genomic Context Tiers",
        choices = c(
          "Tier 1: Early Shared CDS (5'-most shared coding exons)" = "1",
          "Tier 2: 5' UTR / Upstream Promoter"                     = "2",
          "Tier 3: Late Shared CDS (>50% CDS position)"            = "3",
          "Tier 4: Isoform-Specific CDS"                           = "4",
          "Tier 5: 3' UTR / Intron"                                = "5"
        ),
        selected = c("1", "2", "3", "4")
      ),
      fluidRow(
        column(6,
          sliderInput("min_on_target", "Min On-Target Score", min = 0, max = 100, value = 35, step = 5)
        ),
        column(6,
          sliderInput("min_specificity", "Min CFD Specificity", min = 0, max = 100, value = 50, step = 5)
        )
      ),
      fluidRow(
        column(6,
          selectInput(
            "offtarget_mode",
            "Off-Target Scan Scope",
            choices = c(
              "Fast (Chromosome + Seed)" = "fast",
              "Locus Window (Ultrafast)" = "locus",
              "Full Genome (PDict)"      = "full"
            ),
            selected = "fast"
          )
        ),
        column(6,
          numericInput("upstream_promoter_bp", "Promoter Window (bp)", value = 500, min = 100, max = 5000, step = 100)
        )
      ),
      checkboxInput("exclude_tttt", "Exclude U6 Pol III TTTT terminators", value = TRUE),

      hr(style = "margin:10px 0;"),

      # 4. Homology arm provenance (design lives in the Plasmid Designer)
      div(class = "cas9-section-hdr", "4. Homology Arms \u2014 Export Only"),
      helpText(paste0("Homology-arm design has moved to the Plasmid Designer app ",
                      "(port 7790), where arms are sized to the protocol-validated ",
                      "250\u20131000 bp window and assembled into attL donor fragments. ",
                      "This app exports the full genomic arms with each guide ",
                      "(rule R4: arms are kept in genomic orientation); the ",
                      "Plasmid Designer trims them to your chosen length on import.")),
      checkboxInput("orient_5p_to_3p", "Flip minus-strand loci 5' -> 3' in track plot", value = FALSE),

      hr(style = "margin:10px 0;"),

      # 5. Donor Payload (Selectable Marker)
      div(class = "cas9-section-hdr", "5. Donor Payload (Selectable Marker)"),
      textAreaInput(
        "donor_payload_seq",
        "Donor DNA sequence (FASTA or raw)",
        placeholder = "Paste FASTA or raw sequence here...",
        rows = 5
      ),
      fileInput(
        "donor_payload_file",
        "Or upload FASTA",
        accept = c(".fa", ".fasta", ".fa.gz", ".fasta.gz")
      ),
      checkboxInput("donor_payload_orient_minus", "Reverse-complement payload", value = FALSE),
      actionButton(
        "run_donor_btn",
        "Build Donor Construct (5' arm + payload + 3' arm)",
        class = "btn-info btn-lg",
        style = "width:100%; font-weight:700; margin-top:6px;"
      ),

      actionButton(
        "run_design_btn",
        "Run Cas9 gRNA + Homology Arm Design",
        class = "btn-primary btn-lg",
        style = "width:100%; font-weight:700; margin-top:6px;"
      )
    ),

    mainPanel(
      width = 8,
      uiOutput("download_bar_ui"),
      uiOutput("gene_banner_ui"),
      uiOutput("track_plot_card_ui"),
      uiOutput("grna_table_card_ui"),
      uiOutput("cut_site_zoom_card_ui"),
      uiOutput("donor_card_ui")
    )
  )
)

# ============================================================================
# 5. Shiny Server Reactive Logic
# ============================================================================

server <- function(input, output, session) {
  ref_state_bump <- reactiveVal(0L)
  workflow_state <- reactiveVal(NULL)
  disambig_hits  <- reactiveVal(NULL)
  selected_gid   <- reactiveVal(NULL)

  # ---- 1. Species Status & One-Click NCBI RefSeq Download/Indexing ----------
  output$species_status_box <- renderUI({
    ref_state_bump()
    sp <- canonical_species(input$species)
    reg <- tryCatch(get_species_registry(), error = function(e) NULL)
    if (is.null(reg) || !(sp %in% reg$species_key)) return(NULL)
    row <- reg[reg$species_key == sp, , drop = FALSE]
    st  <- row$status[1]
    badge_cls <- if (st == "ready") "cas9-badge-emerald" else if (st == "compressed") "cas9-badge-amber" else "cas9-badge-red"
    badge_txt <- if (st == "ready") "Ready (FASTA + FAI + GTF)" else if (st == "compressed") "Compressed (.gz — Ready on Demand)" else "Needs NCBI RefSeq Download"
    div(style = "font-size:11.5px; color:#475569; margin-top:-4px; margin-bottom:6px;",
      span(class = paste("cas9-badge", badge_cls), badge_txt),
      span(style = "margin-left:6px;", code(row$assembly_accession[1]))
    )
  })

  observeEvent(input$dl_prepare_species_btn, {
    sp <- canonical_species(input$species)
    req(sp)
    withProgress(message = paste0("Preparing D. ", sp, " reference genome..."), value = 0.2, {
      tryCatch({
        reg <- get_species_registry()
        st <- reg$status[reg$species_key == sp][1]
        if (identical(st, "needs_download")) {
          incProgress(0.3, detail = "Downloading genomic.fna.gz + genomic.gtf.gz from NCBI RefSeq...")
          download_species_ref(sp, with_genome = TRUE)
        }
        incProgress(0.3, detail = "Ensuring uncompressed FASTA/GTF and building .fai index...")
        ensure_uncompressed_or_indexed(sp, decompress_gz = TRUE)
        ref_state_bump(ref_state_bump() + 1L)
        updateSelectInput(session, "species", choices = species_choices_labeled(), selected = sp)
        showNotification(paste0("D. ", sp, " genome reference is ready and indexed."), type = "message", duration = 6)
      }, error = function(e) {
        showNotification(paste("Species preparation error:", conditionMessage(e)), type = "error", duration = 12)
      })
    })
  })

  # ---- 2. Presets & Gene Lookup Disambiguation ------------------------------
  observeEvent(input$preset_ebony,  { updateTextInput(session, "gene_query", value = "ebony") })
  observeEvent(input$preset_white,  { updateTextInput(session, "gene_query", value = "white") })
  observeEvent(input$preset_yellow, { updateTextInput(session, "gene_query", value = "yellow") })
  observeEvent(input$preset_abdb,   { updateTextInput(session, "gene_query", value = "Abd-B") })
  observeEvent(input$preset_plexb,  { updateTextInput(session, "gene_query", value = "PlexB") })

  observeEvent(input$search_disambig_btn, {
    q  <- trimws(input$gene_query)
    sp <- canonical_species(input$species)
    if (!nzchar(q)) return()
    hits <- tryCatch(search_genes(q, species = sp, max_hits = 8L), error = function(e) NULL)
    disambig_hits(hits)
  })

  output$gene_disambig_ui <- renderUI({
    hits <- disambig_hits()
    if (is.null(hits) || !is.data.frame(hits) || nrow(hits) == 0L) return(NULL)
    div(style = "background:#f8fafc; border:1px solid #cbd5e1; border-radius:6px; padding:6px 8px; margin-bottom:8px; font-size:11px;",
      div(style = "font-weight:700; color:#334155; margin-bottom:4px;",
          paste0("Matching loci (", nrow(hits), " hits):")),
      tags$table(class = "table table-condensed", style = "margin-bottom:0; font-size:11px;",
        tags$thead(tags$tr(tags$th("Symbol"), tags$th("FBgn / ID"), tags$th("Locus"), tags$th("Match"))),
        tags$tbody(lapply(seq_len(nrow(hits)), function(i) {
          r <- hits[i, ]
          tags$tr(
            tags$td(strong(r$gene_symbol)),
            tags$td(code(ifelse(!is.na(r$flybase_id) && nzchar(r$flybase_id), r$flybase_id, r$ncbi_gene_id))),
            tags$td(sprintf("%s:%d-%d(%s)", r$chromosome, r$start, r$end, r$strand)),
            tags$td(r$match_type)
          )
        }))
      )
    )
  })

  # ---- 3. Homology arm defaults (design now lives in the Plasmid Designer) --
  # The app still extracts the full genomic arms so they can be exported with
  # each guide; the protocol-validated window is 250-1000 bp, so we extract at
  # the 1000 bp ceiling and let the Plasmid Designer trim to the user's length.
  CAS9_ARM_DEFAULT_BP <- 1000L

  # ---- 4. Core Design Execution ---------------------------------------------
  execute_design <- function() {
    q  <- trimws(input$gene_query)
    sp <- canonical_species(input$species)
    if (!nzchar(q)) {
      showNotification("Please enter a gene symbol, FlyBase ID, CG number, or coordinate.", type = "warning")
      return()
    }
    arm5 <- CAS9_ARM_DEFAULT_BP
    arm3 <- CAS9_ARM_DEFAULT_BP
    tiers <- if (length(input$allowed_tiers)) as.integer(input$allowed_tiers) else 1:5

    withProgress(message = paste0("Designing SpCas9 gRNAs & Homology Arms for '", q, "' (D. ", sp, ")..."), value = 0.15, {
      res <- tryCatch({
        incProgress(0.25, detail = "Resolving gene & isoform model from GTF...")
        wf <- run_crispr_design_workflow(
          query                = q,
          species              = sp,
          target_mode          = input$target_mode,
          offtarget_mode       = input$offtarget_mode,
          upstream_promoter_bp = as.integer(input$upstream_promoter_bp),
          allowed_tiers        = tiers,
          min_on_target        = as.numeric(input$min_on_target),
          min_specificity      = as.numeric(input$min_specificity),
          exclude_tttt         = isTRUE(input$exclude_tttt),
          arm_5p_len           = arm5,
          arm_3p_len           = arm3,
          arm_3p_anchor        = "cut_site",
          selected_rank        = 1L,
          seed                 = 42L
        )
        incProgress(0.50, detail = "Completed gRNA ranking and homology arm extraction.")
        wf
      }, error = function(e) {
        showNotification(paste("Design error:", conditionMessage(e)), type = "error", duration = 12)
        NULL
      })

      if (!is.null(res)) {
        workflow_state(res)
        if (!is.null(res$selected_guide) && nrow(res$selected_guide) > 0L) {
          selected_gid(res$selected_guide$guide_id[1])
        }
        for (w in res$warnings) showNotification(w, type = "warning", duration = 8)
      }
    })
  }

  observeEvent(input$run_design_btn, { execute_design() })

  # ---- Donor construct (section 5) ----------------------------------------
  pjd_donor_construct <- reactiveVal(NULL)
  donor_status <- reactiveVal(NULL)

  observeEvent(input$run_donor_btn, {
    wf <- workflow_state()
    az <- live_arm_and_zoom()
    if (is.null(wf) || is.null(az)) {
      donor_status(list(ok = FALSE, msg = "Run the Cas9 design first to generate homology arms."))
      return()
    }
    arm <- az$arm
    if (is.null(arm) || is.null(arm$arm_5p_seq) || is.null(arm$arm_3p_seq)) {
      donor_status(list(ok = FALSE, msg = "Homology arms not available for the selected guide."))
      return()
    }

    # Read payload
    payload <- NULL
    if (!is.null(input$donor_payload_seq) && nzchar(trimws(input$donor_payload_seq))) {
      payload <- tryCatch(read_donor_payload(input$donor_payload_seq, is_path = FALSE),
                          error = function(e) {
                            donor_status(list(ok = FALSE, msg = conditionMessage(e))); NULL
                          })
    } else if (!is.null(input$donor_payload_file) && nzchar(input$donor_payload_file$datapath)) {
      payload <- tryCatch(read_donor_payload(input$donor_payload_file$datapath, is_path = TRUE),
                          error = function(e) {
                            donor_status(list(ok = FALSE, msg = conditionMessage(e))); NULL
                          })
    }
    if (is.null(payload)) return()

    orient <- if (isTRUE(input$donor_payload_orient_minus)) "-" else "+"
    dc <- tryCatch(
      assemble_donor_construct(arm, payload, orient = orient,
                               name = paste0(wf$gene_model$gene_symbol[1], "_donor")),
      error = function(e) {
        donor_status(list(ok = FALSE, msg = conditionMessage(e))); NULL
      }
    )
    if (is.null(dc)) return()

    donor_status(list(ok = TRUE, msg = "Donor construct built."))
    pjd_donor_construct(dc)
  })

  # Auto-run initial default design on startup for immediate interactive view
  observeEvent(TRUE, {
    if (is.null(workflow_state())) {
      tryCatch({
        wf <- run_crispr_design_workflow(
          query          = "ebony",
          species        = "melanogaster",
          target_mode    = "coding_ko",
          offtarget_mode = "fast",
          arm_5p_len     = 1000L,
          arm_3p_len     = 1000L,
          seed           = 42L
        )
        workflow_state(wf)
        if (!is.null(wf$selected_guide) && nrow(wf$selected_guide) > 0L) {
          selected_gid(wf$selected_guide$guide_id[1])
        }
      }, error = function(e) NULL)
    }
  }, once = TRUE)

  # ---- 5. Dynamic Reactive Filtering & Live Homology Arm Update -------------
  live_filtered_df <- reactive({
    wf <- workflow_state()
    if (is.null(wf) || is.null(wf$all_candidates)) return(NULL)
    tiers <- if (length(input$allowed_tiers)) as.integer(input$allowed_tiers) else 1:5
    reranked <- rank_grna_candidates(wf$all_candidates, target_mode = input$target_mode)
    filt <- filter_grna_candidates(
      candidates_df   = reranked,
      allowed_tiers   = tiers,
      min_on_target   = as.numeric(input$min_on_target),
      min_specificity = as.numeric(input$min_specificity),
      exclude_tttt    = isTRUE(input$exclude_tttt)
    )
    if (nrow(filt) == 0L && nrow(reranked) > 0L) {
      filt <- reranked
      filt$filtered_rank <- seq_len(nrow(filt))
    }
    filt
  })

  live_selected_row <- reactive({
    df <- live_filtered_df()
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    gid <- selected_gid()
    if (!is.null(gid) && nzchar(gid) && (gid %in% df$guide_id)) {
      return(df[match(gid, df$guide_id), , drop = FALSE])
    }
    df[1L, , drop = FALSE]
  })

  live_arm_and_zoom <- reactive({
    wf  <- workflow_state()
    row <- live_selected_row()
    if (is.null(wf) || is.null(row) || nrow(row) == 0L) return(NULL)
    arm5 <- CAS9_ARM_DEFAULT_BP
    arm3 <- CAS9_ARM_DEFAULT_BP
    arm <- safe_design_homology_arms(
      guide_row     = row,
      species       = wf$species,
      arm_5p_len    = arm5,
      arm_3p_len    = arm3,
      arm_3p_anchor = "cut_site",
      gene_model    = wf$gene_model
    )
    zm <- compute_cut_site_zoom_view(arm, flank_bp = 30L)
    list(arm = arm, zoom = zm)
  })

  # Sync selection from dropdown or DT table click
  observeEvent(input$selected_guide_picker, {
    if (!is.null(input$selected_guide_picker) && nzchar(input$selected_guide_picker)) {
      selected_gid(input$selected_guide_picker)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$grna_dt_rows_selected, {
    df <- live_filtered_df()
    idx <- input$grna_dt_rows_selected
    if (!is.null(df) && length(idx) == 1L && idx >= 1L && idx <= nrow(df)) {
      gid <- df$guide_id[idx]
      selected_gid(gid)
      updateSelectInput(session, "selected_guide_picker", selected = gid)
    }
  }, ignoreInit = TRUE)

  # ---- 6. Main Panel Rendered Cards & Outputs -------------------------------
  output$download_bar_ui <- renderUI({
    wf <- workflow_state()
    if (is.null(wf)) return(NULL)
    div(class = "cas9-dl-bar",
      div(
        strong(style = "color:#0369a1; font-size:13px;", "Reproducible Export Bundle: "),
        span(style = "font-size:12px; color:#475569;",
             "All files include exact genomic coordinates, scoring versions, and RNG seed.")
      ),
      div(style = "display:flex; flex-wrap:wrap; gap:6px;",
        downloadButton("dl_zip",   "Download All (ZIP)",       class = "btn-primary btn-sm"),
        downloadButton("dl_fasta", "Homology Arms (.fasta)",   class = "btn-default btn-sm"),
        downloadButton("dl_csv",   "Ranked gRNAs (.csv)",      class = "btn-default btn-sm"),
        downloadButton("dl_json",  "Design Summary (.json)",   class = "btn-default btn-sm"),
        downloadButton("dl_log",   "Reproducibility Log (.txt)", class = "btn-default btn-sm"),
        downloadButton("dl_pd_bundle", "Export guides for Plasmid Designer",
                       class = "btn-success btn-sm")
      )
    )
  })

  output$gene_banner_ui <- renderUI({
    wf <- workflow_state()
    if (is.null(wf)) return(NULL)
    build_gene_summary_banner(wf$gene_model, wf$all_candidates, live_filtered_df(), live_selected_row())
  })

  output$track_plot_card_ui <- renderUI({
    wf <- workflow_state()
    if (is.null(wf)) return(NULL)
    n_tx <- if (is.data.frame(wf$gene_model$transcripts)) nrow(wf$gene_model$transcripts) else length(wf$gene_model$exons_by_tx)
    plot_h <- max(280L, min(640L, 140L + 32L * (n_tx + 4L)))
    div(class = "cas9-card",
      div(style = "display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; margin-bottom:6px;",
        h4(style = "margin:0; font-weight:700; color:#0f172a;",
           "Isoform Architecture, Candidate Cut Sites & Flanking Homology Arms"),
        div(style = "font-size:11.5px; color:#475569; display:flex; gap:10px; flex-wrap:wrap;",
          span(style = "color:#0ea5e9; font-weight:600;", "■ Isoform CDS"),
          span(style = "color:#94a3b8; font-weight:600;", "■ UTR"),
          span(style = "color:#10b981; font-weight:600;", "■ Shared CDS"),
          span(style = "color:#8b5cf6; font-weight:600;", "■ 5' UTR / Promoter"),
          span(style = "color:#dc2626; font-weight:700;", "┆ Selected Cut (3 bp 5' of PAM)"),
          span(style = "color:#2563eb; font-weight:600;", "■ 5' Arm"),
          span(style = "color:#d97706; font-weight:600;", "■ 3' Arm")
        )
      ),
      plotOutput("isoform_track_plot", height = paste0(plot_h, "px"))
    )
  })

  output$isoform_track_plot <- renderPlot({
    wf <- workflow_state()
    req(wf)
    az <- live_arm_and_zoom()
    plot_crispr_isoform_tracks(
      gene_model          = wf$gene_model,
      candidates_df       = live_filtered_df(),
      selected_guide      = live_selected_row(),
      arm_design          = if (!is.null(az)) az$arm else NULL,
      orient_5p_to_3p     = isTRUE(input$orient_5p_to_3p),
      show_all_candidates = TRUE
    )
  }, res = 96)

  output$grna_table_card_ui <- renderUI({
    df <- live_filtered_df()
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    sel_row <- live_selected_row()
    cur_gid <- if (!is.null(sel_row) && nrow(sel_row) > 0L) sel_row$guide_id[1] else df$guide_id[1]
    top_df  <- utils::head(df, 100L)
    picker_choices <- stats::setNames(
      top_df$guide_id,
      sprintf("#%d | %s | %s (%s) | Cut %s:%d | %s | Score %.1f | CFD %.1f",
              top_df$rank, top_df$guide_id, top_df$spacer_seq, top_df$pam_seq,
              top_df$chr, top_df$cut_site, top_df$context_label,
              top_df$composite_score, top_df$cfd_specificity_score)
    )

    div(class = "cas9-card",
      div(style = "display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; margin-bottom:8px; gap:10px;",
        div(
          h4(style = "margin:0; font-weight:700; color:#0f172a;",
             paste0("Ranked SpCas9 gRNA Candidates (", nrow(df), " candidates)")),
          span(style = "font-size:12px; color:#64748b;",
               "Click any row in the table or use the selector to inspect its cut site and 5'/3' homology arms.")
        ),
        div(style = "min-width:340px;",
          selectInput("selected_guide_picker", "Active gRNA for Homology Arms & Zoom:",
                      choices = picker_choices, selected = cur_gid, width = "100%")
        )
      ),
      if (requireNamespace("DT", quietly = TRUE)) {
        DT::DTOutput("grna_dt")
      } else {
        tableOutput("grna_fallback_table")
      }
    )
  })

  format_display_table <- function(df) {
    if (is.null(df) || nrow(df) == 0L) return(data.frame())
    flags <- ifelse(df$u6_TTTT_flag, "TTTT;", "")
    flags <- paste0(flags, ifelse(df$homopolymer_flag, "Homopoly;", ""))
    flags <- paste0(flags, ifelse(df$splice_junction_flag, "SpliceJunc;", ""))
    flags <- ifelse(nzchar(flags), sub(";$", "", flags), "OK")
    data.frame(
      Rank        = df$rank,
      Guide_ID    = df$guide_id,
      Spacer_20nt = df$spacer_seq,
      PAM         = df$pam_seq,
      Strand      = df$strand,
      Cut_Site    = df$cut_site,
      Tier        = paste0("T", df$context_tier, ": ", df$context_label),
      Isoforms    = sprintf("%d/%d (%.0f%%)", df$n_isoforms_hit, df$n_coding_isoforms, 100 * df$isoform_fraction),
      CDS_Pos     = ifelse(is.na(df$cds_rel_pos), "-", sprintf("%.2f", df$cds_rel_pos)),
      Composite   = round(df$composite_score, 1),
      RS2_Azimuth = round(df$rs2_azimuth_score, 1),
      CRISPRscan  = round(df$crisprscan_score, 1),
      CFD_Spec    = round(df$cfd_specificity_score, 1),
      MIT_Spec    = round(df$mit_specificity_score, 1),
      OffTargets  = sprintf("%d/%d/%d/%d", df$n_mm0, df$n_mm1, df$n_mm2, df$n_mm3),
      GC_pct      = round(df$gc_percent, 0),
      Flags       = flags,
      stringsAsFactors = FALSE
    )
  }

  if (requireNamespace("DT", quietly = TRUE)) {
    output$grna_dt <- DT::renderDT({
      df <- live_filtered_df()
      req(df)
      disp <- format_display_table(df)
      sel_row <- live_selected_row()
      sel_idx <- if (!is.null(sel_row) && nrow(sel_row) > 0L) match(sel_row$guide_id[1], df$guide_id) else 1L
      if (is.na(sel_idx)) sel_idx <- 1L
      DT::datatable(
        disp,
        selection = list(mode = "single", selected = sel_idx, target = "row"),
        rownames  = FALSE,
        options   = list(pageLength = 10, scrollX = TRUE, dom = "ftip")
      )
    }, server = FALSE)
  }

  output$grna_fallback_table <- renderTable({
    df <- live_filtered_df()
    req(df)
    utils::head(format_display_table(df), 25L)
  }, striped = TRUE, bordered = TRUE, spacing = "xs")

  output$cut_site_zoom_card_ui <- renderUI({
    az <- live_arm_and_zoom()
    if (is.null(az)) return(NULL)
    build_cut_site_zoom_ui(az$arm, az$zoom, live_selected_row())
  })

  output$donor_card_ui <- renderUI({
    dc <- pjd_donor_construct()
    st <- donor_status()
    if (is.null(dc) && is.null(st)) return(NULL)
    tagList(
      div(class = "cas9-card",
          h4("Donor Construct (5' arm + payload + 3' arm)", style = "margin-top:0;"),
          if (!is.null(st) && !isTRUE(st$ok)) {
            div(class = "cas9-badge cas9-badge-red", st$msg)
          } else {
            div(style = "font-size:12px;", "Payload: ", dc$qc$len, " bp | GC: ", round(dc$qc$gc_pct, 1), "% | Max homopolymer: ", dc$qc$max_homopolymer, " | Inserted at donor positions ", dc$payload_start, " to ", dc$payload_end)
            div(style = "font-size:12px; color:#475569; margin-top:4px;",
                "Junction left: ", dc$junction_left, " | Junction right: ", dc$junction_right)
            div(style = "font-size:12px; margin-top:6px;",
                paste(if (isTRUE(dc$qc$ok)) "Donor QC: PASS" else "Donor QC: REVIEW",
                      if (length(dc$qc$warnings)) paste("Warnings:", paste(dc$qc$warnings, collapse = "; "))
                      else ""))
          }
      )
    )
  })
  # ---- 7. Reproducible Download Handlers ------------------------------------
  file_prefix <- reactive({
    wf <- workflow_state()
    sym <- if (!is.null(wf) && !is.null(wf$gene_model$gene_symbol)) wf$gene_model$gene_symbol[1] else "locus"
    sp  <- if (!is.null(wf)) wf$species else "melanogaster"
    gsub("[^A-Za-z0-9_\\-]", "_", paste0("Cas9_", sym, "_", sp))
  })

  output$dl_fasta <- downloadHandler(
    filename = function() paste0(file_prefix(), "_homology_arms.fasta"),
    content  = function(file) {
      az <- live_arm_and_zoom()
      req(az)
      safe_export_fasta(az$arm, file = file, include_plus_strand = TRUE)
    }
  )

  output$dl_csv <- downloadHandler(
    filename = function() paste0(file_prefix(), "_gRNA_candidates.csv"),
    content  = function(file) {
      df <- live_filtered_df()
      az <- live_arm_and_zoom()
      req(df)
      safe_export_csv(df, arm_design = if (!is.null(az)) az$arm else NULL, file = file)
    }
  )

  output$dl_json <- downloadHandler(
    filename = function() paste0(file_prefix(), "_design_summary.json"),
    content  = function(file) {
      wf <- workflow_state()
      df <- live_filtered_df()
      az <- live_arm_and_zoom()
      req(wf, df)
      safe_export_json(df, arm_design = if (!is.null(az)) az$arm else NULL, metadata = wf$metadata, file = file)
    }
  )

  output$dl_log <- downloadHandler(
    filename = function() paste0(file_prefix(), "_CRISPR_design_log.txt"),
    content  = function(file) {
      wf <- workflow_state()
      df <- live_filtered_df()
      az <- live_arm_and_zoom()
      req(wf, df)
      safe_export_log(df, arm_design = if (!is.null(az)) az$arm else NULL, metadata = wf$metadata, file = file)
    }
  )

  output$dl_pd_bundle <- downloadHandler(
    filename = function() {
      az <- live_arm_and_zoom()
      gid <- if (!is.null(az) && !is.null(az$arm)) az$arm$guide_id else "selected"
      paste0(file_prefix(), "_", gid, "_guides_for_plasmid_designer.csv")
    },
    content = function(file) {
      wf <- workflow_state()
      df <- live_filtered_df()
      az <- live_arm_and_zoom()
      validate(need(!is.null(wf),  "Run the design workflow first."),
               need(!is.null(df) && nrow(df) > 0L, "No ranked guides available."),
               need(!is.null(az) && !is.null(az$arm),
                    "Design homology arms for a selected guide first."))
      write_plasmid_designer_bundle(
        grna_df     = df,
        arm_design  = az$arm,
        metadata    = wf$metadata,
        gene_symbol = wf$gene_model$gene_symbol,
        app_version = CAS9_APP_VERSION,
        file        = file
      )
    }
  )

  output$dl_zip <- downloadHandler(
    filename = function() paste0(file_prefix(), "_CRISPR_bundle.zip"),
    content  = function(file) {
      wf <- workflow_state()
      df <- live_filtered_df()
      az <- live_arm_and_zoom()
      req(wf, df)
      plot_fn <- function() {
        plot_crispr_isoform_tracks(
          gene_model          = wf$gene_model,
          candidates_df       = df,
          selected_guide      = live_selected_row(),
          arm_design          = if (!is.null(az)) az$arm else NULL,
          orient_5p_to_3p     = isTRUE(input$orient_5p_to_3p),
          show_all_candidates = TRUE
        )
      }
      safe_export_zip(
        zip_path   = file,
        grna_df    = df,
        arm_design = if (!is.null(az)) az$arm else NULL,
        gene_model = wf$gene_model,
        metadata   = wf$metadata,
        plot_fn    = plot_fn
      )
    }
  )

  # ---- Donor payload exports (section 5) ------------------------------------
  output$dl_donor_fasta <- downloadHandler(
    filename = function() paste0(file_prefix(), "_donor_payload.fasta"),
    content = function(file) {
      dc <- pjd_donor_construct()
      req(dc)
      write_donor_fasta(dc, name = file_prefix(), file = file)
    }
  )

  output$dl_donor_gb <- downloadHandler(
    filename = function() paste0(file_prefix(), "_donor_payload.gb"),
    content = function(file) {
      dc <- pjd_donor_construct()
      req(dc)
      write_donor_genbank(dc, name = file_prefix(), file = file)
    }
  )

  output$dl_donor_csv <- downloadHandler(
    filename = function() paste0(file_prefix(), "_donor_payload.csv"),
    content = function(file) {
      dc <- pjd_donor_construct()
      req(dc)
      out <- data.frame(
        construct = file_prefix(),
        arm_5p_length = nchar(dc$left_arm),
        arm_5p_seq = dc$left_arm,
        arm_5p_coords = paste0("chr1:", 1L, "-", nchar(dc$left_arm)),
        payload_len = nchar(dc$payload),
        payload_seq = dc$payload,
        payload_start = dc$payload_start,
        payload_end = dc$payload_end,
        payload_orient = dc$orient,
        qc_len = dc$qc$len,
        qc_gc_pct = round(dc$qc$gc_pct, 1),
        qc_max_homopolymer = dc$qc$max_homopolymer,
        sequence_provenance = dc$sequence_provenance,
        stringsAsFactors = FALSE
      )
      utils::write.csv(out, file = file, row.names = FALSE)
    }
  )

  output$dl_donor_json <- downloadHandler(
    filename = function() paste0(file_prefix(), "_donor_payload.json"),
    content = function(file) {
      dc <- pjd_donor_construct()
      req(dc)
      payload <- list(
        schema = "donor_payload/1.0",
        app_version = CAS9_APP_VERSION,
        name = file_prefix(),
        left_arm = list(
          length = nchar(dc$left_arm),
          sequence = dc$left_arm
        ),
        payload = list(
          length = nchar(dc$payload),
          sequence = dc$payload,
          start = dc$payload_start,
          end = dc$payload_end,
          orient = dc$orient,
          sequence_provenance = dc$sequence_provenance
        ),
        right_arm = list(
          length = nchar(dc$right_arm),
          sequence = dc$right_arm
        ),
        qc = list(
          ok = dc$qc$ok,
          warnings = dc$qc$warnings,
          len = dc$qc$len,
          gc_pct = round(dc$qc$gc_pct, 2),
          max_homopolymer = dc$qc$max_homopolymer
        ),
        sequence_rule_r6_share = dc$seq_rule_r6_share,
        seq_rule_synthesis_ok = dc$qc$ok,
        sequence_provenance = dc$sequence_provenance
      )
      json_txt <- jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE)
      writeLines(json_txt, file)
    }
  )

  output$dl_donor_zip <- downloadHandler(
    filename = function() paste0(file_prefix(), "_donor_zip.zip"),
    content = function(file) {
      dc <- pjd_donor_construct()
      req(dc)
      tmp <- tempdir()
      fa <- file.path(tmp, "donor_payload.fasta")
      gb <- file.path(tmp, "donor_payload.gb")
      js <- file.path(tmp, "donor_payload.json")
      write_donor_fasta(dc, name = file_prefix(), file = fa)
      write_donor_genbank(dc, name = file_prefix(), file = gb)
      writeLines(jsonlite::toJSON(list(
        schema = "donor_payload/1.0",
        app_version = CAS9_APP_VERSION,
        left_arm = list(length = nchar(dc$left_arm), sequence = dc$left_arm),
        payload = list(
          length = nchar(dc$payload), sequence = dc$payload,
          start = dc$payload_start, end = dc$payload_end,
          orient = dc$orient,
          sequence_provenance = dc$sequence_provenance
        ),
        right_arm = list(length = nchar(dc$right_arm), sequence = dc$right_arm),
        qc = list(
          ok = dc$qc$ok, warnings = dc$qc$warnings,
          len = dc$qc$len, gc_pct = round(dc$qc$gc_pct, 2),
          max_homopolymer = dc$qc$max_homopolymer
        ),
        sequence_rule_r6_share = dc$seq_rule_r6_share,
        sequence_provenance = dc$sequence_provenance
      ), pretty = TRUE, auto_unbox = TRUE), js)
      utils::zip(zipfile = file, files = c(fa, gb, js), flags = "-q")
    }
  )
}

# ============================================================================
# 6. Headless CLI Batch Mode & Interactive Shiny Entry Point
#    (Modeled on HCR_probe_design_v41.R: skipped when source()d in tests)
# ============================================================================

parse_cli_flag <- function(args, flag, default = NULL) {
  eq_pat <- paste0("^", flag, "=")
  hit_eq <- grep(eq_pat, args, value = TRUE)
  if (length(hit_eq)) return(sub(eq_pat, "", hit_eq[1]))
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1L])
  default
}

if (sys.nframe() == 0L) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  gene_arg <- parse_cli_flag(cli_args, "--gene", NULL)
  if (is.null(gene_arg) && length(cli_args) > 0L && !grepl("^--", cli_args[1])) {
    gene_arg <- cli_args[1]
  }

  if (!is.null(gene_arg) && nzchar(gene_arg)) {
    sp_arg     <- parse_cli_flag(cli_args, "--species",   "melanogaster")
    mode_arg   <- parse_cli_flag(cli_args, "--mode",      "coding_ko")
    ot_arg     <- parse_cli_flag(cli_args, "--offtarget", "fast")
    arm5_arg   <- as.integer(parse_cli_flag(cli_args, "--arm5", "1000"))
    arm3_arg   <- as.integer(parse_cli_flag(cli_args, "--arm3", as.character(arm5_arg)))
    anchor_arg <- parse_cli_flag(cli_args, "--anchor",    "cut_site")
    rank_arg   <- as.integer(parse_cli_flag(cli_args, "--rank", "1"))
    seed_arg   <- as.integer(parse_cli_flag(cli_args, "--seed", "42"))
    out_arg    <- parse_cli_flag(cli_args, "--out",       paste0("crispr_out_", gsub("[^A-Za-z0-9_\\-]", "_", gene_arg)))

    cat("\n============================================================\n")
    cat("  CRISPR Cas9 gRNA & Homology Arm Designer (Headless CLI)\n")
    cat("============================================================\n")
    cat("  Query   :", gene_arg, "\n")
    cat("  Species :", sp_arg, "\n")
    cat("  Mode    :", mode_arg, "| Off-target:", ot_arg, "\n")
    cat("  Arms    : 5' =", arm5_arg, "bp | 3' =", arm3_arg, "bp (anchor:", anchor_arg, ")\n")

    wf <- run_crispr_design_workflow(
      query          = gene_arg,
      species        = sp_arg,
      target_mode    = mode_arg,
      offtarget_mode = ot_arg,
      arm_5p_len     = arm5_arg,
      arm_3p_len     = arm3_arg,
      arm_3p_anchor  = anchor_arg,
      selected_rank  = rank_arg,
      seed           = seed_arg
    )

    bundle <- write_crispr_output_bundle(out_arg, wf)
    cat("  Candidates Found :", nrow(wf$filtered_candidates), "\n")
    if (!is.null(wf$arm_design)) {
      cat("  Selected Guide   :", wf$arm_design$guide_id, "(", wf$arm_design$spacer_seq, wf$arm_design$pam_seq, ")\n")
      cat("  Blunt Cut Site   :", wf$arm_design$chr, ":", wf$arm_design$cut_site, "(strand", wf$arm_design$strand, ")\n")
      cat("  5' Homology Arm  :", wf$arm_design$chr, ":", wf$arm_design$arm_5p_start, "-", wf$arm_design$arm_5p_end, "\n")
      cat("  3' Homology Arm  :", wf$arm_design$chr, ":", wf$arm_design$arm_3p_start, "-", wf$arm_design$arm_3p_end, "\n")
    }
    cat("  Output Folder    :", bundle$out_dir, "\n")
    cat("============================================================\n\n")
    quit(save = "no", status = 0)
  }

  cat("\n============================================================\n")
  cat("  CRISPR Cas9 gRNA & Homology Arm Designer (v", CAS9_APP_VERSION, ")\n", sep = "")
  cat("  Open your browser to: http://127.0.0.1:", PORT, "\n", sep = "")
  cat("  Reference folder    : ", get_reference_dir(), "\n", sep = "")
  cat("============================================================\n\n")
  shiny::shinyApp(ui, server, options = list(port = PORT, launch.browser = TRUE))
}
