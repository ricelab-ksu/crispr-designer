library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "homology_arms.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "genome_registry.R"))
source(file.path(project_root, "R", "gene_lookup.R"))
source(file.path(project_root, "R", "grna_engine.R"))
source(file.path(project_root, "R", "homology_arms.R"))

test_that("Exact 1-based SpCas9 cut-site and 5'/3' arm geometry holds on '+' and '-' strands across ebony, white, yellow, Abd-B, and PlexB", {
  genes <- c("ebony", "white", "yellow", "Abd-B", "PlexB")

  for (g in genes) {
    grnas <- design_grnas(
      g,
      species = "melanogaster",
      target_mode = "coding_ko",
      offtarget_mode = "locus",
      max_offtarget_candidates = 10L,
      seed = 42L
    )
    expect_gt(nrow(grnas), 10L)
    gm <- attr(grnas, "gene_model")

    plus_rows  <- grnas[grnas$strand == "+", , drop = FALSE]
    minus_rows <- grnas[grnas$strand == "-", , drop = FALSE]
    expect_gt(nrow(plus_rows), 0L)
    expect_gt(nrow(minus_rows), 0L)

    # Test '+' strand candidate
    p_row <- plus_rows[1L, , drop = FALSE]
    arm_p <- design_homology_arms(p_row, species = "melanogaster", arm_5p_len = 1000L, arm_3p_len = 1000L, gene_model = gm)

    expect_equal(arm_p$strand, "+")
    expect_equal(arm_p$cut_5p_genomic, p_row$pam_start - 4L)
    expect_equal(arm_p$cut_3p_genomic, p_row$pam_start - 3L)
    expect_equal(arm_p$cut_site,       p_row$pam_start - 3L)
    expect_equal(arm_p$arm_5p_start,   arm_p$cut_5p_genomic - 1000L + 1L)
    expect_equal(arm_p$arm_5p_end,     arm_p$cut_5p_genomic)
    expect_equal(arm_p$arm_3p_start,   arm_p$cut_3p_genomic)
    expect_equal(arm_p$arm_3p_end,     arm_p$cut_3p_genomic + 1000L - 1L)
    expect_equal(arm_p$arm_5p_end + 1L, arm_p$arm_3p_start)
    expect_equal(nchar(arm_p$arm_5p_seq), 1000L)
    expect_equal(nchar(arm_p$arm_3p_seq), 1000L)
    expect_equal(arm_p$pam_offset_in_arm_3p, 4L)
    expect_equal(arm_p$arm_3p_pam_start, p_row$pam_start)
    expect_equal(substr(arm_p$arm_5p_seq, 1000L - 16L, 1000L), substr(p_row$spacer_seq, 1L, 17L))
    expect_equal(substr(arm_p$arm_3p_seq, 1L, 3L), substr(p_row$spacer_seq, 18L, 20L))
    expect_equal(substr(arm_p$arm_3p_seq, 4L, 6L), p_row$pam_seq)

    cont_plus <- extract_genomic_subseq("melanogaster", arm_p$accession, arm_p$arm_5p_start, arm_p$arm_3p_end, strand = "+")
    expect_identical(paste0(arm_p$arm_5p_seq, arm_p$arm_3p_seq), cont_plus)

    # Test '-' strand candidate
    m_row <- minus_rows[1L, , drop = FALSE]
    arm_m <- design_homology_arms(m_row, species = "melanogaster", arm_5p_len = 1000L, arm_3p_len = 1000L, gene_model = gm)

    expect_equal(arm_m$strand, "-")
    expect_equal(arm_m$cut_5p_genomic, m_row$pam_end + 4L)
    expect_equal(arm_m$cut_3p_genomic, m_row$pam_end + 3L)
    expect_equal(arm_m$cut_site,       m_row$pam_end + 3L)
    expect_equal(arm_m$arm_5p_start,   arm_m$cut_5p_genomic)
    expect_equal(arm_m$arm_5p_end,     arm_m$cut_5p_genomic + 1000L - 1L)
    expect_equal(arm_m$arm_3p_start,   arm_m$cut_3p_genomic - 1000L + 1L)
    expect_equal(arm_m$arm_3p_end,     arm_m$cut_3p_genomic)
    expect_equal(arm_m$arm_3p_end + 1L, arm_m$arm_5p_start)
    expect_equal(nchar(arm_m$arm_5p_seq), 1000L)
    expect_equal(nchar(arm_m$arm_3p_seq), 1000L)
    expect_equal(arm_m$pam_offset_in_arm_3p, 4L)
    expect_equal(substr(arm_m$arm_5p_seq, 1000L - 16L, 1000L), substr(m_row$spacer_seq, 1L, 17L))
    expect_equal(substr(arm_m$arm_3p_seq, 1L, 3L), substr(m_row$spacer_seq, 18L, 20L))
    expect_equal(substr(arm_m$arm_3p_seq, 4L, 6L), m_row$pam_seq)

    cont_minus <- extract_genomic_subseq("melanogaster", arm_m$accession, arm_m$arm_3p_start, arm_m$arm_5p_end, strand = "-")
    expect_identical(paste0(arm_m$arm_5p_seq, arm_m$arm_3p_seq), cont_minus)

    cont_plus_ref <- extract_genomic_subseq("melanogaster", arm_m$accession, arm_m$genomic_left_arm_start, arm_m$genomic_right_arm_end, strand = "+")
    expect_identical(paste0(arm_m$genomic_left_arm_seq_plus, arm_m$genomic_right_arm_seq_plus), cont_plus_ref)
  }
})

test_that("design_homology_arms supports non-melanogaster species (d_simulans) and symmetric/asymmetric slider lengths (1000-2000 bp)", {
  grnas_sim <- design_grnas("yellow", species = "d_simulans", target_mode = "coding_ko",
                            offtarget_mode = "locus", max_offtarget_candidates = 10L, seed = 42L)
  expect_gt(nrow(grnas_sim), 10L)

  length_pairs <- list(
    c(1000L, 1000L),
    c(1500L, 1500L),
    c(1200L, 1800L),
    c(2000L, 2000L)
  )

  for (lp in length_pairs) {
    L <- lp[1L]
    R <- lp[2L]
    arm_sim <- design_homology_arms(grnas_sim[1L, , drop = FALSE], species = "d_simulans",
                                    arm_5p_len = L, arm_3p_len = R)
    expect_equal(arm_sim$species, "simulans")
    expect_equal(nchar(arm_sim$arm_5p_seq), L)
    expect_equal(nchar(arm_sim$arm_3p_seq), R)
    expect_equal(arm_sim$arm_5p_len, L)
    expect_equal(arm_sim$arm_3p_len, R)

    if (arm_sim$strand == "+") {
      cont <- extract_genomic_subseq("simulans", arm_sim$accession, arm_sim$arm_5p_start, arm_sim$arm_3p_end, strand = "+")
    } else {
      cont <- extract_genomic_subseq("simulans", arm_sim$accession, arm_sim$arm_3p_start, arm_sim$arm_5p_end, strand = "-")
    }
    expect_identical(paste0(arm_sim$arm_5p_seq, arm_sim$arm_3p_seq), cont)
  }

  # Also verify arm_3p_anchor = "pam_start" option
  arm_pam <- design_homology_arms(grnas_sim[1L, , drop = FALSE], species = "d_simulans",
                                  arm_5p_len = 1000L, arm_3p_len = 1000L, arm_3p_anchor = "pam_start")
  expect_equal(arm_pam$pam_offset_in_arm_3p, 1L)
  expect_equal(substr(arm_pam$arm_3p_seq, 1L, 3L), grnas_sim$pam_seq[1L])
})

test_that("QC flags detect N tracts, low-complexity repeats, extreme GC, chromosome boundary truncation, and PAM shielding", {
  near_left_params <- list(
    species           = "melanogaster",
    chr               = "X",
    strand            = "+",
    protospacer_start = 200L,
    protospacer_end   = 219L,
    pam_start         = 220L,
    pam_end           = 222L
  )
  arm_trunc5 <- design_homology_arms(near_left_params, arm_5p_len = 1000L, arm_3p_len = 1000L)
  expect_true(arm_trunc5$boundary_truncated_5p)
  expect_false(arm_trunc5$boundary_truncated_3p)
  expect_true(arm_trunc5$boundary_truncated)
  expect_equal(arm_trunc5$arm_5p_start, 1L)
  expect_equal(arm_trunc5$arm_5p_len, 216L)
  expect_equal(arm_trunc5$boundary_distance_bp, 0L)
  expect_false(arm_trunc5$qc_pass)

  x_len <- get_chr_length("melanogaster", "X")
  near_right_params <- list(
    species           = "melanogaster",
    chr               = "X",
    strand            = "+",
    protospacer_start = x_len - 300L,
    protospacer_end   = x_len - 281L,
    pam_start         = x_len - 280L,
    pam_end           = x_len - 278L
  )
  arm_trunc3 <- design_homology_arms(near_right_params, arm_5p_len = 1000L, arm_3p_len = 1000L)
  expect_false(arm_trunc3$boundary_truncated_5p)
  expect_true(arm_trunc3$boundary_truncated_3p)
  expect_equal(arm_trunc3$arm_3p_end, x_len)

  qc_n <- .analyze_single_arm_qc(paste0(paste(rep("ACGT", 100), collapse = ""), "NNNNNNN", paste(rep("ACGT", 100), collapse = "")), 807L)
  expect_equal(qc_n$n_count, 7L)
  expect_equal(qc_n$max_n_run, 7L)
  expect_true(qc_n$has_n_tract)

  qc_hp <- .analyze_single_arm_qc(paste0(paste(rep("ACGT", 100), collapse = ""), "AAAAAAAAAA", paste(rep("ACGT", 100), collapse = "")), 810L)
  expect_gte(qc_hp$max_homopolymer, 10L)
  expect_true(qc_hp$repeat_warning)

  qc_gc <- .analyze_single_arm_qc(paste(rep("GGCC", 250), collapse = ""), 1000L)
  expect_equal(qc_gc$gc_pct, 100.0)
  expect_true(qc_gc$gc_warning)

  grnas_y <- design_grnas("yellow", species = "melanogaster", target_mode = "coding_ko",
                          offtarget_mode = "locus", max_offtarget_candidates = 10L, seed = 42L)
  gm_y <- attr(grnas_y, "gene_model")

  arm_cargo <- design_homology_arms(grnas_y[1L, , drop = FALSE], donor_type = "cargo_insertion", gene_model = gm_y)
  expect_true(arm_cargo$cargo_disrupts_target)
  expect_true(arm_cargo$scarless_recleavage_risk)

  arm_scar <- design_homology_arms(grnas_y[1L, , drop = FALSE], donor_type = "scarless_point_edit", gene_model = gm_y)
  expect_true(arm_scar$recleavage_risk)
  expect_type(arm_scar$pam_shielding_mutation, "list")
  expect_true(nzchar(arm_scar$pam_shielding_summary))
  expect_equal(nchar(arm_scar$shielded_arm_3p_seq), arm_scar$arm_3p_len)
  expect_false(identical(arm_scar$shielded_arm_3p_seq, arm_scar$arm_3p_seq))
  expect_false(.has_spacer_pam_match(paste0(arm_scar$arm_5p_seq, arm_scar$shielded_arm_3p_seq), arm_scar$spacer_seq))
})

test_that("build_cut_site_zoom and FASTA, CSV, JSON, Log, and ZIP exports work cleanly", {
  grnas_w <- design_grnas("white", species = "melanogaster", target_mode = "coding_ko",
                          offtarget_mode = "locus", max_offtarget_candidates = 10L, seed = 42L)
  gm_w <- attr(grnas_w, "gene_model")

  plus_row  <- grnas_w[grnas_w$strand == "+", , drop = FALSE][1L, , drop = FALSE]
  minus_row <- grnas_w[grnas_w$strand == "-", , drop = FALSE][1L, , drop = FALSE]

  arm_p <- design_homology_arms(plus_row,  arm_5p_len = 1000L, arm_3p_len = 1000L, gene_model = gm_w)
  arm_m <- design_homology_arms(minus_row, arm_5p_len = 1200L, arm_3p_len = 1400L, gene_model = gm_w)

  for (arm in list(arm_p, arm_m)) {
    zoom <- build_cut_site_zoom(arm, flank_bp = 30L)
    expect_equal(nchar(zoom$top_strand_5to3), 60L)
    expect_equal(nchar(zoom$bottom_strand_3to5), 60L)
    expect_equal(length(zoom$genomic_positions), 60L)
    expect_equal(zoom$cut_index, 30L)
    expect_equal(zoom$genomic_positions[zoom$cut_index], arm$cut_left)
    expect_equal(zoom$genomic_positions[zoom$cut_index + 1L], arm$cut_right)
    expect_match(zoom$formatted_ascii, "\\|")
    expect_equal(substr(zoom$target_strand_5to3, 30L - 16L, 30L + 3L), arm$spacer_seq)
    expect_equal(substr(zoom$target_strand_5to3, 30L + 4L, 30L + 6L), arm$pam_seq)
  }

  tmp_dir <- tempfile("test_exports_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fa_path  <- file.path(tmp_dir, "arms.fasta")
  csv_path <- file.path(tmp_dir, "design.csv")
  js_path  <- file.path(tmp_dir, "design.json")
  log_path <- file.path(tmp_dir, "design.log")
  zip_path <- file.path(tmp_dir, "bundle.zip")

  fa_txt  <- export_homology_arms_fasta(arm_m, file = fa_path, include_plus_strand = TRUE)
  csv_txt <- export_design_csv(head(grnas_w, 10), arm_design = arm_m, file = csv_path)
  js_txt  <- export_design_json(head(grnas_w, 5), arm_design = arm_m, file = js_path)
  log_txt <- export_design_log(head(grnas_w, 10), arm_design = arm_m, file = log_path)
  zip_out <- export_design_bundle_zip(zip_path, head(grnas_w, 10), arm_design = arm_m, gene_model = gm_w)

  expect_true(file.exists(fa_path))
  expect_match(fa_txt, ">w\\|.*\\|arm_5p\\|")
  expect_match(fa_txt, ">w\\|.*\\|arm_3p\\|")
  expect_match(fa_txt, "genomic_left_arm_plus_strand")

  expect_true(file.exists(csv_path))
  df_read <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  expect_equal(nrow(df_read), 10L)
  expect_true(all(c("arm_5p_start", "arm_5p_end", "arm_3p_start", "arm_3p_end", "arm_5p_seq", "arm_3p_seq") %in% names(df_read)))

  expect_true(file.exists(js_path))
  expect_match(js_txt, "\"homology_arms\"")
  expect_match(js_txt, "\"cut_site_zoom\"")

  expect_true(file.exists(log_path))
  expect_match(log_txt, "Selected Guide & SpCas9 Blunt Cut-Site Geometry")
  expect_match(log_txt, "Homology Arms & Sequence Quality Control")

  expect_true(file.exists(zip_out))
  expect_gt(file.info(zip_out)$size, 1000)
  zip_contents <- utils::unzip(zip_out, list = TRUE)$Name
  expect_true(all(c("grna_candidates.csv", "homology_arms.fasta", "design_bundle.json", "design_reproducibility.log") %in% zip_contents))
})

