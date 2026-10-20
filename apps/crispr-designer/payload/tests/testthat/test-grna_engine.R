library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "grna_engine.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "genome_registry.R"))
source(file.path(project_root, "R", "gene_lookup.R"))
source(file.path(project_root, "R", "grna_engine.R"))

test_that("enumerate_spcas9_sites computes exact + and - strand cut-site geometry and context windows", {
  syn_seq <- paste0(
    "ACGTACGTAC",
    "GAGTCCGAGCAGAAGAAGAA", "GGG",
    "ATCGATCGAT",
    "CCA", "TTCTTCTTCTGCTCGGACTC",
    "ACGTACGTACGTAC"
  )
  expect_equal(nchar(syn_seq), 80L)

  cands <- enumerate_spcas9_sites(
    chr = "synthetic_chr",
    region_start = 1001L,
    dna_seq = syn_seq,
    gene_symbol = "syn"
  )
  expect_true(nrow(cands) >= 2L)

  plus_hit <- cands[cands$strand == "+" & cands$spacer_seq == "GAGTCCGAGCAGAAGAAGAA", ]
  expect_equal(nrow(plus_hit), 1L)
  expect_equal(plus_hit$protospacer_start, 1011L)
  expect_equal(plus_hit$protospacer_end,   1030L)
  expect_equal(plus_hit$pam_start,         1031L)
  expect_equal(plus_hit$pam_end,           1033L)
  expect_equal(plus_hit$pam_seq,           "GGG")
  expect_equal(plus_hit$cut_site,          plus_hit$pam_start - 3L) # 1028
  expect_equal(plus_hit$cut_5p_genomic,    plus_hit$pam_start - 4L) # 1027
  expect_equal(plus_hit$cut_3p_genomic,    plus_hit$pam_start - 3L) # 1028
  expect_equal(plus_hit$cut_left,          1027L)
  expect_equal(plus_hit$cut_right,         1028L)
  expect_equal(nchar(plus_hit$context_30mer), 30L)
  expect_equal(nchar(plus_hit$context_35mer), 35L)
  expect_equal(substr(plus_hit$context_30mer, 5L, 24L), "GAGTCCGAGCAGAAGAAGAA")
  expect_equal(substr(plus_hit$context_30mer, 25L, 27L), "GGG")

  minus_hit <- cands[cands$strand == "-" & cands$spacer_seq == "GAGTCCGAGCAGAAGAAGAA", ]
  expect_equal(nrow(minus_hit), 1L)
  expect_equal(minus_hit$pam_start,         1044L)
  expect_equal(minus_hit$pam_end,           1046L)
  expect_equal(minus_hit$protospacer_start, 1047L)
  expect_equal(minus_hit$protospacer_end,   1066L)
  expect_equal(minus_hit$pam_seq,           "TGG")
  expect_equal(minus_hit$cut_site,          minus_hit$pam_end + 3L) # 1049
  expect_equal(minus_hit$cut_5p_genomic,    minus_hit$pam_end + 4L) # 1050 (5' side on '-' strand)
  expect_equal(minus_hit$cut_3p_genomic,    minus_hit$pam_end + 3L) # 1049 (3' side on '-' strand)
  expect_equal(minus_hit$cut_left,          1049L)
  expect_equal(minus_hit$cut_right,         1050L)
  expect_equal(nchar(minus_hit$context_30mer), 30L)
  expect_equal(nchar(minus_hit$context_35mer), 35L)
  expect_equal(substr(minus_hit$context_30mer, 5L, 24L), "GAGTCCGAGCAGAAGAAGAA")
  expect_equal(substr(minus_hit$context_30mer, 25L, 27L), "TGG")
})

test_that("Doench 2016 Rule Set 2 / Azimuth and CRISPRscan distinguish high- and low-efficacy guides", {
  high_30 <- "ACAGGCTGCGGCATGCCGACAAGGCGGACT"
  low_30  <- "TTTTAAAATTTTAAAATTTTCCTCGGGGAA"

  s_rs2 <- score_rs2_azimuth(c(high_30, low_30))
  expect_length(s_rs2, 2L)
  expect_true(all(s_rs2 >= 0 & s_rs2 <= 100))
  expect_gt(s_rs2[1], 75.0)
  expect_lt(s_rs2[2], 25.0)
  expect_gt(s_rs2[1] - s_rs2[2], 50.0)

  high_35 <- paste0("AC", high_30, "ACT")
  low_35  <- paste0("TT", low_30,  "TTT")
  s_cs <- score_crisprscan(c(high_35, low_35))
  expect_length(s_cs, 2L)
  expect_gt(s_cs[1], 70.0)
  expect_lt(s_cs[2], 30.0)

  comp_df <- score_on_target(c(high_30, low_30), c(high_35, low_35))
  expect_equal(comp_df$u6_TTTT_flag, c(FALSE, TRUE))
  expect_gt(comp_df$on_target_score[1], 75.0)
  expect_lt(comp_df$on_target_score[2], 15.0)
})

test_that("Doench 2016 CFD and Hsu 2013 MIT pairwise and aggregated off-target math is exact", {
  sp <- "GAGTCCGAGCAGAAGAAGAA"

  expect_equal(compute_cfd_pair(sp, sp, "AGG"), 1.0)
  expect_equal(compute_mit_pair(sp, sp, "AGG"), 1.0)

  expect_equal(compute_cfd_pair(sp, sp, "AAG"), 0.259259, tolerance = 1e-6)
  expect_equal(compute_cfd_pair(sp, sp, "ACG"), 0.107142, tolerance = 1e-6)
  expect_equal(compute_cfd_pair(sp, sp, "ATG"), 0.069444, tolerance = 1e-6)
  expect_equal(compute_cfd_pair(sp, sp, "ATT"), 0.0)

  mm_pos1  <- "AAGTCCGAGCAGAAGAAGAA"
  mm_pos20 <- "GAGTCCGAGCAGAAGAAGAG"

  cfd_1  <- compute_cfd_pair(sp, mm_pos1,  "AGG")
  cfd_20 <- compute_cfd_pair(sp, mm_pos20, "AGG")
  expect_gt(cfd_1, 0.75)
  expect_equal(cfd_20, 0.0)

  mit_1  <- compute_mit_pair(sp, mm_pos1,  "AGG")
  mit_20 <- compute_mit_pair(sp, mm_pos20, "AGG")
  expect_equal(mit_1, 1.0)
  expect_equal(mit_20, 1.0 - 0.583, tolerance = 1e-6)

  mm_adj <- "GAGTCCGAGCAGAAGAAGCG"
  mm_spr <- "GAGTCTGAGCAGAAGAAGAG"
  mit_adj <- compute_mit_pair(sp, mm_adj, "AGG")
  mit_spr <- compute_mit_pair(sp, mm_spr, "AGG")
  expect_gt(mit_spr, mit_adj)
  expect_lt(mit_spr, mit_20)

  expect_equal(aggregate_specificity_score(numeric(0)), 100.0)
  expect_equal(aggregate_specificity_score(c(1.0)), 50.0)

  # Test score_off_targets on synthetic reference with 1 exact match + 1 1-mm off-target
  ref_seq <- paste0(
    "ACGTACGTAC", sp, "AGG", "ACGTACGTAC",
    "TTTTGGGGCC", mm_pos1, "AGG", "ACGTACGTAC"
  )
  cands <- enumerate_spcas9_sites(chr = "syn", region_start = 1L, dna_seq = substr(ref_seq, 1L, 43L))
  cands <- score_off_targets(cands, reference_seq = ref_seq)
  hit <- cands[cands$spacer_seq == sp, ]
  expect_equal(hit$n_mm0, 1L)
  expect_equal(hit$n_mm1, 1L)
  expect_lt(hit$cfd_specificity_score, 60.0)
  expect_equal(hit$mit_specificity_score, 50.0)
})

test_that("classify_genomic_context and rank_grna_candidates prioritize coding_ko vs promoter_5utr", {
  # Test on Drosophila melanogaster yellow (y) gene
  res_ko <- design_grnas("yellow", species = "melanogaster", target_mode = "coding_ko",
                         offtarget_mode = "locus", max_offtarget_candidates = 30L)
  res_pr <- design_grnas("yellow", species = "melanogaster", target_mode = "promoter_5utr",
                         offtarget_mode = "locus", max_offtarget_candidates = 30L)

  expect_s3_class(res_ko, "data.frame")
  expect_gt(nrow(res_ko), 100L)
  expect_equal(res_ko$context_tier[1], 1L)
  expect_equal(res_ko$context_label[1], "early_shared_cds")
  expect_false(res_ko$u6_TTTT_flag[1])

  expect_equal(res_pr$context_tier[1], 2L)
  expect_equal(res_pr$context_label[1], "5p_utr_promoter")
  expect_true(res_pr$in_5p_utr[1] || res_pr$in_promoter[1])
})

test_that("Positive controls on Drosophila melanogaster ebony, white, yellow, and Abd-B are 100% deterministic", {
  for (g in c("ebony", "white", "yellow", "Abd-B")) {
    r1 <- design_grnas(g, species = "melanogaster", target_mode = "coding_ko",
                       offtarget_mode = "fast", max_offtarget_candidates = 25L, seed = 42L)
    r2 <- design_grnas(g, species = "melanogaster", target_mode = "coding_ko",
                       offtarget_mode = "fast", max_offtarget_candidates = 25L, seed = 42L)
    expect_identical(r1, r2, info = paste("Determinism check for", g))

    expect_gt(nrow(r1), 50L)
    top5 <- head(r1, 5)
    expect_true(all(top5$context_tier == 1L), info = g)
    expect_true(all(top5$context_label == "early_shared_cds"), info = g)
    expect_true(all(!top5$u6_TTTT_flag), info = g)
    expect_true(all(top5$n_mm0 == 1L), info = g)
    expect_true(all(top5$on_target_score >= 75.0), info = g)
    expect_true(all(top5$cfd_specificity_score >= 70.0), info = g)
    expect_true(all(top5$mit_specificity_score >= 70.0), info = g)

    # Verify cut-site geometry on every top candidate
    for (k in seq_len(nrow(top5))) {
      if (top5$strand[k] == "+") {
        expect_equal(top5$cut_site[k], top5$pam_start[k] - 3L)
        expect_equal(top5$cut_5p_genomic[k], top5$pam_start[k] - 4L)
        expect_equal(top5$cut_3p_genomic[k], top5$pam_start[k] - 3L)
      } else {
        expect_equal(top5$cut_site[k], top5$pam_end[k] + 3L)
        expect_equal(top5$cut_5p_genomic[k], top5$pam_end[k] + 4L)
        expect_equal(top5$cut_3p_genomic[k], top5$pam_end[k] + 3L)
      }
      expect_equal(top5$cut_right[k] - top5$cut_left[k], 1L)
    }

    meta <- attr(r1, "metadata")
    expect_type(meta, "list")
    expect_equal(meta$engine_version, GRNA_ENGINE_VERSION)
    expect_equal(meta$species, "melanogaster")
    expect_match(meta$genome_assembly, "GCF_000001215")
  }
})

test_that("design_grnas supports non-melanogaster Drosophila species (d_simulans)", {
  res_sim <- design_grnas("yellow", species = "d_simulans", target_mode = "coding_ko",
                          offtarget_mode = "locus", max_offtarget_candidates = 15L, seed = 42L)
  expect_s3_class(res_sim, "data.frame")
  expect_gt(nrow(res_sim), 20L)
  expect_equal(res_sim$context_tier[1], 1L)
  expect_equal(res_sim$context_label[1], "early_shared_cds")
  meta_sim <- attr(res_sim, "metadata")
  expect_equal(meta_sim$species, "simulans")
})


