library(testthat)

project_root <- "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design"
source(file.path(project_root, "R", "gene_lookup.R"))

test_that("load_gene_index loads D. melanogaster and non-melanogaster species fast and with standard schema", {
  expected_cols <- c(
    "gene_symbol", "flybase_id", "cg_number", "ncbi_gene_id",
    "chromosome", "start", "end", "strand", "biotype",
    "synonyms", "product"
  )

  idx_dmel <- load_gene_index("melanogaster")
  expect_true(is.data.frame(idx_dmel))
  expect_gt(nrow(idx_dmel), 15000)
  expect_true(all(expected_cols %in% names(idx_dmel)))

  idx_sim <- load_gene_index("d_simulans")
  expect_true(is.data.frame(idx_sim))
  expect_gt(nrow(idx_sim), 12000)
  expect_true(all(expected_cols %in% names(idx_sim)))

  idx_eug <- load_gene_index("eugracilis")
  expect_true(is.data.frame(idx_eug))
  expect_gt(nrow(idx_eug), 12000)
  expect_true(all(expected_cols %in% names(idx_eug)))
})

test_that("search_genes resolves symbols, FBgn IDs, CG numbers, NCBI IDs, synonyms, and coordinates in D. melanogaster", {
  std_cols <- c(
    "gene_symbol", "flybase_id", "cg_number", "ncbi_gene_id",
    "chromosome", "start", "end", "strand", "biotype",
    "synonyms", "product", "match_type"
  )

  # 1. Gene symbols (case-insensitive)
  s_abdb <- search_genes("Abd-B", "melanogaster")
  expect_identical(names(s_abdb), std_cols)
  expect_equal(s_abdb$gene_symbol[1], "Abd-B")
  expect_equal(s_abdb$flybase_id[1], "FBgn0000015")
  expect_equal(s_abdb$cg_number[1], "CG11648")
  expect_equal(s_abdb$match_type[1], "gene_symbol")

  s_abdb_ci <- search_genes("abd-b", "melanogaster")
  expect_equal(s_abdb_ci$gene_symbol[1], "Abd-B")
  expect_equal(s_abdb_ci$flybase_id[1], "FBgn0000015")

  s_e <- search_genes("e", "melanogaster")
  expect_equal(s_e$gene_symbol[1], "e")
  expect_equal(s_e$flybase_id[1], "FBgn0000527")
  expect_equal(s_e$cg_number[1], "CG3331")

  s_w <- search_genes("w", "melanogaster")
  expect_equal(s_w$gene_symbol[1], "w")
  expect_equal(s_w$flybase_id[1], "FBgn0003996")
  expect_equal(s_w$cg_number[1], "CG2759")

  s_y <- search_genes("y", "melanogaster")
  expect_equal(s_y$gene_symbol[1], "y")
  expect_equal(s_y$flybase_id[1], "FBgn0004034")
  expect_equal(s_y$cg_number[1], "CG3757")

  s_plexb <- search_genes("PlexB", "melanogaster")
  expect_equal(s_plexb$gene_symbol[1], "PlexB")
  expect_equal(s_plexb$flybase_id[1], "FBgn0025740")
  expect_equal(s_plexb$cg_number[1], "CG17245")

  s_dsx <- search_genes("dsx", "melanogaster")
  expect_equal(s_dsx$gene_symbol[1], "dsx")
  expect_equal(s_dsx$flybase_id[1], "FBgn0000504")

  # 2. FlyBase IDs
  expect_equal(search_genes("FBgn0000015", "melanogaster")$gene_symbol[1], "Abd-B")
  expect_equal(search_genes("FBgn0000527", "melanogaster")$gene_symbol[1], "e")
  expect_equal(search_genes("FBgn0003996", "melanogaster")$gene_symbol[1], "w")
  expect_equal(search_genes("FBgn0004053", "melanogaster")$gene_symbol[1], "zen")
  expect_equal(search_genes("FBgn0025740", "melanogaster")$gene_symbol[1], "PlexB")
  expect_equal(search_genes("fbgn0025740", "melanogaster")$match_type[1], "flybase_id")

  # 3. CG numbers
  expect_equal(search_genes("CG11648", "melanogaster")$gene_symbol[1], "Abd-B")
  expect_equal(search_genes("CG3331", "melanogaster")$gene_symbol[1], "e")
  expect_equal(search_genes("CG2759", "melanogaster")$gene_symbol[1], "w")
  expect_equal(search_genes("CG3757", "melanogaster")$gene_symbol[1], "y")
  expect_equal(search_genes("CG17245", "melanogaster")$gene_symbol[1], "PlexB")
  expect_equal(search_genes("cg17245", "melanogaster")$match_type[1], "cg_number")

  # 4. NCBI Gene ID / locus ID
  s_ncbi1 <- search_genes("43766", "melanogaster")
  expect_equal(s_ncbi1$gene_symbol[1], "PlexB")
  expect_equal(s_ncbi1$match_type[1], "ncbi_gene_id")

  s_ncbi2 <- search_genes("Dmel_CG17245", "melanogaster")
  expect_equal(s_ncbi2$gene_symbol[1], "PlexB")
  expect_equal(s_ncbi2$match_type[1], "ncbi_gene_id")

  # 5. Synonyms / aliases
  expect_equal(search_genes("ebony", "melanogaster")$gene_symbol[1], "e")
  expect_equal(search_genes("white", "melanogaster")$gene_symbol[1], "w")
  expect_equal(search_genes("yellow", "melanogaster")$gene_symbol[1], "y")
  expect_equal(search_genes("plex B", "melanogaster")$gene_symbol[1], "PlexB")

  # 6. Genomic coordinates ('chr:start-end', 'chr:start..end:strand')
  s_coord1 <- search_genes("3R:16800000-16815000", "melanogaster")
  expect_true("abd-A" %in% s_coord1$gene_symbol)
  expect_equal(s_coord1$match_type[1], "coordinates")

  s_coord2 <- search_genes("NT_033777.3:16800000-16815000", "melanogaster")
  expect_true("abd-A" %in% s_coord2$gene_symbol)

  s_coord3 <- search_genes("X:2790000-2797000:-", "melanogaster")
  expect_equal(s_coord3$gene_symbol[1], "w")
  expect_equal(s_coord3$flybase_id[1], "FBgn0003996")
})

test_that("lookup_gene returns exact + strand feature geometry for yellow (y / CG3757)", {
  t0 <- Sys.time()
  g_y <- lookup_gene("yellow", species = "melanogaster", promoter_bp = 500L)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  expect_lt(elapsed, 1.0)

  expect_equal(g_y$symbol, "y")
  expect_equal(g_y$flybase_id, "FBgn0004034")
  expect_equal(g_y$cg_number, "CG3757")
  expect_equal(g_y$chr, "NC_004354.4")
  expect_equal(g_y$start, 356509L)
  expect_equal(g_y$end, 361245L)
  expect_equal(g_y$strand, "+")

  expect_equal(nrow(g_y$transcripts), 1L)
  expect_equal(g_y$transcripts$transcript_id[1], "NM_057444.3")
  expect_equal(g_y$transcripts$n_exons[1], 2L)
  expect_true("NM_057444.3" %in% names(g_y$exons_by_tx))

  ex <- g_y$exons_by_tx[["NM_057444.3"]]
  expect_equal(names(ex), c("chr", "start", "end", "strand", "transcript_id"))
  expect_equal(ex$start, c(356509L, 359617L))
  expect_equal(ex$end,   c(356917L, 361245L))

  cds <- g_y$cds_by_tx[["NM_057444.3"]]
  expect_equal(names(cds), c("chr", "start", "end", "strand", "phase", "transcript_id"))
  expect_equal(cds$start, c(356680L, 359617L))
  expect_equal(cds$end,   c(356917L, 361001L))
  expect_equal(nrow(g_y$shared_cds), 2L)
  expect_equal(g_y$shared_cds$start, c(356680L, 359617L))

  expect_equal(nrow(g_y$five_prime_utr), 1L)
  expect_equal(g_y$five_prime_utr$start[1], 356509L)
  expect_equal(g_y$five_prime_utr$end[1],   356679L)

  expect_equal(nrow(g_y$three_prime_utr), 1L)
  expect_equal(g_y$three_prime_utr$start[1], 361002L)
  expect_equal(g_y$three_prime_utr$end[1],   361245L)

  expect_equal(nrow(g_y$promoter_regions), 1L)
  expect_equal(g_y$promoter_regions$tss[1],   356509L)
  expect_equal(g_y$promoter_regions$start[1], 356009L)
  expect_equal(g_y$promoter_regions$end[1],   356508L)

  expect_equal(nrow(g_y$first_coding_exons), 1L)
  expect_equal(g_y$first_coding_exons$start[1], 356680L)
  expect_equal(g_y$first_coding_exons$end[1],   356917L)
})

test_that("lookup_gene returns exact - strand feature geometry for ebony, white, Abd-B, and PlexB", {
  g_e <- lookup_gene("ebony", species = "melanogaster", promoter_bp = 500L)
  expect_equal(g_e$symbol, "e")
  expect_equal(g_e$flybase_id, "FBgn0000527")
  expect_equal(g_e$cg_number, "CG3331")
  expect_equal(g_e$chr, "NT_033777.3")
  expect_equal(g_e$strand, "-")
  expect_equal(nrow(g_e$transcripts), 1L)
  expect_equal(g_e$transcripts$n_exons[1], 7L)

  expect_equal(g_e$promoter_regions$tss[1],   21237177L)
  expect_equal(g_e$promoter_regions$start[1], 21237178L)
  expect_equal(g_e$promoter_regions$end[1],   21237677L)

  expect_equal(g_e$first_coding_exons$start[1], 21231910L)
  expect_equal(g_e$first_coding_exons$end[1],   21233017L)

  expect_true(all(g_e$five_prime_utr$start > 21233017L))
  expect_true(21237177L %in% g_e$five_prime_utr$end)
  expect_true(all(g_e$three_prime_utr$end < 21230045L))

  g_w <- lookup_gene("white", species = "melanogaster")
  expect_equal(g_w$symbol, "w")
  expect_equal(g_w$flybase_id, "FBgn0003996")
  expect_equal(g_w$strand, "-")
  expect_equal(g_w$first_coding_exons$start[1], 2796237L)
  expect_equal(g_w$first_coding_exons$end[1],   2796308L)
  expect_equal(g_w$promoter_regions$start[1],   2796467L)
  expect_equal(g_w$promoter_regions$end[1],     2796966L)

  t0 <- Sys.time()
  g_abdb <- lookup_gene("Abd-B", species = "melanogaster")
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 1.0)
  expect_equal(g_abdb$symbol, "Abd-B")
  expect_equal(g_abdb$flybase_id, "FBgn0000015")
  expect_equal(nrow(g_abdb$transcripts), 10L)
  expect_equal(nrow(g_abdb$shared_cds), 4L)
  expect_equal(attr(g_abdb$shared_cds, "mode"), "intersection")
  expect_true(all(diff(g_abdb$shared_cds$start) < 0))

  g_plexb <- lookup_gene("FBgn0025740", species = "melanogaster")
  expect_equal(g_plexb$symbol, "PlexB")
  expect_equal(g_plexb$cg_number, "CG17245")
  expect_equal(nrow(g_plexb$transcripts), 2L)
  expect_equal(nrow(g_plexb$shared_cds), 5L)

  g_coord <- lookup_gene("X:2790000-2797000:-", species = "melanogaster")
  expect_equal(g_coord$symbol, "w")
  expect_equal(nrow(g_coord$shared_cds), 6L)
})

test_that("search_genes and lookup_gene work on non-melanogaster species (d_simulans and d_eugracilis)", {
  s_sim_abdb <- search_genes("homeobox protein abdominal-B", species = "simulans")
  expect_gt(nrow(s_sim_abdb), 0L)
  expect_equal(s_sim_abdb$gene_symbol[1], "LOC6728102")

  g_sim_abdb <- lookup_gene("LOC6728102", species = "simulans")
  expect_equal(g_sim_abdb$symbol, "LOC6728102")
  expect_equal(g_sim_abdb$chr, "NC_052523.2")
  expect_equal(g_sim_abdb$strand, "+")
  expect_gt(nrow(g_sim_abdb$transcripts), 1L)
  expect_gt(nrow(g_sim_abdb$shared_cds), 0L)
  expect_gt(nrow(g_sim_abdb$five_prime_utr), 0L)
  expect_gt(nrow(g_sim_abdb$promoter_regions), 0L)
  expect_gt(nrow(g_sim_abdb$first_coding_exons), 0L)

  g_sim_w <- lookup_gene("protein white", species = "simulans")
  expect_equal(g_sim_w$symbol, "LOC27207232")
  expect_equal(g_sim_w$chr, "NC_052525.2")
  expect_equal(g_sim_w$strand, "-")
  expect_gt(nrow(g_sim_w$shared_cds), 0L)
  expect_equal(g_sim_w$promoter_regions$start[1], g_sim_w$end + 1L)

  s_eug_abdb <- search_genes("homeobox protein abdominal-B", species = "eugracilis")
  expect_gt(nrow(s_eug_abdb), 0L)
  expect_equal(s_eug_abdb$gene_symbol[1], "LOC108103749")

  g_eug_abdb <- lookup_gene("LOC108103749", species = "eugracilis")
  expect_equal(g_eug_abdb$symbol, "LOC108103749")
  expect_equal(g_eug_abdb$strand, "-")
  expect_gt(nrow(g_eug_abdb$shared_cds), 0L)
  expect_gt(nrow(g_eug_abdb$promoter_regions), 0L)
})


