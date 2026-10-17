# ============================================================================
# tests/testthat/test-shiny_app.R — End-to-End Tests for Cas9_sgRNA_design.R
# ============================================================================

library(testthat)

project_root <- "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design"
app_path     <- file.path(project_root, "Cas9_sgRNA_design.R")
wrapper_path <- file.path(project_root, "app.R")

test_that("Cas9_sgRNA_design.R sources cleanly and defines ui, server, workflow, and plot functions", {
  expect_true(file.exists(app_path))
  expect_true(file.exists(wrapper_path))

  source(app_path, local = FALSE)

  expect_true(exists("ui", envir = .GlobalEnv))
  expect_true(exists("server", envir = .GlobalEnv))
  expect_true(is.function(server))
  expect_true(inherits(ui, "shiny.tag.list") || inherits(ui, "shiny.tag"))
  expect_true(exists("run_crispr_design_workflow", mode = "function"))
  expect_true(exists("plot_crispr_isoform_tracks", mode = "function"))
  expect_true(exists("write_crispr_output_bundle", mode = "function"))
})

test_that("run_crispr_design_workflow executes end-to-end on real D. melanogaster genes and coordinates", {
  source(app_path, local = FALSE)

  wf_ebony <- run_crispr_design_workflow(
    query          = "ebony",
    species        = "melanogaster",
    target_mode    = "coding_ko",
    offtarget_mode = "locus",
    arm_5p_len     = 1200L,
    arm_3p_len     = 1500L,
    arm_3p_anchor  = "cut_site",
    selected_rank  = 1L,
    seed           = 42L
  )

  expect_equal(wf_ebony$status, "ok")
  expect_equal(wf_ebony$species, "melanogaster")
  expect_equal(wf_ebony$gene_model$gene_symbol, "e")
  expect_true(is.data.frame(wf_ebony$all_candidates))
  expect_gt(nrow(wf_ebony$all_candidates), 50L)
  expect_true(is.data.frame(wf_ebony$filtered_candidates))
  expect_gt(nrow(wf_ebony$filtered_candidates), 10L)

  # Verify selected guide and homology arm lengths & exact cut-site adjacency
  arm <- wf_ebony$arm_design
  expect_false(is.null(arm))
  expect_equal(arm$arm_5p_length, 1200L)
  expect_equal(arm$arm_3p_length, 1500L)
  expect_equal(nchar(arm$arm_5p_seq), 1200L)
  expect_equal(nchar(arm$arm_3p_seq), 1500L)

  # Verify cut-site zoom view
  zm <- wf_ebony$cut_site_zoom
  expect_false(is.null(zm))
  expect_equal(nchar(zm$spacer_5p_17), 17L)
  expect_equal(nchar(zm$spacer_3p_3), 3L)
  expect_equal(paste0(zm$spacer_5p_17, zm$spacer_3p_3), arm$spacer_seq)
  expect_true(grepl("\\|", zm$target_5to3_text))
})

test_that("plot_crispr_isoform_tracks renders to PDF without error in both genomic and 5'->3' orientations", {
  source(app_path, local = FALSE)

  wf <- run_crispr_design_workflow(
    query          = "white",
    species        = "melanogaster",
    target_mode    = "coding_ko",
    offtarget_mode = "locus",
    arm_5p_len     = 1000L,
    arm_3p_len     = 1000L,
    seed           = 42L
  )

  tmp_pdf1 <- tempfile(fileext = ".pdf")
  tmp_pdf2 <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(tmp_pdf1, tmp_pdf2)), add = TRUE)

  grDevices::pdf(tmp_pdf1, width = 10, height = 5)
  res1 <- plot_crispr_isoform_tracks(
    gene_model          = wf$gene_model,
    candidates_df       = wf$filtered_candidates,
    selected_guide      = wf$selected_guide,
    arm_design          = wf$arm_design,
    orient_5p_to_3p     = FALSE,
    show_all_candidates = TRUE
  )
  grDevices::dev.off()
  expect_true(isTRUE(res1))
  expect_true(file.exists(tmp_pdf1))
  expect_gt(file.info(tmp_pdf1)$size, 1000)

  grDevices::pdf(tmp_pdf2, width = 10, height = 5)
  res2 <- plot_crispr_isoform_tracks(
    gene_model          = wf$gene_model,
    candidates_df       = wf$filtered_candidates,
    selected_guide      = wf$selected_guide,
    arm_design          = wf$arm_design,
    orient_5p_to_3p     = TRUE,
    show_all_candidates = TRUE
  )
  grDevices::dev.off()
  expect_true(isTRUE(res2))
  expect_true(file.exists(tmp_pdf2))
  expect_gt(file.info(tmp_pdf2)$size, 1000)
})

test_that("write_crispr_output_bundle writes complete reproducible ZIP, FASTA, CSV, JSON, PDF, and log files", {
  source(app_path, local = FALSE)

  wf <- run_crispr_design_workflow(
    query          = "yellow",
    species        = "melanogaster",
    target_mode    = "promoter_5utr",
    offtarget_mode = "locus",
    arm_5p_len     = 1100L,
    arm_3p_len     = 1400L,
    seed           = 42L
  )

  out_dir <- tempfile("crispr_bundle_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  bundle <- write_crispr_output_bundle(out_dir, wf)

  expect_true(file.exists(bundle$csv_path))
  expect_true(file.exists(bundle$fasta_path))
  expect_true(file.exists(bundle$json_path))
  expect_true(file.exists(bundle$log_path))
  expect_true(file.exists(bundle$pdf_path))
  expect_true(file.exists(bundle$zip_path))

  expect_gt(file.info(bundle$csv_path)$size,   500)
  expect_gt(file.info(bundle$fasta_path)$size, 2500)
  expect_gt(file.info(bundle$json_path)$size,  500)
  expect_gt(file.info(bundle$log_path)$size,   200)
  expect_gt(file.info(bundle$zip_path)$size,   1000)

  zip_contents <- utils::unzip(bundle$zip_path, list = TRUE)$Name
  expect_true(all(c("gRNA_candidates.csv", "homology_arms.fasta", "design_summary.json", "CRISPR_design_log.txt") %in% zip_contents))
})
