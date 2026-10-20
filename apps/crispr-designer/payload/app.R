#!/usr/bin/env Rscript
# ============================================================================
# app.R — Shiny Directory Entry Point Wrapper for Cas9_sgRNA_design.R
# ============================================================================
# Allows launching via `shiny::runApp(".")` or RStudio / Electron runner.R.
# ============================================================================

app_script <- file.path(getwd(), "Cas9_sgRNA_design.R")
if (!file.exists(app_script)) {
  app_script <- file.path(dirname(getwd()), "Cas9_sgRNA_design.R")
}
if (!file.exists(app_script)) {
  stop("Cannot locate Cas9_sgRNA_design.R")
}
source(app_script, local = FALSE)
shiny::shinyApp(ui = ui, server = server)
