#!/usr/bin/env Rscript
# ============================================================================
# app.R — Shiny Directory Entry Point Wrapper for Plasmid_Designer.R
# ============================================================================
# Allows launching via `shiny::runApp(".")` or RStudio / Electron runner.R.
# ============================================================================

app_script <- file.path(getwd(), "Plasmid_Designer.R")
if (!file.exists(app_script)) {
  app_script <- "/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer/Plasmid_Designer.R"
}
source(app_script, local = FALSE)
shiny::shinyApp(ui = ui, server = server)
