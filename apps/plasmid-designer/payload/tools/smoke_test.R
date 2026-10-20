#!/usr/bin/env Rscript
# Smoke test: build the app and drive the server headlessly.
suppressPackageStartupMessages(library(shiny))
options(shiny.launch.browser = FALSE)
setwd("/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer")
src <- readLines("Plasmid_Designer.R")
tail_start <- grep("^\\.args <- commandArgs", src)[1]
eval(parse(text = paste(src[1:(tail_start - 1)], collapse = "\n")))
app <- shiny::shinyApp(ui, server)
cat("shinyApp object built OK\n")
shiny::testServer(app, {
  session$setInputs(pjat_backbone = "pJAT74", insert_strategy = "linker",
                    donor_source = "marker", donor_marker = "3xP3_dsRed",
                    donor_orient = "+", arm_source = "bundle", donor_arm_len = 500,
                    guide_source = "paste", spacer_paste = "GAGTCTGACCTGACTGACTG",
                    insert_orient = "+", map_view = "plasmidplot")
  cat("backbone:", pjd_backbone()$name, "\n")
  cat("payload bp:", nchar(pjd_donor_payload()$sequence), "\n")
  cat("gRNA site:", pjd_grna_site_label(pjd_backbone()$u6_promoters, FALSE), "\n")
  cat("registry backbones:", length(backbone_pjat_adapter()), "\n")
  cat("payload resolver label:", pjd_donor_payload()$label, "\n")
  cat("TEST SERVER OK\n")
})
