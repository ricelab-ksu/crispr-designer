#!/usr/bin/env Rscript
# ============================================================================
# tools/r-closure.R — report the R package closure the desktop apps need.
#
# Prints one package name per line: the apps' direct dependencies plus their
# recursive Depends/Imports/LinkingTo closure, minus base and recommended
# packages (those are copied wholesale from the framework).
#
# Usage: Rscript tools/r-closure.R
# ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

need <- c(
  # Plasmid Designer
  "shiny", "yaml", "jsonlite", "httr", "plasmidplot", "png", "readr",
  # CRISPR Designer
  "Biostrings", "DT", "data.table",
  # Shared infrastructure
  "httpuv", "testthat"
)

avail <- rownames(installed.packages())
base <- rownames(installed.packages(priority = "base"))
rec  <- rownames(installed.packages(priority = "recommended"))

# tools::package_dependencies() only knows about CRAN, so Bioconductor
# packages (Biostrings and its BiocGenerics/IRanges/S4Vectors/XVector
# dependencies) must be added by hand and expanded with the local package
# database.
bioc_seed <- c("BiocGenerics", "IRanges", "S4Vectors", "XVector",
               "GenomeInfoDb", "UCSC.utils")

deps <- tools::package_dependencies(
  need,
  which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE
)

# Expand the Bioconductor seeds using the installed package database.
# package_dependencies() returns the dependencies of the seeds, not the seeds
# themselves, so the seeds are added back in explicitly.
bioc_present <- intersect(bioc_seed, avail)
bioc_deps <- tools::package_dependencies(
  bioc_present,
  db = installed.packages(),
  which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE
)

all_needed <- unique(c(need, unlist(deps), bioc_present, unlist(bioc_deps)))
closure <- setdiff(intersect(all_needed, avail), c(base, rec))

missing <- setdiff(need, avail)
if (length(missing)) {
  warning("Required packages are not installed: ", paste(missing, collapse = ", "))
}

cat(sort(closure), sep = "\n")
