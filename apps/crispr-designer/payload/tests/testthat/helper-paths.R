# ============================================================================
# tests/testthat/helper-paths.R — portable path resolution for the test suite.
#
# The tests used to point at a fixed path on an external volume:
#
#   /Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design
#
# That path does not exist anywhere else, so the project root is now resolved
# from the location of the test files themselves. Each test file then uses
# `project_root` from this helper instead of hard-coding a path.
# ============================================================================

# The directory that contains R/ — i.e. the application root.
project_root <- local({
  # tests/testthat/ → tests/ → project root
  candidates <- c(
    dirname(dirname(dirname(normalizePath(
      testthat::test_path("helper-paths.R"), mustWork = FALSE)))),
    getwd(),
    dirname(getwd())
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "R", "genome_registry.R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  normalizePath(candidates[1], mustWork = FALSE)
})

# The Janelia pJAT plasmid collection.
pjats_dir <- local({
  if (exists("get_app_pjat_dir", mode = "function")) {
    d <- get_app_pjat_dir("CRISPR")
    if (nzchar(d) && dir.exists(d)) return(d)
  }
  for (cand in c(
    file.path(project_root, "janelia_atalanta_plasmids"),
    file.path(project_root, "reference_files", "janelia_atalanta_plasmids"),
    file.path(dirname(project_root), "reference_files", "janelia_atalanta_plasmids")
  )) {
    if (dir.exists(cand)) return(normalizePath(cand, mustWork = FALSE))
  }
  file.path(project_root, "janelia_atalanta_plasmids")
})

# The genome store, resolved the same way the application resolves it.
genome_store_dir <- local({
  if (exists("get_app_reference_dir", mode = "function")) {
    return(get_app_reference_dir("CRISPR"))
  }
  # The shared store written by HCR Probe Designer.
  for (cand in c(
    Sys.getenv("CRISPR_REFERENCE_DIR", unset = ""),
    Sys.getenv("HCR_REFERENCE_DIR", unset = ""),
    file.path(path.expand("~"), "Library", "Application Support",
              "hcr-probe-designer", "reference_files"),
    file.path(path.expand("~"), "Downloads", "reference_files"),
    file.path(project_root, "reference_files"),
    file.path(dirname(project_root), "reference_files")
  )) {
    if (nzchar(cand) && dir.exists(cand)) return(normalizePath(cand, mustWork = FALSE))
  }
  file.path(project_root, "reference_files")
})

# The sibling Plasmid Designer application, used by the cross-application
# bundle round-trip test. It lives next to this project in the repository.
plasmid_designer_root <- local({
  for (cand in c(
    file.path(dirname(project_root), "plasmid-designer", "payload"),
    file.path(dirname(project_root), "Plasmid_Designer"),
    Sys.getenv("PLASMID_DESIGNER_ROOT", unset = "")
  )) {
    if (nzchar(cand) && file.exists(file.path(cand, "R", "guide_bundle.R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  file.path(dirname(project_root), "plasmid-designer", "payload")
})
