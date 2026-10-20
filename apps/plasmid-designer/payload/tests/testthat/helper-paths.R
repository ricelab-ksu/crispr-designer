# ============================================================================
# tests/testthat/helper-paths.R — portable path resolution for the test suite.
#
# The tests used to point at a fixed path on an external volume:
#
#   /Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer
#
# That path does not exist anywhere else, so the project root is now resolved
# from the location of the test files themselves. Each test file then uses
# `project_root` from this helper instead of hard-coding a path.
# ============================================================================

# The directory that contains R/ — i.e. the application root.
# `sys.frame(1)$ofile` is NULL when testthat sources via its own wrapper, so
# walk up from the working directory until R/backbone_registry.R is found.
project_root <- local({
  candidates <- c(getwd(), dirname(getwd()), dirname(dirname(getwd())))
  d <- getwd()
  for (i in seq_len(6)) {
    candidates <- c(candidates, d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  for (cand in unique(candidates)) {
    if (nzchar(cand) && file.exists(file.path(cand, "R", "backbone_registry.R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  normalizePath(candidates[1], mustWork = FALSE)
})

# Make the portable path helpers available before anything that uses them.
# app_paths.R sits at apps/plasmid-designer/payload/R/app_paths.R.
for (p in unique(c(
  file.path(project_root, "R", "app_paths.R"),
  file.path(dirname(project_root), "shared", "r", "app_paths.R")
))) {
  if (file.exists(p)) {
    tryCatch(source(p, local = FALSE), error = function(e) NULL)
    break
  }
}

# The Janelia pJAT plasmid collection.
pjats_dir <- local({
  if (exists("get_app_pjat_dir", mode = "function")) {
    d <- get_app_pjat_dir("PLASMID")
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

# The sibling CRISPR Designer application, used by the cross-application
# bundle round-trip test.
crispr_designer_root <- local({
  for (cand in c(
    file.path(dirname(project_root), "crispr-designer", "payload"),
    file.path(dirname(project_root), "Cas9_sgRNA_design"),
    Sys.getenv("CRISPR_DESIGNER_ROOT", unset = "")
  )) {
    if (nzchar(cand) && file.exists(file.path(cand, "R", "grna_engine.R"))) {
      return(normalizePath(cand, mustWork = FALSE))
    }
  }
  file.path(dirname(project_root), "crispr-designer", "payload")
})
