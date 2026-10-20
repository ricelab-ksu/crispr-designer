# ============================================================================
# R/app_paths.R — portable reference-data resolution.
#
# The applications used to point at a fixed path on an external volume:
#
#   /Volumes/Backup_Plus/projects/bioinformatics/reference_files/…
#
# That path does not exist on any other machine, so the desktop builds
# resolve the reference stores at runtime instead. Resolution order:
#
#   1. an explicit environment variable (set by the Electron main process),
#   2. the shared genome store shipped with HCR Probe Designer,
#   3. a store inside the application payload,
#   4. a store next to the working directory.
#
# Nothing here writes to disk; the resolver only reports where things are.
# ============================================================================

#' Locate and load R/app_paths.R.
#'
#' The helper modules are sourced both from the application root and from the
#' test suite, where the working directory is tests/testthat. The most
#' reliable anchor is the path of the file doing the sourcing, which R
#' exposes as `sys.frame(1)$ofile`; the working directory is used as a
#' fallback.
#'
#' @param env_prefix Prefix for the app-specific reference variables.
#' @return TRUE when the helpers are available, FALSE otherwise.
load_app_paths <- function(env_prefix = "CRISPR") {
  if (exists("app_root", mode = "function")) return(TRUE)

  roots <- character(0)

  # The file that called us, if any.
  for (i in seq_along(sys.frames())) {
    of <- tryCatch(sys.frames()[[i]]$ofile, error = function(e) NULL)
    if (!is.null(of) && nzchar(of)) {
      roots <- c(roots, dirname(normalizePath(of, mustWork = FALSE)))
    }
  }

  # The working directory and its parents.
  for (start in c(getwd(), dirname(getwd()))) {
    d <- start
    for (i in seq_len(6)) {
      roots <- c(roots, d)
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    roots <- c(roots, dirname(sub("^--file=", "", file_arg[1])))
  }
  roots <- c(roots, Sys.getenv("APP_ROOT", unset = ""))

  for (r in unique(roots)) {
    if (!nzchar(r) || is.na(r)) next
    p <- file.path(r, "R", "app_paths.R")
    if (file.exists(p)) {
      tryCatch(source(p, local = FALSE), error = function(e) NULL)
      break
    }
  }
  exists("app_root", mode = "function")
}

#' The directory of the file that is currently being sourced.
#'
#' @return A path, or NA when the caller is not a sourced file.
sourced_file_dir <- function() {
  for (i in seq_along(sys.frames())) {
    of <- tryCatch(sys.frames()[[i]]$ofile, error = function(e) NULL)
    if (!is.null(of) && nzchar(of)) {
      return(dirname(normalizePath(of, mustWork = FALSE)))
    }
  }
  NA_character_
}

#' Resolve the application root directory.
#'
#' The application root is the directory that contains the `R/` folder. It is
#' discovered by walking up from the working directory and from the path of
#' the script that is currently running, so it works both in a source
#' checkout and inside a packaged application.
#'
#' @return The application root, or the working directory if none is found.
app_root <- function() {
  candidates <- unique(c(
    getwd(),
    dirname(getwd()),
    dirname(dirname(getwd())),
    # The script currently being run, e.g. .../payload/Cas9_sgRNA_design.R
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
                             value = TRUE)[1]),
    Sys.getenv("APP_ROOT", unset = "")
  ))
  for (cand in candidates) {
    if (is.na(cand) || !nzchar(cand)) next
    if (dir.exists(file.path(cand, "R"))) return(normalizePath(cand, mustWork = FALSE))
  }
  getwd()
}

#' Candidate directories for the per-species genome store.
#'
#' @param env_prefix Prefix for the app-specific override variable,
#'   e.g. "CRISPR" yields "CRISPR_REFERENCE_DIR".
#' @return character vector of candidate paths, most specific first.
genome_store_candidates <- function(env_prefix = "CRISPR") {
  cands <- c()

  # 1. Explicit overrides.
  for (v in c(paste0(env_prefix, "_REFERENCE_DIR"), "HCR_REFERENCE_DIR")) {
    val <- Sys.getenv(v, unset = "")
    if (nzchar(val)) cands <- c(cands, val)
  }

  # 2. The shared store written by HCR Probe Designer.
  hcr <- file.path(path.expand("~"), "Library", "Application Support",
                   "hcr-probe-designer", "reference_files")
  cands <- c(cands, hcr, file.path(path.expand("~"), "Downloads", "reference_files"))

  # 3. Stores inside the app payload / next to the working directory.
  cands <- c(cands,
             file.path(getwd(), "reference_files"),
             file.path(getwd(), "reference_files", "drosophila_species_genomes"),
             file.path(dirname(getwd()), "reference_files"),
             file.path(dirname(dirname(getwd())), "reference_files"))

  unique(cands)
}

#' True when `dir` holds at least one usable per-species genome folder.
#'
#' @param dir Candidate directory.
#' @return logical.
is_genome_store <- function(dir) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(FALSE)
  subs <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
  subs <- subs[grepl("^d_[a-z0-9_]+$", subs, ignore.case = TRUE)]
  for (s in subs) {
    files <- list.files(file.path(dir, s), ignore.case = TRUE)
    if (any(grepl("^genomic\\.fna(\\.gz)?$", files))) return(TRUE)
  }
  FALSE
}

#' Resolve the genome store directory.
#'
#' @param env_prefix Prefix for the app-specific override variable.
#' @return A path. The first candidate is returned when nothing exists yet, so
#'   callers always receive a usable string; `is_genome_store()` can be used
#'   to test whether it is actually populated.
get_app_reference_dir <- function(env_prefix = "CRISPR") {
  cands <- genome_store_candidates(env_prefix)
  for (c in cands) {
    if (is_genome_store(c)) return(normalizePath(c, mustWork = FALSE))
  }
  if (length(cands)) cands[1] else ""
}

#' Candidate directories for the Janelia pJAT plasmid collection.
#'
#' @param env_prefix Prefix for the app-specific override variable.
#' @return character vector of candidate paths.
pjat_store_candidates <- function(env_prefix = "CRISPR") {
  cands <- c()
  for (v in c(paste0(env_prefix, "_PJAT_DIR"), "PJAT_REFERENCE_DIR")) {
    val <- Sys.getenv(v, unset = "")
    if (nzchar(val)) cands <- c(cands, val)
  }
  root <- app_root()
  cands <- c(cands,
             # Bundled inside the payload (packaged builds stage it here).
             file.path(root, "janelia_atalanta_plasmids"),
             # Development checkout: shared/reference/ next to apps/.
             file.path(dirname(root), "shared", "reference", "janelia_atalanta_plasmids"),
             file.path(root, "reference_files", "janelia_atalanta_plasmids"),
             file.path(dirname(getwd()), "janelia_atalanta_plasmids"),
             file.path(dirname(dirname(getwd())), "janelia_atalanta_plasmids"),
             # User locations.
             file.path(path.expand("~" ), "Downloads", "janelia_atalanta_plasmids"),
             file.path(path.expand("~" ), "Library", "Application Support",
                       "hcr-probe-designer", "janelia_atalanta_plasmids"))
  unique(cands)
}

#' True when `dir` holds the pJAT GenBank collection.
#'
#' @param dir Candidate directory.
#' @return logical.
is_pjat_store <- function(dir) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(FALSE)
  files <- list.files(dir, pattern = "\\.(gbk|gb|gbff|genbank)$", ignore.case = TRUE)
  any(grepl("^pJAT", files, ignore.case = TRUE))
}

#' Resolve the pJAT plasmid collection directory.
#'
#' @param env_prefix Prefix for the app-specific override variable.
#' @return A path.
get_app_pjat_dir <- function(env_prefix = "CRISPR") {
  cands <- pjat_store_candidates(env_prefix)
  for (c in cands) {
    if (is_pjat_store(c)) return(normalizePath(c, mustWork = FALSE))
  }
  if (length(cands)) cands[1] else ""
}
