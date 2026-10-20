#!/usr/bin/env Rscript
# ============================================================================
# Shared Electron sidecar bootstrap runner.
#
# Executed by the bundled R runtime that the Electron main process spawns. It:
#   1. Resolves the directory holding the app's R payload.
#   2. Points R's package library at the bundled `library/` tree (if present).
#   3. Reports the reference directories the main process resolved.
#   4. Sources the app script (which must define `ui` and `server`).
#   5. Picks a free TCP port and writes `APP_PORT=<n>` to stdout so Electron
#      can attach once the Shiny server is listening.
#
# Usage (from Electron):
#   R --no-save --no-restore --no-environ --no-site-file --no-init-file \
#     --slave -f runner.R --args <app_dir>
#
# The app script name is passed through the APP_SCRIPT environment variable,
# which the main process sets per app.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
app_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."
if (!dir.exists(app_dir)) app_dir <- getwd()

# ---------------------------------------------------------------------------
# 1. Library path: prefer a bundled library next to the app payload.
#    Layout: <resources>/r/{R.framework|R-x.y}/...  and <resources>/r/library
# ---------------------------------------------------------------------------
bundled_lib <- file.path(dirname(app_dir), "r", "library")
if (dir.exists(bundled_lib)) {
  .libPaths(c(bundled_lib, .libPaths()))
}
cat("APP_LIBS=", paste(.libPaths(), collapse = ";"), "\n", sep = "")

# ---------------------------------------------------------------------------
# 2. Reference data directories, resolved by the Electron main process and
#    passed in through the environment. Each app uses its own prefix.
# ---------------------------------------------------------------------------
for (ev in c("PLASMID_REFERENCE_DIR", "PLASMID_PJAT_DIR",
             "CRISPR_REFERENCE_DIR", "CRISPR_PJAT_DIR",
             "HCR_REFERENCE_DIR", "PJAT_REFERENCE_DIR")) {
  val <- Sys.getenv(ev, unset = "")
  if (nzchar(val)) cat("APP_REF ", ev, "=", val, "\n", sep = "")
}

# ---------------------------------------------------------------------------
# 3. Source the application script. It must define ui/server in the global env.
# ---------------------------------------------------------------------------
app_file <- file.path(app_dir, Sys.getenv("APP_SCRIPT", unset = ""))
if (!nzchar(Sys.getenv("APP_SCRIPT", unset = "")) || !file.exists(app_file)) {
  stop("App script not found: ", app_file, "\nExpected in: ", app_dir)
}
if (getRversion() < "4.0") {
  stop("This app requires R >= 4.0 (found ", getRversion(), ")")
}
source(app_file, local = FALSE, chdir = TRUE)

if (!exists("ui", envir = globalenv()) || !exists("server", envir = globalenv())) {
  stop("App script did not define ui/server objects as expected")
}

# ---------------------------------------------------------------------------
# 4. Choose a free port and launch Shiny headlessly.
# ---------------------------------------------------------------------------
if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The 'shiny' package is not installed in the bundled library ",
       "(looked in ", paste(.libPaths(), collapse = ", "), "). ",
       "Re-run the bundling step: `npm run build:mac`.")
}

port <- NULL
for (i in seq_len(25)) {
  p <- try(httpuv::randomPort(30000L, 60000L), silent = TRUE)
  if (inherits(p, "try-error") || is.null(p) || is.na(p) || length(p) != 1L) next
  port <- as.integer(p)
  break
}
if (is.null(port)) {
  stop("Could not find a free TCP port to host the Shiny server")
}

shiny_app <- shiny::shinyApp(ui = ui, server = server)

# Handshake line that Electron watches for
cat("APP_PORT=", port, "\n", sep = "")
flush(stdout())

shiny::runApp(
  shiny_app,
  host = "127.0.0.1",
  port = port,
  launch.browser = FALSE
)
