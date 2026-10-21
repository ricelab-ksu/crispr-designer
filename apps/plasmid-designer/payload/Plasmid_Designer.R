#!/usr/bin/env Rscript
# ============================================================================
# Plasmid_Designer.R — Janelia Atalanta (pJAT) Vector Designer
# ============================================================================
# Standalone Shiny GUI + Headless CLI for assembling pJAT expression vectors:
#   * pick a pJAT backbone (attR coordinates, U6 count, Addgene/DGRC source)
#   * import Cas9 guides (cas9_guide_bundle/1.0) or paste a 20-nt spacer
#   * insert a gRNA cassette / marker donor between homology arms
#   * design attL donor fragments (250-1000 bp arms, rule R2) via pjat_assemble_fragments()
#   * export GenBank / FASTA / JSON / ZIP
#
# Modeled on the Cas9 sgRNA Designer (Cas9_sgRNA_design.R): same .cas9-*
# CSS skin, same numbered sidebar sections, same download bar.
#
# Port: 7790 (Cas9 app uses 7789 so both can run side by side).
# ============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("shiny", quietly = TRUE)) library(shiny)

# ---- Portable path helpers (must be loaded before anything else) -----------
# R/app_paths.R provides app_root() plus the runtime resolvers for the genome
# and pJAT reference stores. It is loaded first because every other module
# depends on it.

.pjd_source_app_paths <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  roots <- c(getwd(), dirname(getwd()), dirname(dirname(getwd())))
  if (length(file_arg)) roots <- c(dirname(sub("^--file=", "", file_arg[1])), roots)
  for (r in unique(roots)) {
    p <- file.path(r, "R", "app_paths.R")
    if (file.exists(p)) {
      tryCatch(source(p, local = FALSE), error = function(e) NULL)
      break
    }
  }
}
.pjd_source_app_paths()
if (exists("load_app_paths", mode = "function")) load_app_paths("PLASMID")


})

PJD_APP_VERSION <- "1.0.0"
PORT <- as.integer(Sys.getenv("PJD_PORT", unset = "7790"))
if (is.na(PORT) || PORT <= 0L) PORT <- 7790L

# ---- Resolve project root --------------------------------------------------

.resolve_pjd_app_root <- function() {
  candidates <- unique(c(
    getwd(),
    app_root()
  ))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (file.exists(script_path)) {
      candidates <- unique(c(dirname(normalizePath(script_path, mustWork = FALSE)), candidates))
    }
  }
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "R"))) return(normalizePath(cand))
  }
  candidates[1]
}

PJD_APP_ROOT <- .resolve_pjd_app_root()

# ---- Source pJAT modules ---------------------------------------------------

.source_pjd_modules <- function(root = PJD_APP_ROOT) {
  mod_files <- c("addgene_maps.R", "pjat_plasmids.R", "gateway_attl.R",
                 "marker_sequences.R", "sample_inputs.R", "plasmid_map.R",
                 "guide_bundle.R", "plasmid_assembly.R",
                 "backbone_registry.R", "plasmid_assembly_slots.R",
                 "arm_resolution.R", "design_schematic.R", "onboarding.R")
  for (mf in mod_files) {
    p <- file.path(root, "R", mf)
    if (file.exists(p)) {
      tryCatch(source(p, local = FALSE), error = function(e) {
        warning("Could not source ", p, ": ", conditionMessage(e))
      })
    }
  }
}

.source_pjd_modules()

# ---- Static resources -------------------------------------------------------
# Serve the downloaded Addgene vector-map PNGs explicitly. Relying on Shiny's
# implicit `www/` lookup is not enough here: the app is usually started as
# `Rscript Plasmid_Designer.R` (shinyApp() with no appDir), where `www/` is not
# registered automatically.
if (exists("addgene_maps_dir", mode = "function") &&
    dir.exists(addgene_maps_dir()) &&
    requireNamespace("shiny", quietly = TRUE)) {
  shiny::addResourcePath(ADDGENE_MAPS_DIRNAME, addgene_maps_dir())
}

# Interactive vector-map assets (SVG interaction layer + styles).
PJD_ASSET_DIR <- file.path(PJD_APP_ROOT, "www", PJD_MAP_ASSET_DIR)
if (dir.exists(PJD_ASSET_DIR) && requireNamespace("shiny", quietly = TRUE)) {
  shiny::addResourcePath(PJD_MAP_ASSET_DIR, PJD_ASSET_DIR)
}

# ---- Marker catalog (metadata only, no sequences) --------------------------

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0L && all(!is.na(a))) a else b
}

.read_marker_catalog <- function() {
  yml <- file.path(PJD_APP_ROOT, "donor_templates", "marker_cassettes.yml")
  if (!file.exists(yml)) return(data.frame(id = character(0), label = character(0),
                                           description = character(0), selection = character(0),
                                           source = character(0), sequence_status = character(0),
                                           fetch_from = character(0),
                                           expected_features = I(list()),
                                           sequence_provenance = character(0),
                                           sequence_file = character(0),
                                           length_bp = integer(0),
                                           stringsAsFactors = FALSE))
  y <- yaml::read_yaml(yml)
  mk <- lapply(y$markers, function(m) {
    data.frame(
      id = m$id %||% NA_character_,
      label = m$label %||% NA_character_,
      description = m$description %||% NA_character_,
      selection = m$selection %||% NA_character_,
      source = m$source %||% NA_character_,
      sequence_status = m$sequence_status %||% NA_character_,
      fetch_from = m$fetch_from %||% NA_character_,
      expected_features = I(list(m$expected_features %||% character(0))),
      sequence_provenance = m$sequence_provenance %||% NA_character_,
      sequence_file = m$sequence_file %||% NA_character_,
      length_bp = suppressWarnings(as.integer(m$length_bp %||% NA_integer_)),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, mk)
}

PJD_MARKER_CATALOG <- .read_marker_catalog()

read_yaml_marker_catalog <- function() PJD_MARKER_CATALOG

# ============================================================================
# CSS (same skin as the Cas9 Designer, .cas9-* classes)
# ============================================================================

PJD_CUSTOM_CSS <- "
body { background-color: #f8fafc; color: #0f172a; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
.well { background-color: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; box-shadow: 0 1px 2px rgba(15,23,42,0.04); }
.cas9-card { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 16px; margin-bottom: 14px; box-shadow: 0 1px 2px rgba(15,23,42,0.04); }
.cas9-banner { border-left: 4px solid #0ea5e9; }
.cas9-section-hdr { font-size: 13px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em; color: #0284c7; margin-top: 8px; margin-bottom: 8px; }
.cas9-badge { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px; font-weight: 600; }
.cas9-badge-sky { background: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd; }
.cas9-badge-emerald { background: #d1fae5; color: #047857; border: 1px solid #a7f3d0; }
.cas9-badge-amber { background: #fef3c7; color: #b45309; border: 1px solid #fde68a; }
.cas9-badge-red { background: #fee2e2; color: #b91c1c; border: 1px solid #fecaca; }
.cas9-badge-slate { background: #f1f5f9; color: #334155; border: 1px solid #cbd5e1; }
.cas9-dl-bar { background: #f0f9ff; border: 1px solid #bae6fd; border-radius: 8px; padding: 10px 14px; margin-bottom: 14px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 8px; }
.pjd-backbone-panel { border-left: 4px solid #0ea5e9; }
.pjd-design-card-wrap { border-left: 4px solid #db2777; }
.pjd-design-label { font-size: 12px; font-weight: 700; fill: #0f172a; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
.pjd-design-sub { font-size: 10.5px; fill: #475569; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
.pjd-backbone-panel .cas9-section-hdr { margin-top: 0; margin-bottom: 6px; }
.pjd-backbone-title { display: flex; justify-content: space-between; align-items: baseline; gap: 12px; flex-wrap: wrap; margin-bottom: 8px; }
.pjd-backbone-name { font-size: 18px; font-weight: 800; color: #0f172a; letter-spacing: 0.01em; }
.pjd-backbone-purpose { font-size: 13px; color: #334155; }
.pjd-badges { display: flex; gap: 6px; flex-wrap: wrap; margin-bottom: 8px; }
.pjd-att-coords { font-size: 11px; color: #64748b; margin: 4px 0 10px 0; }
.pjd-note-warn { font-size: 12px; color: #475569; background: #fffbeb; border: 1px solid #fde68a; border-radius: 6px; padding: 8px 10px; margin-bottom: 10px; }
.pjd-map-large { width: 100%; border: 1px solid #e2e8f0; border-radius: 8px; background: #ffffff; display: block; margin: 6px auto; padding: 6px; }
.pjd-map-caption { font-size: 11px; color: #64748b; text-align: center; margin-top: 2px; }
.pjd-link { font-size: 12px; font-weight: 600; color: #0284c7; text-decoration: none; }
.pjd-link:hover { color: #0369a1; text-decoration: underline; }
.pjd-map-fallback { font-size: 12px; color: #475569; }
.pjd-map-fallback code { background: #f1f5f9; border: 1px solid #e2e8f0; border-radius: 4px; padding: 1px 5px; font-size: 11px; }
.pjd-map-view-row .form-group { margin-bottom: 8px; }
.pjd-map-view-row .control-label { font-size: 12px; font-weight: 600; color: #334155; }
.pjd-map-view-row .radio-inline { font-size: 12px; }
"

# ============================================================================
# UI builder helpers
# ============================================================================

# ---- Unified backbone panel (main screen) -----------------------------------
# All plasmid identity, purpose, attR coordinates, supplier warning and vector
# map live in one card in the main panel; the sidebar holds only the picker.

# ---- Addgene design-purpose card + vector map -------------------------------

.pjd_addgene_link <- function(cat, label = "View on Addgene") {
  tags$a(class = "pjd-link", href = addgene_page_url(cat$addgene_id),
         target = "_blank", rel = "noopener noreferrer", label)
}

.pjd_supplier_badge <- function(cat) {
  if (!is.na(cat$supplier) && cat$supplier == "DGRC_1632") {
    span(class = "cas9-badge cas9-badge-amber", "DGRC 1632")
  } else if (!is.na(cat$supplier) && !identical(cat$supplier, "Addgene")) {
    span(class = "cas9-badge cas9-badge-slate", cat$supplier)
  } else {
    NULL   # "Addgene" would duplicate the "Addgene <id>" badge
  }
}

#' Static Addgene/SnapGene map block with zoom.
#'
#' Renders the downloaded SnapGene PNG inside a `.pjd-png-view` container.
#' `zoom_pct` scales the image width (100 = fit the card, up to 300 = 3x).
#' All arguments are plain values (no `input$`, no `reactive()`): the server's
#' `renderUI` reads the Shiny inputs and passes them in.
#' @keywords internal
.pjd_addgene_map_block <- function(cat, gb = NULL, zoom_pct = 100) {
  src <- if (isTRUE(cat$map_available)) addgene_map_src(cat$addgene_id) else NULL
  if (is.null(src)) {
    return(div(class = "pjd-map-fallback",
               "No local vector map yet for this backbone. Build the image cache with ",
               tags$code("Rscript tools/fetch_addgene_maps.R"),
               ", or open the map on ",
               .pjd_addgene_link(cat, paste0("Addgene #", cat$addgene_id, " \u2197")),
               "."))
  }
  zoom_pct <- suppressWarnings(as.numeric(zoom_pct[[1L]]))
  if (length(zoom_pct) == 0L || is.na(zoom_pct)) zoom_pct <- 100
  zoom_pct <- max(25, min(300, zoom_pct))

  seq_id_txt <- ""
  if (!is.null(cat$map_sequence_id) && length(cat$map_sequence_id) == 1L &&
      !is.na(cat$map_sequence_id)) {
    seq_id_txt <- paste0("  \u00b7  sequence ID ", cat$map_sequence_id)
  }

  tagList(
    .pjd_map_css_link(),
    tags$div(class = "pjd-png-view",
      tags$div(class = "pjd-png-container",
        tags$img(src = src, class = "pjd-map-large",
                 alt = paste("Vector map for", cat$name, "(Addgene)"),
                 style = sprintf("width:%.0f%%;", zoom_pct))),
      tags$div(class = "pjd-png-toolbar",
        actionButton("png_zoom_out", "\u2212", class = "pjd-view-btn"),
        actionButton("png_zoom_in", "+", class = "pjd-view-btn"),
        actionButton("png_zoom_reset", "Reset", class = "pjd-view-btn"),
        tags$span(class = "pjd-png-zoom-val", paste0(round(zoom_pct), "%"))),
      tags$div(class = "pjd-png-caption",
               "Vector map: Addgene / SnapGene \u2014 click the map to open the source record.")),
    div(class = "pjd-map-caption",
        .pjd_addgene_link(cat, paste0("Addgene #", cat$addgene_id, " \u2197")),
        seq_id_txt)
  )
}

#' Shared plasmidplot style with font sizes scaled up for the app's display.
#'
#' plasmidplot's default "angular" preset uses 9 pt labels / 14 pt title, which
#' reads small inside the ~750 px map container. Bumping every text size by
#' ~40% keeps the layout identical while making labels legible at a glance.
#' @keywords internal
.pjd_plasmidplot_style <- function() {
  plasmidplot::pp_style("angular",
                        label_fontsize   = 18,
                        site_fontsize    = 16,
                        tick_fontsize    = 13,
                        title_fontsize   = 24,
                        subtitle_fontsize = 17)
}

#' Shared stylesheet link for the PNG/plasmidplot map containers.
#' @keywords internal
.pjd_map_css_link <- function() {
  shiny::singleton(shiny::tags$head(
    shiny::tags$link(rel = "stylesheet",
                     href = file.path(PJD_MAP_ASSET_DIR, "pjd_map.css"))
  ))
}

#' Build a plasmidplot plasmid object from a pjd map record.
#'
#' Accepts the record returned by `pjd_map_record_from_genbank()` or
#' `pjd_map_record_from_assembly()` (title, size_bp, features with
#' start/end/strand/label) and rebuilds it as a plasmidplot plasmid with one
#' marker arc per feature. Strand "+" -> clockwise arrow, "-" -> counter-clockwise.
#' @keywords internal
.pjd_plasmidplot_from_record <- function(rec) {
  stopifnot(requireNamespace("plasmidplot", quietly = TRUE))
  p <- plasmidplot::plasmid(name = rec$title %||% "plasmid",
                            length = as.integer(rec$size_bp))
  feats <- rec$features
  if (is.null(feats) || nrow(feats) == 0L) return(p)
  for (i in seq_len(nrow(feats))) {
    f <- feats[i, , drop = FALSE]
    arrow <- if (identical(f$strand, "-")) "start"
             else if (identical(f$strand, "+")) "end"
             else "none"
    p <- plasmidplot::pp_marker(p, start = f$start, end = f$end,
                                label = f$label, arrow = arrow)
  }
  p
}

#' plasmidplot map placeholder (server renders the grid plot).
#'
#' Returns a `plotOutput` slot sized to the same footprint as the Addgene PNG
#' for this backbone (aspect ratio taken from the cached manifest dimensions),
#' or a fallback card with install instructions when plasmidplot is missing.
#' The actual grid drawing happens in `output$plasmidplot_map` (server), which
#' reads the selected backbone's GenBank file. Keeping this helper side-effect
#' free means the whole app still sources and runs when plasmidplot is missing.
#' @keywords internal
.pjd_plasmidplot_map_block <- function(cat) {
  if (!requireNamespace("plasmidplot", quietly = TRUE)) {
    return(div(class = "pjd-map-fallback",
               "plasmidplot is not installed. Install it once with ",
               tags$code("install.packages(\"plasmidplot\")"),
               " to enable this view; the other views keep working."))
  }
  # Match the Addgene PNG footprint: same aspect ratio, same width. Falls back
  # to a square-ish box when the manifest has no dimensions for this plasmid.
  dims <- tryCatch({
    maps <- addgene_maps_read()
    m <- maps[maps$addgene_id == cat$addgene_id, , drop = FALSE]
    if (nrow(m) > 0L && !is.na(m$width_px[[1L]]) && !is.na(m$height_px[[1L]]) &&
        m$width_px[[1L]] > 0 && m$height_px[[1L]] > 0) {
      c(m$width_px[[1L]], m$height_px[[1L]])
    } else c(NA_real_, NA_real_)
  }, error = function(e) c(NA_real_, NA_real_))
  if (any(is.na(dims))) dims <- c(907, 753)

  tagList(
    .pjd_map_css_link(),
    tags$div(class = "pjd-png-view",
      tags$div(style = sprintf("width:100%%; aspect-ratio:%d/%d;",
                               as.integer(dims[[1L]]), as.integer(dims[[2L]])),
        plotOutput("plasmidplot_map", width = "100%", height = "100%"))),
    div(class = "pjd-map-caption",
        paste0("plasmidplot (grid) rendering of ", cat$name,
               " \u2014 parsed straight from the bundled GenBank file.")),
    div(class = "pjd-map-caption",
        .pjd_addgene_link(cat, paste0("Addgene #", cat$addgene_id, " \u2197")))
  )
}

#' Unified backbone panel: identity + purpose + attR coordinates + vector map.
#'
#' Replaces the former sidebar "Backbone" / "Design purpose" cards and the
#' separate main-panel "Backbone" / "Vector map" cards — everything about the
#' selected plasmid is shown once, in the main screen.
.pjd_backbone_panel <- function(cat, gb = NULL, map_view = "plasmidplot",
                                zoom = NULL) {
  src <- if (isTRUE(cat$map_available)) addgene_map_src(cat$addgene_id) else NULL

  att_rows <- c(
    sprintf("attR3: %d..%d (%s)", cat$attR3_start, cat$attR3_end, cat$attR3_strand),
    sprintf("attR4: %d..%d (%s)", cat$attR4_start, cat$attR4_end, cat$attR4_strand),
    sprintf("attR5: %d..%d (%s)", cat$attR5_start, cat$attR5_end, cat$attR5_strand),
    sprintf("attR6: %d..%d (%s)", cat$attR6_start, cat$attR6_end, cat$attR6_strand)
  )

  div(class = "cas9-card pjd-backbone-panel",
      div(class = "cas9-section-hdr", "Selected backbone"),
      div(class = "pjd-backbone-title",
          span(class = "pjd-backbone-name", cat$name),
          span(class = "pjd-backbone-purpose", cat$purpose)),
      if (!is.null(cat$description) && !is.na(cat$description) &&
          nzchar(as.character(cat$description))) {
        div(class = "pjd-backbone-description",
            style = "font-size:12px; color:#334155; margin:6px 0 4px 0; line-height:1.45;",
            as.character(cat$description))
      },
      div(class = "pjd-badges",
          span(class = "cas9-badge cas9-badge-sky",
               paste0("Addgene ", cat$addgene_id)),
          .pjd_supplier_badge(cat),
          span(class = "cas9-badge cas9-badge-slate",
               paste0(format(cat$size_bp, big.mark = ","), " bp")),
          span(class = "cas9-badge cas9-badge-slate",
               paste0(cat$u6_promoters, " U6")),
          span(class = if (isTRUE(cat$has_attp)) "cas9-badge cas9-badge-emerald"
                      else "cas9-badge cas9-badge-slate",
               paste0("attP: ", if (isTRUE(cat$has_attp)) "yes" else "no"))
      ),
      div(class = "pjd-att-coords", paste(att_rows, collapse = "   \u00b7   ")),
      if (!is.na(cat$supplier) && cat$supplier == "DGRC_1632") {
        div(class = "pjd-note-warn",
            span(class = "cas9-badge cas9-badge-amber", "Discontinued at Addgene"),
            paste0(" pJAT55 is discontinued at Addgene; archival stock is DGRC 1632. ",
                   "Treat the DGRC glycerol stock as the archival copy."))
      },
      if (identical(map_view, "png")) {
        .pjd_addgene_map_block(cat, gb = gb, zoom_pct = zoom)
      } else {
        .pjd_plasmidplot_map_block(cat)
      }
  )
}

.pjd_pjat_picker_ui <- function(catalog) {
  choices <- setNames(catalog$name, paste0(catalog$name, "  (",
                                           catalog$size_bp, " bp, ",
                                           catalog$u6_promoters, " U6)"))
  # Only the picker stays in the sidebar; all plasmid details are rendered in the
  # consolidated main-panel card (output$backbone_panel_ui).
  tagList(
    selectInput("pjat_backbone", "pJAT backbone plasmid", choices = choices,
                selected = "pJAT74"),
    # Show the selected plasmid's purpose right under the picker so the choice
    # is informed before the card scrolls into view.
    uiOutput("pjat_picker_purpose")
  )
}

.pjd_guide_import_ui <- function() {
  tagList(
    radioButtons("guide_source", "Guide RNA source",
                 choices = c("Import Cas9 Designer bundle (cas9_guide_bundle/1.0)" = "bundle",
                             "Paste 20-nt spacer" = "paste"),
                 selected = "bundle"),
    conditionalPanel(
      condition = "input.guide_source == 'bundle'",
      fileInput("guide_bundle_csv", "Choose Cas9 guides CSV",
                accept = c(".csv")),
      checkboxInput("guide_bundle_strict", "Strict validation (reject PAM-included spacers)", value = TRUE),
      actionButton("load_sample_guide", "Load sample guide bundle (w gene)",
                   class = "btn-info btn-xs"),
      helpText(paste0("The Cas9 Designer exports 'Cas9_<gene>_<species>_guides_for_plasmid_designer.csv' ",
                      "with spacer_seq, arm_5p_seq, arm_3p_seq and a .json sidecar."))
    ),
    conditionalPanel(
      condition = "input.guide_source == 'paste'",
      textInput("spacer_paste", "20-nt spacer (PAM excluded)",
                placeholder = "e.g. GAGTCTGACCTGACTGACTG"),
      actionButton("spacer_validate_btn", "Validate spacer", class = "btn-info btn-xs")
    )
  )
}

.pjd_assembly_ui <- function() {
  tagList(
    radioButtons("insert_strategy", "Insertion strategy",
                 choices = c("Linker window (attR4_end+1 .. attR5_start-1)" = "linker",
                             "Restriction site (Type IIS)" = "enzyme",
                             "Custom coordinates" = "custom"),
                 selected = "linker"),
    conditionalPanel(
      condition = "input.insert_strategy == 'enzyme'",
      selectInput("insert_enzyme", "Type IIS enzyme",
                  choices = c("BbsI (GAAGAC)" = "BbsI",
                              "BsaI (GGTCTC)" = "BsaI",
                              "BsmBI (CGTCTC)" = "BsmBI"),
                  selected = "BbsI"),
      uiOutput("insert_enzyme_sites_ui"),
      helpText(paste0("Sites are listed for your review. No site is used automatically — ",
                      "confirm the exact cut position before assembling."))
    ),
    conditionalPanel(
      condition = "input.insert_strategy == 'custom'",
      fluidRow(
        column(6, numericInput("insert_start", "Insert start (1-based)", value = NA, min = 1)),
        column(6, numericInput("insert_end", "Insert end (1-based)", value = NA, min = 1))
      )
    ),
    radioButtons("insert_orient", "Insert orientation",
                 choices = c("Plus (+)" = "+", "Minus (-)" = "-"),
                 selected = "+"),
    checkboxInput("payload_destroy_target", "Payload destroys the gRNA target site (rule R6)", value = TRUE)
  )
}

.pjd_donor_ui <- function(catalog) {
  tagList(
    div(class = "cas9-section-hdr", "5. Donor Payload (Selectable Marker)"),
    # ONE payload, three ways to supply it.  Previously the donor textarea and
    # the marker dropdown were two parallel doors to the same single insert,
    # which was confusing: the marker button merely copied marker FASTA into
    # the textarea.  A source chooser makes the relationship explicit.
    radioButtons("donor_source", "Donor payload source",
                 choices = c("Bundled marker library" = "marker",
                             "Paste sequence" = "paste",
                             "Upload FASTA" = "upload"),
                 selected = "marker"),
    conditionalPanel(
      condition = "input.donor_source == 'marker'",
      selectInput("donor_marker", "Marker (bundled sequence, provenance tracked)",
                  choices = c("No marker (scarless edit)" = "none",
                              "3xP3-DsRed (eye marker)" = "3xP3_dsRed",
                              "ie1-DsRed (abdomen/mouthparts)" = "ie1_dsRed",
                              "Neomycin/G418 (neo)" = "neo",
                              "Puromycin (puro)" = "puro",
                              "Hygromycin (hygro)" = "hygro",
                              "Blasticidin (blast)" = "blast",
                              "URA3 (yeast)" = "URA3",
                              "LEU2 (yeast)" = "LEU2",
                              "KanMX (yeast G418)" = "KanMX"),
                  selected = "3xP3_dsRed"),
      uiOutput("donor_marker_info")
    ),
    conditionalPanel(
      condition = "input.donor_source == 'paste'",
      textAreaInput("donor_seq", "Donor DNA sequence (FASTA or raw)",
                    placeholder = "Paste FASTA or raw sequence here...",
                    rows = 6),
      actionButton("load_sample_donor", "Load sample donor FASTA (3xP3-DsRed)",
                   class = "btn-info btn-xs")
    ),
    conditionalPanel(
      condition = "input.donor_source == 'upload'",
      fileInput("donor_file", "Upload FASTA",
                accept = c(".fa", ".fasta", ".fa.gz", ".fasta.gz"))
    ),
    radioButtons("donor_orient", "Payload orientation",
                 choices = c("Plus (+)" = "+", "Minus (-)" = "-"),
                 selected = "+"),
    helpText(paste0("The payload is the sequence inserted at the chosen site. ",
                    "Orientation '+' inserts it as given; '-' inserts its ",
                    "reverse complement. A marker cassette and a donor sequence ",
                    "are the same thing \u2014 the library is just a convenient ",
                    "source for it.")),
    actionButton("run_assembly_btn", "Assemble Vector",
                 class = "btn-primary btn-lg", style = "width:100%; font-weight:700; margin-top:6px;")
  )
}

# ============================================================================
# ui
# ============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML(PJD_CUSTOM_CSS))),
  div(style = "padding: 12px 4px 6px 4px; border-bottom: 1px solid #e2e8f0; margin-bottom: 14px; display:flex; justify-content:space-between; align-items:baseline; flex-wrap:wrap;",
    div(
      span(style = "font-size:24px; font-weight:800; color:#0f172a;", "Plasmid Designer"),
      span(class = "cas9-badge cas9-badge-sky", style = "margin-left:10px;", paste0("v", PJD_APP_VERSION, " · pJAT Vector Design"))
    ),
    div(style = "font-size:12px; color:#64748b; display:flex; align-items:center; gap:8px;",
      "Janelia Atalanta (pJAT) · Gateway attL donor fragments · 250–1000 bp validated arms · SpCas9 20-nt + NGG",
      actionButton("stop_app_btn", "Stop app", icon = icon("power-off"),
                   class = "btn-danger btn-xs", title = "Stop the Plasmid Designer server")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      # 1. Backbone
      div(class = "cas9-section-hdr", "1. pJAT Backbone"),
      .pjd_pjat_picker_ui(pjat_load_catalog()),

      hr(style = "margin:10px 0;"),

      # 2. Guides
      div(class = "cas9-section-hdr", "2. Guide RNA"),
      .pjd_guide_import_ui(),

      hr(style = "margin:10px 0;"),

      # 3. Assembly strategy
      div(class = "cas9-section-hdr", "3. Vector Assembly"),
      .pjd_assembly_ui(),

      hr(style = "margin:10px 0;"),

      # 4. attL donor fragments
      div(class = "cas9-section-hdr", "4. Homology Arms (250\u20131000 bp)"),
      radioButtons("arm_source", "Homology arm source",
                   choices = c("From imported guide bundle (genomic arms \u2014 recommended)" = "bundle",
                               "Plasmid flanks around the insertion site (fallback)" = "flanks"),
                   selected = "bundle"),
      sliderInput("donor_arm_len", "Homology arm length (bp, each side)",
                  min = 250, max = 1000, value = 500, step = 50),
      helpText(paste0("Arms must be genomic sequence flanking the Cas9 cut site ",
                      "(Atalanta rules R1/R4). Importing a Cas9 Designer bundle ",
                      "supplies them automatically; the plasmid-flank fallback is ",
                      "for pasted spacers where no genomic arms exist. Arms longer ",
                      "than 1000 bp are trimmed from the outer end so the bases at ",
                      "the cut site are preserved.")),
      checkboxInput("donor_use_attl_templates",
                    "Use published Janelia Atalanta attL3\u20136 + tRNA templates (recommended)",
                    value = TRUE),
      helpText(paste0("Donor fragments carry attL3/attL4 (fragment 1) and attL5/attL6 ",
                      "(fragment 2); the pJAT backbone carries attR3\u2013attR6. LR Clonase ",
                      "converts attL \u00d7 attR \u2192 attB, so attL sites are absent from ",
                      "finished pJAT plasmids. Sequences: Stern et al. 2024 protocol.")),
      conditionalPanel(
        condition = "!input.donor_use_attl_templates",
        checkboxInput("donor_supply_attl", "I will supply my own attL templates", value = FALSE),
        conditionalPanel(
          condition = "input.donor_supply_attl",
          textInput("attl3", "attL3", placeholder = "e.g. GGGGACAAGTTTGTACAAAAAAGCAGGCT"),
          textInput("attl4", "attL4", placeholder = "e.g. GGGGACCACTTTGTACAAGAAAGCTGGGT"),
          textInput("attl5", "attL5", placeholder = "e.g. GGGGACAGCTTTCTTGTACAAAGT"),
          textInput("attl6", "attL6", placeholder = "e.g. GGGGACAACTTTGTATAGAAAAGTT")
        )
      ),

      hr(style = "margin:10px 0;"),

      # 5. Donor payload
      .pjd_donor_ui(pjat_load_catalog())
    ),

    mainPanel(
      width = 8,
      if (exists("pjd_onboarding_card", mode = "function")) pjd_onboarding_card("plasmid") else NULL,
      uiOutput("download_bar_ui"),
      div(class = "pjd-map-view-row",
          radioButtons("map_view", "Vector map view",
                       choices = c("plasmidplot" = "plasmidplot",
                                   "Addgene PNG" = "png"),
                       selected = "plasmidplot", inline = TRUE)),
      uiOutput("backbone_panel_ui"),
      uiOutput("design_schematic_ui"),
      uiOutput("guide_info_ui"),
      uiOutput("assembly_card_ui"),
      uiOutput("assembly_map_ui"),
      uiOutput("donor_card_ui"),
      uiOutput("sequence_viewer_ui")
    )
  )
)
# ============================================================================
# server
# ============================================================================

server <- function(input, output, session) {

  # ---- 1. Backbone reactive + info card -----------------------------------
  pjd_backbone <- reactive({
    catalog <- pjat_load_catalog()
    cat <- catalog[catalog$name == input$pjat_backbone, , drop = FALSE]
    if (nrow(cat) == 0L) stop("Unknown pJAT plasmid: ", input$pjat_backbone)
    cat[1L, , drop = FALSE]
  })

  # Parsed GenBank for the selected backbone (drives the interactive map).
  pjd_backbone_gb <- reactive({
    b <- pjd_backbone()
    suppressWarnings(read_pjat_genbank(file.path(DEFAULT_PJAT_DIR, b$source_file)))
  })

  # PNG zoom state (server side; only the PNG view reads this).
  png_zoom <- reactiveVal(100)

  observeEvent(input$png_zoom_in, { png_zoom(min(300, png_zoom() + 25)) }, ignoreInit = TRUE)
  observeEvent(input$png_zoom_out, { png_zoom(max(25, png_zoom() - 25)) }, ignoreInit = TRUE)
  observeEvent(input$png_zoom_reset, { png_zoom(100) }, ignoreInit = TRUE)

  # One consolidated card in the main panel (sidebar keeps only the picker).
  output$backbone_panel_ui <- renderUI({
    mv <- input$map_view %||% "plasmidplot"
    .pjd_backbone_panel(pjd_backbone(), gb = pjd_backbone_gb(),
                        map_view = mv, zoom = png_zoom())
  })

  # Purpose + description under the sidebar picker, so the plasmid's role is
  # visible while choosing (answers "what is this plasmid for?").
  output$pjat_picker_purpose <- renderUI({
    b <- tryCatch(pjd_backbone(), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    div(style = "font-size:11px; color:#475569; margin:-6px 0 8px 0; line-height:1.4;",
        strong(paste0(b$name, ": ")), as.character(b$purpose %||% ""))
  })

  # plasmidplot grid rendering (main map view; needs the CRAN package).
  output$plasmidplot_map <- renderPlot({
    req(requireNamespace("plasmidplot", quietly = TRUE))
    b <- pjd_backbone()
    gb_path <- file.path(DEFAULT_PJAT_DIR, b$source_file)
    req(file.exists(gb_path))
    p <- plasmidplot::read_genbank(gb_path)
    p$name <- b$name   # show the catalog name (pJAT74), not the SnapGene locus
    plot(p, style = .pjd_plasmidplot_style())
  })

  # ---- 2b. Live construct-design schematic ---------------------------------
  # Renders the planned construct live as inputs change, so direction, the
  # landing slot and an orientation flip are visible before assembling.
  output$design_schematic_ui <- renderUI({
    b <- tryCatch(pjd_backbone(), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    payload <- tryCatch(pjd_donor_payload(), error = function(e) NULL)
    slot <- switch(input$insert_strategy %||% "linker",
      linker = "linker window (attR4\u2013attR5)",
      enzyme = paste0("Type IIS site: ", sub(" .*", "", input$insert_enzyme %||% "BbsI")),
      custom = "custom coordinates",
      "insertion site")
    has_trna <- any(grepl("trna|tRNA",
                          paste(b$source_file, b$purpose %||% "", collapse = " ")))
    tagList(
      div(class = "cas9-card pjd-design-card-wrap",
          div(class = "cas9-section-hdr", "Construct design (live)"),
          pjd_design_schematic(
            backbone_name = b$name,
            spacer = tryCatch(pjd_selected_spacer(), error = function(e) NULL),
            payload_seq = if (is.null(payload)) NULL else payload$sequence,
            payload_label = if (is.null(payload)) "payload" else payload$label,
            orient = input$donor_orient %||% "+",
            slot_label = slot,
            attr_left = "attR4", attr_right = "attR5",
            grna_site = pjd_grna_site_label(b$u6_promoters, has_trna),
            assembled = !is.null(pjd_assembly())
          ),
          div(class = "pjd-map-caption",
              paste0("Backbone ", b$name, " \u00b7 payload lands at ", slot,
                     " \u00b7 guide site: ",
                     pjd_grna_site_label(b$u6_promoters, has_trna)))
      )
    )
  })

  # Stop button: kill this app's R session so the browser goes blank and the
  # user can relaunch cleanly (code changes only take effect on restart).
  observeEvent(input$stop_app_btn, {
    showNotification("Stopping Plasmid Designer…", type = "message", duration = 2)
    Sys.sleep(0.5)
    stopApp()
  }, ignoreInit = TRUE)

  # ---- 2. Guide import ----------------------------------------------------
  pjd_guide <- reactiveVal(NULL)

  observeEvent(input$guide_bundle_csv, {
    req(input$guide_bundle_csv$datapath)
    tryCatch({
      b <- read_guide_bundle(input$guide_bundle_csv$datapath,
                             strict = isTRUE(input$guide_bundle_strict))
      pjd_guide(b)
      showNotification(sprintf("Imported %d guide(s) from Cas9 bundle (%s)",
                               b$bundle_meta$n_guides, b$bundle_meta$species), type = "message")
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  spacer_valid <- reactiveVal(NULL)

  observeEvent(input$load_sample_guide, {
    tryCatch({
      b <- pjd_load_sample_guide_bundle()
      pjd_guide(b)
      showNotification(
        sprintf("Loaded sample bundle: %d guide(s) \u00b7 spacer %s",
                b$bundle_meta$n_guides, b$guides$spacer_seq[1L]),
        type = "message", duration = 6)
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  observeEvent(input$spacer_validate_btn, {
    sp <- trimws(input$spacer_paste %||% "")
    ok <- tryCatch(pjat_validate_spacer(sp), error = function(e) FALSE)
    spacer_valid(list(spacer = if (isTRUE(ok)) sp else NULL, ok = isTRUE(ok)))
    showNotification(if (isTRUE(ok)) "Spacer valid (20 nt, PAM excluded)"
                     else "Invalid spacer: must be exactly 20 nt, PAM excluded",
                     type = if (isTRUE(ok)) "message" else "error", duration = 5)
  })

  pjd_selected_spacer <- reactive({
    g <- pjd_guide()
    if (!is.null(g) && nrow(g$guides) > 0L) {
      g$guides$spacer_seq[1L]
    } else if (!is.null(spacer_valid()) && isTRUE(spacer_valid()$ok)) {
      spacer_valid()$spacer
    } else {
      NULL
    }
  })

  output$guide_info_ui <- renderUI({
    g <- pjd_guide()
    sp <- pjd_selected_spacer()
    if (is.null(g) && is.null(sp)) return(NULL)
    tagList(
      div(class = "cas9-card",
          h4("Guide RNA", style = "margin-top:0;"),
          if (!is.null(g)) {
            div(style = "font-size:12px;",
                paste("Imported", nrow(g$guides), "guide(s) from",
                      basename(g$bundle_meta$imported_from),
                      "| provenance:", g$provenance))
          },
          div(style = "font-size:12px; color:#475569; margin-top:4px;",
              paste("Spacer for assembly:", if (is.null(sp)) "<none>" else sp))
      )
    )
  })

  # ---- 3. Enzyme site list ------------------------------------------------
  output$insert_enzyme_sites <- renderPrint({
    req(input$insert_strategy == "enzyme")
    b <- pjd_backbone()
    gb <- suppressWarnings(read_pjat_genbank(
      file.path(DEFAULT_PJAT_DIR, b$source_file)))
    enz <- sub(" .*", "", input$insert_enzyme)
    sites <- find_cloning_sites(gb$sequence, enzymes = CLONING_ENZYMES[enz])
    if (nrow(sites) == 0L) {
      cat(sprintf("No %s site found in backbone\n", enz))
      return(invisible())
    }
    print(sites[, c("position", "cut_pos", "overhang", "enzyme")], row.names = FALSE)
  })

  # ---- 4. Assembly --------------------------------------------------------
  pjd_assembly <- reactiveVal(NULL)

  observeEvent(input$run_assembly_btn, {
    b <- pjd_backbone()
    gb <- suppressWarnings(read_pjat_genbank(
      file.path(DEFAULT_PJAT_DIR, b$source_file)))
    sp <- pjd_selected_spacer()
    if (is.null(sp)) {
      showNotification("No spacer/guide available. Import a bundle or paste a spacer.",
                       type = "error", duration = 5)
      return()
    }

    strategy <- input$insert_strategy
    if (identical(strategy, "enzyme")) {
      enz <- sub(" .*", "", input$insert_enzyme)
      sites <- find_cloning_sites(gb$sequence,
                                  enzymes = CLONING_ENZYMES[enz])
      hits <- sites[!grepl(" \\(rc\\)", sites$site), ]
      if (nrow(hits) == 0L) {
        showNotification("No restriction site found", type = "error", duration = 5)
        return()
      }
      cut <- hits$position[1L]
      start <- cut; end <- cut - 1L   # insert-only at cut site
    } else if (identical(strategy, "custom")) {
      start <- input$insert_start; end <- input$insert_end
      if (is.na(start) || is.na(end) || end < start) {
        showNotification("Invalid custom coordinates", type = "error", duration = 5)
        return()
      }
    } else {
      start <- NULL; end <- NULL
    }

    # Build insert from the single resolved payload (section 5 source chooser).
    insert <- NULL
    insert_label <- "user insert"
    payload <- pjd_donor_payload()
    if (!is.null(payload)) {
      insert <- payload$sequence
      insert_label <- payload$label
    }

    if (is.null(insert) && identical(input$insert_strategy, "linker") &&
        is.null(input$guide_source)) {
      insert <- assemble_gRNA_cassette(sp)
      insert_label <- "gRNA cassette"
    }

    asm <- tryCatch(
      assemble_expression_vector(gb, insert = insert, spacer = sp,
                                 strategy = strategy,
                                 start = start, end = end,
                                 orient = input$insert_orient),
      error = function(e) { showNotification(conditionMessage(e), type = "error",
                                             duration = 10); NULL }
    )
    if (is.null(asm)) return()
    asm$insert_label <- insert_label

    # attL donor fragments (rule R2: 250-1000 bp arms).
    #
    # Arms MUST be genomic sequence flanking the Cas9 cut site (rules R1/R4), so
    # the imported cas9_guide_bundle's arms are the default source.  Plasmid
    # flanks are only a fallback for pasted spacers, where no genomic arms exist.
    arm_len <- as.integer(input$donor_arm_len)
    flank_arms <- tryCatch({
      seq_n <- nchar(asm$sequence)
      list(left = substr(asm$sequence, max(1L, asm$insert_start - arm_len),
                         asm$insert_start - 1L),
           right = substr(asm$sequence, asm$insert_end + 1L,
                          min(seq_n, asm$insert_end + arm_len)))
    }, error = function(e) list(left = NULL, right = NULL))

    arms <- tryCatch(
      resolve_homology_arms(
        bundle   = pjd_guide(),
        guide_id = if (!is.null(pjd_guide())) pjd_guide()$guides$guide_id[1] else NULL,
        flank_5p = flank_arms$left,
        flank_3p = flank_arms$right,
        source   = if (identical(input$arm_source, "flanks")) "flanks" else "auto",
        trim_to  = arm_len
      ),
      error = function(e) {
        showNotification(paste("Arm resolution failed:", conditionMessage(e)),
                         type = "error", duration = 10); NULL
      }
    )
    if (is.null(arms)) return()
    left_arm  <- arms$left_arm
    right_arm <- arms$right_arm
    if (length(arms$messages)) {
      showNotification(paste(arms$messages, collapse = "; "),
                       type = "message", duration = 8)
    }

    attl <- if (isTRUE(input$donor_use_attl_templates)) {
      # Published Janelia Atalanta templates (attL3-6 + tRNA), bundled locally.
      tryCatch(
        pjat_attl_template_sequences(),
        error = function(e) {
          showNotification(paste("attL templates unavailable:",
                                 conditionMessage(e)),
                           type = "warning", duration = 10)
          list(attl3 = NULL, attl4 = NULL, attl5 = NULL, attl6 = NULL,
               trna = NULL)
        })
    } else if (isTRUE(input$donor_supply_attl)) {
      list(attl3 = trimws(input$attl3 %||% ""), attl4 = trimws(input$attl4 %||% ""),
           attl5 = trimws(input$attl5 %||% ""), attl6 = trimws(input$attl6 %||% ""),
           trna = NULL)
    } else {
      list(attl3 = NULL, attl4 = NULL, attl5 = NULL, attl6 = NULL, trna = NULL)
    }

    frags <- tryCatch(
      pjat_assemble_fragments(left_arm, right_arm, sp,
                              attl3 = attl$attl3, attl4 = attl$attl4,
                              attl5 = attl$attl5, attl6 = attl$attl6,
                              trna = attl$trna),
      error = function(e) { showNotification(conditionMessage(e), type = "error",
                                             duration = 10); NULL }
    )

    pjd_assembly(list(backbone = b, genbank = gb, assembly = asm,
                      left_arm = left_arm, right_arm = right_arm,
                      fragments = frags, spacer = sp,
                      orient = input$insert_orient,
                      donor_provenance = if (is.null(insert)) "none"
                      else paste0("user-supplied (", nchar(insert), " bp)")))
    showNotification("Vector assembled.", type = "message", duration = 3)
  })

  # ---- 4b. Unified donor payload (single source of truth) -------------------
  # Resolves ONE payload from the chosen source so the assembly never has to
  # decide between a textarea, an upload and a marker library.
  pjd_donor_payload <- reactive({
    src <- input$donor_source %||% "marker"
    if (identical(src, "marker")) {
      id <- input$donor_marker %||% "none"
      if (identical(id, "none")) return(NULL)
      if (!pjd_marker_sequence_available(id)) return(NULL)
      m <- tryCatch(pjd_get_marker_sequence(id), error = function(e) NULL)
      if (is.null(m)) return(NULL)
      return(list(sequence = toupper(gsub("[^ACGTN]", "", m$sequence)),
                  label = m$label %||% id,
                  provenance = paste0("bundled marker library: ", id,
                                      " (", m$length_bp, " bp)")))
    }
    if (identical(src, "paste")) {
      txt <- trimws(input$donor_seq %||% "")
      if (!nzchar(txt)) return(NULL)
      seq <- toupper(gsub("[^ACGTN]", "", txt))
      if (!nzchar(seq)) return(NULL)
      return(list(sequence = seq, label = "pasted donor sequence",
                  provenance = paste0("user-pasted (", nchar(seq), " bp)")))
    }
    # upload
    f <- input$donor_file
    if (is.null(f) || !nzchar(f$datapath %||% "")) return(NULL)
    fa <- readLines(f$datapath, warn = FALSE)
    fa <- fa[!grepl("^>", fa)]
    seq <- toupper(gsub("[^ACGTN]", "", paste(fa, collapse = "")))
    if (!nzchar(seq)) return(NULL)
    list(sequence = seq, label = basename(f$name),
         provenance = paste0("uploaded FASTA ", basename(f$name),
                             " (", nchar(seq), " bp)"))
  })

  # Marker info card (provenance + bundled status) for the library source.
  output$donor_marker_info <- renderUI({
    if (!identical(input$donor_source %||% "marker", "marker")) return(NULL)
    mk <- read_yaml_marker_catalog()
    row <- mk[mk$id == input$donor_marker, , drop = FALSE]
    if (nrow(row) == 0L || identical(input$donor_marker, "none")) return(NULL)
    bundled <- identical(row$sequence_status[1], "bundled") &&
      !is.na(row$length_bp[1])
    div(class = "cas9-card", style = "border-left-color:#0ea5e9;",
        h4("Marker info", style = "margin-top:0;"),
        div(style = "font-size:12px;",
            paste(row$description[1], "| Selection:", row$selection[1])),
        div(style = "font-size:11px; color:#64748b; margin-top:4px;",
            paste("Source:", row$source[1],
                  "| Sequence status:", row$sequence_status[1])),
        div(style = "font-size:11px; margin-top:4px;",
            if (bundled) {
              span(class = "cas9-badge cas9-badge-emerald",
                   paste0("sequence bundled \u00b7 ",
                          format(row$length_bp[1], big.mark = ","), " bp"))
            } else {
              span(class = "cas9-badge cas9-badge-amber", "sequence not bundled")
            }),
        div(style = "font-size:11px; color:#475569; margin-top:4px;",
            paste("Provenance:", if (is.na(row$sequence_provenance[1]) ||
                                        row$sequence_provenance[1] == "unset")
                  "unset \u2014 will block export" else row$sequence_provenance[1]))
    )
  })

  # ---- 5. Sample inputs ----------------------------------------------------
  # (The marker info card and the payload resolver live in section 4b above;
  #  this handler only fills the paste box from the bundled sample FASTA.)
  observeEvent(input$load_sample_donor, {
    tryCatch({
      d <- pjd_load_sample_donor_fasta()
      updateTextAreaInput(session, "donor_seq", value = d$sequence)
      showNotification(
        sprintf("Loaded sample donor: %d bp \u00b7 %s", d$length_bp, d$provenance),
        type = "message", duration = 8)
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  # ---- 6. Assembly summary card -------------------------------------------
  output$assembly_card_ui <- renderUI({
    a <- pjd_assembly()
    if (is.null(a)) return(NULL)
    f <- a$fragments
    tagList(
      div(class = "cas9-card",
          h4("Assembly summary", style = "margin-top:0;"),
          div(style = "font-size:12px;",
              paste("Backbone:", a$backbone$name, "(", format(a$backbone$size_bp,
                      big.mark = ","), "bp)  |  Insert:", a$assembly$insert_len, "bp  |",
                    "Final:", a$assembly$final_len, "bp")),
          div(style = "font-size:12px; color:#475569; margin-top:4px;",
              a$assembly$assembly_summary),
          div(style = "font-size:12px; margin-top:6px;",
              strong("attL donor fragments:  "),
              if (is.null(f)) span(class = "cas9-badge cas9-badge-red", "not built")
              else tagList(
                span(class = if (identical(f$status, "READY")) "cas9-badge cas9-badge-emerald" else "cas9-badge cas9-badge-amber",
                     f$status),
                if (!is.null(f) && !isTRUE(f$checks$attl_sites_supplied)) {
                  span(class = "cas9-badge cas9-badge-amber", "attL placeholders only")
                } else NULL
              )),
          div(style = "font-size:11px; color:#64748b; margin-top:4px;",
              paste("Junction left:", a$assembly$junction_left, "|",
                    "Junction right:", a$assembly$junction_right))
      )
    )
  })

  # ---- 7. Donor card ------------------------------------------------------
  # ---- 6b. Plasmidplot map of the assembled construct ----------------------
  output$assembly_map_ui <- renderUI({
    a <- pjd_assembly()
    if (is.null(a)) return(NULL)
    if (!requireNamespace("plasmidplot", quietly = TRUE)) return(NULL)
    rec <- tryCatch(
      pjd_map_record_from_assembly(a$assembly, backbone = a$backbone$name,
                                   backbone_gb = a$genbank),
      error = function(e) NULL)
    if (is.null(rec)) return(NULL)
    div(class = "cas9-card",
        h4("Assembled vector map", style = "margin-top:0;"),
        tags$div(style = "width:100%; aspect-ratio:1/1;",
                 plotOutput("assembly_plot", width = "100%", height = "100%")),
        div(class = "pjd-map-caption",
            paste0("Backbone ", a$backbone$name, " + ",
                   a$assembly$insert_len, " bp insert = ",
                   format(rec$size_bp, big.mark = ","), " bp final")))
  })

  output$assembly_plot <- renderPlot({
    a <- pjd_assembly()
    req(a)
    req(requireNamespace("plasmidplot", quietly = TRUE))
    rec <- pjd_map_record_from_assembly(a$assembly, backbone = a$backbone$name)
    plot(.pjd_plasmidplot_from_record(rec), style = .pjd_plasmidplot_style())
  })

  output$donor_card_ui <- renderUI({
    a <- pjd_assembly()
    if (is.null(a)) return(NULL)
    qc <- a$assembly$qc
    # pjat_synthesis_qc() reports `n`, `gc_frac`; older wrappers used len/gc_pct.
    qc_len <- if (!is.null(qc$len)) qc$len else qc$n
    gc_pct <- if (!is.null(qc$gc_pct)) qc$gc_pct else
      if (!is.null(qc$gc_frac)) 100 * qc$gc_frac else NA_real_
    tagList(
      div(class = "cas9-card",
          h4("Insert QC", style = "margin-top:0;"),
          div(style = "font-size:12px;",
              paste("Length:", qc_len, "bp | GC:",
                    if (is.na(gc_pct)) "n/a" else paste0(round(gc_pct, 1), "%"),
                    "| Max homopolymer:", qc$max_homopolymer)),
          div(style = "font-size:12px; color:#475569; margin-top:4px;",
              paste(if (isTRUE(qc$ok)) "QC: PASS" else "QC: REVIEW",
                    if (length(qc$warnings)) paste("Warnings:", paste(qc$warnings, collapse = "; "))
                    else "")
      ))
    )
  })

  # ---- 8. Sequence viewer -------------------------------------------------
  output$sequence_viewer_ui <- renderUI({
    a <- pjd_assembly()
    if (is.null(a)) return(NULL)
    tagList(
      div(class = "cas9-card",
          h4("Plasmid sequence (insert highlighted region)", style = "margin-top:0;"),
          pre(style = "background:#f1f5f9; padding:10px; border-radius:4px; font-size:10px; overflow:auto; max-height:400px;",
              substr(a$assembly$sequence,
                     max(1L, a$assembly$insert_start - 50L),
                     min(nchar(a$assembly$sequence), a$assembly$insert_end + 50L)))
      )
    )
  })

  # ---- 9. Download bar ----------------------------------------------------
  output$download_bar_ui <- renderUI({
    a <- pjd_assembly()
    if (is.null(a)) return(NULL)
    div(class = "cas9-dl-bar",
        div(style = "font-size:12px; color:#0f172a; font-weight:600;",
            "Export assembly"),
        div(style = "display:flex; gap:6px; flex-wrap:wrap;",
            downloadButton("dl_plasmid_gb", "plasmid.gb", class = "btn-xs"),
            downloadButton("dl_plasmid_fasta", "plasmid.fasta", class = "btn-xs"),
            downloadButton("dl_attl_fasta", "attL_fragments.fasta", class = "btn-xs"),
            downloadButton("dl_assembly_json", "assembly.json", class = "btn-xs"),
            downloadButton("dl_assembly_zip", "bundle.zip", class = "btn-xs")
        )
    )
  })

  file_prefix <- reactive({
    a <- pjd_assembly()
    sym <- if (!is.null(a) && !is.null(a$backbone)) a$backbone$name else "pJAT"
    gsub("[^A-Za-z0-9_\\-]", "_", paste0("Plasmid_", sym, "_", Sys.Date()))
  })

  output$dl_plasmid_gb <- downloadHandler(
    filename = function() paste0(file_prefix(), "_plasmid.gb"),
    content = function(file) {
      a <- pjd_assembly(); req(a)
      feats <- a$assembly$features
      if (is.null(feats) || nrow(feats) == 0L) {
        feats <- data.frame(type = "insert", start = a$assembly$insert_start,
                            end = a$assembly$insert_end, strand = a$assembly$orient,
                            label = a$assembly$insert_label %||% "user_insert", stringsAsFactors = FALSE)
      }
      write_plasmid_genbank(a$assembly$sequence, feats,
                            name = paste0(a$backbone$name, "_assembly"),
                            topology = "circular", file = file)
    }
  )

  output$dl_plasmid_fasta <- downloadHandler(
    filename = function() paste0(file_prefix(), "_plasmid.fasta"),
    content = function(file) {
      a <- pjd_assembly(); req(a)
      seq <- a$assembly$sequence
      hdr <- sprintf(">%s|%s|circular|len=%dbp|insert=%dbp",
                     a$backbone$name, a$assembly$strategy,
                     nchar(seq), a$assembly$insert_len)
      wrapped <- paste0(strwrap(seq, width = 80, simplify = FALSE), collapse = "\n")
      writeLines(c(hdr, wrapped), con = file)
    }
  )

  output$dl_attl_fasta <- downloadHandler(
    filename = function() paste0(file_prefix(), "_attL_fragments.fasta"),
    content = function(file) {
      a <- pjd_assembly(); req(a)
      f <- a$fragments; req(f)
      lines <- c(
        sprintf(">fragment1_attL3_arm_trna_attL4|status=%s|len=%dbp",
                f$status, nchar(f$fragments$fragment1_attL3_arm_trna_attL4)),
        f$fragments$fragment1_attL3_arm_trna_attL4,
        sprintf(">fragment2_attL5_arm_attL6|status=%s|len=%dbp",
                f$status, nchar(f$fragments$fragment2_attL5_arm_attL6)),
        f$fragments$fragment2_attL5_arm_attL6
      )
      writeLines(lines, con = file)
    }
  )

  output$dl_assembly_json <- downloadHandler(
    filename = function() paste0(file_prefix(), "_assembly.json"),
    content = function(file) {
      a <- pjd_assembly(); req(a)
      out <- list(
        schema = "plasmid_assembly/1.0",
        app_version = PJD_APP_VERSION,
        backbone = as.list(a$backbone),
        spacer = a$spacer,
        insert_len = a$assembly$insert_len,
        final_len = a$assembly$final_len,
        strategy = a$assembly$strategy,
        orient = a$assembly$orient,
        junction_left = a$assembly$junction_left,
        junction_right = a$assembly$junction_right,
        qc = a$assembly$qc,
        donor_provenance = a$donor_provenance,
        attl_status = if (is.null(a$fragments)) NULL else a$fragments$status
      )
      json_txt <- jsonlite::toJSON(out, pretty = TRUE, auto_unbox = TRUE)
      writeLines(json_txt, con = file)
    }
  )

  output$dl_assembly_zip <- downloadHandler(
    filename = function() paste0(file_prefix(), "_bundle.zip"),
    content = function(file) {
      a <- pjd_assembly(); req(a)
      tmp <- tempdir()
      gb <- file.path(tmp, "plasmid.gb")
      fa <- file.path(tmp, "plasmid.fasta")
      af <- file.path(tmp, "attL_fragments.fasta")
      js <- file.path(tmp, "assembly.json")
      write_plasmid_genbank(a$assembly$sequence, a$assembly$features,
                            name = paste0(a$backbone$name, "_assembly"),
                            topology = "circular", file = gb)
      seq <- a$assembly$sequence
      hdr <- sprintf(">%s|circular|len=%dbp", a$backbone$name, nchar(seq))
      writeLines(c(hdr, paste0(strwrap(seq, width = 80, simplify = FALSE), collapse = "\n")), fa)
      f <- a$fragments
      if (!is.null(f)) {
        writeLines(c(
          sprintf(">fragment1|len=%dbp", nchar(f$fragments$fragment1_attL3_arm_trna_attL4)),
          f$fragments$fragment1_attL3_arm_trna_attL4,
          sprintf(">fragment2|len=%dbp", nchar(f$fragments$fragment2_attL5_arm_attL6)),
          f$fragments$fragment2_attL5_arm_attL6
        ), af)
      }
      out <- list(schema = "plasmid_assembly/1.0", backbone = as.list(a$backbone),
                  spacer = a$spacer, insert_len = a$assembly$insert_len,
                  final_len = a$assembly$final_len, strategy = a$assembly$strategy,
                  orient = a$assembly$orient, qc = a$assembly$qc,
                  donor_provenance = a$donor_provenance)
      writeLines(jsonlite::toJSON(out, pretty = TRUE, auto_unbox = TRUE), js)
      utils::zip(zipfile = file, files = c(gb, fa, af, js), flags = "-q")
    }
  )
}

# ============================================================================
# Headless CLI + interactive Shiny entry point
# ============================================================================

.parse_pjd_cli_flag <- function(args, flag, default = NULL) {
  eq_pat <- paste0("^", flag, "=")
  hit_eq <- grep(eq_pat, args, value = TRUE)
  if (length(hit_eq)) return(sub(eq_pat, "", hit_eq[1]))
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1L])
  default
}

.run_pjd_cli <- function(args) {
  gene <- .parse_pjd_cli_flag(args, "--gene")
  species <- .parse_pjd_cli_flag(args, "--species", "melanogaster")
  out <- .parse_pjd_cli_flag(args, "--out", file.path(getwd(), "pjd_out"))
  if (is.null(gene)) {
    stop("CLI mode requires --gene (and optionally --species, --out).")
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  message("Plasmid Designer CLI: gene=", gene, " species=", species, " out=", out)
  invisible(NULL)
}

.pjd_port_in_use <- function(port) {
  con <- tryCatch(
    socketConnection("127.0.0.1", port = port, open = "r+",
                     blocking = TRUE, timeout = 1),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(con)) return(FALSE)
  close(con)
  TRUE
}

# Standalone-execution dispatch (Rscript / `R -f`). When this file is sourced
# by the Electron runner (shared/electron/runner.R), sys.nframe() > 0 and none
# of this runs: sourcing only defines ui/server, and the runner hosts the app
# itself on a random port (printing APP_PORT= for the port handshake).
if (sys.nframe() == 0L) {
  .args <- commandArgs(trailingOnly = TRUE)
  if (any(grepl("^(--help|-h)$", .args))) {
    cat("Usage: Rscript Plasmid_Designer.R [--gene <symbol>] [--species <species>] [--out <dir>]\n")
    invisible(NULL)
  } else if (length(.args) > 0L && !any(grepl("^--file=", .args))) {
    .run_pjd_cli(.args)
  } else if (.pjd_port_in_use(PORT)) {
    # A server is already listening (e.g. the app was double-clicked twice):
    # just bring the user to it instead of failing with "port in use".
    message("Plasmid Designer is already running on port ", PORT,
            " — opening http://127.0.0.1:", PORT, " in your browser.")
    utils::browseURL(sprintf("http://127.0.0.1:%d", PORT))
  } else {
    shiny::shinyApp(ui, server,
                    options = list(port = PORT, launch.browser = TRUE))
  }
}
