# ============================================================================
# R/design_schematic.R — Live construct-design schematic (main window)
# ============================================================================
# Shows the donor / guide / marker BEFORE assembly so the user can see the
# planned construct and verify direction; nothing surprising should appear only
# after clicking Assemble Vector.  Built on the same dependency-free SVG
# approach as plasmid_map.R.
#
# The schematic is a block diagram, not a coordinate map: it answers
#   * what is the spacer?
#   * where does the payload land (slot) and which way round is it?
#   * what does the payload look like after an orientation flip (+/-)?
#   * which attR/attL sites flank the payload?
#   * does the backbone have a gRNA insertion site (U6 / tRNA array)?
# ============================================================================

.pjd_design_block <- function(label, sub = NULL, fill = "#0ea5e9", width = 150,
                              min_width = 90) {
  w <- max(min_width, width)
  shiny::tags$g(
    shiny::tags$rect(x = 0, y = 0, width = w, height = 44, rx = 5,
                     fill = fill, `fill-opacity` = "0.14",
                     stroke = fill, `stroke-width` = "1.4"),
    shiny::tags$text(x = w / 2, y = 19, `text-anchor` = "middle",
                     class = "pjd-design-label", label),
    if (!is.null(sub)) {
      shiny::tags$text(x = w / 2, y = 34, `text-anchor` = "middle",
                       class = "pjd-design-sub", sub)
    }
  )
}

# Reverse-complement helper (matches insert_cassette's semantics).
.pjd_design_revcomp <- function(seq) {
  if (is.null(seq) || !nzchar(seq)) return("")
  comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  paste0(rev(comp[strsplit(toupper(seq), "")[[1L]]]), collapse = "")
}

#' Build the live construct-design schematic.
#'
#' @param backbone_name Character(1) backbone label.
#' @param spacer Character(1) 20-nt spacer, or NULL.
#' @param payload_seq Character(1) donor payload, or NULL.
#' @param payload_label Character(1) label for the payload block.
#' @param orient "+" or "-" payload orientation.
#' @param slot_label Character(1) where the payload lands.
#' @param attr_left,attr_right Character(1) flanking recombination sites.
#' @param grna_site Character(1) backbone gRNA insertion-site descriptor.
#' @param assembled Logical: has the vector been assembled yet?
#' @return shiny tagList (SVG).
#' @export
pjd_design_schematic <- function(backbone_name = "backbone",
                                 spacer = NULL, payload_seq = NULL,
                                 payload_label = "payload", orient = "+",
                                 slot_label = "linker window",
                                 attr_left = "attR4", attr_right = "attR5",
                                 grna_site = "unknown",
                                 assembled = FALSE) {
  orient <- if (identical(orient, "-")) "-" else "+"
  payload_len <- if (is.null(payload_seq)) 0L else nchar(payload_seq)
  payload_shown <- if (is.null(payload_seq)) "" else
    paste0(substr(payload_seq, 1L, 10L), "\u2026",
           substr(payload_seq, max(1L, payload_len - 9L), payload_len))
  payload_final <- if (is.null(payload_seq)) "" else
    if (identical(orient, "-")) .pjd_design_revcomp(payload_seq) else payload_seq
  final_shown <- if (!nzchar(payload_final)) "" else
    paste0(substr(payload_final, 1L, 10L), "\u2026",
           substr(payload_final, max(1L, nchar(payload_final) - 9L),
                  nchar(payload_final)))

  # Horizontal layout: blocks placed left to right with fixed gaps.
  x <- 10
  nodes <- list()
  add <- function(node, width) {
    nodes[[length(nodes) + 1L]] <<- shiny::tags$g(
      transform = sprintf("translate(%d,0)", x), node)
    x <<- x + width + 14
  }

  add(.pjd_design_block(backbone_name, "destination vector",
                        fill = "#64748b", width = 170), 170)
  add(.pjd_design_block(attr_left, "attR site", fill = "#b45309",
                        width = 100), 100)
  add(.pjd_design_block(payload_label, paste0(payload_len, " bp payload"),
                        fill = "#db2777", width = 150), 150)
  add(.pjd_design_block(slot_label, "landing slot", fill = "#7c3aed",
                        width = 150), 150)
  add(.pjd_design_block(attr_right, "attR site", fill = "#b45309",
                        width = 100), 100)
  add(.pjd_design_block("payload orientation", orient,
                        fill = if (orient == "+") "#059669" else "#dc2626",
                        width = 150), 150)

  total_w <- x + 10

  shiny::tags$div(
    class = "pjd-design-card",
    shiny::tags$svg(
      xmlns = "http://www.w3.org/2000/svg",
      viewBox = sprintf("0 0 %d 130", total_w),
      class = "pjd-design-svg", style = "width:100%; height:130px;",
      shiny::tags$line(x1 = 10, y1 = 78, x2 = total_w - 10, y2 = 78,
                       stroke = "#cbd5e1", `stroke-width` = "2"),
      nodes,
      shiny::tags$text(x = 10, y = 108, class = "pjd-design-sub",
                       paste0("gRNA insertion site: ", grna_site)),
      if (!is.null(spacer) && nzchar(spacer)) {
        shiny::tags$text(x = 10, y = 124, class = "pjd-design-sub",
                         paste0("spacer (20 nt): ", spacer))
      } else NULL
    ),
    if (nzchar(payload_shown)) {
      shiny::tags$div(
        class = "pjd-design-readout",
        style = "font-size:11px; color:#475569; margin-top:4px; line-height:1.5;",
        shiny::tags$div(paste0("payload as supplied (", orient, "): ", payload_shown)),
        if (identical(orient, "-")) {
          shiny::tags$div(style = "color:#b91c1c; font-weight:600;",
                          paste0("inserted as (reverse complement): ", final_shown))
        } else {
          shiny::tags$div(style = "color:#047857;",
                          "inserted as supplied (no reverse complement)")
        }
      )
    } else NULL,
    if (assembled) {
      shiny::tags$div(style = "font-size:11px; color:#0369a1; margin-top:4px;",
                      "Assembled \u2014 the vector map below shows the final construct.")
    } else NULL
  )
}


#' Describe a backbone's gRNA insertion site for the design schematic.
#'
#' Answers "does this vector have somewhere for the guide to go?" — typically a
#' U6 promoter (optionally a tRNA array).  Derived from the catalog's
#' `u6_promoters` count and the registry's grna_array slots.
#'
#' @param u6_count Integer number of U6 promoters (may be NA).
#' @param has_trna Logical: does the backbone carry a tRNA array feature?
#' @return Character(1) human-readable descriptor.
#' @export
pjd_grna_site_label <- function(u6_count = NA_integer_, has_trna = FALSE) {
  if (length(u6_count) == 0L || is.na(u6_count)) {
    return("unknown (backbone not characterised)")
  }
  if (u6_count <= 0L) {
    return("none \u2014 add a U6 promoter (empty backbone, e.g. pJAT77)")
  }
  lbl <- sprintf("%d \u00d7 U6-3 promoter%s", u6_count,
                 if (u6_count > 1L) "s" else "")
  if (isTRUE(has_trna)) lbl <- paste0(lbl, " (tRNA array ready)")
  lbl
}

