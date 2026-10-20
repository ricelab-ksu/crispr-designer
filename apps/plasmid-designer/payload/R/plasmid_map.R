# ============================================================================
# R/plasmid_map.R — Interactive circular/linear plasmid map (no dependencies)
# ============================================================================
# Renders a Benchling/SnapGene-style plasmid map as inline SVG, driven entirely
# by the feature table from read_pjat_genbank() or pjd_assembly(). Everything is
# local: no CDN, no JS framework, and the map still works with JavaScript
# disabled (native SVG <title> tooltips).
#
# Architecture approved in project_documentation/attl_and_map_spec.md (Part B):
#   * ~4 KB vanilla JS (www/pjd_assets/pjd_map.js) for hover / zoom / pan / click
#   * native SVG fallback for tooltips
#   * the static Addgene PNG remains available as a secondary reference view
# ============================================================================

PJD_MAP_VERSION <- "1.0.0"

PJD_MAP_ASSET_DIR <- "pjd_assets"

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a) && length(a) == 1L && !is.na(a)) a else b
}

# Feature-class palette (stable across backbones and assemblies).
PJD_MAP_COLORS <- c(
  attR       = "#b45309",   # Gateway recombination sites
  promoter   = "#0284c7",   # U6 / Pol III / viral promoters
  cds        = "#059669",   # coding sequences / markers
  resistance = "#dc2626",   # antibiotic resistance
  negative   = "#64748b",   # ccdB / SacB negative selection
  origin     = "#7c3aed",   # ori / rep_origin
  primer     = "#94a3b8",
  insert     = "#db2777",   # user payload
  other      = "#2563eb"
)

#' Map a GenBank feature to a colour class.
#' @export
pjd_map_feature_class <- function(key, label, product = "") {
  txt <- tolower(paste(key, label, product))
  if (grepl("attr[3-6]|^attr", txt))                     return("attR")
  if (grepl("ccdb|sacb", txt))                           return("negative")
  if (grepl("ampr|cmr|smr|kanr|hygr|puro|bsd|aad|aph|resistance|blasticidin", txt))
                                                         return("resistance")
  if (grepl("promoter|^u6|du6|ie1|3xp3|act5c|hsp70|terminator", txt))
                                                         return("promoter")
  if (grepl("^cds|egfp|dsred|ecfp|marker|transposase|integrase", txt))
                                                         return("cds")
  if (grepl("origin|^ori", txt))                         return("origin")
  if (grepl("primer", txt))                              return("primer")
  if (grepl("insert|payload", txt))                       return("insert")
  "other"
}

#' Normalise a parsed GenBank record into a plottable map record.
#'
#' @param gb list(sequence, size_bp, features) from read_pjat_genbank().
#' @param title Optional map title.
#' @return list(title, size_bp, sequence, features) where features has
#'   type, start, end, strand, label, class, color, row.
#' @export
pjd_map_record_from_genbank <- function(gb, title = NULL) {
  if (is.null(gb) || is.null(gb$sequence) || !nzchar(gb$sequence)) return(NULL)
  size <- nchar(gb$sequence)

  feats <- gb$features
  if (is.null(feats) || nrow(feats) == 0L) {
    feats <- data.frame(key = character(0), start = integer(0),
                        end = integer(0), strand = character(0),
                        label = character(0), product = character(0),
                        stringsAsFactors = FALSE)
  }
  if (!"product" %in% names(feats)) feats$product <- ""
  feats <- feats[feats$key != "source", , drop = FALSE]

  cls <- mapply(pjd_map_feature_class, feats$key, feats$label, feats$product,
                USE.NAMES = FALSE)

  labels <- as.character(feats$label)
  empty <- is.na(labels) | !nzchar(labels)
  labels[empty] <- as.character(feats$key[empty])

  out <- data.frame(
    type   = as.character(feats$key),
    start  = as.integer(feats$start),
    end    = as.integer(feats$end),
    strand = ifelse(is.na(feats$strand) | !nzchar(feats$strand), "+",
                    as.character(feats$strand)),
    label  = labels,
    class  = cls,
    color  = unname(PJD_MAP_COLORS[cls]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$start) & !is.na(out$end), , drop = FALSE]

  # Circular features crossing the origin are split so no arc is drawn the
  # long way round.
  split <- vector("list", nrow(out))
  k <- 0L
  for (i in seq_len(nrow(out))) {
    r <- out[i, , drop = FALSE]
    k <- k + 1L
    if (r$end < r$start) {
      a <- r; a$end <- size
      b <- r; b$start <- 1L
      split[[k]] <- rbind(a, b)
    } else {
      split[[k]] <- r
    }
  }
  out <- if (k > 0L) do.call(rbind, split) else out
  out <- out[out$start >= 1L & out$end <= size & out$end >= out$start, , drop = FALSE]
  out <- out[order(out$start), , drop = FALSE]
  out$row <- if (nrow(out)) (seq_len(nrow(out)) - 1L) %% 3L else integer(0)
  rownames(out) <- NULL

  list(title = title %||% (gb$name %||% "plasmid"),
       size_bp = size, sequence = gb$sequence, features = out)
}

#' Build a map record from an assembled vector (backbone + user insert).
#'
#' Includes backbone features (promoters, CDSs, ori, attR sites, etc.) shifted
#' around the insertion window, plus the insert features with a descriptive
#' label (e.g. "Neomycin/G418" instead of the generic "user insert").
#'
#' @param asm Assembly object from assemble_expression_vector().
#' @param backbone Optional backbone name for the map title.
#' @param backbone_gb Optional parsed GenBank for the backbone (to include its
#'   features on the assembled map). Features overlapping the insertion window
#'   are replaced by the insert; features after the window are shifted by the
#'   net length change.
#' @export
pjd_map_record_from_assembly <- function(asm, backbone = NULL, backbone_gb = NULL) {
  if (is.null(asm) || is.null(asm$sequence)) return(NULL)

  size <- nchar(asm$sequence)
  win_start <- asm$insert_start %||% 1L
  win_end   <- asm$insert_end   %||% size
  net_shift <- asm$insert_len - (win_end - win_start + 1L)

  # ---- Insert features -------------------------------------------------------
  feats <- asm$features
  if (is.null(feats) || nrow(feats) == 0L) {
    feats <- data.frame(
      key = "misc_feature",
      start = win_start,
      end = win_end,
      strand = asm$orient %||% "+",
      label = asm$insert_label %||% "user insert",
      product = "",
      stringsAsFactors = FALSE
    )
  } else {
    if (!"product" %in% names(feats)) feats$product <- ""
    if (!"key" %in% names(feats)) feats$key <- feats$type
    if (!"label" %in% names(feats)) feats$label <- feats$key
  }

  # ---- Backbone features (shifted around insertion window) ------------------
  if (!is.null(backbone_gb) && !is.null(backbone_gb$features) &&
      nrow(backbone_gb$features) > 0L) {
    b_feats <- backbone_gb$features
    b_feats <- b_feats[b_feats$key != "source", , drop = FALSE]
    if (nrow(b_feats) > 0L) {
      if (!"product" %in% names(b_feats)) b_feats$product <- ""

      # Keep features entirely before or entirely after the insertion window;
      # features overlapping the window are replaced by the insert.
      keep_before <- b_feats$end < win_start
      keep_after  <- b_feats$start > win_end
      b_keep <- b_feats[keep_before | keep_after, , drop = FALSE]

      if (nrow(b_keep) > 0L) {
        # Shift features that were after the window by the net length change.
        after_idx <- b_keep$start > win_end
        b_keep$start[after_idx] <- b_keep$start[after_idx] + net_shift
        b_keep$end[after_idx]   <- b_keep$end[after_idx]   + net_shift
        # Clamp to valid range and drop any that became invalid.
        b_keep$start <- pmax(1L, b_keep$start)
        b_keep$end   <- pmin(size, b_keep$end)
        b_keep <- b_keep[b_keep$end >= b_keep$start, , drop = FALSE]
      }

      if (nrow(b_keep) > 0L) {
        common <- intersect(names(feats), names(b_keep))
        feats <- rbind(feats[, common, drop = FALSE],
                       b_keep[, common, drop = FALSE])
      }
    }
  }

  gb <- list(sequence = asm$sequence,
             size_bp = size,
             features = feats,
             name = if (!is.null(backbone)) paste(backbone, "assembly") else "assembly")
  rec <- pjd_map_record_from_genbank(gb)
  if (!is.null(rec) && !is.null(asm$insert_start) && !is.null(asm$insert_end)) {
    rec$insert <- c(start = asm$insert_start, end = asm$insert_end)
  }
  rec
}


# ---- SVG geometry ----------------------------------------------------------

.pjd_map_point <- function(pos, size, radius, cx, cy) {
  a <- 2 * pi * (pos - 1) / size - pi / 2
  c(x = cx + radius * cos(a), y = cy + radius * sin(a))
}

.pjd_map_arc_path <- function(start, end, size, radius, cx, cy) {
  p0 <- .pjd_map_point(start, size, radius, cx, cy)
  p1 <- .pjd_map_point(end, size, radius, cx, cy)
  span <- 2 * pi * (end - start + 1) / size
  large <- if (span > pi) 1L else 0L
  sprintf("M %.2f %.2f A %.2f %.2f 0 %d 1 %.2f %.2f",
          p0[["x"]], p0[["y"]], radius, radius, large, p1[["x"]], p1[["y"]])
}

.pjd_map_legend <- function(feats) {
  present <- unique(feats$class)
  present <- present[present %in% names(PJD_MAP_COLORS)]
  if (!length(present)) return(NULL)
  shiny::div(
    class = "pjd-map-legend",
    lapply(present, function(cl) {
      shiny::span(class = "pjd-map-legend-item",
                  shiny::span(class = "pjd-map-swatch",
                              style = paste0("background:", PJD_MAP_COLORS[[cl]], ";")),
                  cl)
    })
  )
}

.pjd_map_tooltip_attrs <- function(f) {
  list(
    `data-pjd-map` = "feature",
    `data-label` = f$label,
    `data-type` = f$type,
    `data-start` = f$start,
    `data-end` = f$end,
    `data-strand` = f$strand,
    `data-length` = f$end - f$start + 1L,
    `data-color` = f$color
  )
}

.pjd_map_circular_svg <- function(record, height) {
  size <- record$size_bp
  cx <- 230; cy <- 230
  r_backbone <- 176
  feats <- record$features

  ticks <- unique(round(seq(0, size, length.out = 9)))
  # 0 and size are the same point on a circle: keep one label.
  if (length(ticks) > 1L && ticks[length(ticks)] == size) ticks <- head(ticks, -1L)
  tick_nodes <- lapply(ticks, function(pos) {
    p  <- .pjd_map_point(max(1, pos + 1), size, r_backbone, cx, cy)
    p2 <- .pjd_map_point(max(1, pos + 1), size, r_backbone + 10, cx, cy)
    lab <- .pjd_map_point(max(1, pos + 1), size, r_backbone + 26, cx, cy)
    shiny::tagList(
      shiny::tags$line(x1 = p[["x"]], y1 = p[["y"]], x2 = p2[["x"]], y2 = p2[["y"]],
                       stroke = "#cbd5e1", `stroke-width` = 1),
      shiny::tags$text(x = lab[["x"]], y = lab[["y"]], `text-anchor` = "middle",
                       `dominant-baseline` = "middle", class = "pjd-map-tick",
                       paste0(round(pos / 1000, 1), " kb"))
    )
  })

  feature_nodes <- lapply(seq_len(nrow(feats)), function(i) {
    f <- feats[i, , drop = FALSE]
    radius <- 118 + (f$row[[1L]] * 16)
    path <- .pjd_map_arc_path(f$start, f$end, size, radius, cx, cy)
    do.call(shiny::tags$path, c(
      list(d = path, fill = "none", stroke = f$color,
           `stroke-width` = 7, `stroke-linecap` = "round",
           `stroke-opacity` = 0.9),
      .pjd_map_tooltip_attrs(f),
      list(shiny::tags$title(sprintf("%s [%s] %d..%d (%s, %d bp)",
                                     f$label, f$type, f$start, f$end,
                                     f$strand, f$end - f$start + 1L)))
    ))
  })

  label_nodes <- lapply(seq_len(nrow(feats)), function(i) {
    f <- feats[i, , drop = FALSE]
    if ((f$end - f$start + 1L) < size / 60) return(NULL)
    radius <- 118 + (f$row[[1L]] * 16)
    p <- .pjd_map_point((f$start + f$end) / 2, size, radius + 12, cx, cy)
    shiny::tags$text(x = p[["x"]], y = p[["y"]], class = "pjd-map-label",
                     `text-anchor` = "middle", `dominant-baseline` = "middle",
                     f$label)
  })

  shiny::tags$svg(
    xmlns = "http://www.w3.org/2000/svg", viewBox = "0 0 460 460",
    class = "pjd-map-svg", style = paste0("height:", height, ";"),
    `data-pjd-map` = "circular", `data-size` = size,
    shiny::tags$circle(cx = cx, cy = cy, r = r_backbone, fill = "#ffffff",
                       stroke = "#94a3b8", `stroke-width` = 3),
    shiny::tags$circle(cx = cx, cy = cy, r = 96, fill = "none",
                       stroke = "#e2e8f0", `stroke-width` = 1),
    tick_nodes, feature_nodes, label_nodes,
    shiny::tags$text(x = cx, y = cy - 6, `text-anchor` = "middle",
                     class = "pjd-map-title", record$title),
    shiny::tags$text(x = cx, y = cy + 14, `text-anchor` = "middle",
                     class = "pjd-map-subtitle",
                     paste0(format(size, big.mark = ","), " bp"))
  )
}


.pjd_map_linear_svg <- function(record, height) {
  size <- record$size_bp
  feats <- record$features
  x0 <- 40; x1 <- 860
  pos_x <- function(p) x0 + (p - 1) / max(1, size - 1) * (x1 - x0)

  ticks <- unique(round(seq(0, size, length.out = 11)))
  tick_nodes <- lapply(ticks, function(pos) {
    x <- pos_x(max(1, pos + 1))
    shiny::tagList(
      shiny::tags$line(x1 = x, y1 = 100, x2 = x, y2 = 108, stroke = "#cbd5e1"),
      shiny::tags$text(x = x, y = 122, `text-anchor` = "middle",
                       class = "pjd-map-tick",
                       paste0(round(pos / 1000, 1), " kb"))
    )
  })

  feature_nodes <- lapply(seq_len(nrow(feats)), function(i) {
    f <- feats[i, , drop = FALSE]
    xa <- pos_x(f$start); xb <- pos_x(f$end)
    w <- max(2, xb - xa)
    y <- if (identical(f$strand, "-")) 132 else 76
    do.call(shiny::tags$rect, c(
      list(x = xa, y = y, width = w, height = 18, rx = 3,
           fill = f$color, `fill-opacity` = 0.85),
      .pjd_map_tooltip_attrs(f),
      list(shiny::tags$title(sprintf("%s [%s] %d..%d (%s, %d bp)",
                                     f$label, f$type, f$start, f$end,
                                     f$strand, f$end - f$start + 1L)))
    ))
  })

  shiny::tags$svg(
    xmlns = "http://www.w3.org/2000/svg", viewBox = "0 0 900 200",
    class = "pjd-map-svg", style = paste0("height:", height, ";"),
    `data-pjd-map` = "linear", `data-size` = size,
    shiny::tags$line(x1 = x0, y1 = 104, x2 = x1, y2 = 104,
                     stroke = "#94a3b8", `stroke-width` = 2),
    tick_nodes, feature_nodes,
    shiny::tags$text(x = x0, y = 52, class = "pjd-map-label", "5' (+) strand"),
    shiny::tags$text(x = x0, y = 172, class = "pjd-map-label", "3' (-) strand")
  )
}

#' Interactive plasmid map UI (circular or linear).
#'
#' @param record Output of pjd_map_record_from_genbank()/from_assembly().
#' @param mode "circular" or "linear".
#' @param height CSS height for the SVG.
#' @return shiny tagList.
#' @export
.pjd_interactive_map_ui <- function(record, mode = c("circular", "linear"),
                                    height = "520px") {
  mode <- match.arg(mode)
  if (is.null(record) || is.null(record$size_bp) || is.na(record$size_bp) ||
      record$size_bp <= 0) {
    return(shiny::div(class = "cas9-card pjd-map-empty",
                      "No map data for this record."))
  }
  svg <- if (identical(mode, "circular")) {
    .pjd_map_circular_svg(record, height)
  } else {
    .pjd_map_linear_svg(record, height)
  }

  shiny::tagList(
    shiny::singleton(shiny::tags$head(
      shiny::tags$script(src = file.path(PJD_MAP_ASSET_DIR, "pjd_map.js")),
      shiny::tags$link(rel = "stylesheet",
                       href = file.path(PJD_MAP_ASSET_DIR, "pjd_map.css"))
    )),
    shiny::div(class = "pjd-interactive-map",
               `data-pjd-map-mode` = mode,
               shiny::div(class = "pjd-map-toolbar",
                          shiny::tags$button(type = "button", class = "pjd-map-btn",
                                             `data-pjd-action` = "zoom-in", "+"),
                          shiny::tags$button(type = "button", class = "pjd-map-btn",
                                             `data-pjd-action` = "zoom-out", "\u2212"),
                          shiny::tags$button(type = "button", class = "pjd-map-btn",
                                             `data-pjd-action` = "reset", "Reset")),
               .pjd_map_legend(record$features),
               svg,
               shiny::div(class = "pjd-map-tooltip", hidden = "hidden"))
  )
}

