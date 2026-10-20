# ============================================================================
# R/plasmid_assembly_slots.R — Multi-slot cloning engine (backbone_registry/1.0)
# ============================================================================
# Generic cloning into ANY named slot of ANY backbone (not just the pJAT
# attR4..attR5 linker).  Two entry points:
#
#   clone_into_slot(backbone_seq, slot, insert, orient, topology)
#       — insert/replace at a registry slot (grna_array, donor_payload,
#         marker_swap, recombination_site).  Handles slots that span the
#         origin of a circular plasmid (start > end) via wraparound.
#
#   clone_at_restriction_site(backbone_seq, enzyme, site_pos, insert, orient)
#       — cut at a Type IIS site using the real cut offset and insert there,
#         preserving the enzyme recognition site (never auto-selected).
#
# Both reuse insert_cassette() from plasmid_assembly.R so the strand-aware
# reverse-complement and coordinate bookkeeping stay in ONE place.
#
# Requires: plasmid_assembly.R (insert_cassette, CLONING_ENZYMES) sourced first.
# ============================================================================

#' Resolve a slot's [start, end] window, handling circular wraparound.
#'
#' A circular plasmid slot may span the origin, represented as `start > end`.
#'
#' @param slot One-row data.frame from backbone_slot() (needs start, end).
#' @param size_bp Integer backbone length.
#' @param topology "circular" or "linear".
#' @return list(start, end, spans_origin, label).
#' @export
resolve_slot_window <- function(slot, size_bp, topology = c("circular", "linear")) {
  topology <- match.arg(topology)
  stopifnot(is.data.frame(slot), nrow(slot) == 1L,
            all(c("start", "end") %in% names(slot)))
  s <- as.integer(slot$start[1]); e <- as.integer(slot$end[1])
  if (is.na(s) || is.na(e)) {
    stop("slot '", slot$slot_id[1], "' has no coordinates; cannot clone into it.")
  }
  spans <- s > e
  if (spans && !identical(topology, "circular")) {
    stop("slot '", slot$slot_id[1], "' spans the origin (start > end) but the ",
         "backbone is linear; this window is only valid on a circular plasmid.")
  }
  if (s < 1L || e > size_bp) {
    stop("slot '", slot$slot_id[1], "' coordinates (", s, "..", e,
         ") fall outside the backbone (1..", size_bp, ").")
  }
  list(start = s, end = e, spans_origin = spans,
       label = as.character(slot$label[1] %||% slot$slot_id[1]))
}

# Rotate a circular sequence so a wraparound window [start..n]+[1..end] becomes
# a contiguous block at the END of the rotated sequence.  The untouched middle
# [end+1 .. start-1] stays at the front, preserving the full plasmid length.
.paf_rotate_for_window <- function(seq, start, end) {
  n <- nchar(seq)
  middle <- substr(seq, end + 1L, start - 1L)      # untouched backbone
  window <- paste0(substr(seq, start, n), substr(seq, 1L, end))
  paste0(middle, window)
}

#' Clone an insert into a named backbone slot (strand-aware, wraparound-safe).
#'
#' @param backbone_seq Character(1) full backbone sequence.
#' @param slot One-row data.frame from backbone_slot().
#' @param insert Character(1) insert sequence (ACGTN).
#' @param orient "+" (as given) or "-" (reverse-complemented before insertion).
#' @param topology "circular" or "linear".
#' @param replace TRUE to replace the slot window; FALSE to insert before it.
#' @return list: `sequence`, `insert_start`, `insert_end`, `insert_len`,
#'   `slot_id`, `slot_type`, `spans_origin`, `strand`, `window_start`,
#'   `window_end`, `final_len`, `junction_left`, `junction_right`.
#' @export
clone_into_slot <- function(backbone_seq, slot, insert, orient = "+",
                            topology = c("circular", "linear"),
                            replace = TRUE) {
  topology <- match.arg(topology)
  orient <- match.arg(orient, c("+", "-"))
  stopifnot(is.character(backbone_seq), length(backbone_seq) == 1L,
            is.character(insert), length(insert) == 1L, nzchar(insert))
  backbone_seq <- toupper(gsub("[^ACGTN]", "", backbone_seq))
  insert <- toupper(gsub("[^ACGTN]", "", insert))
  n <- nchar(backbone_seq)

  win <- resolve_slot_window(slot, n, topology)

  if (isTRUE(win$spans_origin)) {
    # `rotated` = [untouched middle] ++ [window], so the window sits at the END.
    # Clone there, then rotate back to the original origin.
    rotated <- .paf_rotate_for_window(backbone_seq, win$start, win$end)
    span_len <- (n - win$start + 1L) + win$end      # window length
    middle_len <- n - span_len
    rot_start <- middle_len + 1L
    rot_end <- if (replace) n else middle_len
    asm <- insert_cassette(rotated, rot_start, rot_end, insert, strand = orient)

    # Rotate back: the last `win$end` bases of the ORIGINAL plasmid must return
    # to the front.  In `asm$sequence` they are the final win$end characters.
    back_len <- nchar(asm$sequence)
    tail_seg <- substr(asm$sequence, back_len - win$end + 1L, back_len)
    head_seg <- substr(asm$sequence, 1L, back_len - win$end)
    final <- paste0(tail_seg, head_seg)
    new_start <- if (win$start + asm$insert_len > n) {
      (win$start + asm$insert_len) - n
    } else {
      win$start
    }
    return(list(
      sequence = final,
      insert_start = as.integer(new_start),
      insert_end = as.integer(new_start + asm$insert_len - 1L),
      insert_len = asm$insert_len,
      slot_id = as.character(slot$slot_id[1]),
      slot_type = as.character(slot$slot_type[1]),
      spans_origin = TRUE,
      strand = orient,
      window_start = win$start, window_end = win$end,
      final_len = nchar(final),
      junction_left = asm$junction_left,
      junction_right = asm$junction_right
    ))
  }

  end <- if (replace) win$end else win$start - 1L
  asm <- insert_cassette(backbone_seq, win$start, end, insert, strand = orient)
  c(asm, list(
    slot_id = as.character(slot$slot_id[1]),
    slot_type = as.character(slot$slot_type[1]),
    spans_origin = FALSE,
    window_start = win$start, window_end = win$end
  ))
}


#' Clone an insert at a Type IIS restriction site (real cut offset).
#'
#' The recognition site is preserved; the insert is placed at the enzyme's cut
#' position.  Never auto-selects a site — the caller must pass the position.
#'
#' @param backbone_seq Character(1) full backbone sequence.
#' @param enzyme Enzyme name in CLONING_ENZYMES.
#' @param site_pos 1-based position of the recognition site.
#' @param insert Character(1) insert sequence.
#' @param orient "+" or "-".
#' @return list as clone_into_slot() plus `enzyme`, `cut_pos`, `overhang`.
#' @export
clone_at_restriction_site <- function(backbone_seq, enzyme, site_pos, insert,
                                      orient = "+") {
  orient <- match.arg(orient, c("+", "-"))
  if (!enzyme %in% names(CLONING_ENZYMES)) {
    stop("enzyme must be one of: ", paste(names(CLONING_ENZYMES), collapse = ", "))
  }
  backbone_seq <- toupper(gsub("[^ACGTN]", "", backbone_seq))
  n <- nchar(backbone_seq)
  e <- CLONING_ENZYMES[[enzyme]]
  site_pos <- as.integer(site_pos)
  if (is.na(site_pos) || site_pos < 1L || site_pos + nchar(e$site) - 1L > n) {
    stop("site_pos ", site_pos, " is out of range for ", enzyme,
         " (site is ", nchar(e$site), " bp).")
  }
  found <- substr(backbone_seq, site_pos, site_pos + nchar(e$site) - 1L)
  if (!identical(found, e$site)) {
    stop("no ", enzyme, " recognition site (", e$site, ") at position ",
         site_pos, "; found '", found, "'.")
  }
  cut_pos <- site_pos + nchar(e$site) + as.integer(e$cut_after) - 1L
  if (cut_pos > n) stop("cut position ", cut_pos, " falls outside the backbone.")

  # Insert-only at the cut: replace nothing (end = start - 1).
  asm <- insert_cassette(backbone_seq, cut_pos, cut_pos - 1L, insert,
                         strand = orient)
  overhang <- if (cut_pos + 4L <= n) substr(backbone_seq, cut_pos, cut_pos + 4L) else ""
  c(asm, list(
    slot_id = paste0("restriction_", enzyme, "_", site_pos),
    slot_type = "restriction_site",
    enzyme = enzyme, site_pos = site_pos, cut_pos = cut_pos,
    overhang = overhang,
    window_start = cut_pos, window_end = cut_pos - 1L,
    spans_origin = FALSE
  ))
}

