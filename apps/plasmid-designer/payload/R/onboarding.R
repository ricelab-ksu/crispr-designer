# ============================================================================
# R/onboarding.R — Getting-started help content (both apps)
# ============================================================================
# Plain-language onboarding for new users: the workflow, the vocabulary, and
# worked sequence examples (concatenation and reverse-complement), so the app
# explains itself without a manual.
#
# Everything here is static text; nothing computes on user data.
# ============================================================================

#' Worked example: concatenating sequence blocks.
#'
#' Teaches the arithmetic of building a construct from named pieces.
#'
#' @return Character(1) explanation with literal examples.
#' @export
pjd_help_concat_example <- function() {
  paste0(
    "Building a construct means concatenating sequence blocks end to end.\n",
    "  ATGC + CGGG = ATGCCGGG        (8 bp)\n",
    "  ATGC + CGGG + TTAA = ATGCCGGGTTAA   (12 bp)\n",
    "Nothing is inserted between blocks, and no bases are dropped."
  )
}

#' Worked example: reverse-complementing an insert (the +/- orientation).
#'
#' @return Character(1) explanation with literal examples.
#' @export
pjd_help_orientation_example <- function() {
  paste0(
    "The +/- control sets which strand of the insert reads 5'->3' in the\n",
    "finished plasmid.\n",
    "  + : insert as supplied        ATGC   ->  ATGC\n",
    "  - : reverse complement        ATGC   ->  GCAT\n",
    "Reverse complement = reverse the order AND swap A<->T, C<->G:\n",
    "  ATGC  ->  (complement) TACG  ->  (reverse) GCAT\n",
    "Use - when the element must point the other way (e.g. an inverted\n",
    "cassette or a promoter that must face the opposite direction)."
  )
}

#' Worked example: the insertion window arithmetic.
#'
#' @return Character(1) explanation with literal examples.
#' @export
pjd_help_window_example <- function() {
  paste0(
    "An insertion replaces a window [start..end] of the backbone:\n",
    "  final length = backbone length - window length + insert length\n",
    "Example: 7,986 bp backbone, replace a 158 bp window with a 1,861 bp\n",
    "marker -> 7,986 - 158 + 1,861 = 9,689 bp.\n",
    "To INSERT without deleting anything, set end = start - 1 (a pure\n",
    "insertion at that position)."
  )
}

#' Worked example: homology arms and the attL donor fragments.
#'
#' @return Character(1) explanation.
#' @export
pjd_help_arm_example <- function() {
  paste0(
    "Homology arms are GENOMIC sequence either side of the Cas9 cut site;\n",
    "they direct where the payload integrates (rules R1/R4). Validated\n",
    "window: 250-1000 bp per arm (rule R2).\n",
    "The two synthesized fragments are:\n",
    "  Fragment 1 = attL3 + 5' arm + tRNA-spacer-tRNA + attL4\n",
    "  Fragment 2 = attL5 + 3' arm + attL6\n",
    "LR Clonase recombines attL x attR -> attB, so the finished plasmid\n",
    "carries attB where the attL/attR sites met."
  )
}

#' Vocabulary glossary for new users.
#'
#' @return data.frame(term, meaning).
#' @export
pjd_help_glossary <- function() {
  data.frame(
    term = c("Backbone / destination vector",
             "Slot",
             "Spacer",
             "PAM",
             "Homology arm",
             "Payload",
             "attR / attL",
             "Marker"),
    meaning = c(
      "The plasmid you clone into; supplies the origin of replication, resistance genes and (for pJAT) the attR sites.",
      "A named region of a backbone that can receive an insert (gRNA array, donor payload, marker swap, restriction site).",
      "The 20-nt guide sequence that directs Cas9; always excludes the NGG PAM.",
      "The NGG trinucleotide immediately 3' of the spacer that Cas9 requires.",
      "Genomic sequence flanking the cut site (250-1000 bp) that directs HDR integration.",
      "The sequence inserted at a slot; a marker cassette is one kind of payload.",
      "Gateway recombination sites: the backbone carries attR, the donor fragments carry attL; LR Clonase converts attL x attR to attB.",
      "A selectable or visible cassette (e.g. 3xP3-DsRed, neo) used to find or select the edit."
    ),
    stringsAsFactors = FALSE
  )
}

#' Render the onboarding card (workflow + worked examples + glossary).
#'
#' @param app "plasmid" or "crispr".
#' @return shiny tagList.
#' @export
pjd_onboarding_card <- function(app = c("plasmid", "crispr")) {
  app <- match.arg(app)
  steps <- if (identical(app, "plasmid")) {
    c("Pick a backbone (pJAT catalog, or import any GenBank/FASTA).",
      "Supply the guide: import a Cas9 Designer bundle, or paste a 20-nt spacer.",
      "Choose where the insert lands (linker window, restriction site, custom).",
      "Set homology arms (250-1000 bp) - imported from the guide bundle when available.",
      "Supply the donor payload (marker library, paste, or FASTA upload).",
      "Click Assemble Vector, then export GenBank / FASTA / JSON / ZIP.")
  } else {
    c("Pick a Drosophila species and prepare its reference genome.",
      "Look up your gene (symbol, FBgn, CG number, or coordinates).",
      "Filter guides by context tier and on-target / specificity scores.",
      "Run the design and inspect the ranked guide table and track plot.",
      "Export the guides for the Plasmid Designer (arms included).")
  }
  gl <- pjd_help_glossary()
  shiny::tags$details(
    class = "pjd-onboarding",
    shiny::tags$summary(style = "cursor:pointer; font-weight:700; color:#0284c7;",
                        "How this app works (click to expand)"),
    shiny::tags$div(style = "padding:10px 4px; font-size:12px; line-height:1.6; color:#334155;",
      shiny::tags$ol(lapply(steps, shiny::tags$li)),
      shiny::tags$h5(style = "margin:10px 0 4px 0;", "Worked examples"),
      shiny::tags$pre(style = "font-size:11px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:8px; white-space:pre-wrap;",
                      pjd_help_concat_example()),
      shiny::tags$pre(style = "font-size:11px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:8px; white-space:pre-wrap;",
                      pjd_help_orientation_example()),
      shiny::tags$pre(style = "font-size:11px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:8px; white-space:pre-wrap;",
                      pjd_help_window_example()),
      if (identical(app, "plasmid")) {
        shiny::tags$pre(style = "font-size:11px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:8px; white-space:pre-wrap;",
                        pjd_help_arm_example())
      } else NULL,
      shiny::tags$h5(style = "margin:10px 0 4px 0;", "Glossary"),
      shiny::tags$table(
        style = "font-size:11.5px; width:100%; border-collapse:collapse;",
        shiny::tags$tbody(lapply(seq_len(nrow(gl)), function(i) {
          shiny::tags$tr(
            shiny::tags$td(style = "padding:3px 8px 3px 0; font-weight:600; vertical-align:top; white-space:nowrap;",
                           gl$term[i]),
            shiny::tags$td(style = "padding:3px 0; color:#475569;", gl$meaning[i]))
        }))
      )
    )
  )
}

