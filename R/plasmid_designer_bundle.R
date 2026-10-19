# ============================================================================
# R/plasmid_designer_bundle.R — cas9_guide_bundle/1.0 producer (Plasmid Designer interop)
# ============================================================================
# Emits the interchange bundle consumed by Plasmid_Designer's read_guide_bundle():
#   Cas9_<gene>_<species>_guides_for_plasmid_designer.csv  (+ .json sidecar)
# One row per exported guide. Arms are per selected guide, so only the guide
# that has arm design (arm_design$guide_id) is emitted — never NA-padded and
# never duplicated across rows. Sequence is never fabricated.
# ============================================================================

PLASMID_DESIGNER_BUNDLE_VERSION <- "1.0.0"
GB_SCHEMA     <- "cas9_guide_bundle/1.0"
GB_ARM_POLICY <- "250-1000 bp validated homology arms"

# Columns required by read_guide_bundle(), in contract order.
.gbu_cols <- c(
  "guide_id", "gene_symbol", "species", "assembly", "chr", "accession",
  "strand", "protospacer_start", "protospacer_end", "pam_start", "pam_end",
  "cut_site", "spacer_seq", "pam_seq", "on_target_score", "cfd_specificity",
  "arm_5p_seq", "arm_5p_start", "arm_5p_end", "arm_5p_len",
  "arm_3p_seq", "arm_3p_start", "arm_3p_end", "arm_3p_len",
  "arm_3p_anchor", "pam_shielding_summary", "guide_bundle_schema"
)

.gbu_clean_dna <- function(x) {
  x <- as.character(x)
  x <- toupper(gsub("[^ACGTN]", "", x))
  x[is.na(x)] <- ""
  x
}

.gbu_first <- function(...) {
  for (v in list(...)) {
    if (!is.null(v) && length(v) && !is.na(v[1])) return(as.character(v[1]))
  }
  ""
}

#' Build a cas9_guide_bundle/1.0 from the candidates data.frame + arm design.
#' @return list(csv_text, sidecar (list), guides (data.frame), n_guides)
#' @export
build_plasmid_designer_bundle <- function(grna_df, arm_design, metadata = NULL,
                                          gene_symbol = NULL, app_version = "1.0.0") {
  if (!is.data.frame(grna_df) || nrow(grna_df) == 0L) {
    stop("grna_df must be a non-empty candidates data.frame.")
  }
  if (is.null(arm_design) || !is.list(arm_design)) {
    stop("arm_design is required: design homology arms for the selected guide first.")
  }
  gid <- .gbu_first(arm_design$guide_id)
  if (is.null(gid)) stop("arm_design has no guide_id; cannot select the exported guide.")

  row <- grna_df[grna_df$guide_id == gid, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("selected guide_id '", gid, "' not found (or not unique) in grna_df.")
  }

  spacer <- .gbu_clean_dna(.gbu_first(arm_design$spacer_seq, row$spacer_seq))
  if (!nzchar(spacer)) stop("selected guide '", gid, "' has no spacer_seq.")
  if (nchar(spacer) != 20L) {
    stop("spacer_seq for guide '", gid, "' must be exactly 20 nt (PAM excluded); got ",
         nchar(spacer), " nt.")
  }

  arm5 <- .gbu_clean_dna(.gbu_first(arm_design$arm_5p_seq))
  arm3 <- .gbu_clean_dna(.gbu_first(arm_design$arm_3p_seq))
  if (!nzchar(arm5) || !nzchar(arm3)) {
    stop("selected guide '", gid, "' has empty arm_5p_seq/arm_3p_seq; design arms first.")
  }
  a5_len <- nchar(arm5); a3_len <- nchar(arm3)
  if (a5_len < 250L || a3_len < 250L) {
    stop("selected guide '", gid, "' has arms < 250 bp (5'=", a5_len, ", 3'=", a3_len,
         "); increase arm length to >= 250 bp.")
  }

  # ---- identity / provenance ----------------------------------------------
  species <- .gbu_first(arm_design$species, row$species,
                        if (!is.null(metadata)) metadata$species else NULL,
                        "melanogaster")
  gene    <- .gbu_first(gene_symbol,
                        if (!is.null(metadata)) metadata$gene_symbol else NULL,
                        if (!is.null(metadata)) metadata$query      else NULL,
                        sub("_g[0-9]+$", "", gid))
  assembly <- .gbu_first(row$assembly,
                         if (!is.null(metadata)) metadata$assembly           else NULL,
                         if (!is.null(metadata)) metadata$assembly_accession else NULL,
                         "unknown")
  cfd  <- if ("cfd_specificity_score" %in% names(row)) row$cfd_specificity_score else NA_real_
  on_t <- if ("on_target_score" %in% names(row)) row$on_target_score else NA_real_

  df <- data.frame(
    guide_id            = gid,
    gene_symbol         = gene,
    species             = species,
    assembly            = assembly,
    chr                 = .gbu_first(row$chr, arm_design$chr, "unknown"),
    accession           = .gbu_first(row$accession, arm_design$accession, "unknown"),
    strand              = .gbu_first(row$strand, arm_design$strand, "+"),
    protospacer_start   = as.integer(.gbu_first(row$protospacer_start, arm_design$protospacer_start, 0)),
    protospacer_end     = as.integer(.gbu_first(row$protospacer_end,   arm_design$protospacer_end,   0)),
    pam_start           = as.integer(.gbu_first(row$pam_start, arm_design$pam_start, 0)),
    pam_end             = as.integer(.gbu_first(row$pam_end,   arm_design$pam_end,   0)),
    cut_site            = as.integer(.gbu_first(row$cut_site,  arm_design$cut_site,  0)),
    spacer_seq          = spacer,
    pam_seq             = .gbu_first(row$pam_seq, arm_design$pam_seq, "NGG"),
    on_target_score     = as.numeric(on_t),
    cfd_specificity     = as.numeric(cfd),
    arm_5p_seq          = arm5,
    arm_5p_start        = as.integer(.gbu_first(arm_design$arm_5p_start, 0)),
    arm_5p_end          = as.integer(.gbu_first(arm_design$arm_5p_end,   0)),
    arm_5p_len          = a5_len,
    arm_3p_seq          = arm3,
    arm_3p_start        = as.integer(.gbu_first(arm_design$arm_3p_start, 0)),
    arm_3p_end          = as.integer(.gbu_first(arm_design$arm_3p_end,   0)),
    arm_3p_len          = a3_len,
    arm_3p_anchor       = .gbu_first(arm_design$arm_3p_anchor, "cut_site"),
    pam_shielding_summary = .gbu_first(arm_design$pam_shielding_summary, "none"),
    guide_bundle_schema = GB_SCHEMA,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  df <- df[, .gbu_cols, drop = FALSE]

  sidecar <- list(
    schema      = GB_SCHEMA,
    app_version = app_version,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
    species     = species,
    assembly    = assembly,
    gene_symbol = gene,
    n_guides    = nrow(df),
    arm_policy  = GB_ARM_POLICY
  )

  tc <- textConnection("csv_buf", open = "w", local = TRUE)
  utils::write.csv(df, file = tc, row.names = FALSE)
  close(tc)

  list(csv_text = paste(csv_buf, collapse = "\n"),
       sidecar = sidecar, guides = df, n_guides = nrow(df))
}

#' Write the bundle CSV + JSON sidecar, return their paths.
#' @export
write_plasmid_designer_bundle <- function(grna_df, arm_design, metadata = NULL,
                                          gene_symbol = NULL, app_version = "1.0.0",
                                          file = NULL) {
  b <- build_plasmid_designer_bundle(grna_df, arm_design, metadata,
                                     gene_symbol, app_version)
  if (is.null(file) || !nzchar(file[1])) {
    species_tok <- gsub("[^A-Za-z0-9_\\-]", "_", b$sidecar$species)
    gene_tok    <- gsub("[^A-Za-z0-9_\\-]", "_", b$sidecar$gene_symbol)
    file <- paste0("Cas9_", gene_tok, "_", species_tok, "_guides_for_plasmid_designer.csv")
  }
  json_file <- sub("\\.csv$", ".json", file, ignore.case = TRUE)
  if (identical(json_file, file)) json_file <- paste0(file, ".json")
  writeLines(b$csv_text, con = file)
  writeLines(jsonlite::toJSON(b$sidecar, pretty = TRUE, auto_unbox = TRUE),
             con = json_file)
  invisible(list(csv = file, json = json_file,
                 csv_bytes = file.info(file)$size,
                 json_bytes = file.info(json_file)$size,
                 n_guides = b$n_guides))
}

