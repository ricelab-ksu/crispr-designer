#!/usr/bin/env Rscript
# ============================================================================
# tools/add_pjat_descriptions.R — add curated 2-4 sentence descriptions
# ============================================================================
# Adds a `description` field to every plasmid in donor_templates/pjat_plasmids.yml.
# Sources: the existing one-line `purpose`, plus
# project_documentation/atlanta_plasmid_register.md (Stern, Kim & Behrman 2024).
#
# ADDITIVE ONLY: attR coordinates, sizes and every existing key are preserved.
# Re-running is idempotent (an existing description is replaced).
#
# Usage: Rscript tools/add_pjat_descriptions.R
# ============================================================================

yml <- "/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer/donor_templates/pjat_plasmids.yml"
if (!file.exists(yml)) stop("catalog not found: ", yml)

# Curated descriptions, keyed by plasmid name.  Each states what the plasmid is
# for, its markers, and when to pick it.
DESCRIPTIONS <- list(
  pJAT39 = "Dual-U6 destination vector that marks an HDR integrant with EGFP in flight muscles. Carries two U6-3 promoters so a single plasmid can express two gRNAs, plus the standard Amp/ccdB and Chlor/sacB Gateway selection arms. Choose it when you need two guides and a muscle-visible fluorescent integrant marker.",
  pJAT11 = "attP landing-site destination vector carrying a removable ie1-DsRed marker expressed in the abdomen and mouthparts. Includes one U6-3 promoter for a single gRNA. Use it to establish a phiC31 landing site whose presence is reported by red fluorescence in the abdomen.",
  pJAT31 = "Dual-purpose destination vector that inserts an attP landing site while carrying a removable ie1-DsRed marker, with two U6-3 promoters for two gRNAs. Pick it when you want both an attP site and a two-guide HDR experiment in one plasmid.",
  pJAT7 = "attP landing-site destination vector with a removable 3xP3-DsRed eye marker and one U6-3 promoter. The classic choice for establishing a phiC31 landing site scored by red eyes.",
  pJAT10 = "Single-U6 destination vector that marks an HDR integrant with DsRed in flight muscles. Uses the standard Gateway selection arms and is the single-guide counterpart to pJAT30.",
  pJAT30 = "Dual-U6 destination vector that marks an HDR integrant with DsRed in flight muscles. The second U6-3 promoter allows two gRNAs from one plasmid. Choose it for a two-guide edit with a muscle-visible marker.",
  pJAT77 = "Empty destination backbone with a multiple cloning site and no U6 promoter. Supply your own promoter if you need gRNA expression, or use the MCS and flanking restriction sites to add a U6 cassette. The most flexible starting point when the standard U6 arrangements do not fit.",
  pJAT74 = "Empty destination backbone with a multiple cloning site and a single U6-3 promoter. The default empty backbone for a one-guide experiment where you want to supply the payload yourself.",
  pJAT32 = "Empty destination backbone with a multiple cloning site and two U6-3 promoters. Use it when you need two gRNAs but want to supply the rest of the payload.",
  pJAT17 = "Germline Cas9 source plasmid driven by the nanos promoter, carrying a 3xP3-EYFP eye marker and one U6-3 promoter. Pick it to provide heritable Cas9 activity from the germline.",
  pJAT43 = "Heat-shock-inducible B3 recombinase plasmid carrying an ie1-EGFP marker. Use it to remove marker cassettes flanked by B3 recombination sites after an edit is established.",
  pJAT56 = "phiC31 integrase source driven by the vasa promoter with a vasa 3'UTR, carrying a 3xP3-EGFP eye marker. Choose it for germline phiC31 integrase activity scored by green eyes.",
  pJAT62 = "phiC31 integrase source driven by the vasa promoter with a vasa 3'UTR, carrying an ie1-EGFP marker. The abdomen/mouthparts-marked alternative to pJAT56 for germline phiC31 activity.",
  pJAT55 = "piggyBac transposase source (alpha-tubulin promoter, pBac-Trsp) carrying a 3xP3-ECFP marker. Discontinued at Addgene; the archival stock is DGRC 1632. Use it for piggyBac-mediated transgenesis.",
  pJAT57 = "UAS-Cas9 effector plasmid with a DsRed marker in flight muscles and one U6-3 promoter. Use it to express Cas9 under Gal4/UAS control for spatially restricted editing.",
  pJAT58 = "UAS-driven jGCaMP7s-T2A-tdTomato calcium reporter with a removable 3xP3-DsRed eye marker. Pick it to express a bright calcium sensor and a red reference fluorophore from a single UAS transcript.",
  pJAT71 = "UAS-driven jGCaMP7s-T2A-tdTomato calcium reporter with a removable ie1XP3-DsRed marker. The abdomen/mouthparts-marked alternative to pJAT58 for the same dual sensor/reference cassette.",
  pJAT75 = "Heat-shock-inducible flippase plasmid carrying an ie1-EGFP marker. Use it to excise FRT-flanked cassettes after establishing an edit.",
  pJAT18 = "Inversion plasmid for an N-terminal split-GFP assay in flight muscles. Pairs with a C-terminal partner to report a specific protein-protein interaction or conformation by split-GFP reconstitution.",
  pJAT27 = "Inversion plasmid for a C-terminal split-GFP assay in flight muscles. The C-terminal partner to pJAT18's N-terminal fragment; co-express the two to report a protein interaction or conformation change by split-GFP reconstitution in muscle.",
  pJAT21 = "Inversion plasmid for an N-terminal split-GFP assay in the abdomen and mouthparts. The abdomen/mouthparts-marked counterpart to pJAT18, used when you want the split-GFP readout reported outside the flight muscles.",
  pJAT25 = "Inversion plasmid for a C-terminal split-GFP assay in the abdomen and mouthparts. The abdomen/mouthparts-marked counterpart to pJAT27.",
  pJAT22 = "Inversion plasmid for an N-terminal split-GFP assay scored by the eye. The eye-marked counterpart to pJAT18, useful when you want the split-GFP reconstitution readout to appear as an eye phenotype rather than in the flight muscles.",
  pJAT23 = "Inversion plasmid for a C-terminal split-GFP assay scored by the eye. The eye-marked counterpart to pJAT27; co-express it with an N-terminal partner to report a protein interaction or conformation change as an eye phenotype."
)

lines <- readLines(yml, warn = FALSE)

# Idempotency: drop any existing description lines first so re-running replaces
# them instead of creating a duplicate YAML key.
lines <- lines[!grepl("^    description:", lines)]

# Locate each plasmid block by its name line, then insert a description line
# immediately after that block's `purpose:` line.
name_idx <- grep("^  - name: \"", lines)
if (length(name_idx) == 0L) stop("no plasmid blocks found in ", yml)

out <- character(0)
i <- 1L
added <- 0L
while (i <= length(lines)) {
  out <- c(out, lines[i])
  if (i %in% name_idx) {
    nm <- sub("^  - name: \"([^\"]+)\".*$", "\\1", lines[i])
    # find the purpose line for this block (before the next name line)
    j <- i + 1L
    purpose_idx <- NA_integer_
    while (j <= length(lines) && !(j %in% name_idx)) {
      if (grepl("^    purpose:", lines[j])) { purpose_idx <- j; break }
      j <- j + 1L
    }
    if (!is.na(purpose_idx)) {
      out <- c(out, lines[(i + 1L):purpose_idx])
      i <- purpose_idx
      desc <- DESCRIPTIONS[[nm]]
      if (!is.null(desc)) {
        out <- c(out, sprintf("    description: \"%s\"", gsub("\"", "'", desc)))
        added <- added + 1L
      }
    }
  }
  i <- i + 1L
}

writeLines(out, yml)
cat(sprintf("added/replaced %d description(s) across %d plasmid(s)\n",
            added, length(name_idx)))

