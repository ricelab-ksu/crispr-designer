# Cline Kanban — Delegation Plan

Master prompt: `optimized_prompt.md`. Project root: `/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design/`.
Orchestrator assigns tickets in dependency order and reviews deliverables **only** against acceptance criteria.

## Backlog

### T1 — Existing-Tool Survey & Engine Selection (Phase 1)
- **Owner:** Research subagent (Cline CLI)
- **Scope:** Evaluate `/Volumes/Backup_Plus/projects/bioinformatics/awesome-CRISPR-master/README.md` and GitHub alternatives; rank candidates by usefulness, reproducibility, R/Shiny compatibility, and *Drosophila* support.
- **Deliverable:** Completed `tool_survey.md` with ranked table and architecture decision.
- **Acceptance:** ≥8 candidate tools evaluated and ranked; reuse/adapt/build verdict + license documented.
- **Status:** **Done** (completed during prompt setup; ready for T4 subagent verification of installed Bioconductor/CLI packages).

### T2 — Genome Data Inventory & Normalization (Phase 2.1)
- **Owner:** Data subagent (Cline CLI / Kanban)
- **Scope:** Inventory `/Volumes/Backup_Plus/projects/bioinformatics/reference_files/drosophila_species_genomes` across all 7 species (`d_melanogaster`, `d_simulans`, `d_eugracilis`, `d_erecta`, `d_teissieri`, `d_ananassae`, `d_virilis`); provide an R helper (reusing the NCBI RefSeq lazy-download pattern from `HCR_probe_design_v41.R`) to decompress `.gz` files, generate `.fai` indices, and fetch missing `genomic.fna.gz` + `genomic.gtf.gz` for `d_ananassae` and `d_virilis`.
- **Deliverable:** `genome_inventory.md` (initial inventory completed) + `R/genome_registry.R`.
- **Acceptance:** Helper resolves genomic FASTA + GTF paths for all 5 complete species immediately and cleanly handles/fetches the 2 RNA-only species.
- **Dependencies:** none.

### T3 — Multi-Species Gene-Lookup Index & API (Phase 2.2)
- **Owner:** R dev subagent
- **Scope:** Build fast gene lookup supporting gene symbol, FlyBase FBgn ID, CG number, NCBI gene ID / `LOC` ID, synonym/alias, protein product name, and genomic coordinates (`chr:start-end`). Use `dmel_gene_reference.csv` + `fbgn_map.csv` + `genomic.gtf` for *D. melanogaster*, and auto-generate cached per-species gene lookup tables from `genomic.gtf(.gz)` for the other *Drosophila* species.
- **Deliverable:** `R/gene_lookup.R` + `tests/testthat/test-gene_lookup.R`.
- **Acceptance:** Unit tests pass for symbol, FBgn, CG/alias, and coordinate queries on `d_melanogaster` (e.g., `Abd-B`, `ebony`, `PlexB`) and ≥1 non-melanogaster species (`d_simulans` or `d_eugracilis`), returning exon/CDS/5′UTR/promoter ranges per isoform.
- **Dependencies:** T2.

### T4 — Reproducible gRNA Candidate Discovery & Ranking Engine (Phase 3)
- **Owner:** Algorithm subagent
- **Scope:** Implement reproducible SpCas9 (`20bp + NGG`) target enumeration and ranking by reusing Bioconductor (`crisprDesign`/`crisprScore`/`CRISPRseek`) and/or a deterministic R + `Biostrings`/Bowtie/FlashFry scoring pipeline. Score candidates by: (1) genomic context tier (5′ UTR/promoter & 5′-most shared coding exons prioritized to disrupt protein coding sequence), (2) on-target efficacy (Rule Set 2 / Azimuth or CRISPRscan), and (3) genome-wide off-target specificity (CFD / MIT score against the species' `genomic.fna`).
- **Deliverable:** `R/grna_engine.R` + `tests/testthat/test-grna_engine.R` + scoring documentation.
- **Acceptance:** Deterministic output for fixed inputs; reproduces ≥3 published *Drosophila* guides in top-ranked candidates; logs exact algorithm and genome versions for reproducibility.
- **Dependencies:** T2, T3.

### T5 — Exact Homology-Arm Designer (Phase 4)
- **Owner:** Algorithm subagent
- **Scope:** Given a selected gRNA target site, compute the exact SpCas9 blunt cut site (**3 bp upstream / 5′ of the first nucleotide of the PAM on the target strand**) and extract flanking 5′ (left) and 3′ (right) homology arms of user-selectable length (**1–2 kb slider / options**, symmetric or independent). Flag internal PAM/protospacer re-cleavage risks, `N` tracts, low-complexity repeats, and chromosome boundary truncation.
- **Deliverable:** `R/homology_arms.R` + `tests/testthat/test-homology_arms.R`.
- **Acceptance:** Unit tests mathematically verify 1-based cut-site and arm coordinates on both `+` and `-` strands across multiple genes (`length(arm_5p) == L`, `length(arm_3p) == R`, adjacent across the exact cut site 3 bp 5′ of PAM); exports FASTA + CSV/JSON + reproducibility log.
- **Dependencies:** T4.

### T6 — Shiny GUI & Electron Bundle (Phase 5)
- **Owner:** UI subagent
- **Scope:** Build the Shiny UI + server (`Cas9_sgRNA_design.R`) modeled on `/Volumes/Backup_Plus/HCR_probe_designer/HCR_probe_design_v41.R` (v42.2.0, newer than the bundled `.app` copy v42.1.0; and compatible with `runner.R` / Electron `main.js`). Include: species selector, gene lookup search bar, isoform + cut-site/homology-arm track plot, interactive ranked candidate table (filterable by 5′ UTR/promoter vs. coding exon), 1–2 kb homology-arm slider/options, nucleotide-level cut-site zoom view, and one-click ZIP/FASTA/CSV download.
- **Deliverable:** Runnable Shiny app in `/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design/` + usage README.
- **Acceptance:** App launches cleanly via `Rscript` and completes the full workflow from gene selection to gRNA ranking and 1–2 kb homology-arm export on real *Drosophila* genome data.
- **Dependencies:** T3, T4, T5.

### T7 — End-to-End Validation & Reproducibility Audit (Phase 6)
- **Owner:** QA subagent
- **Scope:** Execute all Phase 6 validation checks: positive controls on ≥3 *Drosophila* genes, `+`/`-` strand cut-site & homology-arm coordinate verification, off-target cross-check, and repeat-run byte-level reproducibility check.
- **Deliverable:** `project_documentation/validation_report.md` + full `testthat` test runner.
- **Acceptance:** All `testthat` suites pass with 0 failures; `validation_report.md` signed off.
- **Dependencies:** T3, T4, T5, T6.

