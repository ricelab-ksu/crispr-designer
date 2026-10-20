# Project Documentation — CRISPR Cas9 gRNA + Homology-Arm Designer (*Drosophila*)

## Purpose
This folder (`/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design/project_documentation/`) is the single source of truth for the project. The orchestrator agent keeps these files current; all implementation work is delegated per `delegation_plan.md`.

## Files

| File | Purpose |
|------|---------|
| `optimized_prompt.md` | The master prompt (Lyra-optimized v2) governing all orchestrator & subagent work. |
| `delegation_plan.md` | Ticketed Kanban backlog (T1–T7) with owners, dependencies, and acceptance criteria. |
| `tool_survey.md` | Ranked evaluation of candidate tools from `awesome-CRISPR` + GitHub and selected reuse architecture. |
| `genome_inventory.md` | Verified inventory of all 7 *Drosophila* species in `reference_files/drosophila_species_genomes`. |
| `status_board.md` | Live ticket status; updated after every subagent review. |

## Key Decisions
1. **Reuse-First Engine Strategy:** Reuse Bioconductor's CRISPR suite (`crisprDesign`/`crisprScore`/`CRISPRseek`) and/or deterministic R + `Biostrings`/Bowtie/FlashFry scoring for gRNA discovery and on/off-target ranking, paired with a deterministic `Biostrings`/`Rsamtools` module for exact 1–2 kb homology-arm extraction anchored 3 bp 5′ of the PAM.
2. **UI & Packaging Reference:** Mirror the single-file Shiny + Electron sidecar architecture at `/Volumes/Backup_Plus/HCR_probe_designer/HCR_probe_design_v41.R` (v42.2.0, newer than the bundled `.app` copy v42.1.0; `runner.R` + `main.js`).
3. **Multi-Species Genome Handling:** 5 species (`d_melanogaster`, `d_simulans`, `d_eugracilis`, `d_erecta`, `d_teissieri`) have full `genomic.fna` + `genomic.gtf` locally; 2 species (`d_ananassae`, `d_virilis`) currently have `rna.fna` only and will use the HCR app's lazy NCBI RefSeq downloader to fetch `genomic.fna.gz` + `genomic.gtf.gz` when needed.

## Reference Paths
- **Project folder:** `/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design`
- **Documentation folder:** `/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design/project_documentation`
- **Genome files:** `/Volumes/Backup_Plus/projects/bioinformatics/reference_files/drosophila_species_genomes`
- **Local `awesome-CRISPR` repo:** `/Volumes/Backup_Plus/projects/bioinformatics/awesome-CRISPR-master/README.md`
- **Styling & packaging reference Shiny app:** `/Volumes/Backup_Plus/HCR_probe_designer/HCR_probe_design_v41.R` (v42.2.0)

