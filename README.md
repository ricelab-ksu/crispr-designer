# CRISPR Cas9 gRNA & Homology Arm Designer for *Drosophila* (`v1.0.0`)

A standalone, reproducible **R / Shiny GUI + Headless CLI** application for designing *Streptococcus pyogenes* Cas9 (**SpCas9**, `20-nt protospacer + NGG PAM`) guide RNAs and flanking **1–2 kb HDR homology arms** across *Drosophila* reference genomes.

Modeled after the lab's **HCR Probe Designer (`HCR_probe_design_v41.R`)** for seamless UX, visual track plotting, NCBI RefSeq on-demand genome management, and Electron sidecar packaging compatibility.

---

## Key Features

1. **Multi-Species *Drosophila* Reference Registry (`R/genome_registry.R`)**
   - Automatically detects all 7 *Drosophila* species folders in `/Volumes/Backup_Plus/projects/bioinformatics/reference_files/drosophila_species_genomes/`:
     - `d_melanogaster` (`GCF_000001215.4_Release_6_plus_ISO1_MT`)
     - `d_simulans` (`GCF_016746395.2_Prin_Dsim_3.1`)
     - `d_eugracilis` (`GCF_018153835.1_ASM1815383v1`)
     - `d_erecta` (`GCF_003286155.1_DereRS2`)
     - `d_teissieri` (`GCF_016746235.2_Prin_Dtei_1.1`)
     - `d_ananassae` (`GCF_017639315.1_ASM1763931v2` — lazy NCBI RefSeq download supported)
     - `d_virilis` (`GCF_003285735.1_DvirRS2` — lazy NCBI RefSeq download supported)
   - Automatic `.fai` FASTA byte-offset indexing for `<10 ms` random genomic subsequence extraction.

2. **Multi-Modal Gene Lookup & Isoform Feature Extraction (`R/gene_lookup.R`)**
   - Look up targets by **Gene Symbol** (`ebony`, `white`, `yellow`, `Abd-B`, `PlexB`), **FlyBase ID** (`FBgn0000527`), **CG Number** (`CG3331`), **NCBI Gene / LOC ID**, **Synonym/Alias**, **Protein Product Name** (for non-*melanogaster* species), or **Genomic Coordinates** (`3R:21231000-21236000`).
   - Extracts all transcript isoforms, per-isoform exons and CDS intervals, **Shared CDS** (coding in every isoform), **5′ UTR**, **3′ UTR**, **Upstream Promoter**, and **First Coding Exons**.

3. **Reproducible SpCas9 gRNA Discovery & Ranking Engine (`R/grna_engine.R`)**
   - Scans both `+` and `-` genomic strands for canonical `20-nt + NGG` SpCas9 sites (plus NAG/NGA off-target evaluation).
   - **Genomic Context Tiering (CHOPCHOP v3 / flyCRISPR style):**
     - **Tier 1 (`early_shared_cds`):** 5′-most shared coding exons (`CDS rel pos <= 0.50` or first coding exon) — highest priority for knockout / coding disruption.
     - **Tier 2 (`5p_utr_promoter`):** 5′ UTR and upstream promoter window — prioritized for regulatory edits or N-terminal knock-ins.
     - **Tier 3 (`late_shared_cds`):** Shared coding exons in the 3′ half of the CDS.
     - **Tier 4 (`isoform_specific_cds`):** Coding exons present in a subset of isoforms.
     - **Tier 5 (`3p_utr` / `intron`):** Non-coding 3′ UTR or intronic regions.
   - **On-Target Efficacy Scoring:** Doench et al. 2016 Rule Set 2 / Azimuth (`0–100`) and Moreno-Mateos et al. 2015 CRISPRscan (`0–100`).
   - **Off-Target Specificity Scoring:** Doench et al. 2016 CFD (`0–100`) and Hsu et al. 2013 MIT (`0–100`) with exact `mm0/mm1/mm2/mm3` hit counts.

4. **Exact SpCas9 Blunt Cut-Site Geometry & 1–2 kb Homology Arm Designer (`R/homology_arms.R`)**
   - SpCas9 generates a blunt double-strand break **3 bp upstream (5′) of the first nucleotide of the NGG PAM on the target strand** (between protospacer positions 17 and 18).
   - User-selectable **1,000–2,000 bp** sliders (default `1,000 bp`, step `100 bp`), symmetric or independent per 5′/3′ arm.
   - Automatic sequence QC (`N` tracts, homopolymers, Shannon entropy, GC%, scaffold boundary truncation) and reading-frame-aware PAM/seed shielding mutation advisor.

---

## Quick Start

### 1. Launch Interactive Shiny GUI

**One-click (recommended):** Double-click `CRISPR Designer.app` in this folder.
- Opens a new Terminal tab, preflights the drive + R environment, then runs `Rscript Cas9_sgRNA_design.R` — your default browser auto-opens to `http://127.0.0.1:7789`.
- Press `Ctrl+C` in that Terminal tab to stop the server when done.
- First time: right-click → **Open** if macOS Gatekeeper blocks it (not code-signed).
- (Under the hood, the `.app` calls `/Users/work/bin/crispr-designer.sh`, modeled on your `scRNA Dashboard.app` pattern.)

**Alternatively:**
```bash
./launch_app.sh
# or directly:
Rscript Cas9_sgRNA_design.R
```

Or via `shiny::runApp(".")` (launches at `http://127.0.0.1:7789`).

### 2. Run in Headless CLI Mode

```bash
Rscript Cas9_sgRNA_design.R --gene ebony --species melanogaster --arm5 1200 --arm3 1500 --out ./crispr_ebony_out
```

Writes `gRNA_candidates.csv`, `homology_arms.fasta`, `design_summary.json`, `isoform_track_plot.pdf`, `CRISPR_design_log.txt`, and `CRISPR_design_bundle.zip`.

### 3. Run Unit & Integration Tests

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
