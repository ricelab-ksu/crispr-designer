# Plasmid Designer

**Janelia Atalanta (pJAT) vector design suite for CRISPR/Cas9 + Gateway cloning.**

A Shiny-based desktop application for assembling pJAT expression vectors in *Drosophila*. Design Cas9 gRNA targets, insert marker-donor cassettes, and export GenBank/FASTA/JSON/ZIP files — all from a clickable one-click mac app.

---

## Quick Start

### One-click (recommended) — macOS

Double-click **`Plasmid Designer.app`** in this folder.

What happens when you double-click:

1. A new Terminal tab opens, preflights the drive + R environment
2. Starts the Shiny app on a local server at **http://127.0.0.1:7790**
3. Your default browser auto-opens to the app
4. Press `Ctrl+C` in that Terminal tab to stop the server when done

**First time (macOS Gatekeeper):** If macOS blocks the app, right-click the `Plasmid Designer.app` → **Open** → **Open** again.

### Alternatively

Run from Terminal:

```bash
# Using the launch script
./launch_app.sh

# Or directly
Rscript Plasmid_Designer.R
```

Or via `shiny::runApp(".")` in RStudio.

---

## Prerequisites

- **macOS 10.15+** (recommended) or Linux
- **R** with the following packages installed:
  - `shiny`
  - `Biostrings`
  - `yaml`
  - `jsonlite`
  - `data.table`
  - `plasmidplot` (optional — enables the 4th "plasmidplot" vector-map view;
    the other three views work without it)

Install missing packages with:

```r
install.packages(c("shiny", "Biostrings", "yaml", "jsonlite", "data.table"))
install.packages("plasmidplot")  # optional 4th map view
```

---

## Launch Modes

### GUI (Shiny)

```bash
Rscript Plasmid_Designer.R
# or
./launch_app.sh
```

Opens the interactive browser-based designer at `http://127.0.0.1:7790`.

### Headless CLI

```bash
Rscript Plasmid_Designer.R --gene even skipped --species melanogaster --arm5 1200 --arm3 1500 --out ./crispr_output
```

Supported global options:

| Flag | Meaning |
|------|---------|
| `--help`, `-h` | Print usage and exit |
| `--gene <SYMBOL>` | Species gene symbol (e.g., `even`, `cargo`) |
| `--species <SPECIES>` | Species name (e.g., `melanogaster`) for cross-species lookups |
| `--arm5 <BP>` | 5′ homology arm endpoint (bp, from CDS start) |
| `--arm3 <BP>` | 3′ homology arm endpoint (bp, from CDS end) |
| `--out <DIR>` | Output directory |

---

## Vector maps (Addgene)

Every detail of the selected plasmid lives in one **Selected backbone** card in the
main panel: identity badges (Addgene ID, supplier, size, U6 count, attP), the
design purpose, the four attR coordinates, any supplier warning (e.g. pJAT55 /
DGRC 1632) and the SnapGene vector map — with a link back to the Addgene record.
The sidebar keeps only the plasmid picker.

Addgene's map URLs contain a per-plasmid sequence ID and UUID, so they cannot be
templated from the Addgene ID. Build (or refresh) the local image cache once:

```bash
# fetch every missing map (idempotent — cached PNGs are reused)
Rscript tools/fetch_addgene_maps.R

# re-download everything, or only a few plasmids
Rscript tools/fetch_addgene_maps.R --force
Rscript tools/fetch_addgene_maps.R --ids 204289,204290

# resolve URLs without downloading anything
Rscript tools/fetch_addgene_maps.R --dry-run
```

| Path | Role |
|---|---|
| `tools/fetch_addgene_maps.R` | Scrapes `https://www.addgene.org/<id>/`, downloads the PNG, writes the manifest |
| `donor_templates/addgene_maps.yml` | Generated manifest: page URL, image URL, sequence ID, dimensions, md5, status |
| `www/addgene_maps/addgene_<id>.png` | Downloaded map images (git-ignored; 24 files, ~6 MB) |

The app degrades gracefully: without the cache it shows the Addgene link and the
fetch command in place of the image. Map images are (c) Addgene / SnapGene and
are displayed locally with attribution and a link to the source record.

---

## Sample data, bundled markers, attL templates

**One-click samples.** Section 2 has *Load sample guide bundle (w gene)*, which
loads `uploads/Cas9_w_melanogaster_w_g0351_guides_for_plasmid_designer.csv`
(spacer `ATTGCAGGGTGACAGCGGAG` + 1,000 bp arms). Section 5 has *Load sample donor
FASTA (3xP3-DsRed)*, which loads the 1,861 bp 3xP3-DsRed cassette extracted from
`pJAT30.gbk` between attR4 and attR5.

**Marker sequences are real now.** `marker_cassettes.yml` carries sequence
provenance for all nine markers:

| Source | Markers | Provenance |
|---|---|---|
| Local pJAT GenBank | `3xP3_dsRed` (1,861 bp), `ie1_dsRed` (2,306 bp) | attR4..attR5 payload of pJAT30 / pJAT11 |
| NCBI NuCCORE | `neo` 795, `puro` 600, `hygro` 906, `blast` 393, `URA3` 804, `LEU2` 1,095, `KanMX` 810 bp | CDS extracted from the cited accession |

Rebuild with `Rscript tools/fetch_marker_sequences.R` (`--force` to re-fetch,
`--offline` for the pJAT-only markers). Select a marker and press *Insert selected
marker sequence* to drop the verified sequence into the donor box.

**Gateway attL templates.** Section 4 defaults to the published Janelia Atalanta
`attL3`, `attL4`, `attL5`, `attL6` and *Drosophila* tRNA sequences
(`donor_templates/gateway_attl_templates.yml`, sourced from Stern et al. 2024,
protocols.io `bp2l6bjokgqe/v3`). With them, `pjat_assemble_fragments()` reports
**READY** instead of **ARM_ONLY**. The manual entry boxes remain available under
*I will supply my own attL templates*. Why the plasmids have no attL: pJAT
backbones are destination vectors carrying attR3–attR6; LR Clonase converts
attL × attR → attB, so attL is consumed during cloning. See
`project_documentation/attl_and_map_spec.md`.

## Vector maps

The *Vector map view* switch above the backbone card offers:

- **plasmidplot** (default) — the bundled GenBank file is parsed by the CRAN
  package [plasmidplot](https://cran.r-project.org/package=plasmidplot) and
  drawn with grid graphics (coloured feature arcs with strand arrows and
  labels). The plot is sized to the same footprint as the Addgene PNG for that
  backbone (matching aspect ratio from the cached manifest dimensions).
- **Addgene PNG** — the downloaded SnapGene image at full resolution, with
  zoom in / zoom out / reset controls (25–300%).

After *Assemble Vector*, an **Assembled vector map** card renders the construct
(backbone + inserted payload) with plasmidplot as well.

---

## Project Structure

```
Plasmid_Designer/
├── Plasmid_Designer.R       # Standalone Shiny GUI + headless CLI entry point
├── app.R                    # shiny::runApp() entry point
├── launch_app.sh            # One-click launch script (Terminal)
├── Plasmid\ Designer.app/   # macOS .app bundle (click-to-activate)
│   └── Contents/
│       └── MacOS/
│           ├── plasmid-designer.sh      # Preflight + app launcher
│           ├── plasmid-designer         # Wrapper script
│           └── plasmid-designer.app     # App bundle stub
├── R/
│   ├── pjat_plasmids.R      # pJAT backbone catalog, GenBank parser, attL fragment design
│   ├── addgene_maps.R       # Addgene/SnapGene map manifest, catalog join, PNG QC
│   ├── gateway_attl.R       # Published attL3-6 + tRNA templates, attR extraction
│   ├── marker_sequences.R   # Offline marker sequence catalog/loader
│   ├── sample_inputs.R      # Bundled sample guide bundle + donor FASTA loaders
│   ├── plasmid_map.R        # Interactive circular/linear SVG vector map
│   ├── guide_bundle.R       # Cas9 guide bundle (1.0) reader/validator
│   │                          (accepts JSON or YAML .json sidecar)
│   └── plasmid_assembly.R   # insert_cassette, assemble_gRNA_cassette, etc.
├── donor_templates/
│   ├── pjat_plasmids.yml    # Generated pJAT catalog (attR coords, Addgene IDs)
│   ├── addgene_maps.yml     # Generated vector-map manifest (URLs, sizes, md5)
│   ├── gateway_attl_templates.yml      # Published attL3-6 + tRNA (verified)
│   ├── marker_cassettes.yml            # Marker library + generated sequence fields
│   └── marker_sequences/*.fasta        # Bundled marker sequences (provenance)
├── tools/
│   ├── fetch_addgene_maps.R      # Scrape + download Addgene vector maps (idempotent)
│   ├── fetch_marker_sequences.R  # Extract local + NCBI marker sequences
│   └── make_icon.R               # Build the .app icon
├── www/
│   ├── addgene_maps/        # Downloaded map PNGs (git-ignored cache)
│   └── pjd_assets/          # Interactive map JS + CSS (vendored, offline)
├── tests/
│   └── testthat/
│       ├── test-addgene_maps.R     # URLs, manifest I/O, join, PNG QC, cache
│       ├── test-gateway_attl.R     # attL templates -> READY fragments, attR extraction
│       ├── test-marker_sequences.R # Catalog, FASTA/md5 integrity, provenance
│       ├── test-plasmid_map.R      # Record normalisation + SVG output
│       ├── test-sample_inputs.R    # Sample guide bundle + donor FASTA
│       ├── test-guide_bundle.R
│       ├── test-pjat_plasmids.R
│       └── test-plasmid_assembly.R
├── project_documentation/   # Architecture, GC rules, workflow docs
└── README.md                # (this file)
```

---

## Tests

```bash
cd /Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer
Rscript -e 'suppressWarnings(library(testthat)); test_dir("tests/testthat")'
```

All tests pass (addgene_maps, gateway_attl, guide_bundle, marker_sequences,
pjat_plasmids, plasmid_assembly, plasmid_map, sample_inputs) — 421 assertions.

---

## Port

**7790** (Cas9 app uses 7789 so both can run side by side).

---
