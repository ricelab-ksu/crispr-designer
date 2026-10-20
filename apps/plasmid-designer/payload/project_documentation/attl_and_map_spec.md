# Gateway `attL` Donor Templates & Interactive Vector Map Specification

**Document ID:** `PJD-SPEC-ATTL-MAP-1.0` | **Date:** 2026-10-17 | **Status:** Decision-Grade Specification  
**Target Workspace:** `/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer`

---

## PART A — Gateway `attL` Donor Templates (`donor_templates/gateway_attl_templates.yml`)

### A.1 Janelia Atalanta Dual-Gateway Donor Design

The Janelia Atalanta (`pJAT`) CRISPR/Cas9 homology-directed repair (HDR) platform (Stern, Kim & Behrman 2024, *bioRxiv* `2023.06.17.545412`, [doi:10.1101/2023.06.17.545412](https://doi.org/10.1101/2023.06.17.545412); Addgene `#204289–204312`, [Addgene Article 28238211](https://www.addgene.org/browse/article/28238211/); Stern 2024, *protocols.io* V.3, [doi:10.17504/protocols.io.bp2l6bjokgqe/v3](https://dx.doi.org/10.17504/protocols.io.bp2l6bjokgqe/v3); US Patent App. [US20240409947A1](https://www.patents-review.com/a/20240409947-vectors-methods-efficient-cloning-homology-directed-repair.html)) enables one-pot Multisite Gateway cloning of two linear dsDNA fragments (typically synthesized by Twist Bioscience with synthesis adapters left attached) into a `pJAT` destination vector:

1. **Synthesized Fragment 1 (`fragment1_attL3_arm_trna_attL4` in `R/pjat_plasmids.R::pjat_assemble_fragments()`, lines 580–642):**
   - Logical contract: `attL3 + 5' homology arm (left_arm, 250–1000 bp) + tRNA–sgRNA–tRNA + attL4`
2. **Synthesized Fragment 2 (`fragment2_attL5_arm_attL6` in `R/pjat_plasmids.R::pjat_assemble_fragments()`, lines 581–649):**
   - Logical contract: `attL5 + 3' homology arm (right_arm, 250–1000 bp) + attL6`

**Physical Orientation & Pol III Splicing Mechanics (`protocols.io` V.3 Step 4.1–4.2 & `Gateway arms.gb`):**
- In the 23 `U6`-containing `pJAT` backbones (`pJAT7–pJAT75`), a *Drosophila* `dU6-3` Pol III promoter sits immediately upstream of `attR3`. In the published GenBank template `Gateway arms.gb` (`LOCUS attL3-tRNA-gRNA` and `2R_5.gb`), the `tRNA–gRNA–tRNA–U6_3'_terminator` block is placed adjacent to `attL3` (`attL3 + 23bp_leader + tRNA:Ala-CGC-1-1 + 20nt_spacer + sgRNA_scaffold + tRNA:Gln + U6_3'_term + Left_Arm + attL4_(reversed)`), so Pol III transcription from backbone `dU6-3` reads across the recombined `attB3` site through the `tRNA–sgRNA–tRNA` array and terminates at `snRNA:U6:96Ac 3'` prior to the homology arm.
- Endogenous eukaryotic tRNA processing enzymes (RNase P and RNase Z) cleave the primary transcript at the 5' and 3' boundaries of the flanking tRNAs (Port & Bullock 2016, *Nat. Methods* 13:852–854, [doi:10.1038/nmeth.3972](https://doi.org/10.1038/nmeth.3972)), releasing a precise 20-nt `sgRNA` without requiring a 5'-`G` nucleotide.
- In dual-`U6` plasmids (`pJAT30`, `pJAT31`, `pJAT32`, `pJAT39`), a second `U6-3` promoter sits in reverse orientation downstream of `attR6`, enabling an optional second `tRNA–sgRNA–tRNA` array on `attL6-tRNA-gRNA_(reversed)`.
- Both `attL3`/`attL4` and `attL5`/`attL6` pairs are synthesized oriented **inwards** toward the homology arm (`attL3` and `attL5` on the `+` strand; `attL4` and `attL6` reverse-complemented on the `-` strand; `protocols.io` V.3 Step 4.1).

---

### A.2 Why the 24 `pJAT` Destination Plasmids Contain `attR3–attR6` (Not `attL`)

All 24 `pJAT` destination plasmids (`Addgene 204289–204312`; `DGRC 1632` for `pJAT55`) contain four **124 bp `attR` sites** (`attR3`, `attR4`, `attR5`, `attR6`) in the cassette order `dU6-3 – attR3 – AmpR – ccdB – attR4 – [PAYLOAD] – attR5 – ChlorR – sacB – attR6` (`aadA` Spectinomycin resistance on the vector backbone).

1. **LR Clonase Reaction (`attL × attR -> attB + attP`):** Bacteriophage $\lambda$ Integrase (`Int`) and Integration Host Factor (`IHF`) in **LR Clonase II** catalyze directional recombination between each `attL` site on the synthesized entry/donor fragments and its cognate `attR` site on the destination plasmid (`attL3 × attR3 -> attB3 + attP3`, `attL4 × attR4 -> attB4 + attP4`, `attL5 × attR5 -> attB5 + attP5`, `attL6 × attR6 -> attB6 + attP6`).
2. **Excision of Dual Counter-Selectable Cassettes:** Recombination replaces `attR3–AmpR–ccdB–attR4` with Fragment 1 (flanked by `attB3`/`attB4`) and replaces `attR5–ChlorR–sacB–attR6` with Fragment 2 (flanked by `attB5`/`attB6`). Plating transformed *E. coli* on **Spectinomycin + 10% Sucrose** selects against both un-recombined `ccdB` (gyrase poison) and `sacB` (lethal levansucrase on sucrose) (Stern et al. 2024; `US20240409947A1`).
3. **Local Register Citation:** As documented in `projects/Cas9_sgRNA_design/project_documentation/atlanta_plasmid_register.md` (lines 32–35 and 117–122):
   > *"attL sites are absent from these plasmids by design. LR recombination consumes attL and leaves attB, so a finished pJAT plasmid contains attR, not attL. The four attL sites needed for the synthesized donor fragments come from the Gateway donor templates (protocols.io). They are never invented by the code: if not supplied, `pjat_assemble_fragments()` returns status `ARM_ONLY` and emits the arms with `[attL3]`-style placeholders."*

---

### A.3 Published `attL3`, `attL4`, `attL5`, `attL6`, and *Drosophila* tRNA Sequences

All sequences below were extracted verbatim from the published GenBank attachment **`Gateway arms.gb`** (`p6xpbnxwf1.gb`, dated `25-MAR-2024`) and cross-verified against **`attL_pieces.gb`** (`fvkgbnxwf2.gb`) and the Twist synthesis example **`2R_5.gb`** (`f3cbbnxwf3.gb`) in Stern DL (2024) *Janelia Atalanta series plasmid cloning V.3* ([doi:10.17504/protocols.io.bp2l6bjokgqe/v3](https://dx.doi.org/10.17504/protocols.io.bp2l6bjokgqe/v3); static attachment mirror: `https://static.yanyin.tech/literature_test/protocol_io_true/protocols.io.bp2l6bjokgqe/p6xpbnxwf1.gb`).

| Key | Length | Strand in Synthesized Fragment | Status | Source URL / Accession / Coordinates |
|---|---:|---|---|---|
| `attl3` | 96 bp | `+` (5' end of Fragment 1) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb` (`p6xpbnxwf1.gb`), `LOCUS attL3-tRNA-gRNA` (`urn.local...62-hiscem7`), pos `1..96` (`/label="attL3"`). Confirmed in `2R_5.gb` (`1..96`). |
| `attl4` | 97 bp | `-` / inward (`attL4_(reversed)`, 3' end of Fragment 1) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb` (`p6xpbnxwf1.gb`), `LOCUS attL4_(reversed)` (`urn.local...63-hiscem7`), pos `1..97` (`complement(1..97)` `/label="attL4"`). Confirmed in `2R_5.gb` (`1087..1183`). |
| `attl5` | 97 bp | `+` (5' end of Fragment 2) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb` (`p6xpbnxwf1.gb`), `LOCUS attL5`, pos `1..97` (`/label="attL5"`). Confirmed in `2R_5.gb` (`LOCUS attL5-2R_5_R_600-attL6`, `1..97`). |
| `attl6` | 97 bp | `-` / inward (`attL6_(reversed)`, 3' end of Fragment 2) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb` (`p6xpbnxwf1.gb`), `LOCUS attL6_(reversed)`, pos `1..97` (`complement(1..97)` `/label="attL6"`). Confirmed in `2R_5.gb` (`694..790`). |
| `trna` (`trna_5p_ala`) | 72 bp | `+` (5' tRNA spacer upstream of 20-nt protospacer) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb`, `LOCUS attL3-tRNA-gRNA`, pos `120..191` (*D. melanogaster* `tRNA:Ala-CGC-1-1`, FlyBase `FBgn0050198` / `FBtr0072445`, RNAcentral `URS0000067406_7227`). |
| `trna_3p_gln` | 58 bp (core) / 78 bp (full) | `+` (3' tRNA downstream of `sgRNA` scaffold) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb`, `LOCUS attL3-tRNA-gRNA`, pos `308..365` (58 bp `/label="tRNA"`) preceded by 20 bp leader at `288..307` (`protocols.io` V.3 Step 4.2 note per Filip Port; Port & Bullock 2016 `tRNA:Gln-CTG`). |
| `u6_leader_5p` | 23 bp | `+` (leader between `attL3` and `tRNA:Ala-CGC-1-1`) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb`, `LOCUS attL3-tRNA-gRNA`, pos `97..119` (`/label="eORI_mE_46"`). |
| `sgrna_scaffold` | 76 bp | `+` (SpCas9 guide scaffold immediately 3' of spacer) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb`, `LOCUS attL3-tRNA-gRNA`, pos `212..287` (`/label="sgRNA scaffold"`). |
| `u6_terminator_3p` | 50 bp | `+` (*D. melanogaster* `snRNA:U6:96Ac 3'` terminator) | **VERIFIED** | `protocols.io` `bp2l6bjokgqe/v3` `Gateway arms.gb`, `LOCUS attL3-tRNA-gRNA`, pos `366..415` (`/label="snRNA:U6:96Ac 3'"`). |

#### Verbatim Nucleotide Sequences (`5' -> 3'` on Synthesized Top Strand)

- **`attl3` (96 bp, `VERIFIED`):** `AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTTTTATAATGCCAACTTTGTATAATAAAGTTG`
- **`attl4` (97 bp, inward/reversed `attL4_(reversed)`, `VERIFIED`):** `CAACTTTTCTATACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG` *(Canonical forward `+` strand: `CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATGCCAACTTTGTATAGAAAAGTTG`)*
- **`attl5` (97 bp, `VERIFIED`):** `CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTATACAAAAGTTG`
- **`attl6` (97 bp, inward/reversed `attL6_(reversed)`, `VERIFIED`):** `CAACTTTTTAATACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG` *(Canonical forward `+` strand: `CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATGCCAACTTTGTATTAAAAAGTTG`)*
- **`trna` / `trna_5p_ala` (*D. melanogaster* `tRNA:Ala-CGC-1-1`, 72 bp, `VERIFIED`):** `GGGGACGTAGCTCAGTGGTAGAGCGCTCGCTTCGCATGTGAGAAGTCCCGGGTTCAAACCCCGGCGTCTCCA`
- **`trna_3p_gln` (78 bp full = 20 bp leader `288..307` + 58 bp core `308..365`, `VERIFIED`):**
  - 20 bp leader (`288..307`): `TAACAAAGCACCAGTGGTCT`
  - 58 bp core (`308..365`): `AGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA`
  - 78 bp full (`288..365`): `TAACAAAGCACCAGTGGTCTAGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA`
- **`u6_leader_5p` (`eORI_mE_46`, 23 bp, `97..119`, `VERIFIED`):** `CTTTGAGTGTGTGTAGACATCAA`
- **`sgrna_scaffold` (SpCas9 scaffold, 76 bp, `212..287`, `VERIFIED`):** `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC`
- **`u6_terminator_3p` (`snRNA:U6:96Ac 3'`, 50 bp, `366..415`, `VERIFIED`):** `TTTTTTGCCTACCTGGAGCCTGAGAGTTGTTCAATAAAATAAAAATGTTT`

---

### A.4 Contract for `donor_templates/gateway_attl_templates.yml` and `pjat_assemble_fragments()`

#### 1. YAML Schema (`donor_templates/gateway_attl_templates.yml`)

```yaml
version: "1.0.0"
schema: "gateway_attl_templates/1.0"
label: "template from published protocol"
attl3: "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTTTTATAATGCCAACTTTGTATAATAAAGTTG"
attl4: "CAACTTTTCTATACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG"
attl5: "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTATACAAAAGTTG"
attl6: "CAACTTTTTAATACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG"
trna: "GGGGACGTAGCTCAGTGGTAGAGCGCTCGCTTCGCATGTGAGAAGTCCCGGGTTCAAACCCCGGCGTCTCCA"
cassette_parts:
  u6_leader_5p: "CTTTGAGTGTGTGTAGACATCAA"
  trna_5p_ala: "GGGGACGTAGCTCAGTGGTAGAGCGCTCGCTTCGCATGTGAGAAGTCCCGGGTTCAAACCCCGGCGTCTCCA"
  sgrna_scaffold: "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
  trna_3p_leader_20bp: "TAACAAAGCACCAGTGGTCT"
  trna_3p_core_58bp: "AGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA"
  trna_3p_full_78bp: "TAACAAAGCACCAGTGGTCTAGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA"
  u6_terminator_3p: "TTTTTTGCCTACCTGGAGCCTGAGAGTTGTTCAATAAAATAAAAATGTTT"
provenance:
  source: "Stern DL (2024) Janelia Atalanta series plasmid cloning V.3 (Gateway arms.gb); Stern, Kim & Behrman (2024) bioRxiv 2023.06.17.545412"
  url: "https://dx.doi.org/10.17504/protocols.io.bp2l6bjokgqe/v3"
  retrieved: "2026-10-17"
  verified: true
```

#### 2. Consumption by `pjat_assemble_fragments()` (`R/pjat_plasmids.R`)

- **Loader `pjat_load_attl_templates(root = PJAT_PROJECT_ROOT)`:** Reads `donor_templates/gateway_attl_templates.yml`, asserts the required keys (`version`, `attl3`, `attl4`, `attl5`, `attl6`, `trna`, and `provenance` with subkeys `source`, `url`, `retrieved`, `verified`), validates `^[ACGT]+$` on each nucleotide string, and returns the parsed list.
- **Consumption in `pjat_assemble_fragments()`:**
  - Signature: `pjat_assemble_fragments(left_arm, right_arm, spacer, attl3 = NULL, attl4 = NULL, attl5 = NULL, attl6 = NULL, trna = NULL)`.
  - When `Plasmid_Designer.R` (or a caller) supplies `attl3`, `attl4`, `attl5`, `attl6`, and `trna` from `pjat_load_attl_templates()` (or user-overridden text inputs):
    - `array` is built as `paste0(toupper(trna), "-", spacer, "-", toupper(trna))` instead of `"[tRNA]-<spacer>-[tRNA]"`.
    - `have_all` (`all(vapply(list(attl3, attl4, attl5, attl6), function(x) !is.null(x) && nzchar(x), logical(1)))`) evaluates to `TRUE`.
    - `status` transitions from `"ARM_ONLY"` to **`"READY"`**, `checks$attl_sites_supplied` becomes `TRUE`, and `fragment_parts$f1`/`f2` contain actual nucleotide sequences instead of `"[missing]"` / `"[attL3]"` placeholders.

---

### A.5 Biology Guardrail: Labeling & User Override Contract

Because `attL` sites are supplied on the synthesized linear fragments or entry clones (and are absent from the `pJAT` destination plasmids), the implementation must enforce:

1. **Mandatory Provenance Label:** When pre-populated from `donor_templates/gateway_attl_templates.yml`, the UI and exports must label the sequences **`"template from published protocol"`** (`protocols.io bp2l6bjokgqe/v3`, Stern et al. 2024).
2. **User Override Controls:** `Plasmid_Designer.R` must provide:
   - A toggle (`input$donor_supply_attl`, default `TRUE`: *"Supply attL3–attL6 & tRNA (template from published protocol)"*).
   - Editable inputs (`input$attl3`, `input$attl4`, `input$attl5`, `input$attl6`, `input$trna`) pre-filled from `gateway_attl_templates.yml`, plus a *"Reset to published template"* action button.
   - If the user edits any sequence, `donor_provenance` / `attl_source` switches from `"template from published protocol"` to `"user override"`. If the user unchecks `input$donor_supply_attl` or clears any `attL` input, `pjat_assemble_fragments()` receives `NULL`/`""` and returns `status = "ARM_ONLY"`.

---

## PART B — Interactive Vector Map Specification

### B.1 Comparative Evaluation of Candidate Approaches (GitHub & npm Facts as of 2026-10-17)

| Criterion | `Lattice-Automation/seqviz` | `stothard-group/cgview-js` (`sciguy/cgview-js`) | `TeselaGen/openVectorEditor` (`@teselagen/ove`) | **Zero-Dependency Inline SVG/HTML + Vanilla JS (R Shiny)** |
|---|---|---|---|---|
| **Repository & Package** | [Lattice-Automation/seqviz](https://github.com/Lattice-Automation/seqviz) (`npm: seqviz@3.10.25`) | [sciguy/cgview-js](https://github.com/sciguy/cgview-js) (`npm: cgview@1.8.2`, redirects from `stothard-group/cgview-js`) | [TeselaGen/openVectorEditor](https://github.com/TeselaGen/openVectorEditor) (archived; moved to `TeselaGen/tg-oss` `@teselagen/ove@0.8.42`) | First-party R (`R/plasmid_map.R`) + local `www/pjd_map.js` |
| **License** | **MIT** ([api.github.com/repos/Lattice-Automation/seqviz](https://api.github.com/repos/Lattice-Automation/seqviz)) | **Apache-2.0** ([api.github.com/repos/sciguy/cgview-js](https://api.github.com/repos/sciguy/cgview-js); Grant & Stothard 2026, *JOSS* [doi:10.21105/joss.09930](https://doi.org/10.21105/joss.09930)) | **MIT** ([api.github.com/repos/TeselaGen/openVectorEditor](https://api.github.com/repos/TeselaGen/openVectorEditor)) | **Project license** (100% first-party R/JS) |
| **Maintenance & Activity** | **Active:** 331 stars, 70 forks, 24 open issues; latest release `v3.10.25` (`2026-09-23`), last push `2026-10-13`. | **Active:** 49 stars, 5 forks, 2 open issues; latest release `v1.8.2` (`2026-05-04`), last push `2026-10-13`. | **Standalone repo ARCHIVED `2024-01-22`** (last push `2023-07-28`, last tag `app-v18.1.0` `2022-08-10`, 202 stars, 53 open issues). Monorepo `@teselagen/ove@0.8.42` active. | **First-party:** maintained directly in `R/` alongside `Plasmid_Designer.R`, zero upstream breakage risk. |
| **Bundle Size** | **229.2 KB min** / **67.5 KB gzip** (`2.32 MB` unpacked on npm, 3 direct deps; [bundlephobia `seqviz@3.10.25`](https://bundlephobia.com/api/size?package=seqviz%403.10.25)). | **288.9 KB min** / **82.3 KB gzip** (`6.64 MB` unpacked on npm, 2 deps `d3@7.0.1` + `svgcanvas@^2.3.0`; [bundlephobia `cgview@1.8.2`](https://bundlephobia.com/api/size?package=cgview%401.8.2)). | **1.83 MB min** / **713.5 KB gzip** (`22.81 MB` unpacked on npm, 47 runtime deps; [bundlephobia `@teselagen/ove@0.8.42`](https://bundlephobia.com/api/size?package=%40teselagen%2Fove%400.8.42)). | **~4 KB local JS/CSS** (`0 KB` third-party libraries); SVG payload ~12–25 KB per render. |
| **Offline Capability** | **Conditional:** fetches `"Roboto Mono:300,400,500"` from external servers on mount via `webfontloader` unless `disableExternalFonts: true` is passed (`README.md`). | **Full offline** once `dist/cgview.min.js` + `dist/cgview.css` + `d3@7.0.1` are vendored locally into `www/`. | **Poor:** heavy React/Redux/BlueprintJS v3 asset graph; UMD bundle requires multiple CSS/icon font files to avoid broken UI offline. | **100% offline by construction:** inline `<svg>` + local `www/pjd_map.js`, system font stack, zero network calls. |
| **Circular + Feature Tracks** | **Yes:** SVG circular (`Circular`) & base-level linear (`Linear`) views, directional annotations, enzyme cut sites. | **Yes:** HTML5 Canvas circular & linear plasmid/genome maps, concentric feature tracks, legends. | **Yes:** Full plasmid editor with circular & linear views, restriction sites, ORFs, and sequence editing. | **Yes:** Parametric `<svg>` circular ring + concentric strand/type tracks and linear horizontal track view with directional arrow arcs. |
| **Hover / Zoom** | Click selection & scroll rotation (`rotateOnScroll`); zoom switches to linear base view rather than smooth circular SVG zoom. | Smooth Canvas zoom/pan, feature hover highlighting, track toggles, and reset controls. | Hover tooltips, click selection, linear zoom; heavy multi-panel editor toolbar. | **Yes:** Floating hover tooltips (`label`, `type`, `start..end`, length, strand), click-to-select, and smooth wheel/button `viewBox` zoom & pan. |
| **R / Shiny Integration Effort** | **Medium:** `seqviz.Viewer(el, props)` UMD bundle embeds React, or requires React 18/19 + custom Shiny message handler + `disableExternalFonts: true`. | **Low–Medium:** Vanilla JS `new CGView.Viewer('#id', config)`; requires vendoring `cgview.min.js` + `cgview.css` + `d3` into `www/` and custom Shiny JS glue. | **Very High:** React 18 + Redux + BlueprintJS v3 CSS clashes with Shiny's Bootstrap 3 DOM/CSS. | **Low (Lowest):** Pure R function returning `shiny::tagList()` rendered via `renderUI()` + ~120 lines of vanilla JS in `www/pjd_map.js`. |

---

### B.2 Recommended Architecture: Primary Approach + Offline Fallback

1. **Primary Approach — Zero-Dependency Parametric Inline SVG/HTML + Scoped Vanilla JS (`R/plasmid_map.R` + `www/pjd_map.js`):**
   - **Rationale:**
     - **Zero build/dependency hazard:** Uses base R + `shiny` directly. Eliminates `openVectorEditor`'s 1.83 MB bundle and archived repo status (`2024-01-22`), avoids `seqviz`'s default `webfontloader` Google Font network request, and avoids `cgview-js`'s 288.9 KB `d3`+Canvas bundle while keeping crisp vector text at all zoom levels.
     - **Unified Circular & Linear Geometry:** `.pjd_interactive_map_ui()` computes SVG polar arc paths (`<path d="M ... A ...">` with directional arrowheads and concentric track offsets for overlapping features) in `mode = "circular"` and horizontal arrow tracks in `mode = "linear"`.
     - **Interactive Hover, Click & Zoom:** `www/pjd_map.js` (~4 KB vanilla JS) provides floating hover cards (`label`, `type`, `start..end`, length, strand), click-to-inspect events (`Shiny.setInputValue`), and smooth mouse-wheel / button (`+`, `-`, `Reset`) `viewBox` zoom and pan.
2. **Offline Fallback — Zero-JS Native SVG `<title>` Tooltips + Local Addgene/SnapGene PNG Cache (`R/addgene_maps.R`):**
   - **Rationale:**
     - Every `<g class="pjd-map-feature">` embeds a native SVG `<title>` tag (`label [type]: start..end (strand, width bp)`), so hover tooltips and vector rendering work with **zero JavaScript**.
     - For unmodified catalog backbones (`pJAT7–pJAT77`), the local Addgene/SnapGene PNG cache (`www/addgene_maps/addgene-plasmid-<id>-map.png` via `R/addgene_maps.R`) remains available as a static reference fallback.

---

### B.3 Exact R & JavaScript Contract

#### 1. Renderer Function Signature (`R/plasmid_map.R`)

```r
.pjd_interactive_map_ui <- function(record,
                                    mode = c("circular", "linear"),
                                    height = "520px")
```

- **`record`:** Normalized plasmid list (Section B.3.2). Returns a `.cas9-card` empty-state box if `record` is `NULL` or `size_bp <= 0`.
- **`mode`:** `match.arg(mode)` — `"circular"` (polar ring with concentric feature tracks, kb tick marks, center title + `size_bp`) or `"linear"` (horizontal ruler `1..size_bp` with `+`/`-` strand tracks).
- **`height`:** CSS height string (e.g., `"520px"`).
- **Return Value:** `shiny::tagList` containing a `singleton(tags$head(tags$script(src = "pjd_assets/pjd_map.js")))` and a `.pjd-interactive-map` container holding toolbar buttons (`Zoom +`, `Zoom -`, `Reset`), legend badges, inline `<svg>`, and a floating `<div class="pjd-map-tooltip" hidden></div>`.

#### 2. Input `record` Shape (`record$features` `data.frame`)

```r
list(
  name     = character(1),   # e.g., "pJAT74" or "pJAT74_assembly"
  size_bp  = integer(1),     # total plasmid length in bp
  topology = character(1),   # "circular" or "linear"
  sequence = character(1),   # full nucleotide sequence ("" allowed)
  features = data.frame(
    type   = character(),    # "attR", "promoter", "CDS", "insert", "spacer", "attP", "rep_origin", ...
    start  = integer(),      # 1-based start coordinate
    end    = integer(),      # 1-based end coordinate
    strand = character(),    # "+" or "-"
    label  = character(),    # feature display name
    color  = character(),    # hex color string (populated via .pjd_feature_color(type, label) if NA)
    stringsAsFactors = FALSE
  )
)
```

**Canonical Color Palette (`.pjd_feature_color(type, label)`):**
- `attR` / `attL` / `attB`: `#dc2626` (red) | `attP`: `#059669` (emerald) | `promoter` (`dU6-3`, `3xP3`, `ie1`, `UAS`, `nanos`, `vasa`): `#0284c7` (sky)
- `insert` (`user_insert`): `#7c3aed` (violet) | `spacer` (`spacer_20nt`): `#d97706` (amber)
- `CDS` / `gene`: `#2563eb` (blue; `#be123c` crimson for `ccdB`/`sacB`) | `rep_origin`: `#64748b` (slate) | default: `#475569`

#### 3. JavaScript Message Contract (`www/pjd_map.js`)

- **Delegated DOM Binding:** Event delegation on `document` handles hover tooltips, wheel/drag `viewBox` zoom/pan, and toolbar button clicks for any `.pjd-interactive-map` rendered via `renderUI()`.
- **Client -> Server Event (`pjd_map_selected_feature`):** Clicking `.pjd-map-feature` calls `Shiny.setInputValue("pjd_map_selected_feature", { label, type, start, end, strand, nonce: Date.now() }, { priority: "event" })`.
- **Server -> Client Handler (`pjd_map_focus`):** `Shiny.addCustomMessageHandler("pjd_map_focus", function(msg) { ... })` highlights and zooms to `{ container_id, start, end, label }`.

#### 4. Dual-Stage Rendering: `read_pjat_genbank()` (Backbone) and `pjd_assembly()` (Assembled Vector)

1. **Backbone Adapter (`.pjd_record_from_genbank(gb, name = NULL)`):**
   - Accepts `gb <- read_pjat_genbank(file)` (`R/pjat_plasmids.R`, line 110).
   - Converts `gb$features` (`key, location, start, end, strand, label, product, qualifiers`) into `features` (`type = key`, `start`, `end`, `strand`, `label = ifelse(nzchar(label), label, ifelse(nzchar(product), product, key))`, `color = .pjd_feature_color(type, label)`), dropping `key == "source"` or `NA` coordinates.
2. **Post-Assembly Adapter (`.pjd_record_from_assembly(a)` where `a <- pjd_assembly()`):**
   - Takes `a$genbank` (`read_pjat_genbank()` record) and `a$assembly` (`assemble_expression_vector()` output, `R/plasmid_assembly.R`, line 185).
   - For replacement window `[w_start, w_end] = [a$assembly$window$start, a$assembly$window$end]` and `delta = a$assembly$insert_len - (w_end - w_start + 1L)`:
     - Retains backbone features with `end < w_start` unchanged;
     - Shifts backbone features with `start > w_end` by `+ delta`;
     - Drops excised features strictly inside `[w_start, w_end]`;
     - Appends `a$assembly$features` (`user_insert` at `insert_start..insert_end` and `spacer_20nt` at `insert_start..(insert_start+19L)`), assigns `color`, and sets `size_bp = a$assembly$final_len`.
   - Both stages pass their resulting `record` directly into `.pjd_interactive_map_ui(record, mode = input$map_mode %||% "circular", height = "520px")`.

---

### B.4 Offline Requirement (Hard Runtime Guardrail)

1. **No CDN at Runtime:** Runtime HTML/JS/CSS must never load scripts, stylesheets, or fonts from external CDNs (`unpkg.com`, `jsdelivr.net`, `cdnjs.cloudflare.com`, `fonts.googleapis.com`).
2. **Local Asset Registration (`www/`):** All static assets (`www/pjd_map.js`, `www/addgene_maps/*.png`, or any vendored bundle under `www/vendor/`) must reside inside `/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer/www/` and be mounted via `shiny::addResourcePath()` in `Plasmid_Designer.R` (matching lines 71–75).




