# Integration Report — CRISPR + Plasmid Designer refactor

**Date:** 2026-10-18 · **Scope:** 7 workstreams across two apps
**Repos:** `/Volumes/Backup_Plus/projects/bioinformatics/projects/{Plasmid_Designer,Cas9_sgRNA_design}`

---

## 1. Test results (integration gate)

| App | Pass | Fail | Warn | Baseline | Δ |
|---|---:|---:|---:|---:|---:|
| Plasmid_Designer | **915** | 0 | 1 | 442 | +473 |
| Cas9_sgRNA_design | **978** | 0 | 0 | 978 | 0 (net; 3 updated for v1.1) |

PD grew from 9 → 16 test files. The single warning is pre-existing
(`plasmidplot` built under R 4.6.1).

New PD test files: `test-backbone_registry.R`, `test-plasmid_assembly_slots.R`,
`test-arm_resolution.R`, `test-guide_bundle_arm_window.R`,
`test-pjat_descriptions.R`, `test-design_schematic.R`, `test-onboarding.R`.

---

## 2. Acceptance criteria (Leader verification checklist)

| # | Gate | Result |
|---|---|---|
| 1 | Baseline captured | ✅ `project_documentation/baseline_tests.md` (PD 442) |
| 2 | pJAT behaviour unchanged | ✅ adapter == `pjat_load_catalog()`; 0 attR mismatches across 24 plasmids |
| 3 | Bundle arms actually used | ✅ `resolve_homology_arms()` is the default path; flanks are fallback only |
| 4 | Arm policy 250–1000 | ✅ no `min=1000,max=2000` sliders remain in either app |
| 5 | Cloning math proven | ✅ length identity, revcomp, wraparound, cut offsets all tested |
| 6 | One payload input | ✅ `use_marker_seq` removed; `donor_source` chooser added |
| 7 | Descriptions present | ✅ 24/24 ≥ 120 chars, all distinct |
| 8 | Live diagram renders | ✅ `pjd_design_schematic()` + orientation-flip tests |
| 9 | Docs complete | ✅ `app_map.md` ×2 (Mermaid + section→server→card maps) |
| 10 | Rules audit | ✅ see §4 |

---

## 3. Defects fixed

| Defect | Before | After |
|---|---|---|
| **Dead bundle arms** | `arm_5p_seq`/`arm_3p_seq` validated then discarded; fragments built from plasmid flanks | `resolve_homology_arms()` uses bundle arms by default (rules R1/R4) |
| **Arm-length mismatch** | CRISPR sliders 1000–2000 bp vs protocol 250–1000 | Both apps 250–1000; over-long arms trimmed from the OUTER end |
| **Duplicate payload doors** | textarea + upload + marker dropdown (3 doors, 1 insert) | one `donor_source` chooser (marker \| paste \| upload) |
| **No live design view** | construct only visible after Assemble | live schematic with slot, orientation, revcomp readout, gRNA site |
| **No plasmid purpose text** | one-line `purpose` only | curated 2–4 sentence `description` for all 24 |
| **pJAT-only cloning** | hard-wired pJAT/attR assumptions | `backbone_registry/1.0` — any GenBank/FASTA + named slots |

---

## 4. Rules audit (R1–R6, from `R/pjat_plasmids.R` header)

| Rule | Requirement | Status |
|---|---|---|
| R1 | Arms precisely flank the cut site | ✅ arms come from the guide bundle's cut-site-anchored design |
| R2 | 250 bp min, 1000 bp max | ✅ enforced in `arm_resolution.R`, `guide_bundle.R`, `pjat_assemble_fragments()` |
| R3 | Spacer 20 nt, PAM excluded | ✅ unchanged; `pjat_validate_spacer()` still enforced |
| R4 | Arms in genomic orientation, never flipped | ✅ bundle arms used as-is; trimming preserves the cut-adjacent end |
| R5 | Synthesis QC (≤10 bp homopolymer, GC 30–70%) | ✅ `pjat_synthesis_qc()` unchanged, still applied |
| R6 | Payload should destroy the target site | ✅ advisory flag preserved (`payload_destroy_target`) |
| Invariant | No fabricated sequence | ✅ trimming copies substrings; nothing is generated |

---

## 5. Interop contract change

`cas9_guide_bundle/1.0` → **`1.1`** (1.0 still accepted on read).

- `arm_policy`: `"250-1000 bp validated homology arms"`
- new sidecar keys: `arm_min_bp`, `arm_max_bp`, `arms_trimmed`, `arms_trim_note`
- producer trims arms > 1000 bp from the outer end (cut-site bases preserved)

Round-trip verified: CRISPR writes 1.1 → PD `read_guide_bundle()` imports with
arms intact.
