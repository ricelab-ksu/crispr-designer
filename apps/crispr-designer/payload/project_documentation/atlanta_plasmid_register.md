# Janelia Atalanta (pJAT) Plasmid Register

Authoritative reference for the pJAT CRISPR/Cas9 HDR platform, used by
`R/pjat_plasmids.R` and the donor-fragment designer.

- **Primary source:** Stern DL, Kim E, Behrman EL. *The Janelia Atalanta plasmids
  provide a simple and efficient CRISPR/Cas9-mediated homology directed repair
  platform for Drosophila.* bioRxiv 2023.06.17.545412.
  doi:[10.1101/2023.06.17.545412](https://doi.org/10.1101/2023.06.17.545412)
- **Addgene article:** https://www.addgene.org/browse/article/28238211/
- **Protocol:** Stern DL. *Janelia Atalanta series plasmid cloning*, protocols.io
  **V.3**. doi:[10.17504/protocols.io.bp2l6bjokgqe/v3](https://dx.doi.org/10.17504/protocols.io.bp2l6bjokgqe/v3)
- **Sequence source:** `/Volumes/Backup_Plus/projects/bioinformatics/reference_files/janelia_atalanta_plasmids/` (24 plasmids, 2026-10-17)
- **Machine-readable catalog:** `donor_templates/pjat_plasmids.yml`
- **Parser & API:** `R/pjat_plasmids.R` · tests: `tests/testthat/test-pjat_plasmids.R`

## Fields common to all plasmids

| Field | Value |
|---|---|
| Vector backbone | Janelia Atalanta cloning plasmid |
| Vector type | CRISPR |
| Selectable markers | ccdB, SacB (double negative selection) |
| Bacterial resistance | Spectinomycin, 50 µg/mL (backbone); Amp and Chlor are in the Gateway arms |
| Growth temperature | 37 °C |
| Growth strain | ccdB Survival |
| Copy number | High Copy |
| Terms | Academic/Nonprofit: UBMTA · Not Available to Industry |
| Depositing lab | David Stern · HHMI, Janelia Research Campus, USA |
| Publication | Stern et al. bioRxiv 2023.06.17.545412 |

Backbone cassette order: `dU6-3 – attR3 – Amp – ccdB – attR4 – [PAYLOAD] –
attR5 – Chlor – sacB – attR6` (pJAT77 omits the U6; pJAT32/30/31/39 add a
second U6 downstream of attR6). The **attR** sites are 124 bp each.

## Catalog

`U6` = number of Drosophila U6-3 promoters available to drive gRNA(s).
`bp` = plasmid size. Source files are Addgene GenBank (`.gbk`) except pJAT55 (DGRC).

| ID | Name | bp | U6 | attP | Purpose |
|---|---|---:|:---:|:---:|---|
| 204289 | pJAT39 | 10,998 | 2 | | Mark HDR integrant with EGFP in flight muscles |
| 204290 | pJAT11 | 10,135 | 1 | ✅ | Insert phiC31 attP site; removable ie1-DsRed (abdomen/mouthparts) |
| 204291 | pJAT31 | 10,941 | 2 | ✅ | Insert attP; ie1-DsRed; second U6 for two gRNAs |
| 204292 | pJAT7 | 9,454 | 1 | ✅ | Insert attP; removable 3XP3-DsRed (eyes) |
| 204293 | pJAT10 | 9,690 | 1 | | Mark HDR integrant with DsRed in flight muscles |
| 204294 | pJAT30 | 10,496 | 2 | | Mark HDR integrant with DsRed; second U6 for two gRNAs |
| 204295 | pJAT77 | 7,560 | **0** | | **Empty backbone**, MCS, no U6 (restriction sites to add flanking U6) |
| 204296 | pJAT74 | 7,986 | **1** | | **Empty backbone**, MCS, one U6 |
| 204297 | pJAT32 | 8,792 | **2** | | **Empty backbone**, MCS, two U6 |
| 204298 | pJAT17 | 15,088 | 1 | | Germline Cas9 (nanos promoter); 3XP3-EYFP (eyes) |
| 204299 | pJAT43 | 13,161 | 1 | | Heat-shock B3 recombinase; ie1-EGFP |
| 204300 | pJAT56 | 13,669 | 1 | | phiC31 (vasa promoter + vasa 3'UTR); 3XP3-EGFP |
| 204301 | pJAT62 | 14,663 | 1 | | phiC31 (vasa promoter + vasa 3'UTR); ie1-EGFP |
| **204302** | **pJAT55** | **13,951** | 1 | | **piggyBac transposase** (alpha-tub-pBac-Trsp); 3XP3-ECFP — **Addgene DISCONTINUED → DGRC 1632** |
| 204303 | pJAT57 | 15,316 | 1 | | UAS-Cas9; DsRed in flight muscles |
| 204304 | pJAT58 | 13,856 | 1 | | jGCaMP7s-T2A-tdTomato (UAS); removable 3XP3-DsRed |
| 204305 | pJAT71 | 14,671 | 1 | | jGCaMP7s-T2A-tdTomato (UAS); removable ie1XP3-DsRed |
| 204306 | pJAT75 | 12,462 | 1 | | Heat-shock flippase; ie1-EGFP |
| 204307 | pJAT18 | 9,620 | 1 | | Inversion: N-terminal split GFP (flight muscles) |
| 204308 | pJAT27 | 10,147 | 1 | | Inversion: C-terminal split GFP (flight muscles) |
| 204309 | pJAT21 | 9,868 | 1 | | Inversion: N-terminal split GFP (abdomen/mouthparts) |
| 204310 | pJAT25 | 10,393 | 1 | | Inversion: C-terminal split GFP (abdomen/mouthparts) |
| 204311 | pJAT22 | 8,987 | 1 | | Inversion: N-terminal split GFP (eyes) |
| 204312 | pJAT23 | 9,512 | 1 | | Inversion: C-terminal split GFP (eyes) |

### Supplier note — pJAT55

`pJAT55` is **discontinued at Addgene** (record + sequence remain at
[Addgene 204302](https://www.addgene.org/204302/), but it cannot be ordered).
The live supplier is the **Drosophila Genomics Resource Center (DGRC)**:

- **RRID:** DGRC_1632 · **Price:** $63.00 · **Donor:** David Stern
- Page: https://dgrc.bio.indiana.edu/product/View?product=1632

## attR coordinates (1-based, all 124 bp)

Coordinates are 1-based inclusive; `(c)` denotes the reverse strand. Verified
from the GenBank files by `R/pjat_plasmids.R`.

| Plasmid | attR3 | attR4 | attR5 | attR6 |
|---|---|---|---|---|
| pJAT7 | 2477..2600 | c(4341..4464) | 6090..6213 | c(8894..9017) |
| pJAT10 | 2477..2600 | c(4341..4464) | 6326..6449 | c(9130..9253) |
| pJAT11 | 2761..2884 | c(4625..4748) | 7055..7178 | c(9859..9982) |
| pJAT17 | 5570..5693 | c(7434..7557) | 14817..14940 | c(2533..2656) |
| pJAT18 | 6621..6744 | c(8483..8606) | 780..903 | c(3584..3707) |
| pJAT21 | 6631..6754 | c(8495..8618) | 790..913 | c(3594..3717) |
| pJAT22 | 7421..7544 | c(298..421) | 1580..1703 | c(4384..4507) |
| pJAT23 | 2478..2601 | c(4342..4465) | 6149..6272 | c(8953..9076) |
| pJAT25 | 2477..2600 | c(4341..4464) | 7029..7152 | c(9833..9956) |
| pJAT27 | 2499..2622 | c(4363..4486) | 6805..6928 | c(9609..9732) |
| pJAT30 | 2761..2884 | c(4625..4748) | 6610..6733 | c(9414..9537) |
| pJAT31 | 5748..5871 | c(7612..7735) | 10042..10165 | c(1905..2028) |
| pJAT32 | 2761..2884 | c(4625..4748) | 4906..5029 | c(7710..7833) |
| pJAT39 | c(5390..5513) | 3526..3649 | c(1039..1162) | 9233..9356 |
| pJAT43 | c(6472..6595) | 4608..4731 | c(12313..12436) | 9509..9632 |
| pJAT55 | 2761..2884 | c(4625..4748) | 10871..10994 | c(13675..13798) |
| pJAT56 | 2476..2599 | c(4340..4463) | 10304..10427 | c(13108..13231) |
| pJAT57 | 2761..2884 | c(4625..4748) | 12236..12359 | c(15040..15163) |
| pJAT58 | 4954..5077 | c(6818..6941) | 12969..13092 | c(1917..2040) |
| pJAT62 | 2725..2848 | c(4589..4712) | 11547..11670 | c(14351..14474) |
| pJAT71 | 12379..12502 | c(14243..14366) | 6538..6661 | c(9342..9465) |
| pJAT74 | 3980..4103 | c(5844..5967) | 6125..6248 | c(943..1066) |
| pJAT75 | 4257..4380 | c(6121..6244) | 10878..11001 | c(1220..1343) |
| pJAT77 | c(1007..1130) | 6703..6826 | c(6422..6545) | 3618..3741 |

Two invariants: **pJAT74 and pJAT32 share identical attR3/attR4 coordinates**
(the two empty backbones differ only downstream), and three plasmids
(**pJAT39, pJAT43, pJAT77**) carry an inverted cassette — attR3/attR4 on
opposite strands from the rest. The parser records strand per site, so this is
never assumed away.

## Notes for the donor-fragment designer

- **attL sites are absent from these plasmids by design.** LR recombination
  consumes attL and leaves attB, so a finished pJAT plasmid contains attR, not
  attL. The four attL sites needed for the *synthesized donor fragments* come
  from the Gateway donor templates (protocols.io). They are **never invented** by
  the code: if not supplied, `pjat_assemble_fragments()` returns status
  `ARM_ONLY` and emits the arms with `[attL3]`-style placeholders.
- **Arms must be supplied in genomic orientation** (rule R4). The app records
  this as an explicit check rather than trusting the caller.
- **The spacer must be exactly 20 nt with the NGG PAM excluded** (rule R3);
  inclusion is a hard error.
- **250 bp is the validated minimum arm length** (rule R2).

## Parsing hazards (handled by `R/pjat_plasmids.R`)

1. **H1 — Geneious HTML-in-notes.** `pJAT55.gbff` embeds `<a href="urn:local:...">`
   fragments inside qualifier values, so `file(1)` reports it as *"HTML document
   text"*. Detection is by `LOCUS`/`FEATURES`/`ORIGIN` content, never by
   extension or MIME type. Its LOCUS name is also quoted with an embedded comma
   (`LOCUS pJAT55{...} 13951 bp`), so size is parsed by regex on `<digits> bp`
   rather than by token position.
2. **H2 — Custom att feature keys.** Addgene annotates attR as `protein_bind`
   with `/label="attR3"`; DGRC annotates it under a bare `attR3` key. The parser
   matches on the feature key *or* its label/qualifiers. Filtering on
   `protein_bind` alone silently drops all four sites from pJAT55.
- Sequence file: `pJAT55.gbff` (Geneious export, local URN accession
  `urn.local.sternd.1ws-gdrtrmo`)
- DGRC lists resistance as "Ampicillin, Chloramphenicol" — **incomplete** (arm
  markers only); the backbone `aadA` confers spectinomycin/streptomycin
  resistance and is present in the sequence.

Because pJAT55 encodes the **piggyBac transposase** needed to excise removable
markers in the scarless-edit workflow, it is the one member that cannot be
re-ordered from Addgene. Treat the DGRC glycerol stock as the archival copy.