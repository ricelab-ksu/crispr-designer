# Genome Inventory — *Drosophila* Reference Files

> Base path: `/Volumes/Backup_Plus/projects/bioinformatics/reference_files/drosophila_species_genomes`
> Inspected & verified on 2026-10-16.

## Summary by Species (7 Species Total)

| Species Folder | RefSeq Assembly / Build | Genomic FASTA (`genomic.fna`) | Annotation (`genomic.gtf`) | Transcript FASTA (`rna.fna`) | Gene Symbol / ID Lookup Tables | CRISPR Readiness |
|----------------|-------------------------|-------------------------------|----------------------------|------------------------------|--------------------------------|------------------|
| `d_melanogaster` | `GCF_000001215.4` (Release 6 plus ISO1 MT) | `GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna` (139 MB) + `.fai` | `genomic.gtf` (177 MB), `genomic_sorted.gtf` (+ `.idx`), `dmel-all-r6.68.gtf` (75 MB) | `rna.fna` (92 MB) + `.fai` | `dmel_gene_reference.csv` (3.4 MB), `fbgn_map.csv` (53 KB), `marker_database_from_gtf.csv` | **READY (Gold Standard)** |
| `d_simulans` | `GCF_016746395.2` (`Prin_Dsim_3.1`) | `GCF_016746395.2_Prin_Dsim_3.1_genomic.fna` (127 MB) | `genomic.gtf` (160 MB) | `rna.fna` (89 MB) + `.gz` | Extracted from `genomic.gtf` (`gene`, `locus_tag`, `db_xref`, `product`) | **READY** (needs `.fai` index on first run) |
| `d_eugracilis` | `GCF_018153835.1` (`ASM1815383v1`) | `GCF_018153835.1_ASM1815383v1_genomic.fna` (159 MB) | `genomic.gtf` (123 MB) | `rna.fna` (67 MB) + `.gz` | Extracted from `genomic.gtf` | **READY** (needs `.fai` index on first run) |
| `d_erecta` | `GCF_003286155.1` (`DereRS2`) | `genomic.fna.gz` (43 MB compressed) | `genomic.gtf.gz` (3.9 MB compressed) | `rna.fna.gz` (14 MB compressed) | Extracted from `genomic.gtf.gz` | **READY (Compressed)** — decompress or bgzip/index on first use |
| `d_teissieri` | `GCF_016746235.2` (`Prin_Dtei_1.1`) | `genomic.fna.gz` (44 MB compressed) | `genomic.gtf.gz` (4.2 MB compressed) | `rna.fna.gz` (12 MB compressed) | Extracted from `genomic.gtf.gz` | **READY (Compressed)** — decompress or bgzip/index on first use |
| `d_ananassae` | `GCF_017639315.1` (`ASM1763931v2`) | **MISSING locally** | **MISSING locally** | `GCF_017639315.1_ASM1763931v2_rna.fna` (89 MB) + `.gz` | Currently transcript headers only | **NEEDS GENOMIC FASTA + GTF** (required for promoters, introns, and 1–2 kb homology arms) |
| `d_virilis` | `GCF_003285735.1` (`DvirRS2`) | **MISSING locally** | **MISSING locally** | `GCF_003285735.1_DvirRS2_rna.fna` (76 MB) + `.gz` | Currently transcript headers only | **NEEDS GENOMIC FASTA + GTF** (required for promoters, introns, and 1–2 kb homology arms) |

## Key Notes for Subagents (T2 & T3)

1. **`dmel_gene_reference.csv` Schema (`d_melanogaster`):**
   Columns: `"search_name","gene_symbol","ncbi_gene_id","cg_number","flybase_id","ncbi_geneid","transcript_ids","chromosome","start","end","strand","biotype","synonyms"`.
   - Provides instant lookup by gene symbol (`PlexB`), FlyBase ID (`FBgn0025740`), CG number (`CG17245`), NCBI gene ID (`43766`), and semicolon-delimited synonyms.
2. **Non-melanogaster Species Gene Lookup & Lazy Download:**
   - As documented in `HCR_probe_design_v41.R`, non-melanogaster RefSeq GTFs index genes by `gene` / `locus_tag` (e.g., `LOC...`) and `product` (protein name). T3 must build or cache a lightweight `gene_reference.csv` per species from `genomic.gtf(.gz)` on first load so gene lookup is fast and uniform across all 7 species.
   - For `d_ananassae` and `d_virilis`, reuse the NCBI RefSeq download helper from `HCR_probe_design_v41.R` so the app can automatically fetch `genomic.fna.gz` and `genomic.gtf.gz` when those species are selected.
