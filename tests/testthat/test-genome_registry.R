library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "genome_registry.R"))) {
  project_root <- "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design"
}
source(file.path(project_root, "R", "genome_registry.R"))

test_that("canonical_species normalizes prefixes and common misspellings", {
  expect_equal(canonical_species("melanogaster"), "melanogaster")
  expect_equal(canonical_species("d_melanogaster"), "melanogaster")
  expect_equal(canonical_species("Drosophila melanogaster"), "melanogaster")
  expect_equal(canonical_species("D. simulans"), "simulans")
  expect_equal(canonical_species("tesseri"), "teissieri")
  expect_equal(canonical_species("d_tesseri"), "teissieri")
  expect_equal(canonical_species("dmel"), "melanogaster")
  expect_equal(canonical_species("dvir"), "virilis")
  expect_true(is.na(canonical_species(NULL)))
  expect_true(is.na(canonical_species(NA)))
})

test_that("get_reference_dir respects environment overrides and defaults cleanly", {
  old_cas9 <- Sys.getenv("CAS9_REFERENCE_DIR", unset = "")
  old_hcr  <- Sys.getenv("HCR_REFERENCE_DIR", unset = "")
  on.exit({
    if (nzchar(old_cas9)) Sys.setenv(CAS9_REFERENCE_DIR = old_cas9) else Sys.unsetenv("CAS9_REFERENCE_DIR")
    if (nzchar(old_hcr))  Sys.setenv(HCR_REFERENCE_DIR = old_hcr)   else Sys.unsetenv("HCR_REFERENCE_DIR")
  }, add = TRUE)

  Sys.unsetenv("CAS9_REFERENCE_DIR")
  Sys.unsetenv("HCR_REFERENCE_DIR")
  expect_equal(
    get_reference_dir(),
    "/Volumes/Backup_Plus/projects/bioinformatics/reference_files/drosophila_species_genomes"
  )

  Sys.setenv(HCR_REFERENCE_DIR = "/tmp/hcr_test_ref")
  expect_equal(get_reference_dir(), "/tmp/hcr_test_ref")

  Sys.setenv(CAS9_REFERENCE_DIR = "/tmp/cas9_test_ref")
  expect_equal(get_reference_dir(), "/tmp/cas9_test_ref")
})

test_that("get_species_registry inventories all 7 Drosophila species accurately", {
  reg <- get_species_registry()
  expect_s3_class(reg, "data.frame")
  expected_species <- c(
    "melanogaster", "simulans", "eugracilis",
    "erecta", "teissieri", "ananassae", "virilis"
  )
  expect_true(all(expected_species %in% reg$species_key))
  expect_equal(nrow(reg), 7L)

  non_na_paths <- c(
    stats::na.omit(reg$genomic_fna_path),
    stats::na.omit(reg$genomic_gtf_path),
    stats::na.omit(reg$rna_fna_path)
  )
  expect_false(any(grepl("^\\._", basename(non_na_paths))))

  for (sp in c("melanogaster", "simulans", "eugracilis")) {
    row <- reg[sp, ]
    expect_equal(row$status, "ready", info = sp)
    expect_true(file.exists(row$genomic_fna_path), info = sp)
    expect_true(file.exists(row$genomic_gtf_path), info = sp)
    expect_true(file.exists(row$rna_fna_path), info = sp)
  }

  for (sp in c("erecta", "teissieri")) {
    row <- reg[sp, ]
    expect_equal(row$status, "compressed", info = sp)
    expect_true(file.exists(row$genomic_fna_path), info = sp)
    expect_true(file.exists(row$genomic_gtf_path), info = sp)
    expect_true(file.exists(row$rna_fna_path), info = sp)
  }

  for (sp in c("ananassae", "virilis")) {
    row <- reg[sp, ]
    expect_true(file.exists(row$rna_fna_path), info = sp)
    if (is.na(row$genomic_fna_path) || is.na(row$genomic_gtf_path)) {
      expect_equal(row$status, "needs_download", info = sp)
    }
  }
})


test_that("chromosome alias mapping resolves X, 2L, 2R, 3L, 3R, 4, Y, M and accessions", {
  cmap <- get_chromosome_map("melanogaster")
  expect_equal(resolve_chromosome("melanogaster", "X"), "NC_004354.4")
  expect_equal(resolve_chromosome("melanogaster", "chrX"), "NC_004354.4")
  expect_equal(resolve_chromosome("melanogaster", "2L"), "NT_033779.5")
  expect_equal(resolve_chromosome("melanogaster", "chr2L"), "NT_033779.5")
  expect_equal(resolve_chromosome("melanogaster", "2R"), "NT_033778.4")
  expect_equal(resolve_chromosome("melanogaster", "3L"), "NT_037436.4")
  expect_equal(resolve_chromosome("melanogaster", "3R"), "NT_033777.3")
  expect_equal(resolve_chromosome("melanogaster", "chr3R"), "NT_033777.3")
  expect_equal(resolve_chromosome("melanogaster", "4"), "NC_004353.4")
  expect_equal(resolve_chromosome("melanogaster", "Y"), "NC_024512.1")
  expect_equal(resolve_chromosome("melanogaster", "M"), "NC_024511.2")
  expect_equal(resolve_chromosome("melanogaster", "MT"), "NC_024511.2")
  expect_equal(resolve_chromosome("melanogaster", "mitochondrion"), "NC_024511.2")

  expect_equal(resolve_chromosome("melanogaster", "NT_033777.3"), "NT_033777.3")
  expect_equal(resolve_chromosome("melanogaster", "NC_004354.4"), "NC_004354.4")
  expect_equal(cmap$accession_to_alias[["NT_033777.3"]], "3R")
  expect_equal(cmap$accession_to_alias[["NC_004354.4"]], "X")

  expect_equal(resolve_chromosome("simulans", "2L"), "NC_052520.2")
  expect_equal(resolve_chromosome("simulans", "3R"), "NC_052523.2")
  expect_equal(resolve_chromosome("simulans", "X"), "NC_052525.2")
})

test_that("ensure_fai_index builds SAMtools-identical .fai index and decompress_gzip_file works", {
  tmp_dir <- tempfile("fai_test_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fa_path <- file.path(tmp_dir, "test_genome.fna")
  seq1 <- paste(rep("ACGTACGTAA", 25), collapse = "")
  seq2 <- paste(rep("TTGGAACCGG", 18), collapse = "")
  lines <- c(
    ">chr2L synthetic chromosome 2L",
    substring(seq1, seq(1, 201, by = 50), seq(50, 250, by = 50)),
    ">chr3R synthetic chromosome 3R",
    substring(seq2, seq(1, 151, by = 50), c(50, 100, 150, 180))
  )
  writeLines(lines, fa_path)

  fai_path <- ensure_fai_index(fa_path, force = TRUE)
  expect_true(file.exists(fai_path))
  fai_df <- read_fai_index(fai_path)
  expect_equal(nrow(fai_df), 2L)
  expect_equal(fai_df$name, c("chr2L", "chr3R"))
  expect_equal(fai_df$length, c(250, 180))
  expect_equal(fai_df$linebases, c(50L, 50L))
  expect_equal(fai_df$linewidth, c(51L, 51L))

  dmel_fna <- find_genome_file("melanogaster")
  dmel_ref_fai <- read_fai_index(paste0(dmel_fna, ".fai"))
  tmp_dmel_fai <- file.path(tmp_dir, "dmel_rebuilt.fai")
  ensure_fai_index(dmel_fna, fai_path = tmp_dmel_fai, force = TRUE)
  rebuilt_df <- read_fai_index(tmp_dmel_fai)
  expect_equal(rebuilt_df$name, dmel_ref_fai$name)
  expect_equal(rebuilt_df$length, dmel_ref_fai$length)
  expect_equal(rebuilt_df$offset, dmel_ref_fai$offset)
  expect_equal(rebuilt_df$linebases, dmel_ref_fai$linebases)
  expect_equal(rebuilt_df$linewidth, dmel_ref_fai$linewidth)

  gz_path <- file.path(tmp_dir, "compressed.fna.gz")
  gz_con <- gzfile(gz_path, "wb")
  writeLines(lines, gz_con)
  close(gz_con)
  decomp_path <- decompress_gzip_file(gz_path, overwrite = TRUE)
  expect_equal(readLines(decomp_path), lines)
})

test_that("extract_genomic_subseq matches Biostrings exactly on + and - strands", {
  sub_plus_alias <- extract_genomic_subseq("melanogaster", "3R", 1000, 2250, strand = "+")
  sub_plus_acc   <- extract_genomic_subseq("melanogaster", "NT_033777.3", 1000, 2250, strand = "+")
  expect_equal(sub_plus_alias, sub_plus_acc)
  expect_equal(nchar(sub_plus_alias), 1251L)

  sub_minus <- extract_genomic_subseq("melanogaster", "3R", 1000, 2250, strand = "-")
  expect_equal(sub_minus, revcomp_dna(sub_plus_alias))
  expect_equal(
    sub_minus,
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(sub_plus_alias)))
  )

  dmel_fna <- find_genome_file("melanogaster")
  con <- file(dmel_fna, "rt")
  first_lines <- readLines(con, n = 100)
  close(con)
  expect_match(first_lines[1], "^>NC_004354\\.4")
  first_seq <- paste(first_lines[2:50], collapse = "")
  s_coord <- 75L
  e_coord <- 1850L
  expected_plus <- toupper(substr(first_seq, s_coord, e_coord))
  expected_minus <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(expected_plus)))

  expect_equal(extract_genomic_subseq("melanogaster", "X", s_coord, e_coord, "+"), expected_plus)
  expect_equal(extract_genomic_subseq("melanogaster", "NC_004354.4", s_coord, e_coord, "-"), expected_minus)

  expect_equal(get_chr_length("melanogaster", "X"), 23542271L)
  expect_equal(get_chr_length("melanogaster", "NC_004354.4"), 23542271L)

  sim_sub <- extract_genomic_subseq("simulans", "2L", 500, 650, strand = "+")
  expect_equal(nchar(sim_sub), 151L)
})

test_that("download_species_ref constructs valid NCBI RefSeq URLs and performs atomic rename", {
  tmp_ref <- tempfile("ref_dl_")
  dir.create(tmp_ref, recursive = TRUE)
  on.exit(unlink(tmp_ref, recursive = TRUE), add = TRUE)

  captured_urls <- character(0)
  mock_dl <- function(url, destfile, mode = "wb", quiet = TRUE) {
    captured_urls <<- c(captured_urls, url)
    writeBin(charToRaw(">mock_seq\nACGT\n"), destfile)
    0L
  }

  msg <- download_species_ref(
    "d_ananassae",
    with_genome = TRUE,
    ref_dir = tmp_ref,
    progress = function(...) {},
    download_fn = mock_dl
  )
  expect_match(msg, "ananassae reference saved")
  expect_length(captured_urls, 4L)
  expect_true(all(grepl("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/", captured_urls, fixed = TRUE)))
  expect_true(file.exists(file.path(tmp_ref, "d_ananassae", "rna.fna.gz")))
  expect_true(file.exists(file.path(tmp_ref, "d_ananassae", "genomic.fna.gz")))
  expect_true(file.exists(file.path(tmp_ref, "d_ananassae", "genomic.gtf.gz")))
})

test_that("ensure_uncompressed_or_indexed handles ready, compressed, and needs_download species", {
  res_mel <- ensure_uncompressed_or_indexed("melanogaster")
  expect_equal(res_mel$status, "ready")
  expect_false(res_mel$is_compressed)
  expect_true(file.exists(res_mel$fai_path))

  res_sim <- ensure_uncompressed_or_indexed("simulans")
  expect_equal(res_sim$status, "ready")
  expect_true(file.exists(res_sim$fai_path))

  res_eug <- ensure_uncompressed_or_indexed("eugracilis")
  expect_equal(res_eug$status, "ready")
  expect_true(file.exists(res_eug$fai_path))

  res_ere <- ensure_uncompressed_or_indexed("erecta", decompress_gz = FALSE)
  expect_equal(res_ere$status, "compressed")
  expect_true(res_ere$is_compressed)
})
