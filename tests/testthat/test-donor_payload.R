# ============================================================================
# tests/testthat/test-donor_payload.R — Tests for donor_payload.R module
# ============================================================================

library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "pjat_plasmids.R"))) {
  project_root <- "/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design"
}

source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "donor_payload.R"))

test_that("read_donor_payload reads FASTA and raw text", {
  # Test raw string input
  txt <- ">mock_record\nACGTACGT\nTTTTGGGG"
  seq <- read_donor_payload(txt, is_path = FALSE)
  expect_equal(seq, paste0("ACGTACGT", "TTTTGGGG"))  # "ACGTACGT" + "TTTTGGGG"
  expect_equal(seq, paste0("ACGTACGT", "TTTTGGGG"))

  # Test file input
  tmp <- tempfile(fileext = ".fa")
  writeLines(c(">test", "ACGT", "TTTT"), con = tmp)
  on.exit(unlink(tmp), add = TRUE)
  seq2 <- read_donor_payload(tmp, is_path = TRUE)
  expect_equal(seq2, paste0("ACGT", "TTTT"))
})

test_that("read_donor_payload rejects empty sequence", {
  expect_error(read_donor_payload("   ", is_path = FALSE), "empty")
  tmp <- tempfile(fileext = ".fa")
  writeLines(">empty", con = tmp)
  on.exit(unlink(tmp), add = TRUE)
  expect_error(read_donor_payload(tmp, is_path = TRUE), "empty")
})

test_that("donor_payload_qc returns correct QC metrics", {
  qc <- donor_payload_qc("ACGTACGT", name = "test")
  expect_equal(qc$len, 8)
  expect_equal(qc$gc_pct, 50)
  expect_equal(qc$max_homopolymer, 1)
  expect_equal(qc$ok, TRUE)
  expect_true(is.data.frame(qc$cloning_sites))
})

test_that("donor_payload_qc flags homopolymers", {
  qc <- donor_payload_qc(paste(rep("A", 11), collapse = ""), name = "test")
  expect_false(qc$ok)
  expect_gt(qc$max_homopolymer, 10)
})

test_that("donor_payload_qc flags internal BbsI/BsaI/BsmBI sites", {
  # BbsI site GAAGAC
  qc <- donor_payload_qc("GATAGAAGACGT", name = "test")
  expect_true(nrow(qc$cloning_sites) > 0)
  expect_true("BbsI" %in% qc$cloning_sites$enzyme)
})

test_that("assemble_donor_construct builds 5'arm + payload + 3'arm", {
  arm <- list(arm_5p_seq = "ACGT", arm_3p_seq = "TGCA")
  dc <- assemble_donor_construct(arm, "AAAA", orient = "+")
  expect_equal(dc$sequence, "ACGTAAAATGCA")
  expect_equal(dc$payload_start, 5)
  expect_equal(dc$payload_end, 8)
  expect_equal(dc$left_arm, "ACGT")
  expect_equal(dc$right_arm, "TGCA")
  expect_equal(dc$sequence_provenance, "user-supplied (paste/upload)")
  expect_equal(dc$seq_rule_r6_share, 25)
})

test_that("assemble_donor_construct reverse-complements minus strand", {
  arm <- list(arm_5p_seq = "ACGT", arm_3p_seq = "TGCA")
  dc <- assemble_donor_construct(arm, "AAAA", orient = "-")
  expect_equal(dc$sequence, "ACGTTTTTTGCA")
  expect_equal(dc$payload, "TTTT")  # reversed
  expect_equal(dc$orient, "-")
})

test_that("assemble_donor_construct handles empty payload", {
  arm <- list(arm_5p_seq = "ACGT", arm_3p_seq = "TGCA")
  expect_error(assemble_donor_construct(arm, "   ", orient = "+"), "empty")
})

test_that("assemble_donor_construct rejects invalid orient", {
  arm <- list(arm_5p_seq = "ACGT", arm_3p_seq = "TGCA")
  expect_error(assemble_donor_construct(arm, "AAAA", orient = "x"),
               "'arg' should be one of")
})

test_that("write_donor_fasta produces correct FASTA", {
  arm <- list(arm_5p_seq = "ACGT", arm_3p_seq = "TGCA")
  dc <- assemble_donor_construct(arm, "AAAA", orient = "+")
  fa <- write_donor_fasta(dc, name = "test")
  expect_true(grepl("^>", fa))
  expect_true(grepl("test", fa))
  expect_true(grepl("ACGTAAAATGCA", fa))
})

test_that("write_donor_fasta includes arm/payload lengths", {
  arm <- list(arm_5p_seq = "ACGT", arm_3p_seq = "TGCA")
  dc <- assemble_donor_construct(arm, "AAAA", orient = "+")
  fa <- write_donor_fasta(dc, name = "test")
  expect_true(grepl("len=12bp", fa))
  expect_true(grepl("arm=5p:4bp\\+payload:4bp\\+3p:4bp", fa))
})

test_that("write_donor_genbank produces valid GenBank with features", {
  arm <- list(arm_5p_seq = "ACGT", arm_3p_seq = "TGCA")
  dc <- assemble_donor_construct(arm, "AAAA", orient = "+")
  gb <- write_donor_genbank(dc, name = "test")
  expect_true(grepl("^LOCUS", gb))
  expect_true(grepl("donor_payload", gb))
  expect_true(grepl("5p_homology_arm", gb))
  expect_true(grepl("3p_homology_arm", gb))
  expect_true(grepl("1..4", gb))       # 5' arm
  expect_true(grepl("5..8", gb))       # payload
  expect_true(grepl("9..12", gb))      # 3' arm
})

test_that("R6 overlap detection identifies payload sharing identity with arms", {
  arm <- list(arm_5p_seq = "ACGTACGTACGTACGTACGT", arm_3p_seq = "AAAA")
  dc <- assemble_donor_construct(arm, "ACGTACGTACGTACGTACGT", orient = "+")
  expect_equal(dc$seq_rule_r6_share, 100)

  arm2 <- list(arm_5p_seq = "AAAAAAAAAAAAAAAAAAAA", arm_3p_seq = "TTTTTTTTTTTTTTTTTTTT")
  dc2 <- assemble_donor_construct(arm2, "AAAAACGTAAAAACGTAAAA", orient = "+")
  # payload shares >50% with 5' arm (8/20 match at start)
  expect_gt(dc2$seq_rule_r6_share, 50)
})

test_that("donor_payload_qc gc_pct computed correctly", {
  # 100% GC
  qc1 <- donor_payload_qc(paste(rep("G", 10), collapse = ""))
  expect_equal(qc1$gc_pct, 100)
  expect_false(qc1$ok)  # GC > 70%

  # 0% GC
  qc2 <- donor_payload_qc(paste(rep("A", 10), collapse = ""))
  expect_equal(qc2$gc_pct, 0)
  expect_false(qc2$ok)  # GC < 30%

  # 50% GC with diverse composition (not low-complexity)
  # Use an acronym-like sequence with exactly 50% GC and no runs > 1
  s <- "ATCGATCGCA"  # 10 chars, 50% GC (C,G,C,G,C), no homopolymer runs
  qc3 <- donor_payload_qc(s, name = "test")
  # GC should be 50% (5/10)
  expect_equal(qc3$gc_pct, 50, tolerance = 0.01)
  # No homopolymer runs
  expect_equal(qc3$max_homopolymer, 1)
})
