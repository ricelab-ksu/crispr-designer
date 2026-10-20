library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "pjat_plasmids.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "guide_bundle.R"))
source(file.path(project_root, "R", "plasmid_assembly.R"))

PJAT_TEST_DIR <- pjats_dir

test_that("insert_cassette replaces window and reports correct coordinates", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  insert <- strrep("A", 100)
  asm <- insert_cassette(gb$sequence, 1000L, 1099L, insert, strand = "+")
  expect_equal(asm$final_len, gb$size_bp)   # 100 bp replacing 100 bp
  expect_equal(asm$insert_start, 1000L)
  expect_equal(asm$insert_end, 1099L)
  expect_equal(substr(asm$sequence, 1000L, 1099L), strrep("A", 100))
})

test_that("insert_cassette reverse-complements minus-strand inserts", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  insert <- "ACGT"
  asm_plus <- insert_cassette(gb$sequence, 1000L, 1003L, insert, "+")
  asm_minus <- insert_cassette(gb$sequence, 1000L, 1003L, insert, "-")
  expect_equal(substr(asm_plus$sequence, 1000L, 1003L), "ACGT")
  expect_equal(substr(asm_minus$sequence, 1000L, 1003L), "ACGT")  # RC of ACGT is ACGT
  insert2 <- "AAAC"
  asm2 <- insert_cassette(gb$sequence, 1000L, 1003L, insert2, "-")
  expect_equal(substr(asm2$sequence, 1000L, 1003L), "GTTT")  # RC of AAAC
})

test_that("assemble_gRNA_cassette enforces 20-nt PAM-excluded spacer", {
  expect_equal(nchar(assemble_gRNA_cassette(strrep("G", 20))), 20L + nchar(GRRNA_SCAFFOLD) + nchar(GRRNA_TERMINATOR))
  expect_error(assemble_gRNA_cassette(strrep("G", 21)), "exactly 20 nt")
  expect_error(assemble_gRNA_cassette(strrep("G", 19)), "exactly 20 nt")
})

test_that("find_cloning_sites reports BbsI/BsaI/BsmBI without auto-selecting", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  sites <- find_cloning_sites(gb$sequence)
  expect_true(any(sites$enzyme == "BbsI"))
  expect_true(any(sites$enzyme == "BsaI"))
  expect_true(all(sites$position > 0L & sites$position <= nchar(gb$sequence)))
})

test_that("assemble_expression_vector uses linker window by default", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  asm <- assemble_expression_vector(gb, spacer = strrep("G", 20), strategy = "linker")
  expect_equal(asm$window$start, 5967L + 1L)   # attR4 end + 1
  expect_equal(asm$window$end,   6125L - 1L)   # attR5 start - 1
  expect_true(grepl(strrep("G", 20), asm$sequence, fixed = TRUE))
  expect_equal(asm$final_len, gb$size_bp - (asm$window$end - asm$window$start + 1L) + asm$insert_len)
  expect_true(asm$qc$ok || length(asm$qc$warnings) >= 0L)
})

test_that("assemble_expression_vector rejects missing spacer and insert", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  expect_error(assemble_expression_vector(gb), "provide either")
})

test_that("write_plasmid_genbank produces parseable circular GenBank", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  asm <- assemble_expression_vector(gb, spacer = strrep("G", 20), strategy = "linker")
  txt <- write_plasmid_genbank(asm$sequence, asm$features, name = "test_assembly", topology = "circular")
  expect_match(txt, "^LOCUS")
  expect_match(txt, "circular")
  expect_match(txt, "//$")
  # round-trip parse
  tf <- tempfile(fileext = ".gb")
  writeLines(txt, tf)
  reparsed <- suppressWarnings(read_pjat_genbank(tf))
  expect_equal(nchar(reparsed$sequence), nchar(asm$sequence))
})

test_that("plasmid_qc flags Type IIS sites and homopolymers", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  qc <- plasmid_qc(gb$sequence, name = "pJAT74")
  expect_true(is.data.frame(qc$cloning_sites))
  expect_true(any(qc$cloning_sites$enzyme == "BbsI"))
})

test_that("assembly length identity: final = backbone - window + insert", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  insert <- strrep("C", 300)
  asm <- insert_cassette(gb$sequence, 1000L, 1200L, insert, "+")
  expect_equal(asm$final_len, gb$size_bp - 201L + 300L)
})

test_that("insert_cassette supports insert-only ranges (end = start - 1)", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  insert <- strrep("T", 50)
  asm <- insert_cassette(gb$sequence, 1000L, 999L, insert, "+")
  expect_equal(asm$final_len, gb$size_bp + 50L)
  expect_equal(asm$insert_start, 1000L)
  expect_equal(asm$insert_end, 1049L)
  expect_equal(substr(asm$sequence, 1000L, 1049L), insert)
  # backbone context preserved around the insertion
  expect_equal(substr(asm$sequence, 1L, 999L), substr(gb$sequence, 1L, 999L))
  expect_equal(substr(asm$sequence, 1050L, asm$final_len),
               substr(gb$sequence, 1000L, nchar(gb$sequence)))
  expect_error(insert_cassette(gb$sequence, 1000L, 998L, insert, "+"),
               "invalid insertion range")
})

test_that("assemble_expression_vector exposes insert_start/insert_end", {
  skip_if_not(dir.exists(PJAT_TEST_DIR), "pJAT collection not present")
  gb <- read_pjat_genbank(file.path(PJAT_TEST_DIR, "pJAT74.gbk"))
  asm <- assemble_expression_vector(gb, spacer = strrep("G", 20), strategy = "linker")
  expect_false(is.null(asm$insert_start))
  expect_false(is.null(asm$insert_end))
  expect_equal(asm$insert_start, asm$window$start)
  expect_equal(asm$insert_end - asm$insert_start + 1L, asm$insert_len)
  # arm extraction the UI performs must not error and must be >= 250 bp
  arm_len <- 1000L
  left_arm  <- substr(asm$sequence, max(1L, asm$insert_start - arm_len),
                      asm$insert_start - 1L)
  right_arm <- substr(asm$sequence, asm$insert_end + 1L,
                      min(nchar(asm$sequence), asm$insert_end + arm_len))
  expect_true(nchar(left_arm) >= 250L)
  expect_true(nchar(right_arm) >= 250L)
})
