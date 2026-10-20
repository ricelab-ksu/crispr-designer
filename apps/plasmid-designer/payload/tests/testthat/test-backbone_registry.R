# ============================================================================
# tests/testthat/test-backbone_registry.R — backbone_registry/1.0 contract
# ============================================================================
# The registry is the abstraction that lets Plasmid Designer accept ANY
# plasmid.  These tests pin the two things that make that safe:
#   1. the pJAT adapter is behavior-preserving (identical to pjat_load_catalog())
#   2. a non-pJAT backbone registers and gets sensible auto-discovered slots
# ============================================================================

source_modules <- function() {
  root <- project_root
  for (f in c("addgene_maps.R", "pjat_plasmids.R", "gateway_attl.R",
              "marker_sequences.R", "sample_inputs.R", "plasmid_map.R",
              "guide_bundle.R", "plasmid_assembly.R", "backbone_registry.R")) {
    p <- file.path(root, "R", f)
    if (file.exists(p)) suppressWarnings(source(p, local = FALSE))
  }
}
source_modules()

test_that("backbone_slot validates slot_type", {
  s <- backbone_slot("s1", "grna_array", "test", 10L, 20L, "+")
  expect_equal(nrow(s), 1L)
  expect_equal(s$slot_type, "grna_array")
  expect_error(backbone_slot("s2", "not_a_type"), "slot_type")
  expect_error(backbone_slot("", "grna_array"), "nzchar")
})

test_that("backbone_entry normalises sequence and enforces size identity", {
  e <- backbone_entry("t1", "test", sequence = "acgtn acgt")
  expect_equal(e$sequence, "ACGTNACGT")
  expect_equal(e$size_bp, 9L)            # size derived from cleaned sequence
  expect_error(backbone_entry("t2", "test", sequence = "ACGT", size_bp = 99L))
  expect_error(backbone_entry("t3", "test", sequence = ""))
})

test_that("backbone_registry rejects unnamed entries and duplicate ids", {
  e1 <- backbone_entry("a", "A", sequence = "ACGT")
  e2 <- backbone_entry("b", "B", sequence = "ACGT")
  reg <- backbone_registry(a = e1, b = e2)
  expect_s3_class(reg, "backbone_registry")
  expect_length(reg, 2L)
  expect_error(backbone_registry(e1, e2), "named")
  expect_error(backbone_registry(a = e1, a = e2), "duplicate")
})

test_that("backbone_get and backbone_slots resolve by id and by name", {
  e <- backbone_entry("pJATx", "pJAT-x", sequence = "ACGTACGT",
                      slots = backbone_slot("g1", "grna_array", "u6", 1L, 4L, "+"))
  reg <- backbone_registry(pJATx = e)
  expect_equal(backbone_get(reg, "pJATx")$name, "pJAT-x")
  expect_equal(backbone_get(reg, "pJAT-x")$id, "pJATx")
  expect_null(backbone_get(reg, "nope"))
  expect_equal(nrow(backbone_slots(reg, "pJATx", "grna_array")), 1L)
  expect_equal(nrow(backbone_slots(reg, "pJATx", "marker_swap")), 0L)
})

test_that("backbone_register_fasta builds a sequence-only entry", {
  fa <- tempfile(fileext = ".fasta")
  writeLines(c(">test plasmid", "ACGTACGTAC", "GTACGTACGT"), fa)
  e <- backbone_register_fasta(fa, id = "myvec", name = "My Vector")
  expect_equal(e$format, "fasta")
  expect_equal(e$size_bp, 20L)
  expect_equal(nrow(e$features), 0L)
  expect_equal(nrow(e$slots), 0L)
  unlink(fa)
})

test_that("backbone_discover_slots does not mistake T7/SP6/lac promoters for gRNA arrays", {
  # Precision guard: a broad 'promoter' match would sweep these in.
  feats <- data.frame(
    key = c("promoter", "promoter", "promoter", "CDS"),
    start = c(10L, 100L, 200L, 300L),
    end = c(40L, 130L, 230L, 600L),
    strand = "+", label = c("T7 promoter", "SP6 promoter", "lac UV5 promoter", "AmpR"),
    product = "", qualifiers = "", stringsAsFactors = FALSE
  )
  seq <- strrep("ACGT", 200L)
  e <- backbone_entry("x", "x", sequence = seq,
                      slots = backbone_discover_slots(feats, seq, "circular"))
  reg <- backbone_registry(x = e)
  expect_equal(nrow(backbone_slots(reg, "x", "grna_array")), 0L)
  expect_equal(nrow(backbone_slots(reg, "x", "marker_swap")), 1L)  # AmpR is real
})

test_that("backbone_discover_slots finds a dU6-3 promoter as a gRNA array", {
  feats <- data.frame(
    key = "promoter", start = 100L, end = 520L, strand = "+",
    label = "dU6-3 promoter", product = "", qualifiers = "",
    stringsAsFactors = FALSE
  )
  slots <- backbone_discover_slots(feats, strrep("ACGT", 200L), "circular")
  expect_equal(sum(slots$slot_type == "grna_array"), 1L)
})

test_that("pJAT adapter reproduces the catalog exactly (behavior-preserving)", {
  root <- project_root
  yml <- file.path(root, "donor_templates", "pjat_plasmids.yml")
  ref_dir <- pjats_dir
  skip_if_not(file.exists(yml), "pJAT catalog YAML not available")
  skip_if_not(dir.exists(ref_dir), "pJAT GenBank reference dir not mounted")

  cat_df <- pjat_load_catalog(yml = yml)
  reg <- backbone_pjat_adapter(yml = yml, ref_dir = ref_dir)

  # same set of plasmids, no more, no fewer
  expect_setequal(names(reg), cat_df$name)
  expect_equal(length(reg), nrow(cat_df))

  # size and every attR coordinate must match the catalog byte-for-byte
  for (i in seq_len(nrow(cat_df))) {
    e <- reg[[cat_df$name[i]]]
    expect_equal(e$size_bp, cat_df$size_bp[i])
    att <- find_pjat_recombination_sites(list(features = e$features,
                                              sequence = e$sequence))
    for (k in 3:6) {
      nm <- paste0("attR", k)
      r <- att[att$att == nm, , drop = FALSE]
      expect_equal(as.integer(r$start[1]), cat_df[[paste0(nm, "_start")]][i])
      expect_equal(as.integer(r$end[1]),   cat_df[[paste0(nm, "_end")]][i])
    }
  }
})

test_that("pJAT adapter discovers one gRNA slot per catalog U6 promoter", {
  root <- project_root
  yml <- file.path(root, "donor_templates", "pjat_plasmids.yml")
  ref_dir <- pjats_dir
  skip_if_not(file.exists(yml) && dir.exists(ref_dir), "pJAT assets not mounted")

  cat_df <- pjat_load_catalog(yml = yml)
  reg <- backbone_pjat_adapter(yml = yml, ref_dir = ref_dir)
  for (i in seq_len(nrow(cat_df))) {
    n_grna <- nrow(backbone_slots(reg, cat_df$name[i], "grna_array"))
    expect_equal(n_grna, cat_df$u6_promoters[i],
                 info = paste("gRNA slots for", cat_df$name[i]))
  }
})

test_that("every pJAT backbone exposes a donor_payload slot (inverted cassettes included)", {
  root <- project_root
  yml <- file.path(root, "donor_templates", "pjat_plasmids.yml")
  ref_dir <- pjats_dir
  skip_if_not(file.exists(yml) && dir.exists(ref_dir), "pJAT assets not mounted")

  reg <- backbone_pjat_adapter(yml = yml, ref_dir = ref_dir)
  for (nm in names(reg)) {
    s <- backbone_slots(reg, nm, "donor_payload")
    expect_gte(nrow(s), 1L)
    # a discovered window must be a real, in-range interval
    expect_true(all(s$start >= 1L & s$end <= reg[[nm]]$size_bp & s$start <= s$end),
                info = nm)
  }
})

