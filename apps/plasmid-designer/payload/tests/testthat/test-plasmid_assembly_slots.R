# ============================================================================
# tests/testthat/test-plasmid_assembly_slots.R — multi-slot cloning engine
# ============================================================================
# Cloning math must be provably correct: length identity, strand-aware
# reverse-complement, circular wraparound, and exact Type IIS cut offsets.
# ============================================================================

library(testthat)

project_root <- normalizePath(
  file.path(dirname(normalizePath(testthat::test_path())), ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(project_root, "R", "plasmid_assembly.R"))) {
  project_root <- project_root
}
source(file.path(project_root, "R", "pjat_plasmids.R"))
source(file.path(project_root, "R", "plasmid_assembly.R"))
source(file.path(project_root, "R", "backbone_registry.R"))
source(file.path(project_root, "R", "plasmid_assembly_slots.R"))

# ---- fixtures ---------------------------------------------------------------

make_slot <- function(id, type, start, end, strand = "+") {
  backbone_slot(id, type, paste("slot", id), start, end, strand)
}

# ---- resolve_slot_window ----------------------------------------------------

test_that("resolve_slot_window reports a normal window", {
  w <- resolve_slot_window(make_slot("s", "donor_payload", 100L, 200L), 1000L,
                           "circular")
  expect_equal(w$start, 100L); expect_equal(w$end, 200L)
  expect_false(w$spans_origin)
})

test_that("resolve_slot_window detects a wraparound window on a circular plasmid", {
  w <- resolve_slot_window(make_slot("s", "donor_payload", 950L, 50L), 1000L,
                           "circular")
  expect_true(w$spans_origin)
})

test_that("resolve_slot_window rejects a wraparound window on a linear plasmid", {
  expect_error(
    resolve_slot_window(make_slot("s", "donor_payload", 950L, 50L), 1000L, "linear"),
    "spans the origin")
})

test_that("resolve_slot_window rejects out-of-range coordinates", {
  expect_error(resolve_slot_window(make_slot("s", "donor_payload", 900L, 1200L),
                                   1000L, "circular"), "outside the backbone")
  expect_error(resolve_slot_window(make_slot("s", "donor_payload", NA_integer_,
                                             NA_integer_), 1000L, "circular"),
               "no coordinates")
})

# ---- clone_into_slot: length identity and content ---------------------------

test_that("clone_into_slot replaces a window: final = backbone - window + insert", {
  backbone <- paste(rep("ACGT", 250L), collapse = "")   # 1000 bp
  slot <- make_slot("donor_payload_linker", "donor_payload", 101L, 200L)
  insert <- strrep("G", 300L)
  r <- clone_into_slot(backbone, slot, insert, orient = "+", topology = "circular")
  expect_equal(nchar(r$sequence), 1000L - 100L + 300L)
  expect_equal(r$insert_len, 300L)
  expect_equal(r$slot_id, "donor_payload_linker")
  expect_equal(r$slot_type, "donor_payload")
  expect_false(r$spans_origin)
  # the insert must be present verbatim at the reported coordinates
  expect_equal(substr(r$sequence, r$insert_start, r$insert_end), insert)
})

test_that("clone_into_slot reverse-complements a minus-orientation insert", {
  backbone <- paste(rep("ACGT", 250L), collapse = "")
  slot <- make_slot("s", "donor_payload", 101L, 200L)
  insert <- "AAAACCCCGGGGTTTT"
  plus <- clone_into_slot(backbone, slot, insert, orient = "+", topology = "circular")
  minus <- clone_into_slot(backbone, slot, insert, orient = "-", topology = "circular")
  expect_equal(substr(plus$sequence, plus$insert_start, plus$insert_end), insert)
  expect_equal(substr(minus$sequence, minus$insert_start, minus$insert_end),
               "AAAACCCCGGGGTTTT" |> chartr("ACGT", "TGCA", x = _) |> strsplit("") |>
                 (\(s) paste(rev(s[[1]]), collapse = ""))())
})

test_that("clone_into_slot insert-only (replace = FALSE) grows without deleting", {
  backbone <- paste(rep("ACGT", 250L), collapse = "")
  slot <- make_slot("s", "grna_array", 401L, 500L)
  insert <- strrep("T", 50L)
  r <- clone_into_slot(backbone, slot, insert, orient = "+",
                       topology = "circular", replace = FALSE)
  expect_equal(nchar(r$sequence), 1050L)
  expect_equal(r$insert_len, 50L)
})

test_that("clone_into_slot handles a circular wraparound window", {
  backbone <- paste(rep("ACGT", 250L), collapse = "")   # 1000 bp
  # window spans the origin: 951..1000 + 1..50  (100 bp total)
  slot <- make_slot("s", "marker_swap", 951L, 50L)
  insert <- strrep("A", 200L)
  r <- clone_into_slot(backbone, slot, insert, orient = "+", topology = "circular")
  expect_true(r$spans_origin)
  # replaced 100 bp with 200 bp -> 1100 bp, circular length preserved contract
  expect_equal(nchar(r$sequence), 1000L - 100L + 200L)
  expect_equal(r$insert_len, 200L)
  # the insert must appear exactly once and verbatim
  expect_equal(length(gregexpr(insert, r$sequence, fixed = TRUE)[[1L]]), 1L)
})

# ---- clone_at_restriction_site ---------------------------------------------

test_that("clone_at_restriction_site inserts at the real BbsI cut offset", {
  # BbsI: GAAGAC, cut_after = 2 -> cut_pos = site_pos + 6 + 2 - 1
  backbone <- paste0(strrep("A", 99L), "GAAGAC", strrep("C", 95L))  # 200 bp
  site_pos <- 100L
  insert <- "TTTTGGGG"
  r <- clone_at_restriction_site(backbone, "BbsI", site_pos, insert, orient = "+")
  expect_equal(r$enzyme, "BbsI")
  expect_equal(r$cut_pos, site_pos + 6L + 2L - 1L)
  expect_equal(nchar(r$sequence), 208L)          # insert-only
  expect_equal(substr(r$sequence, r$insert_start, r$insert_end), insert)
  # the recognition site is preserved (never auto-selected, never destroyed)
  expect_true(grepl("GAAGAC", r$sequence, fixed = TRUE))
})

test_that("clone_at_restriction_site rejects a site that is not there", {
  backbone <- strrep("A", 200L)
  expect_error(clone_at_restriction_site(backbone, "BbsI", 100L, "TTTT"),
               "no BbsI recognition site")
  expect_error(clone_at_restriction_site(backbone, "NotAnEnzyme", 100L, "TTTT"),
               "enzyme must be one of")
})

