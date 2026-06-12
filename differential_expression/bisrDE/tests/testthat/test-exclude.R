# v1.5.0 — sample/group/contrast exclusion mechanisms.

base_ss <- function() {
  data.frame(
    SampleID = paste0("S", 1:8),
    GroupID  = c("A", "A", "B", "B", "C", "C", "D", "D"),
    a_vs_b   = c(1, 1, 0, 0, NA, NA, NA, NA),
    c_vs_d   = c(NA, NA, NA, NA, 1, 1, 0, 0),
    stringsAsFactors = FALSE
  )
}

test_that("filter_samplesheet drops samples by SampleID", {
  ss <- filter_samplesheet(base_ss(), exclude_samples = c("S1", "S5"))
  expect_false(any(c("S1", "S5") %in% ss$SampleID))
  expect_equal(nrow(ss), 6)
})

test_that("filter_samplesheet drops samples by GroupID", {
  ss <- filter_samplesheet(base_ss(), exclude_groups = "C")
  expect_false("C" %in% ss$GroupID)
  expect_equal(nrow(ss), 6)
})

test_that("filter_samplesheet honours a declarative Exclude column", {
  ss0 <- base_ss()
  ss0$Exclude <- c(1, 0, 0, 0, 0, 0, 0, 0)
  ss <- filter_samplesheet(ss0)
  expect_false("S1" %in% ss$SampleID)
  expect_equal(nrow(ss), 7)
})

test_that("filter_samplesheet Exclude accepts TRUE/yes spellings", {
  ss0 <- base_ss()
  ss0$Exclude <- c("yes", "", "TRUE", "", "", "", "", "")
  ss <- filter_samplesheet(ss0)
  expect_false(any(c("S1", "S3") %in% ss$SampleID))
  expect_equal(nrow(ss), 6)
})

test_that("filter_samplesheet Exclude does NOT drop falsy values", {
  ss0 <- base_ss()
  # '0', '', 'no', 'false', NA must all be treated as keep.
  ss0$Exclude <- c("0", "", "no", "false", NA, "FALSE", "n", "0")
  ss <- filter_samplesheet(ss0)
  expect_equal(nrow(ss), 8)  # nothing dropped
})

test_that("filter_samplesheet unions all three exclusion sources", {
  ss0 <- base_ss()
  ss0$Exclude <- c(1, 0, 0, 0, 0, 0, 0, 0)          # S1
  ss <- filter_samplesheet(ss0,
                           exclude_samples = "S2",    # S2
                           exclude_groups  = "D")     # S7, S8
  expect_false(any(c("S1", "S2", "S7", "S8") %in% ss$SampleID))
  expect_equal(nrow(ss), 4)
})

test_that("filter_samplesheet errors if <2 samples remain", {
  # Drop 7 of 8 samples -> only S8 remains (< 2) -> abort.
  expect_error(
    filter_samplesheet(base_ss(), exclude_samples = paste0("S", 1:7)),
    "remain"
  )
})

test_that("filter_samplesheet is a no-op for unknown IDs (nothing dropped)", {
  # cli alerts are not R warning conditions, so we assert on the effect
  # (no rows dropped) rather than expect_warning().
  ss <- filter_samplesheet(base_ss(), exclude_samples = "NOPE")
  expect_equal(nrow(ss), 8)  # nothing actually dropped
})

test_that("parse_contrasts include_contrasts is an allowlist", {
  comparisons <- parse_contrasts(base_ss(), include_contrasts = "a_vs_b")
  expect_length(comparisons, 1)
  expect_equal(comparisons[[1]]$name, "a_vs_b")
})

test_that("parse_contrasts exclude_contrasts is a denylist", {
  comparisons <- parse_contrasts(base_ss(), exclude_contrasts = "a_vs_b")
  expect_length(comparisons, 1)
  expect_equal(comparisons[[1]]$name, "c_vs_d")
})

test_that("parse_contrasts skips Exclude/DisplayName metadata columns", {
  ss <- base_ss()
  # Insert metadata columns between GroupID and the contrasts.
  ss <- data.frame(
    SampleID    = ss$SampleID,
    GroupID     = ss$GroupID,
    Exclude     = 0,
    DisplayName = paste0("disp", 1:8),
    a_vs_b      = ss$a_vs_b,
    c_vs_d      = ss$c_vs_d,
    stringsAsFactors = FALSE
  )
  comparisons <- parse_contrasts(ss)
  nm <- vapply(comparisons, function(x) x$name, character(1))
  expect_setequal(nm, c("a_vs_b", "c_vs_d"))
  expect_false(any(c("Exclude", "DisplayName") %in% nm))
})

test_that("parse_contrasts handles a non-numeric contrast column gracefully", {
  ss <- data.frame(
    SampleID = paste0("S", 1:4),
    GroupID  = c("A", "A", "B", "B"),
    bad      = c("yes", "yes", "no", "no"),  # mis-encoded (not 0/1)
    good     = c(1, 1, 0, 0),
    stringsAsFactors = FALSE
  )
  # Should not error; the mis-encoded column is skipped (no exp/ctrl), the
  # valid one is kept.
  comparisons <- parse_contrasts(ss)
  nm <- vapply(comparisons, function(x) x$name, character(1))
  expect_equal(nm, "good")
})

test_that("parse_contrasts still respects first_contrast_col override", {
  ss <- data.frame(
    SampleID = paste0("S", 1:4),
    Tissue   = c("Liver", "Liver", "Lung", "Lung"),
    GroupID  = c("A", "A", "B", "B"),
    a_vs_b   = c(1, 1, 0, 0),
    stringsAsFactors = FALSE
  )
  comparisons <- parse_contrasts(ss, group_col = "GroupID",
                                 first_contrast_col = 4)
  expect_length(comparisons, 1)
  expect_equal(comparisons[[1]]$exp, "A")
})
