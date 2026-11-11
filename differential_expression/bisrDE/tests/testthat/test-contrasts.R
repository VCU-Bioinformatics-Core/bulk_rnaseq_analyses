test_that("parse_contrasts builds list(name, exp, ctrl) entries from a clean samplesheet", {
  ss <- data.frame(
    SampleID = paste0("S", 1:6),
    GroupID  = c("A", "A", "A", "B", "B", "B"),
    a_vs_b   = c(1, 1, 1, 0, 0, 0),
    stringsAsFactors = FALSE
  )

  comparisons <- parse_contrasts(ss)

  expect_length(comparisons, 1)
  expect_equal(comparisons[[1]]$name, "a_vs_b")
  expect_equal(comparisons[[1]]$exp,  "A")
  expect_equal(comparisons[[1]]$ctrl, "B")
})

test_that("parse_contrasts handles multiple contrasts", {
  ss <- data.frame(
    SampleID = paste0("S", 1:8),
    GroupID  = c("A", "A", "B", "B", "C", "C", "D", "D"),
    a_vs_b   = c(1, 1, 0, 0, NA, NA, NA, NA),
    c_vs_d   = c(NA, NA, NA, NA, 1, 1, 0, 0),
    stringsAsFactors = FALSE
  )

  comparisons <- parse_contrasts(ss)

  expect_length(comparisons, 2)
  expect_equal(comparisons[[1]]$name, "a_vs_b")
  expect_equal(comparisons[[2]]$name, "c_vs_d")
  expect_equal(comparisons[[2]]$exp,  "C")
  expect_equal(comparisons[[2]]$ctrl, "D")
})

test_that("parse_contrasts skips contrasts missing exp or ctrl group", {
  ss <- data.frame(
    SampleID = paste0("S", 1:4),
    GroupID  = c("A", "A", "B", "B"),
    only_exp = c(1, 1, NA, NA),         # no ctrl side
    only_ctrl = c(NA, NA, 0, 0),        # no exp side
    fine     = c(1, 1, 0, 0),           # ok
    stringsAsFactors = FALSE
  )

  comparisons <- parse_contrasts(ss)

  expect_length(comparisons, 1)
  expect_equal(comparisons[[1]]$name, "fine")
})

test_that("parse_contrasts respects first_contrast_col override", {
  ss <- data.frame(
    SampleID = paste0("S", 1:4),
    Tissue   = c("Liver", "Liver", "Lung", "Lung"),
    GroupID  = c("A", "A", "B", "B"),
    a_vs_b   = c(1, 1, 0, 0),
    stringsAsFactors = FALSE
  )

  # GroupID is now in column 3; contrasts start at column 4.
  comparisons <- parse_contrasts(ss, group_col = "GroupID",
                                 first_contrast_col = 4)
  expect_length(comparisons, 1)
  expect_equal(comparisons[[1]]$exp,  "A")
})

test_that("read_samplesheet round-trips a CSV correctly", {
  ss <- data.frame(
    SampleID = paste0("S", 1:4),
    GroupID  = c("A", "A", "B", "B"),
    a_vs_b   = c(1, 1, 0, 0),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(ss, tmp, row.names = FALSE)

  ss_read <- read_samplesheet(tmp)
  expect_equal(ss_read$SampleID, ss$SampleID)
  expect_equal(ss_read$GroupID,  ss$GroupID)
  expect_equal(ss_read$a_vs_b,   ss$a_vs_b)
})
