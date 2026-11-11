# v1.7: QC at a glance — computed verdict on a synthetic two-group matrix.

test_that("qc_at_a_glance reports separation and no outliers for clean groups", {
  set.seed(7)
  n_genes <- 300
  base <- rnorm(n_genes, 8, 2)
  # Group B shifts 60 genes strongly; replicates get small noise.
  shift <- c(rep(3, 60), rep(0, n_genes - 60))
  mk <- function(s) base + s * shift + rnorm(n_genes, 0, 0.2)
  qc_mat <- cbind(A1 = mk(0), A2 = mk(0), A3 = mk(0), B1 = mk(1), B2 = mk(1), B3 = mk(1))
  attr(qc_mat, "transform") <- "vst"
  counts <- round(2^qc_mat)
  sample_info <- data.frame(sample = colnames(qc_mat),
                            condition = rep(c("A", "B"), each = 3))

  q <- qc_at_a_glance(counts, qc_mat, sample_info)
  expect_true(q$separation > 2)
  expect_length(q$outliers, 0)
  expect_equal(q$transform, "vst")
  expect_true(grepl("separate clearly", q$text))
  expect_true(grepl("PC1 explains", q$text))
})

test_that("qc_at_a_glance flags a sample far from its group", {
  set.seed(8)
  n_genes <- 300
  base <- rnorm(n_genes, 8, 2)
  mk <- function(extra = 0) base + rnorm(n_genes, 0, 0.2) + extra
  qc_mat <- cbind(A1 = mk(), A2 = mk(), A3 = mk(), A4 = mk(),
                  B1 = mk(), B2 = mk(), B3 = mk(),
                  B4 = base + rnorm(n_genes, 0, 0.2) + c(rep(6, 100), rep(0, 200)))
  attr(qc_mat, "transform") <- "log2_tmm_cpm"
  counts <- round(2^qc_mat)
  sample_info <- data.frame(sample = colnames(qc_mat),
                            condition = rep(c("A", "B"), each = 4))
  q <- qc_at_a_glance(counts, qc_mat, sample_info)
  expect_true("B4" %in% q$outliers)
  expect_true(grepl("far from every replicate", q$text))
})


test_that("qc_at_a_glance flags a single outlier in a 3 vs 3 design", {
  set.seed(9)
  n_genes <- 2000
  base <- rnorm(n_genes, 8, 2)
  mk <- function() base + rnorm(n_genes, 0, 0.2)
  shift <- c(rep(8, 300), rep(0, n_genes - 300))
  qc_mat <- cbind(A1 = mk(), A2 = mk(), A3 = mk(), B1 = mk(), B2 = mk(), B3 = mk() + shift)
  attr(qc_mat, "transform") <- "vst"
  counts <- round(2^qc_mat)
  sample_info <- data.frame(sample = colnames(qc_mat), condition = rep(c("A", "B"), each = 3))
  q <- qc_at_a_glance(counts, qc_mat, sample_info)
  expect_equal(q$outliers, "B3")
  expect_true(q$outlier_ratio[["B3"]] > 2)
  expect_true(grepl("far from every replicate", q$text))
  expect_equal(q$n_groups, 2)
  expect_equal(length(q$libsize_M_range), 2)
})

test_that("qc_at_a_glance does not flag clean 3 vs 3 replicates", {
  set.seed(10)
  n_genes <- 2000
  base <- rnorm(n_genes, 8, 2)
  shift <- c(rep(3, 100), rep(0, n_genes - 100))
  mk <- function(s) base + s * shift + rnorm(n_genes, 0, 0.3)
  qc_mat <- cbind(A1 = mk(0), A2 = mk(0), A3 = mk(0), B1 = mk(1), B2 = mk(1), B3 = mk(1))
  attr(qc_mat, "transform") <- "vst"
  q <- qc_at_a_glance(round(2^qc_mat), qc_mat,
                      data.frame(sample = colnames(qc_mat), condition = rep(c("A", "B"), each = 3)))
  expect_length(q$outliers, 0)
})

test_that("qc_at_a_glance distinguishes one group from no estimable spread", {
  set.seed(11)
  m <- matrix(rnorm(400, 8, 2), nrow = 100, dimnames = list(NULL, c("S1", "S2", "S3", "S4")))
  attr(m, "transform") <- "vst"
  one <- qc_at_a_glance(round(2^m), m, data.frame(sample = colnames(m), condition = "A"))
  expect_true(grepl("Only one group", one$text))
  singletons <- qc_at_a_glance(round(2^m), m, data.frame(sample = colnames(m), condition = c("A", "B", "C", "D")))
  expect_true(grepl("no within-group spread", singletons$text))
  expect_false(grepl("Only one group", singletons$text))
})
