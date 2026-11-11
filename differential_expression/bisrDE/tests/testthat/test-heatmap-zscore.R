# v1.7: heatmap z-scores are computed on log2(counts + 1), and a gene with zero
# variance across the shown samples gets z = 0 rather than NaN (no jitter).

test_that(".heatmap_matrix z-scores on the log scale and handles zero variance", {
  sample_info <- data.frame(sample = c("S1", "S2", "S3", "S4"),
                            condition = c("A", "A", "B", "B"))
  counts <- rbind(
    big   = c(1000, 100, 10, 1),     # geometric spread: log z-scores are linear in sample index
    small = c(10,   10,   1,  1),
    flat  = c(5, 5, 5, 5)            # zero variance
  )
  colnames(counts) <- sample_info$sample
  results_df <- data.frame(padj = c(0.001, 0.001, 0.001),
                           log2FoldChange = c(3, 3, 3),
                           row.names = rownames(counts))

  hm <- bisrDE:::.heatmap_matrix(results_df, counts, sample_info,
                                 p = 0.05, lfc = 0.58,
                                 exp_name = "A", ctrl_name = "B")
  z <- hm$zscores
  expect_false(any(is.nan(z)))
  expect_equal(unname(z["flat", ]), rep(0, 4))
  # Explicit expectation on the log scale; the linear-scale z-scores of the
  # same row (c(1.49, 0.05, -0.49, -0.55)) would fail this.
  expected <- as.numeric(scale(log2(c(1000, 100, 10, 1) + 1)))
  expect_equal(unname(z["big", ]), expected, tolerance = 1e-6)
  linear <- as.numeric(scale(c(1000, 100, 10, 1)))
  expect_false(isTRUE(all.equal(unname(z["big", ]), linear, tolerance = 0.05)))
})
