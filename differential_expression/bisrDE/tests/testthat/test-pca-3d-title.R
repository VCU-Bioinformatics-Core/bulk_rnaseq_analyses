# The 3D PCA title must report the variance explained by the three plotted
# components, not the sum over every component (which is always 100).

test_that("pca_plotly_3d titles the plot with the PC1-PC3 variance sum", {
  sample_info <- data.frame(
    sample    = paste0("S", 1:8),
    condition = rep(c("A", "B"), each = 4),
    stringsAsFactors = FALSE
  )
  set.seed(11)
  tmm <- matrix(rlnorm(80, 5, 1), nrow = 10, ncol = 8,
                dimnames = list(paste0("g", 1:10), paste0("S", 1:8)))

  var_pct <- bisrDE:::.compute_pca(tmm, rank = 3)$var_pct
  expect_gt(length(var_pct), 3)              # more than three components
  expected_pct <- round(sum(var_pct[1:3]), 1)
  expect_lt(expected_pct, 100)               # the first three do not sum to 100
  expected <- paste0("Variance explained by PC1 + PC2 + PC3 = ",
                     expected_pct, "%")

  p <- pca_plotly_3d(tmm, sample_info)
  title <- plotly::plotly_build(p)$x$layout$title
  if (is.list(title)) title <- title$text

  expect_equal(title, expected)
  shown <- as.numeric(sub("^.*= ([0-9.]+)%$", "\\1", title))
  expect_lt(shown, 100)
})
