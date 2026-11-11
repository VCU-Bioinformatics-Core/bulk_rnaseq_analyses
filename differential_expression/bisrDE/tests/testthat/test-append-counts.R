# .append_counts joins per-sample normalized counts onto a DE table (v1.5.3:
# "normalized counts in the same spreadsheet as the log2FC").

test_that(".append_counts joins per-sample counts with original SampleID headers", {
  df  <- data.frame(SYMBOL_ID = c("g1", "g2", "g3"),
                    log2FoldChange = c(1, -1, 0), stringsAsFactors = FALSE)
  # Matrix columns are make.names()-munged (e.g. "6mon1" -> "X6mon1").
  mat <- matrix(1:9, nrow = 3,
                dimnames = list(c("g1", "g2", "g3"),
                                c("X6mon1", "met1", "X6mon2")))
  orig <- c(X6mon1 = "6mon1", met1 = "met1", X6mon2 = "6mon2")

  out <- bisrDE:::.append_counts(df, df$SYMBOL_ID, mat, orig, "TMM")

  expect_true(all(c("TMM_6mon1", "TMM_met1", "TMM_6mon2") %in% colnames(out)))
  expect_equal(out[["TMM_6mon1"]], c(1L, 2L, 3L)) # aligned by gene, first column
  expect_equal(out[["TMM_met1"]], c(4L, 5L, 6L))
  expect_equal(nrow(out), 3L)
  # original stat columns are preserved
  expect_true(all(c("SYMBOL_ID", "log2FoldChange") %in% colnames(out)))
})

test_that(".append_counts fills NA for genes absent from the matrix", {
  df  <- data.frame(SYMBOL_ID = c("g1", "gX"), stringsAsFactors = FALSE)
  mat <- matrix(c(10, 20), nrow = 1, dimnames = list("g1", c("s1", "s2")))

  out <- bisrDE:::.append_counts(df, df$SYMBOL_ID, mat,
                                 c(s1 = "s1", s2 = "s2"), "DESeq2norm")

  expect_equal(out[["DESeq2norm_s1"]][1], 10)
  expect_true(is.na(out[["DESeq2norm_s1"]][2])) # gX not in the matrix
})
