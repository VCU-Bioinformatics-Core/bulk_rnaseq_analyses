# v1.7: DESeq2 options — independent filtering flag, LFC shrinkage columns and
# the provenance attributes — on a small simulated dataset.

.sim_dds <- function() {
  skip_if_not_installed("DESeq2")
  set.seed(11)
  dds <- DESeq2::makeExampleDESeqDataSet(n = 300, m = 6, betaSD = 1)
  dds$condition <- factor(rep(c("ctrl", "trt"), each = 3))
  dds
}

test_that("perform_deseq2_analysis records options and shrinks with type = normal", {
  dds <- .sim_dds()
  res <- suppressMessages(perform_deseq2_analysis(dds, "trt", "ctrl",
                                                  independent_filtering = FALSE,
                                                  lfc_shrink = "normal"))
  expect_true(all(c("log2FoldChange", "log2FC_shrunken", "lfcSE_shrunken") %in% colnames(res)))
  expect_false(any(is.na(res$padj[res$baseMean > 0 & !is.na(res$pvalue)])))  # no NA padj from filtering
  o <- attr(res, "de_options")
  expect_equal(o$lfc_shrink, "normal")
  expect_false(o$independent_filtering)
  expect_equal(o$reference_level, "ctrl")
  expect_equal(o$coef, "condition_trt_vs_ctrl")
  expect_s4_class(attr(res, "dds"), "DESeqDataSet")
  # Shrinkage pulls toward zero on average, never changes the MLE column.
  expect_true(mean(abs(res$log2FC_shrunken), na.rm = TRUE) <=
                mean(abs(res$log2FoldChange), na.rm = TRUE))
})

test_that("perform_deseq2_analysis with apeglm produces shrunken columns when installed", {
  skip_if_not_installed("apeglm")
  dds <- .sim_dds()
  res <- suppressMessages(perform_deseq2_analysis(dds, "trt", "ctrl", lfc_shrink = "apeglm"))
  expect_equal(attr(res, "de_options")$lfc_shrink, "apeglm")
  expect_true("log2FC_shrunken" %in% colnames(res))
})

test_that("independent filtering = TRUE can introduce NA padj; FALSE does not", {
  dds <- .sim_dds()
  a <- suppressMessages(perform_deseq2_analysis(dds, "trt", "ctrl", independent_filtering = TRUE))
  b <- suppressMessages(perform_deseq2_analysis(dds, "trt", "ctrl", independent_filtering = FALSE))
  expect_true(attr(a, "de_options")$independent_filtering)
  expect_gt(sum(is.na(a$padj)), sum(is.na(b$padj)))   # strict: proves the flag reached results()
  expect_false("log2FC_shrunken" %in% colnames(b))  # default lfc_shrink = "none"
})


test_that(".contrast_coef matches DESeq2's make.names of the whole coefficient", {
  skip_if_not_installed("DESeq2")
  set.seed(12)
  dds <- DESeq2::makeExampleDESeqDataSet(n = 200, m = 6, betaSD = 1)
  dds$condition <- factor(rep(c("0h", "24h+drug"), each = 3), levels = c("0h", "24h+drug"))
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  expect_equal(bisrDE:::.contrast_coef(dds, "24h+drug", "0h"), "condition_24h.drug_vs_0h")
  expect_null(bisrDE:::.contrast_coef(dds, "nope", "0h"))
})

test_that("a shrinkage failure does not lose the comparison", {
  skip_if_not_installed("DESeq2")
  dds <- .sim_dds()
  res <- testthat::with_mocked_bindings(
    suppressMessages(perform_deseq2_analysis(dds, "trt", "ctrl", lfc_shrink = "normal")),
    lfcShrink = function(...) stop("simulated shrinkage failure"),
    .package = "DESeq2"
  )
  expect_true(all(c("log2FoldChange", "padj") %in% colnames(res)))
  expect_false("log2FC_shrunken" %in% colnames(res))
  expect_equal(attr(res, "de_options")$lfc_shrink, "none")
})
