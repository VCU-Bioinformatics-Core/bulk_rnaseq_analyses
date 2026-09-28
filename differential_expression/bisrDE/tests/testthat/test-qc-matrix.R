# v1.7: qc_matrix() returns the blind VST on the dds genes, or log2(TMM + 1)
# restricted to the same genes when VST is unavailable.

test_that("qc_matrix returns VST on the filtered genes with a transform attribute", {
  skip_if_not_installed("DESeq2")
  set.seed(3)
  dds <- DESeq2::makeExampleDESeqDataSet(n = 1500, m = 6)
  tmm <- run_tmm(DESeq2::counts(dds))
  m <- qc_matrix(dds, tmm)
  expect_equal(attr(m, "transform"), "vst")
  expect_equal(dim(m), c(nrow(dds), ncol(dds)))
  expect_setequal(rownames(m), rownames(dds))
})

test_that("qc_matrix falls back to log2(TMM+1) on the dds gene universe", {
  skip_if_not_installed("DESeq2")
  set.seed(4)
  dds <- DESeq2::makeExampleDESeqDataSet(n = 40, m = 4)
  full <- rbind(DESeq2::counts(dds), extra = rep(0, 4))   # a gene the dds no longer has
  rownames(full)[nrow(full)] <- "dropped_gene"
  tmm <- run_tmm(full)
  m <- testthat::with_mocked_bindings(qc_matrix(dds, tmm), run_vst = function(...) NULL)
  expect_equal(attr(m, "transform"), "log2_tmm_cpm")
  expect_setequal(rownames(m), rownames(dds))
  expect_false("dropped_gene" %in% rownames(m))
})
