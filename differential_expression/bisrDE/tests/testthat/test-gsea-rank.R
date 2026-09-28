# v1.7: the GSEA ranking vector is built in one place, ranks by the Wald
# statistic by default, falls back to log2FC when `stat` is absent, dedupes by
# first occurrence, and orders ties deterministically.

test_that(".gsea_rank_vector ranks by stat and is deterministic", {
  df <- data.frame(
    ENSEMBL_ID     = c("g3", "g1", "g2", "g2", "g4", "g5"),
    baseMean       = c(10, 10, 10, 10, 0, 10),
    stat           = c(1, 5, -2, 9, 4, 5),
    log2FoldChange = c(0.1, 2, -1, 3, 1, 2)
  )
  v <- bisrDE:::.gsea_rank_vector(df, "ENSEMBL_ID", "stat")
  expect_equal(attr(v, "metric"), "stat")
  expect_equal(names(v), c("g1", "g5", "g3", "g2"))   # ties (g1, g5) broken by ID; g4 dropped (baseMean 0); g2 first occurrence
  expect_equal(as.numeric(v), c(5, 5, 1, -2))
})

test_that(".gsea_rank_vector falls back to log2FoldChange without a stat column", {
  df <- data.frame(ENTREZID = c("1", "2"), baseMean = c(1, 1),
                   log2FoldChange = c(-1, 2))
  expect_message(v <- bisrDE:::.gsea_rank_vector(df, "ENTREZID", "stat"), "no `stat` column")
  expect_equal(attr(v, "metric"), "log2fc")
  expect_equal(names(v), c("2", "1"))
})

test_that(".gsea_rank_vector drops NA ids and NA metrics", {
  df <- data.frame(ENTREZID = c(NA, "", "7", "8"), baseMean = 1,
                   stat = c(1, 2, NA, 3), log2FoldChange = 1)
  v <- bisrDE:::.gsea_rank_vector(df, "ENTREZID", "stat")
  expect_equal(names(v), "8")
})

test_that(".enrichment_table adds a significance flag and the metric", {
  skip_if_not_installed("clusterProfiler")
  # Minimal gseaResult stand-in: only @result is read.
  obj <- methods::new("gseaResult",
                      result = data.frame(ID = c("a", "b"), p.adjust = c(0.01, 0.2),
                                          NES = c(2, -1)),
                      geneList = c(x = 1), geneSets = list(), organism = "x",
                      setType = "x", keytype = "x", permScores = matrix(0),
                      params = list(), gene2Symbol = character(0), readable = FALSE)
  tb <- bisrDE:::.enrichment_table(obj, padj = 0.05, gsea_rank = "stat")
  expect_equal(tb$significant, c(TRUE, FALSE))
  expect_equal(unique(tb$ranking_metric), "stat")
  sig <- bisrDE:::.significant_sets(obj, 0.05)
  expect_equal(nrow(sig@result), 1)
  expect_null(bisrDE:::.significant_sets(obj, 0.001))
})
