# v1.7: one significance predicate everywhere. A gene exactly at the padj
# cutoff is NOT significant in filter_significant(), the heatmap matrix, or the
# volcano colouring; a gene exactly at the |log2FC| threshold IS.

.fixture <- function() {
  data.frame(
    padj           = c(0.05, 0.049, 0.05, 0.001, NA),
    log2FoldChange = c(2,    2,     0.58, 0.58,  3),
    SYMBOL         = paste0("g", 1:5),
    row.names      = paste0("g", 1:5)
  )
}

test_that("filter_significant excludes padj == cutoff and includes |lfc| == threshold", {
  x <- .fixture()
  keep <- filter_significant(x, padj_threshold = 0.05, lfc_threshold = 0.58)
  expect_setequal(rownames(keep), c("g2", "g4"))
})

test_that("volcano colouring agrees with filter_significant at the boundaries", {
  x <- .fixture()
  p <- generate_volcano(x, "A", "B", p = 0.05, lfc = 0.58, n_labels = 5)
  tags <- p$data$color_tag
  names(tags) <- rownames(p$data)
  expect_equal(as.character(tags[c("g1", "g3", "g5")]),
               rep("Not significant", 3))
  expect_true(all(as.character(tags[c("g2", "g4")]) == "Over expressed"))
})

test_that(".heatmap_matrix uses the same predicate", {
  x <- .fixture()
  sample_info <- data.frame(sample = c("S1", "S2", "S3", "S4"),
                            condition = c("A", "A", "B", "B"))
  set.seed(1)
  counts <- matrix(rpois(20, 50), nrow = 5,
                   dimnames = list(rownames(x), sample_info$sample))
  hm <- bisrDE:::.heatmap_matrix(x, counts, sample_info, p = 0.05, lfc = 0.58,
                                 exp_name = "A", ctrl_name = "B")
  expect_setequal(rownames(hm$zscores), c("g2", "g4"))
})
