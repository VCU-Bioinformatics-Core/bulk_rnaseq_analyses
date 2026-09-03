# v1.7: with a log2FC_shrunken column the volcano places points at the shrunken
# value while colouring by the MLE-based DEG call; use_shrunken = FALSE reverts.

.vfix <- function() {
  data.frame(
    log2FoldChange  = c(3, -3, 0.2),
    log2FC_shrunken = c(1.5, -1.2, 0.1),
    padj            = c(0.001, 0.001, 0.5),
    SYMBOL          = c("A", "B", "C"),
    ENSEMBL_ID      = c("E1", "E2", "E3")
  )
}

test_that("generate_volcano uses the shrunken LFC for x and the MLE for the call", {
  p <- generate_volcano(.vfix(), "X", "Y", p = 0.05, lfc = 0.58)
  expect_equal(p$labels$x, "Log2 Fold-Change (shrunken)")
  built <- ggplot2::ggplot_build(p)
  expect_equal(sort(built$data[[1]]$x), sort(c(1.5, -1.2, 0.1)))
  expect_equal(as.character(p$data$color_tag), c("Over expressed", "Under expressed", "Not significant"))
})

test_that("generate_volcano(use_shrunken = FALSE) keeps the MLE axis", {
  p <- generate_volcano(.vfix(), "X", "Y", use_shrunken = FALSE)
  expect_equal(p$labels$x, "Log2 Fold-Change (FC)")
  expect_equal(sort(ggplot2::ggplot_build(p)$data[[1]]$x), sort(c(3, -3, 0.2)))
})
