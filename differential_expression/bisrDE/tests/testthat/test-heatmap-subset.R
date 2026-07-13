# Phase 1.2 fix verification: the heatmap matrix must contain ONLY samples
# that belong to the current comparison's exp + ctrl groups, NOT every
# sample in the full samplesheet.

build_fixture <- function() {
  # 6 samples across 3 groups: A (2), B (2), C (2). The "A vs B" comparison
  # should yield a 4-column matrix (S1, S2, S3, S4) — never include C.
  sample_info <- data.frame(
    sample    = paste0("S", 1:6),
    condition = c("A", "A", "B", "B", "C", "C"),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  norm_counts <- matrix(
    rpois(60, lambda = 100),
    nrow = 10, ncol = 6,
    dimnames = list(paste0("g", 1:10), paste0("S", 1:6))
  )

  results_df <- data.frame(
    padj           = rep(0.001, 10),
    log2FoldChange = c(rep(2, 5), rep(-2, 5)),
    row.names      = paste0("g", 1:10)
  )

  list(
    sample_info = sample_info,
    norm_counts = norm_counts,
    results_df  = results_df
  )
}

test_that(".heatmap_matrix returns only the comparison's exp+ctrl samples", {
  fx <- build_fixture()

  hm <- bisrDE:::.heatmap_matrix(
    fx$results_df, fx$norm_counts, fx$sample_info,
    p = 0.05, lfc = 0.58,
    exp_name = "A", ctrl_name = "B"
  )

  expect_false(is.null(hm))
  expect_equal(sort(colnames(hm$zscores)), c("S1", "S2", "S3", "S4"))
  expect_false(any(c("S5", "S6") %in% colnames(hm$zscores)))
  expect_equal(nrow(hm$zscores), 10)
})

test_that(".heatmap_matrix returns NULL when filtering yields <2 samples", {
  fx <- build_fixture()
  # Comparison "A vs Z" — Z has no samples in sample_info, so only the 2
  # A-group samples remain. <2 of each group = NULL output.
  hm <- bisrDE:::.heatmap_matrix(
    fx$results_df, fx$norm_counts, fx$sample_info,
    p = 0.05, lfc = 0.58,
    exp_name = "Z", ctrl_name = "Z"
  )
  expect_null(hm)
})

test_that(".heatmap_matrix returns NULL when filtering yields <2 genes", {
  fx <- build_fixture()
  # Tighten thresholds: padj < 1e-100 keeps zero genes
  hm <- bisrDE:::.heatmap_matrix(
    fx$results_df, fx$norm_counts, fx$sample_info,
    p = 1e-100, lfc = 0.58,
    exp_name = "A", ctrl_name = "B"
  )
  expect_null(hm)
})

test_that(".heatmap_matrix top_n trims to the top N genes by ascending padj", {
  fx <- build_fixture()
  # Vary padj across genes so top-N selection is testable
  fx$results_df$padj <- seq(0.001, 0.01, length.out = 10)

  hm <- bisrDE:::.heatmap_matrix(
    fx$results_df, fx$norm_counts, fx$sample_info,
    p = 0.05, lfc = 0.58,
    exp_name = "A", ctrl_name = "B",
    top_n = 5
  )

  expect_equal(nrow(hm$zscores), 5)
  # Top 5 by ascending padj = g1..g5
  expect_equal(sort(rownames(hm$zscores)), paste0("g", 1:5))
})

test_that(".heatmap_matrix uses ENSEMBL_ID column when available", {
  fx <- build_fixture()
  fx$results_df$ENSEMBL_ID <- paste0("ENSG_", 1:10)
  rownames(fx$norm_counts) <- paste0("ENSG_", 1:10)

  hm <- bisrDE:::.heatmap_matrix(
    fx$results_df, fx$norm_counts, fx$sample_info,
    p = 0.05, lfc = 0.58,
    exp_name = "A", ctrl_name = "B"
  )

  expect_true(all(rownames(hm$zscores) %in% paste0("ENSG_", 1:10)))
})

test_that(".heatmap_matrix matches a symbol-keyed count matrix (id_type=symbol)", {
  fx <- build_fixture()
  # Simulate an --id-type symbol run: the count matrix is keyed by gene
  # SYMBOLS, the annotated DE table stores the input under SYMBOL_ID, and
  # ENSEMBL_ID holds *mapped* accessions that do NOT match the count matrix.
  # Before the fix this returned NULL (hardcoded ENSEMBL_ID never matched).
  syms <- paste0("Sym", 1:10)
  rownames(fx$norm_counts) <- syms
  fx$results_df$SYMBOL_ID  <- syms
  fx$results_df$SYMBOL     <- syms
  fx$results_df$ENSEMBL_ID <- paste0("ENSMUSG", 1:10) # mapped, non-matching

  hm <- bisrDE:::.heatmap_matrix(
    fx$results_df, fx$norm_counts, fx$sample_info,
    p = 0.05, lfc = 0.58, exp_name = "A", ctrl_name = "B"
  )

  expect_false(is.null(hm))
  expect_setequal(rownames(hm$zscores), syms)
})

test_that(".match_id_column picks the column overlapping the count matrix", {
  df <- data.frame(
    ENSEMBL_ID = paste0("ENSMUSG", 1:3),
    SYMBOL_ID  = c("Actb", "Gapdh", "H2-M10.1"),
    stringsAsFactors = FALSE
  )
  expect_equal(bisrDE:::.match_id_column(df, c("Actb", "Gapdh", "H2-M10.1")), "SYMBOL_ID")
  expect_equal(bisrDE:::.match_id_column(df, paste0("ENSMUSG", 1:3)), "ENSEMBL_ID")
  expect_null(bisrDE:::.match_id_column(df, c("nope", "nada")))
})
