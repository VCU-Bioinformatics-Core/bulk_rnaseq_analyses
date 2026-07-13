test_that(".id_type_to_keytype maps the three valid inputs", {
  expect_equal(bisrDE:::.id_type_to_keytype("ensembl"), "ENSEMBL")
  expect_equal(bisrDE:::.id_type_to_keytype("entrez"),  "ENTREZID")
  expect_equal(bisrDE:::.id_type_to_keytype("symbol"),  "SYMBOL")
})

test_that(".id_type_to_keytype is case-insensitive", {
  expect_equal(bisrDE:::.id_type_to_keytype("ENSEMBL"), "ENSEMBL")
  expect_equal(bisrDE:::.id_type_to_keytype("Entrez"),  "ENTREZID")
})

test_that(".id_type_to_keytype rejects invalid input", {
  expect_error(bisrDE:::.id_type_to_keytype("refseq"), "Invalid id_type")
  expect_error(bisrDE:::.id_type_to_keytype(""), "Invalid id_type")
})

test_that("annotate_results adds ENTREZID, SYMBOL, GENENAME, and ENSEMBL_ID columns", {
  skip_if_not_installed("org.Hs.eg.db")
  org_hs <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))

  # Real human Ensembl IDs (TP53, BRCA1) plus a deliberately-bogus one to
  # exercise the all.x = TRUE merge.
  ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG00009999999")
  results <- data.frame(
    baseMean       = c(100, 200, 50),
    log2FoldChange = c(1.5, -2.0, 0.3),
    lfcSE          = c(0.2, 0.3, 0.1),
    pvalue         = c(0.001, 0.002, 0.5),
    padj           = c(0.01, 0.02, 0.6),
    row.names      = ids
  )

  annotated <- annotate_results(results, id_type = "ensembl",
                                annotation_db = org_hs)

  expect_true(all(c("ENSEMBL_ID", "ENTREZID", "SYMBOL", "GENENAME") %in% colnames(annotated)))
  expect_true(all(c("baseMean", "log2FoldChange", "padj") %in% colnames(annotated)))
  expect_equal(nrow(annotated), 3)

  # TP53 / BRCA1 should map to their canonical Symbols.
  tp53_row  <- annotated[annotated$ENSEMBL_ID == "ENSG00000141510", ]
  brca1_row <- annotated[annotated$ENSEMBL_ID == "ENSG00000012048", ]
  expect_equal(tp53_row$SYMBOL,  "TP53")
  expect_equal(brca1_row$SYMBOL, "BRCA1")

  # The bogus ID should still appear as a row (left join), with NA Entrez/Symbol.
  bogus_row <- annotated[annotated$ENSEMBL_ID == "ENSG00009999999", ]
  expect_equal(nrow(bogus_row), 1)
  expect_true(is.na(bogus_row$ENTREZID))
})

test_that("annotate_results emits a canonical SYMBOL column for id_type = 'symbol'", {
  skip_if_not_installed("org.Hs.eg.db")
  org_hs <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))

  # Input IDs are gene SYMBOLS. Regression: the annotated table must carry a
  # canonical `SYMBOL` column (not only `SYMBOL_ID`), else the volcano / heatmap
  # code that hardcodes `.data$SYMBOL` errors and aborts the whole comparison.
  results <- data.frame(
    log2FoldChange = c(1.5, -2.0),
    padj           = c(0.01, 0.02),
    row.names      = c("TP53", "BRCA1")
  )
  annotated <- annotate_results(results, id_type = "symbol",
                                annotation_db = org_hs)

  expect_true(all(c("ENSEMBL_ID", "ENTREZID", "SYMBOL", "GENENAME") %in%
                    colnames(annotated)))
  expect_setequal(annotated$SYMBOL, c("TP53", "BRCA1"))
})

test_that("annotate_results errors when annotation_db is unset", {
  results <- data.frame(
    log2FoldChange = 1,
    padj           = 0.01,
    row.names      = "ENSG00000141510"
  )
  expect_error(
    annotate_results(results, id_type = "ensembl", annotation_db = NULL),
    "annotation_db"
  )
})
