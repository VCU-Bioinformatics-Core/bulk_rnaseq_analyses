# read_counts: Ensembl version suffixes are stripped, but gene symbols that
# legitimately contain dots (e.g. mouse "H2-M10.1", "Tex19.1") must NOT be
# collapsed — doing so produced duplicate row names and a hard error on
# salmon.merged.gene_counts.tsv files whose gene_id column holds symbols.

test_that("read_counts strips Ensembl versions but preserves symbol dots", {
  tsv <- tempfile(fileext = ".tsv")
  writeLines(c(
    "gene_id\tgene_name\tS1\tS2",
    "ENSMUSG00000000001.3\tGnai3\t10\t20",
    "ENSMUSG00000000002.15\tPbsn\t5\t6",
    "H2-M10.1\tH2-M10.1\t1\t2",
    "H2-M10.2\tH2-M10.2\t3\t4",
    "Tex19.1\tTex19.1\t7\t8"
  ), tsv)

  x <- read_counts(tsv)

  # Ensembl accessions: version suffix removed.
  expect_true(all(c("ENSMUSG00000000001", "ENSMUSG00000000002") %in% rownames(x)))
  # Symbols with dots: kept intact, distinct, no duplicate row names.
  expect_true(all(c("H2-M10.1", "H2-M10.2", "Tex19.1") %in% rownames(x)))
  expect_false(any(duplicated(rownames(x))))
  # Non-numeric gene_name column dropped; only the two sample columns remain.
  expect_equal(ncol(x), 2L)
})

test_that("read_counts loads an all-symbol counts matrix without collapsing", {
  tsv <- tempfile(fileext = ".tsv")
  # Four distinct H2-M10.* symbols that a blanket "\\..*" strip would collapse.
  writeLines(c(
    "gene_id\tgene_name\tA\tB",
    "H2-M10.1\tH2-M10.1\t1\t2",
    "H2-M10.2\tH2-M10.2\t3\t4",
    "H2-M10.3\tH2-M10.3\t5\t6",
    "Rn4.5s\tRn4.5s\t7\t8"
  ), tsv)

  expect_silent(x <- read_counts(tsv))
  expect_equal(nrow(x), 4L)
  expect_setequal(rownames(x), c("H2-M10.1", "H2-M10.2", "H2-M10.3", "Rn4.5s"))
})
