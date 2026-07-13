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

test_that("read_counts drops Ensembl _PAR_Y pseudo-autosomal duplicates", {
  tsv <- tempfile(fileext = ".tsv")
  writeLines(c(
    "gene_id\tS1\tS2",
    "ENSG00000002586.20\t10\t20",
    "ENSG00000002586.20_PAR_Y\t5\t6", # collides with the chrX copy after strip
    "ENSG00000123456.3\t1\t2"
  ), tsv)

  x <- read_counts(tsv)

  expect_false(any(duplicated(rownames(x))))         # no dup -> no hard error
  expect_true("ENSG00000002586" %in% rownames(x))    # chrX copy kept
  expect_false(any(grepl("_PAR_Y", rownames(x))))    # PAR_Y row dropped
  expect_equal(nrow(x), 2L)
})

test_that("read_counts drops featureCounts numeric annotation columns", {
  tsv <- tempfile(fileext = ".tsv")
  # featureCounts-style: Chr/Start/End/Length are NUMERIC but not samples.
  writeLines(c(
    "Geneid\tChr\tStart\tEnd\tStrand\tLength\tS1\tS2",
    "ENSG00000000001.3\t1\t100\t900\t+\t800\t10\t20",
    "ENSG00000000002.5\t2\t200\t700\t-\t500\t5\t6"
  ), tsv)

  x <- read_counts(tsv)

  expect_setequal(colnames(x), c("S1", "S2")) # annotation cols dropped
  expect_equal(nrow(x), 2L)
  expect_setequal(rownames(x), c("ENSG00000000001", "ENSG00000000002"))
  expect_equal(unname(x[["S1"]]), c(10, 5))
})

test_that("read_counts munges sample column names with make.names (ID matching)", {
  tsv <- tempfile(fileext = ".tsv")
  # Raw salmon headers can contain hyphens; downstream matching uses
  # make.names(SampleID), so read_counts must munge the columns identically.
  writeLines(c(
    "gene_id\tgene_name\tHCT116-P53-minus-1\tHCT116-P53-plus-1",
    "ENSG00000000001.2\tActb\t10\t20"
  ), tsv)

  x <- read_counts(tsv)

  expect_equal(colnames(x), c("HCT116.P53.minus.1", "HCT116.P53.plus.1"))
})

test_that("read_counts drops the RSEM transcript_id(s) column", {
  tsv <- tempfile(fileext = ".tsv")
  writeLines(c(
    "gene_id\ttranscript_id(s)\tS1\tS2",
    "ENSG00000000001.3\tENST1,ENST2\t10\t20",
    "ENSG00000000002.5\tENST3\t5\t6"
  ), tsv)

  x <- read_counts(tsv)

  expect_setequal(colnames(x), c("S1", "S2"))
  expect_equal(nrow(x), 2L)
})

test_that("read_counts errors clearly when there are no sample columns", {
  tsv <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tgene_name", "ENSG1\tActb"), tsv)
  expect_error(read_counts(tsv), "no numeric sample columns")
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
