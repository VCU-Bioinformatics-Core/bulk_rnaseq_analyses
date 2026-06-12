# v1.5.0 — optional DisplayName plot labels. The matrix join key
# (make.names(SampleID)) must stay intact; only the visible labels change.

test_that(".display_lookup is identity when no display column", {
  si <- data.frame(
    sample    = c("5479_18", "S2", "S3"),
    condition = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  lk <- bisrDE:::.display_lookup(si)
  munged <- make.names(si$sample)            # "X5479_18", "S2", "S3"
  expect_equal(unname(lk[munged]), munged)   # identity fallback
  expect_equal(names(lk), munged)
})

test_that(".display_lookup uses display column when present", {
  si <- data.frame(
    sample    = c("SRR1", "SRR2", "SRR3"),
    condition = c("pos", "pos", "neg"),
    display   = c("HPV+ #1", "HPV+ #2", "HPV- #1"),
    stringsAsFactors = FALSE
  )
  lk <- bisrDE:::.display_lookup(si)
  expect_equal(unname(lk["SRR1"]), "HPV+ #1")
  expect_equal(unname(lk["SRR3"]), "HPV- #1")
})

test_that(".display_lookup falls back to munged ID per blank display cell", {
  si <- data.frame(
    sample    = c("SRR1", "SRR2", "SRR3"),
    condition = c("pos", "pos", "neg"),
    display   = c("HPV+ #1", "", NA),
    stringsAsFactors = FALSE
  )
  lk <- bisrDE:::.display_lookup(si)
  expect_equal(unname(lk["SRR1"]), "HPV+ #1")
  expect_equal(unname(lk["SRR2"]), "SRR2")  # blank -> munged
  expect_equal(unname(lk["SRR3"]), "SRR3")  # NA    -> munged
})

test_that(".heatmap_matrix keeps munged join key with NO display column", {
  sample_info <- data.frame(
    sample    = paste0("S", 1:6),
    condition = c("A", "A", "B", "B", "C", "C"),
    stringsAsFactors = FALSE
  )
  set.seed(1)
  nc <- matrix(rpois(60, 100), nrow = 10, ncol = 6,
               dimnames = list(paste0("g", 1:10), paste0("S", 1:6)))
  rdf <- data.frame(padj = rep(0.001, 10),
                    log2FoldChange = c(rep(2, 5), rep(-2, 5)),
                    row.names = paste0("g", 1:10))

  hm <- bisrDE:::.heatmap_matrix(rdf, nc, sample_info,
                                 exp_name = "A", ctrl_name = "B")
  # Default behaviour unchanged: columns are the munged SampleIDs.
  expect_setequal(colnames(hm$zscores), c("S1", "S2", "S3", "S4"))
})

test_that(".heatmap_matrix relabels columns when display column present", {
  sample_info <- data.frame(
    sample    = paste0("S", 1:6),
    condition = c("A", "A", "B", "B", "C", "C"),
    display   = paste0("Sample-", 1:6),
    stringsAsFactors = FALSE
  )
  set.seed(2)
  nc <- matrix(rpois(60, 100), nrow = 10, ncol = 6,
               dimnames = list(paste0("g", 1:10), paste0("S", 1:6)))
  rdf <- data.frame(padj = rep(0.001, 10),
                    log2FoldChange = c(rep(2, 5), rep(-2, 5)),
                    row.names = paste0("g", 1:10))

  hm <- bisrDE:::.heatmap_matrix(rdf, nc, sample_info,
                                 exp_name = "A", ctrl_name = "B")
  # Labels are display names; the A+B samples are S1..S4 -> Sample-1..4.
  expect_setequal(colnames(hm$zscores),
                  c("Sample-1", "Sample-2", "Sample-3", "Sample-4"))
})

test_that("pca_static label aesthetic falls back to sample with no display", {
  sample_info <- data.frame(
    sample    = paste0("S", 1:4),
    condition = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  set.seed(3)
  tmm <- matrix(rpois(40, 100), nrow = 10, ncol = 4,
                dimnames = list(paste0("g", 1:10), paste0("S", 1:4)))
  p <- pca_static(tmm, sample_info)
  expect_s3_class(p, "ggplot")
  expect_setequal(p$data$Display, paste0("S", 1:4))
})

test_that("pca_static uses display labels when present", {
  sample_info <- data.frame(
    sample    = paste0("S", 1:4),
    condition = c("A", "A", "B", "B"),
    display   = c("Ctrl-1", "Ctrl-2", "Trt-1", "Trt-2"),
    stringsAsFactors = FALSE
  )
  set.seed(4)
  tmm <- matrix(rpois(40, 100), nrow = 10, ncol = 4,
                dimnames = list(paste0("g", 1:10), paste0("S", 1:4)))
  p <- pca_static(tmm, sample_info)
  expect_setequal(p$data$Display, c("Ctrl-1", "Ctrl-2", "Trt-1", "Trt-2"))
})
