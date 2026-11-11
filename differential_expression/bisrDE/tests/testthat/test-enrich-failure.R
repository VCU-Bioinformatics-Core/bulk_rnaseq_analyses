# A backend that throws must be told apart from one with nothing to report.
# Both return NULL; the failure also signals a `bisrde_enrich_failure`
# condition that run_analysis() records in `enrich_status`, and the report's
# summary table turns the pair into "no result", "unreachable" or "failed".

.kegg_input <- function() {
  data.frame(ENTREZID = as.character(1:5), baseMean = 10,
             stat = c(3, 2, 1, -1, -2), log2FoldChange = c(1, 1, 0, -1, -1))
}

# Run process_kegg_gsea() with gseKEGG throwing `msg`; return the value and
# every failure condition signalled on the way.
.kegg_with_error <- function(msg) {
  failures <- list()
  testthat::local_mocked_bindings(
    gseKEGG = function(...) stop(msg, call. = FALSE),
    .package = "clusterProfiler"
  )
  value <- suppressMessages(withCallingHandlers(
    process_kegg_gsea(.kegg_input(), annotation = "mouse",
                      annotation_db = "unused"),
    bisrde_enrich_failure = function(cnd) failures[[length(failures) + 1L]] <<- cnd
  ))
  list(value = value, failures = failures)
}

test_that("a KEGG download error returns NULL and signals a network failure", {
  skip_if_not_installed("clusterProfiler")
  out <- .kegg_with_error("Failed to download KEGG data. Wrong 'species' or the network is unreachable")
  expect_null(out$value)
  expect_length(out$failures, 1L)
  f <- out$failures[[1]]
  expect_s3_class(f, "bisrde_enrich_failure")
  expect_equal(f$backend, "KEGG")
  expect_true(f$network)
  expect_match(conditionMessage(f), "Failed to download KEGG data")
})

test_that("a non-network KEGG error is flagged network = FALSE", {
  skip_if_not_installed("clusterProfiler")
  out <- .kegg_with_error("subscript out of bounds")
  expect_null(out$value)
  expect_length(out$failures, 1L)
  expect_equal(out$failures[[1]]$backend, "KEGG")
  expect_false(out$failures[[1]]$network)
})

test_that("with no handler installed a failure still just returns NULL", {
  skip_if_not_installed("clusterProfiler")
  testthat::local_mocked_bindings(
    gseKEGG = function(...) stop("timeout", call. = FALSE),
    .package = "clusterProfiler"
  )
  expect_null(suppressMessages(
    process_kegg_gsea(.kegg_input(), annotation = "mouse", annotation_db = "unused")
  ))
})

test_that("the failure helper names its backend and classifies the message", {
  f <- function(backend, msg) {
    got <- NULL
    suppressMessages(withCallingHandlers(
      bisrDE:::.enrich_failure(backend, simpleError(msg)),
      bisrde_enrich_failure = function(cnd) got <<- cnd
    ))
    got
  }
  expect_true(f("Hallmark", "Timeout was reached [zenodo.org]")$network)
  expect_true(f("KEGG", "cannot open URL 'https://rest.kegg.jp/x'")$network)
  expect_true(f("Hallmark", "Could not resolve host: x.org")$network)
  expect_false(f("KEGG", "object 'x' not found")$network)
  expect_equal(f("Reactome", "boom")$backend, "Reactome")
  # GO and Reactome read local annotation packages: never a network fault.
  expect_false(f("GO", "cannot open the connection")$network)
  expect_false(f("Reactome", "cannot open URL 'https://x.org'")$network)
})

# Each exported backend function must route its errors through the helper
# with its own label. The first call into the backend is mocked to throw.
.backend_input <- function() {
  data.frame(ENSEMBL_ID = paste0("G", 1:5), ENTREZID = as.character(1:5),
             baseMean = 10, stat = c(3, 2, 1, -1, -2),
             log2FoldChange = c(1, 1, 0, -1, -1))
}

.failure_of <- function(expr) {
  failures <- list()
  value <- suppressMessages(withCallingHandlers(
    expr,
    bisrde_enrich_failure = function(cnd) failures[[length(failures) + 1L]] <<- cnd
  ))
  list(value = value, failures = failures)
}

test_that("GO, Reactome and Hallmark route errors to the failure helper", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("ReactomePA")
  skip_if_not_installed("msigdbr")
  boom <- function(...) stop("boom", call. = FALSE)
  x <- .backend_input()

  testthat::local_mocked_bindings(gseGO = boom, .package = "clusterProfiler")
  out <- .failure_of(process_gsea(x, annotation_db = "unused"))
  expect_null(out$value)
  expect_equal(vapply(out$failures, `[[`, "", "backend"), "GO")

  testthat::local_mocked_bindings(gsePathway = boom, .package = "ReactomePA")
  out <- .failure_of(process_reactome_gsea(x, annotation = "mouse",
                                           annotation_db = "unused"))
  expect_null(out$value)
  expect_equal(vapply(out$failures, `[[`, "", "backend"), "Reactome")

  testthat::local_mocked_bindings(msigdbr = boom, .package = "msigdbr")
  out <- .failure_of(process_msigdb_hallmark(x, annotation = "mouse",
                                             annotation_db = "unused"))
  expect_null(out$value)
  expect_equal(vapply(out$failures, `[[`, "", "backend"), "Hallmark")
})

test_that("a backend that returns a result signals no failure", {
  skip_if_not_installed("clusterProfiler")
  obj <- methods::new("gseaResult",
                      result = data.frame(ID = "a", p.adjust = 0.01, NES = 2),
                      geneList = c(x = 1), geneSets = list(), organism = "x",
                      setType = "x", keytype = "x", permScores = matrix(0),
                      params = list(), gene2Symbol = character(0),
                      readable = FALSE)
  testthat::local_mocked_bindings(
    gseGO = function(...) obj,
    setReadable = function(x, ...) x,
    .package = "clusterProfiler"
  )
  out <- .failure_of(process_gsea(.backend_input(), annotation_db = "unused"))
  expect_identical(out$value, obj)
  expect_length(out$failures, 0L)
})

test_that("enrichment_cell counts significant sets when a result exists", {
  skip_if_not_installed("clusterProfiler")
  obj <- methods::new("gseaResult",
                      result = data.frame(ID = c("a", "b", "c"),
                                          p.adjust = c(0.01, 0.2, NA),
                                          NES = c(2, -1, 1)),
                      geneList = c(x = 1), geneSets = list(), organism = "x",
                      setType = "x", keytype = "x", permScores = matrix(0),
                      params = list(), gene2Symbol = character(0),
                      readable = FALSE)
  expect_equal(enrichment_cell(obj, list(state = "ok"), 0.05), "1")
  expect_equal(enrichment_cell(obj, NULL, 0.05), "1")
  expect_equal(enrichment_cell(obj, list(state = "ok"), 0.5), "2")
})

test_that("enrichment_cell words each NULL state", {
  expect_equal(enrichment_cell(NULL, list(state = "empty"), 0.05), "no result")
  expect_equal(enrichment_cell(NULL, NULL, 0.05), "no result")
  expect_equal(
    enrichment_cell(NULL, list(state = "failed", message = "m", network = TRUE), 0.05),
    "unreachable")
  expect_equal(
    enrichment_cell(NULL, list(state = "failed", message = "m", network = FALSE), 0.05),
    "failed")
})

test_that("run_analysis records a failed KEGG backend next to empty ones", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Mm.eg.db")
  counts_file <- file.path("..", "..", "..", "assets", "example_counts.tsv")
  sheet_file  <- file.path("..", "..", "..", "assets", "example_samplesheet.csv")
  skip_if_not(file.exists(counts_file) && file.exists(sheet_file),
              "bundled example fixture not found")

  counts <- suppressMessages(read_counts(counts_file))
  gene_meta <- attr(counts, "gene_meta")
  sheet  <- suppressMessages(read_samplesheet(sheet_file))
  comps  <- suppressMessages(parse_contrasts(sheet))
  countsdf <- counts |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),
                                \(x) as.integer(round(x))))
  countsdf <- suppressMessages(align_counts_to_samplesheet(countsdf, sheet))
  tmm <- run_tmm(countsdf)
  sample_info <- data.frame(sample = sheet$SampleID, condition = sheet$GroupID,
                            stringsAsFactors = FALSE)
  sample_info$display <- bisrDE:::.auto_display(sample_info$condition)
  dds <- suppressWarnings(suppressMessages(DESeq2::DESeqDataSetFromMatrix(
    countData = countsdf, colData = sample_info, design = ~condition)))
  outdir <- tempfile("bisrde_enrich_")
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  out_dirs <- suppressMessages(setup_directories(outdir))
  org_db <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))

  # KEGG goes through the real process_kegg_gsea() with only the download
  # failing; the other three backends return NULL without a failure.
  testthat::local_mocked_bindings(
    gseKEGG = function(...) stop("Failed to download KEGG data", call. = FALSE),
    .package = "clusterProfiler"
  )
  testthat::local_mocked_bindings(
    process_gsea            = function(...) NULL,
    process_reactome_gsea   = function(...) NULL,
    process_msigdb_hallmark = function(...) NULL
  )
  res <- suppressWarnings(suppressMessages(run_analysis(
    comps[[1]], dds, tmm, sample_info, out_dirs,
    annotation = "mouse", annotation_db = org_db, gene_meta = gene_meta
  )))

  expect_false(is.null(res))
  st <- res$enrich_status
  expect_named(st, c("GO", "KEGG", "Reactome", "Hallmark"))
  expect_equal(st$KEGG$state, "failed")
  expect_true(st$KEGG$network)
  expect_match(st$KEGG$message, "Failed to download KEGG data")
  expect_equal(st$GO$state, "empty")
  expect_equal(st$Reactome$state, "empty")
  expect_equal(st$Hallmark$state, "empty")
  expect_null(res$kegg)
})

test_that(".enrich_status_entry maps a result, NULL and a failure to states", {
  expect_equal(bisrDE:::.enrich_status_entry("obj")$state, "ok")
  expect_equal(bisrDE:::.enrich_status_entry(NULL)$state, "empty")
  failure <- structure(class = c("bisrde_enrich_failure", "condition"),
                       list(message = "m", call = NULL, backend = "KEGG",
                            network = TRUE))
  e <- bisrDE:::.enrich_status_entry(NULL, failure)
  expect_equal(e$state, "failed")
  expect_equal(e$message, "m")
  expect_true(e$network)
})

test_that("the run summary JSON carries enrichment_status, and tolerates old bundles", {
  deseq <- data.frame(padj = c(0.01, 0.5), log2FoldChange = c(2, 0.1))
  status <- list(
    GO       = list(state = "ok", message = NULL, network = FALSE),
    KEGG     = list(state = "failed", message = "Failed to download KEGG data",
                    network = TRUE),
    Reactome = list(state = "empty", message = NULL, network = FALSE),
    Hallmark = list(state = "empty", message = NULL, network = FALSE)
  )
  results <- list(list(deseq = deseq, enrich_status = status),
                  list(deseq = deseq))           # bundle without enrich_status
  comparisons <- list(list(name = "a_vs_b"), list(name = "c_vs_d"))
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)

  suppressMessages(bisrDE:::.write_report_json(
    path = path, runid = "r", brs_ticket = "", annotation = "mouse",
    id_type = "ensembl", counts_path = "c", samplesheet_path = "s",
    outdir = "o", started = Sys.time(), finished = Sys.time(),
    n_genes_input = 2L, n_genes_filtered = 2L, n_samples = 4L,
    n_samplesheet_rows = 4L, comparisons = comparisons, results = results,
    rds_path = "x.rds", session_log = "log"
  ))
  js <- jsonlite::read_json(path)

  es <- js$comparisons[[1]]$enrichment_status
  expect_equal(es$KEGG$state, "failed")
  expect_equal(es$KEGG$message, "Failed to download KEGG data")
  expect_equal(es$GO$state, "ok")
  expect_null(es$GO$message)
  expect_equal(es$Reactome$state, "empty")
  expect_null(js$comparisons[[2]]$enrichment_status)
  expect_true(js$comparisons[[2]]$ok)
})
