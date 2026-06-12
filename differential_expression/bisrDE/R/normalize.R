# Normalization wrappers: TMM (edgeR) and VST (DESeq2).
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R.

#' TMM-normalize a count matrix
#'
#' @description Wrap edgeR's TMM normalization (DGEList -> calcNormFactors ->
#'   cpm) into a single call. Returns the cpm matrix; the intermediate
#'   `DGEList` is discarded since downstream code only consumes the cpm
#'   output.
#'
#' @param counts A count matrix or data frame (Genes x Samples) of
#'   integer-coercible counts. Typically the output of
#'   [align_counts_to_samplesheet()] cast to integer.
#' @return A numeric matrix of TMM-normalized counts-per-million (Genes x
#'   Samples), preserving rownames and colnames of the input.
#' @details Persisting the result to disk is the caller's responsibility
#'   (the parent `de.R` writes a `normalizedCounts_tmm<DATE>.csv`; the
#'   package leaves that as a pipeline-orchestrator concern). If
#'   normalization-factor inspection becomes a need later, extend the
#'   signature with a `return_dge = FALSE` flag rather than changing the
#'   default return shape.
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @export
run_tmm <- function(counts) {
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge)
}


#' Variance-stabilized transform with a robust fallback
#'
#' @description Wrap `DESeq2::vst()` with a `tryCatch` that falls back to
#'   `DESeq2::varianceStabilizingTransformation()` when the fast vst path
#'   fails. The fast path subsamples 1000 genes for the dispersion fit and
#'   raises an error on small/sparse datasets; the fallback is slower but
#'   always-correct.
#'
#' @param dds A `DESeqDataSet`. Typically post-pre-filter, pre-`DESeq()` —
#'   `blind = TRUE` is the default to keep the transform design-agnostic
#'   (appropriate for sample-level QC across the whole experiment, not for
#'   any one contrast).
#' @param blind Passed through to `DESeq2::vst()` /
#'   `DESeq2::varianceStabilizingTransformation()`. Default `TRUE`.
#' @return A `DESeqTransform` object on success, or `NULL` if both the fast
#'   path and the fallback fail. Callers should branch on `is.null()` and
#'   skip downstream vst-dependent steps rather than abort the run.
#'
#' @importFrom DESeq2 vst varianceStabilizingTransformation
#' @export
run_vst <- function(dds, blind = TRUE) {
  tryCatch(
    DESeq2::vst(dds, blind = blind),
    error = function(e) {
      cli::cli_alert_warning(
        "run_vst: vst() failed ({e$message}); falling back to varianceStabilizingTransformation()"
      )
      tryCatch(
        DESeq2::varianceStabilizingTransformation(dds, blind = blind),
        error = function(e2) {
          cli::cli_alert_warning(
            "run_vst: VST fallback also failed ({e2$message}); returning NULL"
          )
          NULL
        }
      )
    }
  )
}
