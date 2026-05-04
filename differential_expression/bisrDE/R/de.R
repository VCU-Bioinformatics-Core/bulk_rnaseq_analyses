# Differential expression: DESeq2 driver + significance filter.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.2 migration spec.

#' Run DESeq2 for one contrast
#'
#' @description Run `DESeq()` and extract `results()` for a single
#'   experimental-vs-control contrast on a pre-built `DESeqDataSet`.
#'   Re-levels `dds$condition` so `ctrl` is the reference, then queries
#'   results with `cooksCutoff = TRUE` and `independentFiltering = FALSE`
#'   (matches the parent pipeline's behavior).
#'
#' @param dds A `DESeqDataSet` whose `condition` factor includes both
#'   `exp` and `ctrl` as levels. Pre-filtering and design specification
#'   should already be done by the caller.
#' @param exp Experimental group label (must be a level of
#'   `dds$condition`).
#' @param ctrl Control group label (must be a level of `dds$condition`).
#' @return A data frame of DESeq2 results (one row per gene) with the
#'   standard columns: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`,
#'   `pvalue`, `padj`. Rownames are gene IDs from `dds`.
#' @details Errors loudly via `cli::cli_alert_danger` + `stop()` if either
#'   `exp` or `ctrl` is absent from `levels(dds$condition)` — that case
#'   typically indicates a samplesheet/counts mismatch upstream.
#'
#' @importFrom DESeq2 DESeq results
#' @importFrom stats relevel
#' @export
perform_deseq2_analysis <- function(dds, exp, ctrl) {
  tryCatch(
    {
      cli::cli_alert_info("DESeq2: {exp} vs {ctrl}")
      cli::cli_inform("Condition levels: {.val {levels(dds$condition)}}")

      if (!(exp %in% levels(dds$condition))) {
        stop(paste("Experimental group", exp, "not found in condition levels"))
      }
      if (!(ctrl %in% levels(dds$condition))) {
        stop(paste("Control group", ctrl, "not found in condition levels"))
      }

      dds$condition <- stats::relevel(dds$condition, ref = ctrl)

      cli::cli_alert_info("Running DESeq()...")
      dds <- DESeq2::DESeq(dds)

      cli::cli_alert_info("Extracting contrast results...")
      res <- DESeq2::results(
        dds,
        contrast              = c("condition", exp, ctrl),
        cooksCutoff           = TRUE,
        independentFiltering  = FALSE
      )

      cli::cli_alert_success("DESeq2 done; {nrow(res)} genes in result table")
      as.data.frame(res)
    },
    error = function(e) {
      cli::cli_alert_danger("DESeq2 error: {e$message}")
      stop(e)
    }
  )
}


#' Filter DESeq2 results to significant genes
#'
#' @description Subset a DESeq2 (or annotated) results data frame to rows
#'   whose adjusted p-value and absolute log2 fold change pass the BISR
#'   default thresholds. Rows with `NA` `padj` (e.g. genes filtered out by
#'   DESeq2's Cook's-cutoff or low-count rules) are dropped.
#'
#' @param results A data frame containing at least `padj` and
#'   `log2FoldChange` columns. Accepts either raw [perform_deseq2_analysis()]
#'   output or its annotated counterpart.
#' @param padj_threshold Maximum adjusted p-value for retention (default
#'   `0.05`).
#' @param lfc_threshold Minimum absolute log2 fold change for retention
#'   (default `0.58`, i.e. ~1.5x fold change).
#' @return A data frame containing only rows passing both thresholds, in
#'   the original order.
#' @details Adoption sites for this helper are migrated in later phases:
#'   `generate_volcano()` and `run_analysis()` currently inline the same
#'   predicate. Until those move, this function is a no-op for the
#'   parent pipeline (which still goes through the inline path).
#'
#' @export
filter_significant <- function(results,
                               padj_threshold = 0.05,
                               lfc_threshold  = 0.58) {
  keep <- !is.na(results$padj) &
    results$padj <= padj_threshold &
    abs(results$log2FoldChange) >= lfc_threshold
  results[keep, , drop = FALSE]
}
