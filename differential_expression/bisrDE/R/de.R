# Differential expression: DESeq2 driver + significance filter.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.2 migration spec.

#' Run DESeq2 for one contrast
#'
#' @description Run `DESeq()` and extract `results()` for a single
#'   experimental-vs-control contrast on a pre-built `DESeqDataSet`.
#'   Re-levels `dds$condition` so `ctrl` is the reference, then queries
#'   results with `cooksCutoff = TRUE` and `independentFiltering =
#'   independent_filtering` (default `FALSE`, the historical BISR setting;
#'   DESeq2's own default is `TRUE`). Optionally appends a shrunken
#'   log2 fold change (`DESeq2::lfcShrink`).
#'
#' @param dds A `DESeqDataSet` whose `condition` factor includes both
#'   `exp` and `ctrl` as levels. Pre-filtering and design specification
#'   should already be done by the caller.
#' @param exp Experimental group label (must be a level of
#'   `dds$condition`).
#' @param ctrl Control group label (must be a level of `dds$condition`).
#' @param independent_filtering Passed to `DESeq2::results(independentFiltering =)`.
#'   Default `FALSE` keeps every tested gene in the table (no padj set to
#'   `NA` by the mean-count filter). `TRUE` is DESeq2's default and gains
#'   power on low-count genes at the cost of `NA` padj for the filtered ones.
#' @param lfc_shrink One of `"apeglm"`, `"normal"`, `"none"`. When not
#'   `"none"`, `DESeq2::lfcShrink()` is run and two columns are appended:
#'   `log2FC_shrunken` and `lfcSE_shrunken`. `"apeglm"` needs the apeglm
#'   package; if it is not installed the call falls back to `"normal"` with a
#'   warning. The MLE `log2FoldChange` column is never altered, so DEG calls
#'   are unaffected. Default `"none"`.
#' @return A data frame of DESeq2 results (one row per gene) with the
#'   standard columns: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`,
#'   `pvalue`, `padj` (plus `log2FC_shrunken` / `lfcSE_shrunken` when
#'   shrinkage is on). Rownames are gene IDs from `dds`. Two attributes carry
#'   provenance for the caller: `"dds"` (the fitted `DESeqDataSet`) and
#'   `"de_options"` (a list with `independent_filtering`, `lfc_shrink` as
#'   actually used, `design`, `reference_level`, `coef`).
#' @details Errors loudly via `cli::cli_alert_danger` + `stop()` if either
#'   `exp` or `ctrl` is absent from `levels(dds$condition)` — that case
#'   typically indicates a samplesheet/counts mismatch upstream.
#'
#' @importFrom DESeq2 DESeq results lfcShrink resultsNames
#' @importFrom stats relevel
#' @export
perform_deseq2_analysis <- function(dds, exp, ctrl,
                                    independent_filtering = FALSE,
                                    lfc_shrink = c("none", "apeglm", "normal")) {
  lfc_shrink <- match.arg(lfc_shrink)
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

      cli::cli_alert_info(
        "Extracting contrast results (independentFiltering = {independent_filtering})..."
      )
      res <- DESeq2::results(
        dds,
        contrast              = c("condition", exp, ctrl),
        cooksCutoff           = TRUE,
        independentFiltering  = independent_filtering
      )
      out <- as.data.frame(res)

      # ---- optional LFC shrinkage (MLE column is left untouched) ----
      coef_name <- .contrast_coef(dds, exp, ctrl)
      shrink_used <- "none"
      if (lfc_shrink != "none") {
        method <- lfc_shrink
        if (method == "apeglm" && !requireNamespace("apeglm", quietly = TRUE)) {
          cli::cli_alert_warning(
            "lfc_shrink = 'apeglm' requested but the apeglm package is not installed; using type = 'normal'"
          )
          method <- "normal"
        }
        if (method == "apeglm" && is.null(coef_name)) {
          cli::cli_alert_warning(
            "lfc_shrink = 'apeglm' needs a model coefficient for {exp} vs {ctrl} and none was found in resultsNames(); using type = 'normal'"
          )
          method <- "normal"
        }
        cli::cli_alert_info("Shrinking log2 fold changes (type = {.val {method}})...")
        # Shrinkage decorates a table that already exists; a failure here must
        # not take the comparison down with it.
        shr <- tryCatch(
          if (method == "apeglm") {
            DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
          } else {
            DESeq2::lfcShrink(dds, contrast = c("condition", exp, ctrl),
                              type = "normal", quiet = TRUE)
          },
          error = function(e) {
            cli::cli_alert_warning(
              "lfcShrink (type = {.val {method}}) failed: {conditionMessage(e)}; continuing without shrunken columns"
            )
            NULL
          }
        )
        if (!is.null(shr)) {
          shr <- as.data.frame(shr)[rownames(out), , drop = FALSE]
          out$log2FC_shrunken <- shr$log2FoldChange
          out$lfcSE_shrunken  <- shr$lfcSE
          shrink_used <- method
        }
      }

      attr(out, "dds") <- dds
      attr(out, "de_options") <- list(
        independent_filtering = independent_filtering,
        lfc_shrink            = shrink_used,
        design                = paste(deparse(DESeq2::design(dds)), collapse = ""),
        reference_level       = ctrl,
        coef                  = coef_name
      )

      cli::cli_alert_success("DESeq2 done; {nrow(out)} genes in result table")
      out
    },
    error = function(e) {
      cli::cli_alert_danger("DESeq2 error: {e$message}")
      stop(e)
    }
  )
}


#' Name of the DESeq2 model coefficient for `exp` vs `ctrl`
#'
#' @description After `relevel(ref = ctrl)` the coefficient DESeq2 names
#'   `condition_<exp>_vs_<ctrl>` estimates the contrast directly, which is
#'   what `lfcShrink(type = "apeglm")` needs (it only accepts `coef =`).
#'   DESeq2 may `make.names()` the level labels, so both spellings are tried.
#'
#' @param dds A `DESeqDataSet` after `DESeq()`.
#' @param exp,ctrl Group labels.
#' @return The coefficient name, or `NULL` if none matches.
#' @keywords internal
.contrast_coef <- function(dds, exp, ctrl) {
  rn <- DESeq2::resultsNames(dds)
  # DESeq2 applies make.names() to the WHOLE coefficient string (see
  # renameModelMatrixColumns), not to each label, so "24h+drug" vs "0h" becomes
  # "condition_24h.drug_vs_0h", never "condition_X24h.drug_vs_X0h".
  cands <- unique(c(
    paste0("condition_", exp, "_vs_", ctrl),
    make.names(paste0("condition_", exp, "_vs_", ctrl))
  ))
  hit <- cands[cands %in% rn]
  if (length(hit) == 0) NULL else hit[[1]]
}


#' Filter DESeq2 results to significant genes
#'
#' @description Subset a DESeq2 (or annotated) results data frame to rows
#'   whose adjusted p-value and absolute log2 fold change pass the BISR
#'   default thresholds. Rows with `NA` `padj` (e.g. genes filtered out by
#'   DESeq2's Cook's-cutoff or low-count rules) are dropped. The predicate is
#'   `padj < padj_threshold & abs(log2FoldChange) >= lfc_threshold`, the same
#'   one the volcano, heatmaps, DEG counts and report use, so a gene exactly
#'   at the padj cutoff is NOT significant anywhere in the pipeline.
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
    results$padj < padj_threshold &
    abs(results$log2FoldChange) >= lfc_threshold
  results[keep, , drop = FALSE]
}
