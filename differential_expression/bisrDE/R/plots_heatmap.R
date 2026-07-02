# Heatmap of DE genes (z-score normalized).
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.4 spec.

#' Heatmap of differentially expressed genes
#'
#' @description Render a z-score heatmap of DE-significant genes restricted
#'   to the samples that belong to the current comparison's experimental
#'   and control groups. Has two modes via `top_n`: `NULL` (the default)
#'   keeps every gene that passes thresholds with no row labels;
#'   `top_n = 100` (or any integer) takes the top-N genes by ascending
#'   `padj` and labels rows by gene `SYMBOL`.
#'
#' @param results_df Data frame containing (annotated) DESeq2 results. If a
#'   `SYMBOL` column is present and `top_n` is set, gene symbols are used
#'   as row labels.
#' @param normalized_counts Genes x Samples numeric matrix of normalized
#'   counts (typically the output of [run_tmm()]).
#' @param sample_info Data frame with at least two columns: `sample`
#'   (canonical IDs matching the `make.names()` form of
#'   `colnames(normalized_counts)`) and `condition` (group label).
#' @param p P-value threshold for significance. Default `0.05`.
#' @param lfc Log2 fold change threshold. Default `0.58` (≈1.5x).
#' @param exp_name Experimental group name (or vector of names).
#' @param ctrl_name Control group name (or vector of names).
#' @param top_n If non-NULL, take the top N genes by ascending `padj` after
#'   significance filtering, and label rows by gene `SYMBOL` when
#'   available. If `NULL` (default), include every gene that passes
#'   thresholds with no row labels.
#' @return Invisibly `NULL` if filtering yields fewer than 2 samples or 2
#'   genes; otherwise the side-effect is the rendered heatmap (printed to
#'   the active graphics device).
#' @details Restricts the count matrix to the comparison's samples (fixes
#'   the prior bug where heatmap columns showed every sample in the full
#'   samplesheet). Z-score normalizes (via `t(scale(t(values)))`) and
#'   renders with a blue-white-red palette via `gplots::heatmap.2`.
#'
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette
#' @export
generate_heatmap <- function(results_df, normalized_counts, sample_info,
                             p = 0.05, lfc = 0.58, exp_name, ctrl_name,
                             top_n = NULL) {
  hm <- .heatmap_matrix(results_df, normalized_counts, sample_info,
                        p, lfc, exp_name, ctrl_name, top_n)
  if (is.null(hm)) return(invisible(NULL))

  gplots::heatmap.2(
    hm$zscores,
    col          = grDevices::colorRampPalette(.diverging_stops())(20),
    density.info = "none",
    dendrogram   = "both",
    Colv         = TRUE,
    trace        = "none",
    margins      = c(10, 10),
    labRow       = hm$row_labels,
    cexRow       = 0.7,
    cexCol       = 0.9
  )
}


#' Build the z-score matrix and row labels for `generate_heatmap`
#'
#' @description Internal helper extracted from [generate_heatmap()] for
#'   testability. Performs the significance filter, restricts to the
#'   comparison's exp+ctrl samples (per the Phase 1.2 sample-subset fix),
#'   and z-score normalizes.
#'
#' @inheritParams generate_heatmap
#' @return `NULL` if filtering yields fewer than 2 samples or 2 genes;
#'   otherwise a list with `zscores` (the z-scored matrix, Genes x
#'   comparison_samples) and `row_labels` (gene-symbol labels for the
#'   top-N variant, or `NA` for the all-significant variant).
#' @keywords internal
.heatmap_matrix <- function(results_df, normalized_counts, sample_info,
                            p = 0.05, lfc = 0.58, exp_name, ctrl_name,
                            top_n = NULL) {
  filtered_data <- results_df |>
    dplyr::filter(!is.na(.data$padj) & .data$padj < p &
                    abs(.data$log2FoldChange) >= lfc) |>
    dplyr::arrange(.data$padj)

  if (!is.null(top_n) && nrow(filtered_data) > top_n) {
    filtered_data <- filtered_data |> dplyr::slice_head(n = top_n)
  }

  gene_ids <- if ("ENSEMBL_ID" %in% colnames(filtered_data)) {
    as.character(filtered_data$ENSEMBL_ID)
  } else {
    rownames(filtered_data)
  }

  comparison_samples <- sample_info$sample[
    sample_info$condition %in% c(exp_name, ctrl_name)
  ]
  comparison_samples <- make.names(as.character(comparison_samples))
  comparison_samples <- intersect(comparison_samples, colnames(normalized_counts))

  if (length(comparison_samples) < 2 || nrow(filtered_data) < 2) {
    cli::cli_alert_warning(".heatmap_matrix: fewer than 2 samples or genes after filtering; skipping")
    return(NULL)
  }

  gene_ids <- intersect(gene_ids, rownames(normalized_counts))
  if (length(gene_ids) < 2) {
    cli::cli_alert_warning(".heatmap_matrix: fewer than 2 genes match the count matrix; skipping")
    return(NULL)
  }

  values <- normalized_counts[gene_ids, comparison_samples, drop = FALSE] |>
    as.matrix() |>
    jitter(factor = 1, amount = 0.00001)

  # Relabel sample columns with display names (DisplayName column when
  # present; munged ID otherwise, so default output is unchanged). The
  # column SUBSET above used the munged join key; only the visible labels
  # change here.
  disp <- .display_lookup(sample_info)
  new_cols <- unname(disp[colnames(values)])
  new_cols[is.na(new_cols)] <- colnames(values)[is.na(new_cols)]
  colnames(values) <- new_cols

  zscores <- t(scale(t(values)))

  row_labels <- NA
  if (!is.null(top_n) && "SYMBOL" %in% colnames(filtered_data)) {
    sym_keys <- if ("ENSEMBL_ID" %in% colnames(filtered_data)) {
      as.character(filtered_data$ENSEMBL_ID)
    } else {
      rownames(filtered_data)
    }
    sym_lookup <- stats::setNames(filtered_data$SYMBOL, sym_keys)
    sym_for_rows <- sym_lookup[gene_ids]
    row_labels <- ifelse(is.na(sym_for_rows) | sym_for_rows == "",
                         gene_ids, sym_for_rows)
  }

  list(zscores = zscores, row_labels = row_labels)
}
