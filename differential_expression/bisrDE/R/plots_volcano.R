# Volcano plot of DESeq2 results.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.4 spec.

#' Volcano plot of differential expression results
#'
#' @description Render a volcano plot from an annotated DESeq2 results
#'   data frame. Color, y-axis, threshold line, and label-selection all use
#'   the same significance field (`sig`, default `"padj"`) so the dotted
#'   threshold line aligns with the grey/colored boundary.
#'
#' @param data Annotated DESeq2 results data frame (output of
#'   [annotate_results()]). Must contain `log2FoldChange`, `ENSEMBL_ID`,
#'   `SYMBOL`, and the `sig` column.
#' @param exp_name Experimental group name (rendered in the title).
#' @param ctrl_name Control group name (rendered in the title).
#' @param p P-value threshold for significance. Default `0.05`.
#' @param lfc Log2 fold change threshold. Default `0.58` (≈1.5x).
#' @param sig Column name for significance values. Default `"padj"`. The
#'   same column is used for the colored category, the y-axis transform,
#'   the threshold line, and the top-genes-label selection.
#' @return A `ggplot` object.
#' @details Top 20 up- and down-regulated genes (by ascending `sig`) are
#'   labeled with their `SYMBOL`. Color tags: "Over expressed" (firebrick),
#'   "Under expressed" (steelblue), "Not significant" (grey70).
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggplot2 geom_point geom_hline geom_vline
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_minimal theme element_blank labs
#' @importFrom ggrepel geom_label_repel
#' @export
generate_volcano <- function(data, exp_name, ctrl_name, p = 0.05, lfc = 0.58,
                             sig = "padj") {
  labeled_dat <- data |>
    dplyr::mutate(
      color_tag = dplyr::case_when(
        !is.na(.data[[sig]]) & .data[[sig]] < p & .data$log2FoldChange <= -lfc ~ "Under expressed",
        !is.na(.data[[sig]]) & .data[[sig]] < p & .data$log2FoldChange >=  lfc ~ "Over expressed",
        TRUE ~ "Not significant"
      ),
      color_tag = factor(.data$color_tag,
                         levels = c("Over expressed", "Under expressed", "Not significant"))
    )

  top_genes_up <- labeled_dat |>
    dplyr::filter(!is.na(.data[[sig]]) & .data[[sig]] < p &
                    .data$log2FoldChange >= lfc) |>
    dplyr::arrange(.data[[sig]]) |>
    dplyr::slice_head(n = 20) |>
    dplyr::pull("ENSEMBL_ID")

  top_genes_dn <- labeled_dat |>
    dplyr::filter(!is.na(.data[[sig]]) & .data[[sig]] < p &
                    .data$log2FoldChange <= -lfc) |>
    dplyr::arrange(.data[[sig]]) |>
    dplyr::slice_head(n = 20) |>
    dplyr::pull("ENSEMBL_ID")

  labeled_dat |>
    dplyr::mutate(
      highlight = ifelse(.data$ENSEMBL_ID %in% c(top_genes_up, top_genes_dn),
                         .data$SYMBOL, NA_character_)
    ) |>
    ggplot(aes(
      x = .data$log2FoldChange,
      y = -log10(.data[[sig]]),
      color = .data$color_tag
    )) +
    geom_point(size = 2.5, alpha = 0.6) +
    theme_minimal(base_size = 14) +
    ggrepel::geom_label_repel(
      aes(label = .data$highlight),
      max.overlaps = Inf, show.legend = FALSE,
      size = 4, segment.color = "grey50"
    ) +
    scale_color_manual(values = c(
      "Over expressed"  = "#D55E00",
      "Under expressed" = "#0072B2",
      "Not significant" = "grey70"
    )) +
    geom_hline(yintercept = -log10(p), col = "red", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc), col = "black", linetype = 2) +
    theme(legend.title = element_blank()) +
    labs(
      x = "Log2 Fold-Change (FC)",
      y = paste0("-log10( ", sig, " )"),
      title = paste0("Differentially expressed genes - ",
                     exp_name, " vs. ", ctrl_name)
    )
}
