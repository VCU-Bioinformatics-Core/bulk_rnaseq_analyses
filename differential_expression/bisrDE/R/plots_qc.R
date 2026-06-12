# Sample-level QC plots: correlation heatmap, vst distance heatmap,
# library-size + detected-genes barplot, hclust dendrogram + log-CPM
# density. These act on the full TMM-normalized count matrix + the
# samplesheet (NOT per-comparison subsets) and are intended for the
# report's "Sample Exploration" section.
#
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.4 spec. Notable refactor:
# `qc_vst_dist_heatmap` now delegates to [run_vst()] (migrated in Phase
# 5.2) instead of inlining the vst+fallback pattern.

#' Sample x Sample correlation heatmap
#'
#' @description Spearman (default) correlation heatmap on a chosen gene
#'   set, with a top color bar marking each sample's `condition`.
#'   Samples that correlate poorly with their group are immediately
#'   visible.
#'
#' @param normalized_counts Genes x Samples numeric matrix (typically
#'   [run_tmm()] output). Column names must match the `make.names()`
#'   form of `sample_info$sample`.
#' @param sample_info Data frame with `sample` (canonical IDs) and
#'   `condition` (group label).
#' @param sig_genes Vector of gene IDs (Ensembl) defining the rows used
#'   for the correlation. Typically the union of all DE genes across
#'   comparisons; pass `rownames(normalized_counts)` for an "all genes"
#'   variant.
#' @param fig_path Output PNG path.
#' @param method Correlation method. Default `"spearman"` (robust to
#'   outliers and scale; matches DESeq2-vignette convention).
#' @return Invisibly returns `fig_path` on success, or `NULL` if the
#'   inputs are too small to produce a useful heatmap.
#'
#' @importFrom grDevices png dev.off hcl.colors
#' @export
qc_correlation_heatmap <- function(normalized_counts, sample_info, sig_genes,
                                   fig_path, method = "spearman") {
  if (is.null(sig_genes) || length(sig_genes) < 2) {
    cli::cli_alert_warning("qc_correlation_heatmap: <2 input genes; skipping")
    return(invisible(NULL))
  }

  sig_genes <- intersect(as.character(sig_genes), rownames(normalized_counts))
  if (length(sig_genes) < 2) {
    cli::cli_alert_warning("qc_correlation_heatmap: <2 matching genes after intersect; skipping")
    return(invisible(NULL))
  }

  mat <- normalized_counts[sig_genes, , drop = FALSE]
  if (ncol(mat) < 2) {
    cli::cli_alert_warning("qc_correlation_heatmap: <2 samples; skipping")
    return(invisible(NULL))
  }

  cor_mat <- stats::cor(mat, method = method)

  ss_munged <- make.names(as.character(sample_info$sample))
  group_lookup <- stats::setNames(as.character(sample_info$condition), ss_munged)
  groups <- group_lookup[colnames(cor_mat)]
  groups[is.na(groups)] <- "Unknown"

  # Display labels (DisplayName column when present; munged ID otherwise).
  disp <- .display_lookup(sample_info)
  row_labs <- unname(ifelse(is.na(disp[rownames(cor_mat)]),
                            rownames(cor_mat), disp[rownames(cor_mat)]))
  col_labs <- unname(ifelse(is.na(disp[colnames(cor_mat)]),
                            colnames(cor_mat), disp[colnames(cor_mat)]))

  uniq_groups <- unique(groups)
  group_palette <- stats::setNames(
    grDevices::hcl.colors(length(uniq_groups), palette = "Dark 3"),
    uniq_groups
  )

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Group                = groups,
    col                  = list(Group = group_palette),
    show_legend          = TRUE,
    annotation_name_side = "left"
  )

  hm <- ComplexHeatmap::Heatmap(
    cor_mat,
    name           = sprintf("%s\nrho", tools::toTitleCase(method)),
    top_annotation = ha,
    col            = circlize::colorRamp2(
      breaks = c(min(cor_mat, 0), stats::median(cor_mat), 1),
      colors = c("blue", "white", "firebrick")
    ),
    row_labels           = row_labs,
    column_labels        = col_labs,
    show_row_names       = TRUE,
    show_column_names    = TRUE,
    row_names_gp         = grid::gpar(fontsize = 8),
    column_names_gp      = grid::gpar(fontsize = 8),
    column_title         = sprintf(
      "Sample x Sample %s correlation (%d DE genes)",
      method, length(sig_genes)
    ),
    column_title_gp      = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(direction = "vertical")
  )

  # v1.4.0: -30% from prior 1200x1100 per user feedback (correlation HM
  # was crowding the embed-resources HTML).
  grDevices::png(fig_path, width = 840, height = 770, res = 150)
  ComplexHeatmap::draw(hm, merge_legend = TRUE)
  grDevices::dev.off()

  invisible(fig_path)
}


#' Sample x Sample vst distance heatmap
#'
#' @description Standard DESeq2-vignette QC plot: variance-stabilizing
#'   transform on raw counts, then `dist()` on the transposed matrix to
#'   get sample x sample Euclidean distances, then visualize via
#'   `ComplexHeatmap`. Replicates that don't cluster together flag
#'   potential outliers or mis-labelling.
#'
#' @param dds A `DESeqDataSet` (typically post pre-filtering, before
#'   per-comparison `DESeq()` calls). The vst is `blind = TRUE` so the
#'   transform is design-agnostic.
#' @param sample_info Same shape as in [qc_correlation_heatmap()].
#' @param fig_path Output PNG path.
#' @return Invisibly returns `fig_path` on success, or `NULL` if vst
#'   can't be computed (e.g. <2 samples, all-zero genes, `run_vst`
#'   fallback fails).
#' @details Delegates the vst computation to [run_vst()], which wraps
#'   `DESeq2::vst()` with a `tryCatch` fallback to
#'   `DESeq2::varianceStabilizingTransformation()`. The fast vst path
#'   subsamples 1000 genes and fails on small/sparse datasets; the
#'   fallback is slower but always-correct.
#'
#' @importFrom grDevices png dev.off hcl.colors
#' @export
qc_vst_dist_heatmap <- function(dds, sample_info, fig_path) {
  vst_obj <- run_vst(dds, blind = TRUE)
  if (is.null(vst_obj)) return(invisible(NULL))

  vst_mat <- SummarizedExperiment::assay(vst_obj)
  if (ncol(vst_mat) < 2) {
    cli::cli_alert_warning("qc_vst_dist_heatmap: <2 samples; skipping")
    return(invisible(NULL))
  }

  dist_mat <- as.matrix(stats::dist(t(vst_mat)))

  ss_munged <- make.names(as.character(sample_info$sample))
  group_lookup <- stats::setNames(as.character(sample_info$condition), ss_munged)
  groups <- group_lookup[colnames(dist_mat)]
  groups[is.na(groups)] <- "Unknown"

  # Display labels (DisplayName column when present; munged ID otherwise).
  disp <- .display_lookup(sample_info)
  row_labs <- unname(ifelse(is.na(disp[rownames(dist_mat)]),
                            rownames(dist_mat), disp[rownames(dist_mat)]))
  col_labs <- unname(ifelse(is.na(disp[colnames(dist_mat)]),
                            colnames(dist_mat), disp[colnames(dist_mat)]))

  uniq_groups <- unique(groups)
  group_palette <- stats::setNames(
    grDevices::hcl.colors(length(uniq_groups), palette = "Dark 3"),
    uniq_groups
  )

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Group                = groups,
    col                  = list(Group = group_palette),
    show_legend          = TRUE,
    annotation_name_side = "left"
  )

  hm <- ComplexHeatmap::Heatmap(
    dist_mat,
    name           = "vst dist",
    top_annotation = ha,
    col            = circlize::colorRamp2(
      breaks = c(0, stats::median(dist_mat), max(dist_mat)),
      colors = c("midnightblue", "lightblue", "white")
    ),
    row_labels           = row_labs,
    column_labels        = col_labs,
    show_row_names       = TRUE,
    show_column_names    = TRUE,
    row_names_gp         = grid::gpar(fontsize = 8),
    column_names_gp      = grid::gpar(fontsize = 8),
    column_title         = sprintf(
      "Sample x Sample Euclidean distance (vst-transformed, %d genes)",
      nrow(vst_mat)
    ),
    column_title_gp      = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(direction = "vertical")
  )

  # v1.4.0: -30% from prior 1200x1100 per user feedback.
  grDevices::png(fig_path, width = 840, height = 770, res = 150)
  ComplexHeatmap::draw(hm, merge_legend = TRUE)
  grDevices::dev.off()

  invisible(fig_path)
}


#' Per-sample library size + detected-genes barplot
#'
#' @description Two-panel ggplot stacked vertically via `patchwork`:
#'   library size (in millions of reads) on top, detected genes
#'   (`colSums(counts > 0)`) on bottom. Bars colored by `condition`.
#'   Quick first-look QC for shallow sequencing or low-diversity
#'   libraries.
#'
#' @param counts Genes x Samples raw count matrix (integer or numeric).
#' @param sample_info Data frame with `sample` and `condition` columns.
#' @param fig_path Output PNG path.
#' @return Invisibly returns `fig_path`.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_x_discrete labs theme_minimal theme element_text ggsave
#' @importFrom grDevices hcl.colors
#' @export
qc_libsize_detected_barplot <- function(counts, sample_info, fig_path) {
  if (ncol(counts) < 1) {
    cli::cli_alert_warning("qc_libsize_detected_barplot: empty count matrix; skipping")
    return(invisible(NULL))
  }

  ss_munged <- make.names(as.character(sample_info$sample))
  group_lookup <- stats::setNames(as.character(sample_info$condition), ss_munged)

  # Display-label mapper for the x axis (munged ID -> DisplayName / munged).
  disp <- .display_lookup(sample_info)
  disp_labeller <- function(b) {
    out <- unname(disp[b])
    out[is.na(out)] <- b[is.na(out)]
    out
  }

  per_sample <- data.frame(
    sample    = colnames(counts),
    libsize_M = colSums(counts) / 1e6,
    detected  = colSums(counts > 0),
    stringsAsFactors = FALSE
  )
  per_sample$group <- group_lookup[per_sample$sample]
  per_sample$group[is.na(per_sample$group)] <- "Unknown"
  per_sample$sample <- factor(per_sample$sample, levels = per_sample$sample)

  group_palette <- stats::setNames(
    grDevices::hcl.colors(length(unique(per_sample$group)), palette = "Dark 3"),
    unique(per_sample$group)
  )

  p_lib <- ggplot(per_sample, aes(x = .data$sample, y = .data$libsize_M,
                                  fill = .data$group)) +
    geom_col() +
    scale_fill_manual(values = group_palette) +
    scale_x_discrete(labels = disp_labeller) +
    labs(x = NULL, y = "Library size (M reads)",
         title = "Library size per sample", fill = "Group") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "right")

  p_det <- ggplot(per_sample, aes(x = .data$sample, y = .data$detected,
                                  fill = .data$group)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = group_palette) +
    scale_x_discrete(labels = disp_labeller) +
    labs(x = "Sample", y = "Detected genes (counts > 0)",
         title = "Detected genes per sample") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

  combined <- patchwork::wrap_plots(p_lib, p_det, ncol = 1, heights = c(1, 1)) +
    patchwork::plot_annotation(theme = ggplot2::theme(
      plot.title = element_text(face = "bold")
    ))

  # v1.4.0: -30% from prior 11 x 8 in per user feedback.
  ggsave(fig_path, combined, width = 7.7, height = 5.6, dpi = 150, bg = "white")
  invisible(fig_path)
}


#' Hierarchical clustering dendrogram + log-CPM density per sample
#'
#' @description Two side-by-side panels via `patchwork`:
#'   - Left: Ward dendrogram on `dist(t(log2(tmm + 1)))`. Sample labels
#'     are colored by group so visually-clustered replicates can be
#'     sanity-checked at a glance.
#'   - Right: per-sample density curve of `log2(CPM + 1)`, colored by
#'     group. Distributions should overlap; an outlier curve flags a
#'     library prep / depth / contamination problem.
#'
#' @param tmm TMM-normalized counts (Genes x Samples), e.g. [run_tmm()]
#'   output.
#' @param sample_info Same shape as elsewhere (`sample` + `condition`).
#' @param fig_path Output PNG path.
#' @param dist_method Distance metric for `dist()`. Default `"euclidean"`.
#' @param clust_method Linkage for `hclust()`. Default `"ward.D2"`.
#' @return Invisibly returns `fig_path`.
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggplot2 geom_segment geom_text geom_density
#' @importFrom ggplot2 scale_color_manual coord_cartesian
#' @importFrom ggplot2 labs theme_minimal theme element_blank ggsave
#' @importFrom grDevices hcl.colors
#' @export
qc_hclust_density <- function(tmm, sample_info, fig_path,
                              dist_method  = "euclidean",
                              clust_method = "ward.D2") {
  if (ncol(tmm) < 2) {
    cli::cli_alert_warning("qc_hclust_density: <2 samples; skipping")
    return(invisible(NULL))
  }

  log_tmm <- log2(tmm + 1)

  ss_munged <- make.names(as.character(sample_info$sample))
  group_lookup <- stats::setNames(as.character(sample_info$condition), ss_munged)

  # ---- 1) Hclust dendrogram (via ggdendro) ----
  d  <- stats::dist(t(log_tmm), method = dist_method)
  hc <- stats::hclust(d, method = clust_method)
  ddata <- ggdendro::dendro_data(hc)

  leaf_labs <- ggdendro::label(ddata)
  leaf_labs$group <- group_lookup[as.character(leaf_labs$label)]
  leaf_labs$group[is.na(leaf_labs$group)] <- "Unknown"

  # Display labels (DisplayName when present; munged dendrogram leaf
  # label otherwise). The dendrogram structure stays keyed by the munged
  # IDs; only the visible text changes.
  disp <- .display_lookup(sample_info)
  leaf_labs$display <- unname(disp[as.character(leaf_labs$label)])
  leaf_labs$display[is.na(leaf_labs$display)] <-
    as.character(leaf_labs$label)[is.na(leaf_labs$display)]

  uniq_groups <- unique(leaf_labs$group)
  group_palette <- stats::setNames(
    grDevices::hcl.colors(length(uniq_groups), palette = "Dark 3"),
    uniq_groups
  )

  segs <- ggdendro::segment(ddata)
  y_max <- max(segs$y, na.rm = TRUE)

  p_dend <- ggplot() +
    geom_segment(
      data = segs,
      aes(x = .data$x, y = .data$y,
          xend = .data$xend, yend = .data$yend),
      color = "grey30"
    ) +
    geom_text(
      data = leaf_labs,
      aes(x = .data$x, y = .data$y,
          label = .data$display, color = .data$group),
      hjust = 1, angle = 90, vjust = 0.5, size = 3,
      nudge_y = -y_max * 0.02
    ) +
    scale_color_manual(values = group_palette) +
    coord_cartesian(clip = "off") +
    labs(
      title = sprintf("Hierarchical clustering (%s, %s)", clust_method, dist_method),
      x     = NULL,
      y     = "Distance",
      color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin        = grid::unit(c(5, 5, 30, 5), "pt")
    )

  # ---- 2) Per-sample density of log2(CPM + 1) ----
  long_df <- data.frame(
    sample   = rep(colnames(log_tmm), each = nrow(log_tmm)),
    log2_cpm = as.vector(log_tmm),
    stringsAsFactors = FALSE
  )
  long_df$group <- group_lookup[long_df$sample]
  long_df$group[is.na(long_df$group)] <- "Unknown"

  p_den <- ggplot(long_df, aes(x = .data$log2_cpm,
                               color = .data$group,
                               group = .data$sample)) +
    geom_density(alpha = 0.6) +
    scale_color_manual(values = group_palette) +
    labs(
      title = "log2(CPM + 1) density per sample",
      x     = "log2(CPM + 1)",
      y     = "Density",
      color = "Group"
    ) +
    theme_minimal(base_size = 13)

  combined <- patchwork::wrap_plots(p_dend, p_den, ncol = 2, widths = c(1.2, 1)) +
    patchwork::plot_layout(guides = "collect")

  # v1.4.0: -30% from prior 14 x 7 in per user feedback.
  ggsave(fig_path, combined, width = 9.8, height = 4.9, dpi = 150, bg = "white")
  invisible(fig_path)
}
