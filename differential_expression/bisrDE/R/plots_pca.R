# PCA plots (static ggplot, interactive 2D plotly, interactive 3D plotly).
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.4 spec.
#
# These functions extract logic that lived inline in the parent's main
# block (de.R lines 1740-1850 in the post-Phase-4 layout). Each plotter
# does its own prcomp internally so callers don't need to thread state.

#' Compute prcomp on a TMM matrix (internal)
#'
#' @description Drop all-zero rows, transpose so samples are rows, run
#'   `stats::prcomp` and capture the percent-variance vector.
#'
#' @param tmm Genes x Samples numeric matrix (typically [run_tmm()]).
#' @param rank Optional integer passed to `prcomp(rank. = ...)`. The
#'   3D plot uses `rank = 3`; the 2D path leaves it `NULL`.
#' @return List with elements `prcomp` (the `prcomp` object) and
#'   `var_pct` (numeric, percent variance per component, rounded to 1 dp).
#' @keywords internal
.compute_pca <- function(tmm, rank = NULL) {
  data_reduced_zero <- as.data.frame(tmm) |>
    dplyr::filter(dplyr::if_any(tidyselect::where(is.numeric)))

  pca_matrix <- data_reduced_zero |>
    as.matrix() |>
    t()

  pca <- if (is.null(rank)) {
    stats::prcomp(pca_matrix)
  } else {
    stats::prcomp(pca_matrix, rank. = rank)
  }

  pca_var <- pca$sdev^2
  var_pct <- round(pca_var / sum(pca_var) * 100, 1)

  list(prcomp = pca, var_pct = var_pct)
}


#' Static PCA scatter (PC1 vs PC2)
#'
#' @description Render a 2D PCA scatter plot via ggplot2. Each point is a
#'   sample, labeled by its `sample` ID and colored by `condition`.
#'
#' @param tmm Genes x Samples numeric matrix (e.g. [run_tmm()] output).
#' @param sample_info Data frame with `sample` and `condition` columns.
#'   Row order must match the column order of `tmm`.
#' @return A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot aes geom_text xlab ylab theme_bw
#' @export
pca_static <- function(tmm, sample_info) {
  pca <- .compute_pca(tmm)

  pca_df <- data.frame(
    Sample = rownames(pca$prcomp$x),
    X      = pca$prcomp$x[, 1],
    Y      = pca$prcomp$x[, 2],
    Group  = sample_info$condition,
    stringsAsFactors = FALSE
  )

  ggplot(pca_df, aes(x = .data$X, y = .data$Y,
                     label = .data$Sample, color = .data$Group)) +
    geom_text() +
    xlab(paste0("PC1 - ", pca$var_pct[1], "%")) +
    ylab(paste0("PC2 - ", pca$var_pct[2], "%")) +
    theme_bw()
}


#' Interactive 2D PCA scatter
#'
#' @description Render an interactive 2D PCA scatter via `plotly::plot_ly`.
#'   Each point is a sample, colored by `condition`. Hover reveals the
#'   sample ID.
#'
#' @param tmm Genes x Samples numeric matrix (e.g. [run_tmm()] output).
#' @param sample_info Data frame with `sample` and `condition` columns.
#' @return A `plotly` object.
#'
#' @export
pca_plotly <- function(tmm, sample_info) {
  pca <- .compute_pca(tmm)

  pca_df <- data.frame(
    Sample = rownames(pca$prcomp$x),
    X      = pca$prcomp$x[, 1],
    Y      = pca$prcomp$x[, 2],
    Group  = sample_info$condition,
    stringsAsFactors = FALSE
  )

  fig <- plotly::plot_ly(
    pca_df,
    x      = ~X,
    y      = ~Y,
    color  = ~Group,
    colors = c("steelblue", "firebrick", "olivedrab", "plum"),
    type   = "scatter",
    mode   = "markers",
    size   = 5.7,
    span   = 4.9,
    alpha  = 2.5,
    alpha_stroke = 0.5
  )

  plotly::layout(
    fig,
    legend       = list(title = list(text = "color")),
    plot_bgcolor = "#e5ecf6",
    xaxis = list(
      title         = paste0("PC1 - ", pca$var_pct[1], "%"),
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor     = "#ffff"
    ),
    yaxis = list(
      title         = paste0("PC2 - ", pca$var_pct[2], "%"),
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor     = "#ffff"
    )
  )
}


#' Interactive 3D PCA scatter
#'
#' @description Render an interactive 3D PCA scatter via `plotly::plot_ly`.
#'   Computes a separate `prcomp(rank. = 3)` (matching the parent
#'   pipeline's behavior — the 2D and 3D paths use independent prcomp
#'   calls). PC2 and PC3 axes are sign-flipped for visualization
#'   (preserves the parent pipeline's view orientation).
#'
#' @param tmm Genes x Samples numeric matrix (e.g. [run_tmm()] output).
#' @param sample_info Data frame with `sample` and `condition` columns.
#' @return A `plotly` object with a title showing total explained
#'   variance across the three retained components.
#'
#' @export
pca_plotly_3d <- function(tmm, sample_info) {
  pca <- .compute_pca(tmm, rank = 3)
  prin_comp <- pca$prcomp

  components <- as.data.frame(prin_comp$x)
  components$PC2 <- -components$PC2
  components$PC3 <- -components$PC3
  components$Group <- sample_info$condition

  tot_explained <- summary(prin_comp)[["importance"]]["Proportion of Variance", ]
  tot_explained <- 100 * sum(tot_explained)
  tit <- paste0("Total Explained Variance = ", tot_explained)

  fig <- plotly::plot_ly(
    components,
    x      = ~PC1,
    y      = ~PC2,
    z      = ~PC3,
    color  = ~Group,
    colors = c("steelblue", "firebrick", "olivedrab", "plum")
  )
  fig <- plotly::add_markers(fig, size = 12)

  plotly::layout(
    fig,
    title = tit,
    scene = list(bgcolor = "#e5ecf6")
  )
}
