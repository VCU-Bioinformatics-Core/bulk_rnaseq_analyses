# Dataset-aware QC summary ("QC at a glance") for the report.

#' Compute a short, data-driven QC verdict for the run
#'
#' @description Turns the sample-exploration figures into numbers a reader
#'   can act on: sequencing-depth range, detected-gene range, variance
#'   captured by PC1/PC2, whether the groups separate on PC1-PC2, and which
#'   samples sit unusually far from their own group. The report prints the
#'   returned `text`; the pieces are also stored in the run JSON.
#'
#' @param counts Raw Genes x Samples count matrix (post alignment; the same
#'   one the library-size barplot uses).
#' @param qc_mat Log-scale expression matrix from [qc_matrix()].
#' @param sample_info Data frame with `sample` and `condition` columns.
#' @param transform Label of the transform behind `qc_mat` (`"vst"` or
#'   `"log2_tmm_cpm"`); read from `attr(qc_mat, "transform")` when `NULL`.
#' @return A list: `libsize_M` (named numeric) and `libsize_M_range`,
#'   `depth_fold` (max/min), `detected` (named integer) and
#'   `detected_range`, `pc_var` (numeric, PC1 and PC2 percent variance),
#'   `n_groups`, `separation` (ratio of mean between-group centroid distance
#'   to mean within-group distance on PC1-PC2; `NA` with one group or no
#'   estimable within-group spread), `silhouette` (mean silhouette width on
#'   PC1-PC2 using group labels; `NA` when not computable), `outliers`
#'   (character, see Details), `outlier_ratio` (named numeric, each flagged
#'   sample's nearest-replicate distance over the typical one), `transform`,
#'   `n_genes`, and `text` (one paragraph).
#' @details Separation is descriptive, not a test: a ratio well above 1
#'   means the groups are farther apart than replicates are from each other.
#'   Outliers are scored by the distance to the NEAREST replicate of the same
#'   group on PC1-PC2 (not to the group centroid, which an outlier drags
#'   toward itself so that in a 3-vs-3 design nothing could ever be flagged).
#'   A sample is flagged when that distance is more than 3 robust SDs
#'   (median + 3 MAD over all samples) beyond the typical nearest-replicate
#'   distance AND at least twice the typical distance. In a group of two the
#'   pair simply disagrees and both samples are flagged. Flags are something
#'   to inspect, not a reason to drop a sample.
#' @export
qc_at_a_glance <- function(counts, qc_mat, sample_info, transform = NULL) {
  if (is.null(transform)) transform <- attr(qc_mat, "transform")
  if (is.null(transform)) transform <- "unknown"

  counts <- as.matrix(counts)
  libsize_M <- colSums(counts) / 1e6
  detected  <- colSums(counts > 0)

  grp <- .group_lookup(sample_info)
  m <- as.matrix(qc_mat)
  m <- m[rowSums(m != 0) > 0, , drop = FALSE]
  pca <- stats::prcomp(t(m))
  var_pct <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  pcs <- pca$x[, seq_len(min(2, ncol(pca$x))), drop = FALSE]
  if (ncol(pcs) < 2) pcs <- cbind(pcs, 0)
  groups <- unname(grp[rownames(pcs)])
  groups[is.na(groups)] <- "Unknown"

  centroids <- do.call(rbind, lapply(split(as.data.frame(pcs), groups), colMeans))
  within <- sqrt(rowSums((pcs - centroids[groups, , drop = FALSE])^2))
  n_groups <- nrow(centroids)
  separation <- NA_real_
  silhouette <- NA_real_
  spread_ok  <- FALSE
  if (n_groups >= 2) {
    between <- stats::dist(centroids)
    w <- mean(within[within > 0])
    spread_ok <- is.finite(w) && w > 0
    if (spread_ok) separation <- round(mean(between) / w, 2)
    if (requireNamespace("cluster", quietly = TRUE) && length(unique(groups)) < length(groups)) {
      sil <- tryCatch(cluster::silhouette(as.integer(factor(groups)), stats::dist(pcs)),
                      error = function(e) NULL)
      if (!is.null(sil)) silhouette <- round(mean(sil[, "sil_width"]), 2)
    }
  }

  # Outliers: distance to the nearest replicate of the same group. The group
  # centroid is dragged toward an outlier, so with 3 per group a centroid rule
  # could never flag anything; the nearest-replicate distance is large only
  # for the sample that disagrees with all of its mates.
  d_all <- as.matrix(stats::dist(pcs))
  nn <- vapply(seq_len(nrow(pcs)), function(i) {
    mates <- which(groups == groups[i]); mates <- mates[mates != i]
    if (length(mates) == 0) NA_real_ else min(d_all[i, mates])
  }, numeric(1))
  outliers <- character(0)
  outlier_ratio <- numeric(0)
  med <- stats::median(nn, na.rm = TRUE)
  mad <- stats::mad(nn, na.rm = TRUE)
  if (is.finite(med) && is.finite(mad) && mad > 0 && med > 0) {
    flag <- !is.na(nn) & nn > med + 3 * mad & nn > 2 * med
    outliers <- rownames(pcs)[flag]
    outlier_ratio <- stats::setNames(round(nn[flag] / med, 1), outliers)
  }

  disp <- .display_lookup(sample_info)
  lab <- function(x) { d <- unname(disp[x]); ifelse(is.na(d), x, d) }
  depth_ratio <- round(max(libsize_M) / max(min(libsize_M), 1e-9), 1)
  text <- paste0(
    sprintf("Library sizes range from %.1f to %.1f million reads (%.1f-fold spread%s); ",
            min(libsize_M), max(libsize_M), depth_ratio,
            if (depth_ratio > 3) ", above the 3-fold level at which depth can drive clustering" else ""),
    sprintf("%s to %s genes are detected per sample. ",
            format(min(detected), big.mark = ","), format(max(detected), big.mark = ",")),
    sprintf("PCA on %s (all %s genes): PC1 explains %.1f%% and PC2 %.1f%% of the variance. ",
            switch(transform, vst = "blind variance-stabilised counts",
                   log2_tmm_cpm = "log2(TMM-CPM + 1)", "the QC expression matrix"),
            format(nrow(m), big.mark = ","), var_pct[1], var_pct[2]),
    if (n_groups < 2) "Only one group, so group separation is not assessed. "
    else if (!spread_ok) "Replicates show no within-group spread on PC1-PC2 (singleton groups or identical profiles), so group separation cannot be scored. "
    else
      sprintf("Group centroids are %.1fx farther apart than replicates are from their own centroid on PC1-PC2%s, so the groups %s. ",
              separation,
              if (!is.na(silhouette)) sprintf(" (mean silhouette %.2f)", silhouette) else "",
              if (separation >= 2) "separate clearly" else if (separation >= 1) "separate modestly" else "do not separate on the leading components"),
    if (length(outliers) == 0) "No sample sits unusually far from every replicate of its group."
    else sprintf("Sample%s %s %s far from every replicate of the same group (nearest replicate %s the typical distance): inspect before interpreting, do not drop on this alone.",
                 if (length(outliers) > 1) "s" else "", paste(lab(outliers), collapse = ", "),
                 if (length(outliers) > 1) "sit" else "sits",
                 paste0(paste(outlier_ratio, collapse = "x / "), "x"))
  )

  list(libsize_M = libsize_M, libsize_M_range = round(range(libsize_M), 2), depth_fold = depth_ratio,
       detected = detected, detected_range = as.integer(range(detected)),
       pc_var = var_pct[1:2], n_groups = n_groups,
       separation = separation, silhouette = silhouette,
       outliers = outliers, outlier_ratio = outlier_ratio,
       transform = transform, n_genes = nrow(m), text = text)
}
