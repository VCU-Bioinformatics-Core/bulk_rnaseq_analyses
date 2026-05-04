# IO functions: counts/samplesheet ingestion and contrast parsing.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; these
# package copies become canonical at Phase 5.8 when the parent script
# becomes a thin optparse -> bisrDE::run_pipeline() wrapper.

#' Read a merged gene-counts TSV
#'
#' @description Read an nf-core/rnaseq-style merged gene-counts file into a
#'   numeric count matrix. Strips Ensembl version suffixes from rownames so
#'   downstream annotation joins are stable across `org.*.eg.db` versions.
#'
#' @param path Path to the counts file (TSV with a header and gene-ID first
#'   column; typically `rsem.merged.gene_counts.tsv`).
#' @return A `data.frame` of integer-coercible counts with gene IDs as
#'   rownames and one column per sample. Non-numeric columns from the input
#'   (e.g. `gene_name`) are dropped.
#' @details Equivalent to the inline block at the parent `de.R` data-load
#'   step. The Ensembl-suffix strip (`"\\..*"` on rownames) is intentional —
#'   AnnotationDbi keys do not carry version suffixes, so `ENSG00000123.4`
#'   would silently fail to map. Sample-level QC (alignment to samplesheet
#'   IDs) is the caller's responsibility — see [align_counts_to_samplesheet()].
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr select
#' @importFrom stringr str_remove
#' @importFrom tidyselect where
#' @export
read_counts <- function(path) {
  counts <- data.frame(
    readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE),
    row.names = 1
  ) |>
    dplyr::select(tidyselect::where(is.numeric))

  rownames(counts) <- stringr::str_remove(rownames(counts), "\\..*")
  counts
}


#' Read a comma-delimited samplesheet
#'
#' @description Read the BISR per-run samplesheet into a data frame. The
#'   samplesheet must have at least three columns: `SampleID`, `GroupID`,
#'   then one column per contrast (1 = experimental, 0 = control, NA =
#'   excluded from that contrast).
#'
#' @param path Path to the samplesheet CSV.
#' @return A `data.frame` with `stringsAsFactors = FALSE`, suitable as input
#'   to [parse_contrasts()] and [align_counts_to_samplesheet()].
#' @details Locked to comma-delimited input by design (matches how
#'   `run_analysis.sh` invokes the pipeline). If TSV samplesheets become a
#'   real need, extend the API rather than auto-detecting.
#'
#' @importFrom utils read.delim
#' @export
read_samplesheet <- function(path) {
  utils::read.delim(
    path,
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE
  )
}


#' Parse contrast columns into a comparisons list
#'
#' @description Convert a samplesheet's contrast columns into a list of
#'   `list(name, exp, ctrl)` entries usable by the DE driver. Each contrast
#'   column encodes one comparison: `1` marks samples in the experimental
#'   group, `0` marks samples in the control group, `NA` excludes the sample
#'   from that contrast.
#'
#' @param samplesheet Data frame returned by [read_samplesheet()].
#' @param group_col Column name holding the GroupID labels used to populate
#'   `exp`/`ctrl` in the output (default `"GroupID"`).
#' @param first_contrast_col Index of the first contrast column. Defaults to
#'   `3` to skip the conventional `SampleID`/`GroupID` columns.
#' @return A list. Each element is `list(name = <contrast column name>,
#'   exp = <unique experimental GroupID(s)>, ctrl = <unique control
#'   GroupID(s)>)`. Contrasts with an empty exp or ctrl group are skipped
#'   with a `cli::cli_alert_warning`.
#'
#' @export
parse_contrasts <- function(samplesheet,
                            group_col = "GroupID",
                            first_contrast_col = 3) {
  contrast_cols <- colnames(samplesheet)[first_contrast_col:ncol(samplesheet)]
  cli::cli_inform("Contrast columns found: {.val {contrast_cols}}")

  comparisons <- list()
  group_vec <- samplesheet[[group_col]]

  for (contrast_name in contrast_cols) {
    contrast_data <- samplesheet[[contrast_name]]

    cli::cli_h2("Contrast: {.strong {contrast_name}}")

    exp_group  <- unique(group_vec[contrast_data == 1 & !is.na(contrast_data)])
    ctrl_group <- unique(group_vec[contrast_data == 0 & !is.na(contrast_data)])

    cli::cli_inform(c(
      "*" = "Experimental: {.val {exp_group}}",
      "*" = "Control:      {.val {ctrl_group}}"
    ))

    if (length(exp_group) > 0 && length(ctrl_group) > 0) {
      comparisons[[length(comparisons) + 1]] <- list(
        name = contrast_name,
        exp  = exp_group,
        ctrl = ctrl_group
      )
    } else {
      cli::cli_alert_warning(
        "Skipping {.strong {contrast_name}}: missing exp or ctrl group"
      )
    }
  }

  cli::cli_alert_success("Built {length(comparisons)} comparison{?s}")
  comparisons
}


#' Align a count matrix to a samplesheet
#'
#' @description Subset and reorder the columns of a count matrix to match the
#'   samples declared in the samplesheet. The samplesheet is the source of
#'   truth — if it lists samples missing from the counts, this errors; if
#'   the counts carry extra samples, those are dropped with a `cli` info
#'   message.
#'
#' @param counts A count `data.frame` or matrix (Genes x Samples), typically
#'   from [read_counts()].
#' @param samplesheet Data frame from [read_samplesheet()].
#' @param sample_col Column in `samplesheet` holding sample identifiers
#'   (default `"SampleID"`).
#' @return The input `counts` reordered and subset to follow `samplesheet`
#'   row order, with `make.names()`-munged sample IDs as column names.
#' @details `read_counts()` indirectly applies `make.names()` to sample IDs
#'   (via `data.frame(..., row.names = 1)` on tibble-style input), so this
#'   function applies the same transform to the samplesheet IDs before
#'   intersecting. Without this alignment, `DESeq2::DESeqDataSetFromMatrix`
#'   fails with `ncol(countData) == nrow(colData) is not TRUE` when the
#'   upstream merged-counts file carries samples beyond the scope of the
#'   current DE run.
#'
#' @export
align_counts_to_samplesheet <- function(counts,
                                        samplesheet,
                                        sample_col = "SampleID") {
  ss_samples_munged <- make.names(as.character(samplesheet[[sample_col]]))
  counts_samples    <- colnames(counts)
  missing_in_counts <- setdiff(ss_samples_munged, counts_samples)
  extra_in_counts   <- setdiff(counts_samples, ss_samples_munged)

  if (length(missing_in_counts) > 0) {
    cli::cli_alert_danger(
      "Samplesheet has {length(missing_in_counts)} sample{?s} not in counts: {.val {missing_in_counts}}"
    )
    cli::cli_abort("Counts / samplesheet sample mismatch (missing in counts)")
  }
  if (length(extra_in_counts) > 0) {
    cli::cli_alert_info(
      "Counts has {length(extra_in_counts)} extra sample{?s} not in samplesheet; dropping from analysis"
    )
  }

  counts[, ss_samples_munged, drop = FALSE]
}
