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


#' Drop excluded samples / groups from a samplesheet
#'
#' @description Remove rows from a samplesheet so they are dropped from the
#'   entire analysis (contrasts, normalization, DESeq2, plots). Three
#'   exclusion sources are combined (union):
#'   \enumerate{
#'     \item An optional declarative `Exclude` column — any row whose value
#'           is `1` / `TRUE` / `"yes"` (case-insensitive) is dropped. Good
#'           for a permanent "this sample was contaminated" record kept with
#'           the samplesheet.
#'     \item `exclude_samples` — `SampleID`s to drop (ad-hoc, e.g. a CLI
#'           flag) without editing the canonical samplesheet.
#'     \item `exclude_groups` — every sample whose `GroupID` matches is
#'           dropped.
#'   }
#'
#' @param samplesheet Data frame from [read_samplesheet()].
#' @param exclude_samples Optional character vector of `SampleID`s to drop.
#' @param exclude_groups Optional character vector of `GroupID`s to drop.
#' @param sample_col Sample-ID column name (default `"SampleID"`).
#' @param group_col Group-ID column name (default `"GroupID"`).
#' @return The samplesheet with excluded rows removed (row order otherwise
#'   preserved). Errors via `cli::cli_abort` if fewer than two samples
#'   remain.
#' @details The `SampleID` join key is untouched — exclusion only removes
#'   rows, so downstream `make.names()`-based alignment to the count matrix
#'   is unaffected for the surviving samples.
#'
#' @export
filter_samplesheet <- function(samplesheet,
                               exclude_samples = NULL,
                               exclude_groups  = NULL,
                               sample_col = "SampleID",
                               group_col  = "GroupID") {
  ss <- samplesheet
  n0 <- nrow(ss)
  drop <- rep(FALSE, n0)

  # 1) Declarative Exclude column.
  if ("Exclude" %in% colnames(ss)) {
    ex <- ss[["Exclude"]]
    ex_chr <- trimws(tolower(as.character(ex)))
    ex_drop <- !is.na(ex) & ex_chr %in% c("1", "true", "yes", "y", "t")
    if (any(ex_drop)) {
      cli::cli_alert_info(
        "Exclude column: dropping {sum(ex_drop)} sample{?s}: {.val {as.character(ss[[sample_col]][ex_drop])}}"
      )
    }
    drop <- drop | ex_drop
  }

  # 2) --exclude-samples (by SampleID).
  if (length(exclude_samples) > 0) {
    s_drop <- as.character(ss[[sample_col]]) %in% as.character(exclude_samples)
    matched <- intersect(as.character(exclude_samples),
                         as.character(ss[[sample_col]]))
    missing <- setdiff(as.character(exclude_samples),
                       as.character(ss[[sample_col]]))
    if (length(matched) > 0) {
      cli::cli_alert_info("exclude-samples: dropping {.val {matched}}")
    }
    if (length(missing) > 0) {
      cli::cli_alert_warning(
        "exclude-samples: {length(missing)} ID{?s} not in samplesheet: {.val {missing}}"
      )
    }
    drop <- drop | s_drop
  }

  # 3) --exclude-groups (by GroupID).
  if (length(exclude_groups) > 0) {
    g_drop <- as.character(ss[[group_col]]) %in% as.character(exclude_groups)
    matched <- intersect(as.character(exclude_groups),
                         as.character(ss[[group_col]]))
    missing <- setdiff(as.character(exclude_groups),
                       as.character(ss[[group_col]]))
    if (length(matched) > 0) {
      cli::cli_alert_info("exclude-groups: dropping group{?s} {.val {matched}}")
    }
    if (length(missing) > 0) {
      cli::cli_alert_warning(
        "exclude-groups: {length(missing)} group{?s} not in samplesheet: {.val {missing}}"
      )
    }
    drop <- drop | g_drop
  }

  ss <- ss[!drop, , drop = FALSE]

  if (sum(drop) > 0) {
    cli::cli_alert_success(
      "Sample exclusion: kept {nrow(ss)} / {n0} sample{?s}"
    )
  }
  if (nrow(ss) < 2) {
    cli::cli_abort(
      "After exclusions only {nrow(ss)} sample{?s} remain (need >= 2)."
    )
  }

  ss
}


#' Metadata column names recognised in a samplesheet
#'
#' @description The reserved, non-contrast columns a samplesheet may carry.
#'   Everything from `first_contrast_col` onward that is NOT one of these is
#'   treated as a contrast column. Centralised so [parse_contrasts()] and
#'   [run_pipeline()] agree on what counts as metadata.
#' @keywords internal
.samplesheet_meta_cols <- c("SampleID", "GroupID", "Exclude", "DisplayName")


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
#' @param first_contrast_col Index of the first candidate contrast column.
#'   Defaults to `3` to skip the conventional `SampleID`/`GroupID` columns.
#'   Reserved metadata columns ([.samplesheet_meta_cols]: `SampleID`,
#'   `GroupID`, `Exclude`, `DisplayName`) at or after this index are still
#'   skipped by name, so an `Exclude`/`DisplayName` column inserted between
#'   `GroupID` and the first contrast does not get mistaken for a contrast.
#' @param include_contrasts Optional character vector. When non-NULL, only
#'   these contrast columns are processed (allowlist).
#' @param exclude_contrasts Optional character vector of contrast columns to
#'   skip (denylist). Applied after `include_contrasts`.
#' @return A list. Each element is `list(name = <contrast column name>,
#'   exp = <unique experimental GroupID(s)>, ctrl = <unique control
#'   GroupID(s)>)`. Contrasts with an empty exp or ctrl group are skipped
#'   with a `cli::cli_alert_warning`.
#'
#' @export
parse_contrasts <- function(samplesheet,
                            group_col = "GroupID",
                            first_contrast_col = 3,
                            include_contrasts = NULL,
                            exclude_contrasts = NULL) {
  contrast_cols <- colnames(samplesheet)[first_contrast_col:ncol(samplesheet)]
  # Drop reserved metadata columns by name (preserving order) so an
  # Exclude / DisplayName column sitting before the first contrast is not
  # treated as a contrast.
  contrast_cols <- setdiff(contrast_cols,
                           union(.samplesheet_meta_cols, group_col))

  if (!is.null(include_contrasts)) {
    keep <- intersect(contrast_cols, include_contrasts)
    dropped <- setdiff(contrast_cols, keep)
    if (length(dropped) > 0) {
      cli::cli_alert_info(
        "include-contrasts: keeping {length(keep)} contrast{?s}, skipping {.val {dropped}}"
      )
    }
    missing_req <- setdiff(include_contrasts, contrast_cols)
    if (length(missing_req) > 0) {
      cli::cli_alert_warning(
        "include-contrasts: requested contrast{?s} not in samplesheet: {.val {missing_req}}"
      )
    }
    contrast_cols <- keep
  }
  if (!is.null(exclude_contrasts)) {
    hit <- intersect(contrast_cols, exclude_contrasts)
    if (length(hit) > 0) {
      cli::cli_alert_info("exclude-contrasts: skipping {.val {hit}}")
    }
    contrast_cols <- setdiff(contrast_cols, exclude_contrasts)
  }

  cli::cli_inform("Contrast columns found: {.val {contrast_cols}}")

  comparisons <- list()
  group_vec <- samplesheet[[group_col]]

  for (contrast_name in contrast_cols) {
    contrast_data <- samplesheet[[contrast_name]]

    cli::cli_h2("Contrast: {.strong {contrast_name}}")

    # Diagnostic: only 1 (exp), 0 (ctrl) and blank/NA (exclude) are
    # recognised. Flag anything else (e.g. "yes"/"no") so a mis-encoded
    # column doesn't silently look like "missing exp or ctrl group".
    vals <- contrast_data[!is.na(contrast_data) &
                            nzchar(trimws(as.character(contrast_data)))]
    unexpected <- setdiff(unique(trimws(as.character(vals))), c("0", "1"))
    if (length(unexpected) > 0) {
      cli::cli_alert_warning(
        "Contrast {.strong {contrast_name}} has non-0/1 value{?s} {.val {unexpected}}; only 1 (exp), 0 (ctrl) and blank (exclude) are recognised"
      )
    }

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
