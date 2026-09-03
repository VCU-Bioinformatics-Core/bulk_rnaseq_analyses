# IO functions: counts/samplesheet ingestion and contrast parsing.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; these
# package copies become canonical at Phase 5.8 when the parent script
# becomes a thin optparse -> bisrDE::run_pipeline() wrapper.

#' Read a merged gene-counts TSV
#'
#' @description Read an nf-core/rnaseq-style merged gene-counts file into a
#'   numeric count matrix (gene IDs as rownames, one column per sample). The
#'   first column is the gene ID; known metadata columns (`gene_name`,
#'   `transcript_id(s)`, `Chr`/`Start`/`End`/`Strand`/`Length`, ...) are dropped
#'   by name, so salmon / RSEM / featureCounts-style layouts all work as-is.
#'   Strips Ensembl version suffixes so annotation joins are stable across
#'   `org.*.eg.db` versions.
#'
#' @param path Path to the counts file (TSV with a header and gene-ID first
#'   column; typically `rsem.merged.gene_counts.tsv`).
#' @return A `data.frame` of integer-coercible counts with gene IDs as
#'   rownames and one column per sample. Non-numeric columns from the input
#'   (e.g. `gene_name`) are dropped from the matrix but NOT lost: the
#'   attribute `"gene_meta"` is a data frame with one row per kept gene
#'   (`id` = the rowname, `id_versioned` = the original first-column value,
#'   `gene_name` = the input `gene_name`/`gene_symbol`/`symbol` column or
#'   `NA`), so the annotation step can fall back to the input's gene name
#'   when the OrgDb has none and can report the versioned ID. Subsetting the
#'   data frame drops the attribute, so callers read it right after this
#'   function returns (see [run_pipeline()]).
#' @details The Ensembl version-suffix strip (`"\\..*"`) is applied ONLY to
#'   `^ENS` accessions — AnnotationDbi keys carry no version, so `ENSG00000123.4`
#'   would otherwise fail to map, while gene symbols that legitimately contain
#'   dots (e.g. `H2-M10.1`) are left intact. Ensembl `_PAR_Y` pseudo-autosomal
#'   genes (e.g. `ENSG00000002586.20_PAR_Y`) are dropped first: they are copies
#'   of the chrX gene and would collide with it once the version is stripped
#'   (Ensembl documents that they can be ignored). Sample-level QC (alignment to
#'   samplesheet IDs) is the caller's responsibility — see
#'   [align_counts_to_samplesheet()].
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr select
#' @importFrom stringr str_remove
#' @importFrom tidyselect where
#' @export
read_counts <- function(path) {
  raw <- readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  if (ncol(raw) < 2) {
    cli::cli_abort(c(
      "read_counts: {.path {path}} has fewer than 2 columns.",
      "i" = "Expected a gene-ID column followed by >= 1 sample-count column."
    ))
  }

  # Flexible counts-file handling. The first column is the gene ID (nf-core
  # merged files use `gene_id`; featureCounts uses `Geneid`). Known non-count
  # metadata columns are dropped by NAME (case-insensitive) so the different
  # nf-core organizations all reduce to the same "gene x sample" matrix: salmon
  # (`gene_name`), RSEM (`transcript_id(s)`), and featureCounts-style layouts
  # (`Chr`/`Start`/`End`/`Strand`/`Length`) — the last of which would otherwise
  # leak NUMERIC annotation columns in as if they were samples. Whatever numeric
  # columns remain are the per-sample counts.
  meta_cols <- c(
    "gene_id", "geneid", "gene_name", "gene_symbol", "symbol", "name",
    "transcript_id", "transcript_id(s)", "transcript_ids", "tx_id",
    "chr", "chromosome", "start", "end", "strand", "length",
    "gene_biotype", "biotype", "description"
  )
  sample_cols <- names(raw)[-1][!tolower(names(raw)[-1]) %in% meta_cols]
  # Keep the input's own gene name (salmon/RSEM `gene_name`, or a
  # `gene_symbol`/`symbol` column) for the annotation fallback + the
  # versioned original ID for traceability back to the nf-core row.
  name_col <- names(raw)[tolower(names(raw)) %in% c("gene_name", "gene_symbol", "symbol")][1]
  gene_meta <- data.frame(
    id_versioned = as.character(raw[[1]]),
    gene_name    = if (!is.na(name_col)) as.character(raw[[name_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  counts <- as.data.frame(raw[, sample_cols, drop = FALSE])
  counts <- counts[, vapply(counts, is.numeric, logical(1)), drop = FALSE]
  if (ncol(counts) < 1) {
    cli::cli_abort(c(
      "read_counts: no numeric sample columns found in {.path {path}}.",
      "i" = "After the gene-ID and known metadata columns, no numeric columns remained."
    ))
  }
  # `make.names()` the sample columns so they match `make.names(SampleID)` used
  # by align_counts_to_samplesheet() + the display/group lookups. readr keeps
  # the raw header (e.g. "HCT116-P53-minus-1"); downstream expects the munged
  # "HCT116.P53.minus.1" (this is what data.frame()'s check.names used to do).
  colnames(counts) <- make.names(colnames(counts))
  rownames(counts) <- as.character(raw[[1]])

  # Drop Ensembl PAR_Y pseudo-autosomal duplicates (e.g.
  # "ENSG00000002586.20_PAR_Y"). These are copies of the chrX gene in the
  # pseudo-autosomal region and, once the version suffix is stripped below,
  # collide with the chrX copy -> "duplicate 'row.names'". Ensembl documents
  # that they can be ignored, so we remove them.
  par_y <- grepl("_PAR_Y$", rownames(counts))
  if (any(par_y)) {
    cli::cli_alert_info(
      "Dropping {sum(par_y)} Ensembl _PAR_Y pseudo-autosomal duplicate{?s}"
    )
    counts <- counts[!par_y, , drop = FALSE]
    gene_meta <- gene_meta[!par_y, , drop = FALSE]
  }

  # Strip Ensembl version suffixes (e.g. "ENSMUSG00000000001.3" -> "...001"),
  # but ONLY from Ensembl-style accessions. Gene symbols can carry dots that are
  # part of the name (e.g. "H2-M10.1", "Rn4.5s", "Tex19.1" are *distinct*
  # genes), so a blanket "\\..*" strip would collapse them into duplicate row
  # names and error at the rownames<- assignment.
  is_ens <- grepl("^ENS", rownames(counts))
  rownames(counts)[is_ens] <-
    stringr::str_remove(rownames(counts)[is_ens], "\\..*")
  gene_meta <- data.frame(id = rownames(counts), gene_meta,
                          stringsAsFactors = FALSE)
  attr(counts, "gene_meta") <- gene_meta
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
  ss <- utils::read.delim(
    path,
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  ss <- .canonicalize_samplesheet_cols(ss)

  missing <- setdiff(c("SampleID", "GroupID"), colnames(ss))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Samplesheet is missing required column{?s} {.val {missing}}.",
      "i" = "Columns found: {.val {colnames(ss)}}",
      "i" = "Expected {.val SampleID} then {.val GroupID}, followed by one column per contrast.",
      "x" = "Without {.val GroupID} no contrast can resolve its experimental/control groups."
    ))
  }
  ss
}


#' Samplesheet column aliases (squashed lower-case key -> canonical name)
#'
#' @description Real-world samplesheets spell the ID columns many ways
#'   (`Sample_ID`, `sample id`, `Group_ID`, `condition`, ...). Keys here are
#'   the column name lower-cased with all non-alphanumerics removed.
#' @keywords internal
.samplesheet_aliases <- c(
  sampleid    = "SampleID",
  sample      = "SampleID",
  samplename  = "SampleID",
  samples     = "SampleID",
  groupid     = "GroupID",
  group       = "GroupID",
  groupname   = "GroupID",
  condition   = "GroupID",
  displayname = "DisplayName",
  exclude     = "Exclude"
)


#' Normalize samplesheet column names to the canonical set
#'
#' @description Rename known header variants (see [.samplesheet_aliases]) to
#'   the canonical `SampleID` / `GroupID` / `DisplayName` / `Exclude` names the
#'   pipeline uses. Without this a header like `Sample_ID,Group_ID` parses
#'   fine but yields `samplesheet$GroupID == NULL`, so every contrast silently
#'   finds no groups and the run dies much later with the misleading
#'   "No contrasts to analyse" error.
#'
#'   Contrast columns (anything containing `_vs_`) are never renamed, and a
#'   variant is skipped if its canonical name is already present, so no
#'   duplicate columns can be created.
#'
#' @param df Samplesheet data frame.
#' @return `df` with canonical column names.
#' @keywords internal
.canonicalize_samplesheet_cols <- function(df) {
  nm <- colnames(df)
  squashed <- tolower(gsub("[^A-Za-z0-9]", "", nm))
  renamed <- character(0)

  for (i in seq_along(nm)) {
    if (grepl("_vs_", nm[i], fixed = TRUE)) next          # a contrast column
    key <- squashed[i]
    if (!key %in% names(.samplesheet_aliases)) next
    canon <- .samplesheet_aliases[[key]]
    if (identical(nm[i], canon)) next                     # already canonical
    if (canon %in% nm) next                               # don't duplicate
    renamed <- c(renamed, sprintf("%s -> %s", nm[i], canon))
    nm[i] <- canon
  }

  if (length(renamed) > 0) {
    cli::cli_alert_info("Samplesheet column{?s} normalized: {.val {renamed}}")
  }
  colnames(df) <- nm
  df
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
