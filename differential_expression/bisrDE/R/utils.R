# Small string + file helpers shared by the pipeline orchestrator.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.5 spec.

#' Build a munged-SampleID -> display-label lookup
#'
#' @description Construct a named character vector mapping each sample's
#'   `make.names()`-munged ID (the form used as count-matrix column names
#'   and plot label keys) to the label that should be shown on plots.
#'
#'   When `sample_info` carries a `display` column (populated from the
#'   samplesheet's optional `DisplayName` column), those values are used.
#'   Otherwise the lookup is the identity on the munged ID, so plots render
#'   exactly as they did before display names existed.
#'
#' @param sample_info Data frame with at least a `sample` column; optionally
#'   a `display` column.
#' @return Named character vector: names are `make.names(sample_info$sample)`,
#'   values are the display labels (falling back to the munged ID per cell
#'   when a `display` entry is blank/`NA`).
#' @keywords internal
.display_lookup <- function(sample_info) {
  munged <- make.names(as.character(sample_info$sample))
  disp <- munged
  if (!is.null(sample_info$display)) {
    d <- as.character(sample_info$display)
    has <- !is.na(d) & nzchar(trimws(d))
    disp[has] <- d[has]
  }
  stats::setNames(disp, munged)
}


#' Auto-derive friendly per-sample display labels from group + replicate
#'
#' @description When the samplesheet carries no `DisplayName` column the
#'   pipeline still wants human-friendly plot labels instead of raw SampleIDs
#'   (e.g. SRA accessions). This builds `"<GroupID> <n>"` labels, numbering
#'   replicates `1..n` within each group in samplesheet order (e.g.
#'   `"hpvPositive 1"`, `"hpvPositive 2"`, ...). Used as the default for
#'   `sample_info$display`; an explicit `DisplayName` value overrides it per
#'   sample.
#'
#' @param conditions Character vector of group labels, in samplesheet order.
#' @return Character vector of display labels aligned to `conditions`.
#' @keywords internal
.auto_display <- function(conditions) {
  conditions <- as.character(conditions)
  out <- character(length(conditions))
  for (g in unique(conditions)) {
    idx <- which(conditions == g)
    out[idx] <- paste(g, seq_along(idx))
  }
  out
}


#' Build a munged-SampleID -> condition (group) lookup
#'
#' @description Named vector mapping each sample's `make.names()`-munged ID
#'   to its `condition`/group label. Use this to align group labels to a
#'   matrix's column order (or a `prcomp` rowname order) rather than
#'   relying on `sample_info` row order matching the data column order.
#'
#' @param sample_info Data frame with `sample` and `condition` columns.
#' @return Named character vector: names are `make.names(sample_info$sample)`,
#'   values are the group labels.
#' @keywords internal
.group_lookup <- function(sample_info) {
  munged <- make.names(as.character(sample_info$sample))
  stats::setNames(as.character(sample_info$condition), munged)
}


#' Okabe-Ito colorblind-safe categorical palette
#'
#' @description Return `n` colors from the Okabe-Ito qualitative palette —
#'   the de-facto colorblind-safe categorical set (deuteranopia/protanopia/
#'   tritanopia distinguishable). Used for group/condition colors across all
#'   plots so the scheme is consistent and accessible. Recycles via a ramp
#'   only if `n` exceeds the 8 named colors.
#'
#' @param n Number of colors needed.
#' @return Character vector of `n` hex colors.
#' @keywords internal
.okabe_ito <- function(n) {
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#000000")
  if (n <= length(pal)) pal[seq_len(n)] else grDevices::colorRampPalette(pal)(n)
}


#' Diverging colorblind-safe color stops (low, mid, high)
#'
#' @description Blue -> near-white -> vermillion stops for a diverging
#'   continuous scale (correlation, z-score). Pair with
#'   `circlize::colorRamp2()` or `grDevices::colorRampPalette()`. Avoids the
#'   red/green problem of the previous ad-hoc ramps.
#'
#' @return Character vector of 3 hex colors.
#' @keywords internal
.diverging_stops <- function() c("#0072B2", "#F7F7F7", "#D55E00")


#' Pick the DE-table ID column that matches a count matrix's row names
#'
#' @description A count matrix is keyed by the INPUT gene IDs (ensembl /
#'   symbol / entrez, per `--id-type`), which the annotated DE table stores
#'   under `"<KEYTYPE>_ID"`. Downstream subsetting that hardcodes `ENSEMBL_ID`
#'   silently drops every gene for symbol/entrez inputs. This returns the
#'   candidate ID column of `df` with the most overlap with `target`, so the
#'   right key is used regardless of id_type.
#'
#' @param df Annotated DE data frame.
#' @param target Character vector to match against (e.g. `rownames(counts)`).
#' @return Name of the best-matching column, or `NULL` if none overlaps.
#' @keywords internal
.match_id_column <- function(df, target) {
  candidates <- c("ENSEMBL_ID", "SYMBOL_ID", "ENTREZID_ID", "SYMBOL", "ENTREZID")
  candidates <- candidates[candidates %in% colnames(df)]
  best <- NULL
  best_hits <- 0L
  for (cc in candidates) {
    hits <- length(intersect(as.character(df[[cc]]), target))
    if (hits > best_hits) {
      best_hits <- hits
      best <- cc
    }
  }
  best
}


#' Append per-sample normalized counts to a DE results table
#'
#' @description Join a normalized-count matrix onto an annotated DE data frame
#'   by gene, adding one `<prefix>_<SampleID>` column per sample so the DE
#'   spreadsheet carries expression values next to the log2FC. Matrix columns
#'   (which are `make.names()`-munged SampleIDs) are renamed back to the
#'   original SampleIDs via `orig_ids`.
#'
#' @param df Annotated DE data frame.
#' @param ids Character vector (length `nrow(df)`) of each row's gene ID in the
#'   matrix's row-name namespace (see [.match_id_column()]).
#' @param mat Normalized-count matrix (genes x samples).
#' @param orig_ids Named vector mapping munged SampleID -> original SampleID.
#' @param prefix Column-name prefix, e.g. `"TMM"` or `"DESeq2norm"`.
#' @return `df` with the per-sample count columns appended (`NA` for any gene
#'   absent from `mat`).
#' @keywords internal
.append_counts <- function(df, ids, mat, orig_ids, prefix) {
  sub <- mat[match(ids, rownames(mat)), , drop = FALSE]
  cn  <- unname(orig_ids[colnames(mat)])
  cn[is.na(cn)] <- colnames(mat)[is.na(cn)]
  colnames(sub) <- paste0(prefix, "_", cn)
  cbind(df, as.data.frame(sub, check.names = FALSE, stringsAsFactors = FALSE))
}


#' Build a `<base_dir>/<prefix><name><extension>` file path
#'
#' @description Convenience wrapper around `file.path()` + `paste0()` for
#'   the pipeline's deterministic "prefix + comparison name + extension"
#'   output naming.
#'
#' @param base_dir Base directory path.
#' @param prefix Prefix to add before the name.
#' @param name Main file name.
#' @param extension File extension. Default `".csv"`.
#' @return A complete file path string.
#' @export
create_file_path <- function(base_dir, prefix, name, extension = ".csv") {
  file.path(base_dir, paste0(prefix, name, extension))
}


#' Format an "exp vs. ctrl" comparison label
#'
#' @description Build the canonical `"<prefix><exp> vs. <ctrl>"` string used
#'   in plot titles and console alerts throughout the pipeline.
#'
#' @param exp Experimental group name.
#' @param ctrl Control group name.
#' @param prefix Optional prefix (default `""`).
#' @return The formatted comparison string.
#' @export
create_comparison_name <- function(exp, ctrl, prefix = "") {
  paste0(prefix, exp, " vs. ", ctrl)
}


#' Save a ggplot to disk via ggsave with the BISR defaults
#'
#' @description Thin wrapper around `ggplot2::ggsave` that locks in the
#'   `bg = "white"` background and the BISR default plot dimensions.
#'   Centralized so figure sizing is consistent across the pipeline.
#'
#' @param plot A `ggplot` object.
#' @param filename Output file path.
#' @param width Plot width in inches. Default `12` (was `15` pre-v1.4.0;
#'   reduced ~21% based on report-feedback that 15-inch volcano + dotplots
#'   crowded the embed-resources HTML).
#' @param height Plot height in inches. Default `13.5` (was `17` pre-v1.4.0;
#'   same -21% rationale as `width`).
#' @return Invisibly the path; called for its side effect.
#' @importFrom ggplot2 ggsave
#' @export
save_plot <- function(plot, filename, width = 12, height = 13.5) {
  ggplot2::ggsave(plot, filename = filename, width = width, height = height,
                  bg = "white")
}


#' Export a plotly figure to a self-contained HTML file
#'
#' @description Wraps `htmlwidgets::saveWidget(..., selfcontained = TRUE)`
#'   with input validation and `cli`-style error reporting.
#'
#' @param plotly_obj A `plotly` object.
#' @param file_path Path where the HTML file should be saved.
#' @return Invisibly `NULL`. Errors are caught and logged via
#'   `cli::cli_alert_danger` rather than re-raised, so a failed export does
#'   not abort the pipeline.
#'
#' @importFrom htmlwidgets saveWidget
#' @export
export_plotly_to_html <- function(plotly_obj, file_path) {
  tryCatch(
    {
      if (!inherits(plotly_obj, "plotly")) {
        stop("Invalid plotly object")
      }
      htmlwidgets::saveWidget(plotly_obj, file_path, selfcontained = TRUE)
    },
    error = function(e) {
      cli::cli_alert_danger("export_plotly_to_html: {e$message}")
    }
  )
  invisible(NULL)
}


#' Smallest group size for the low-count pre-filter
#'
#' @description The DESeq2-vignette pre-filter keeps genes with >= 10 reads
#'   in at least `smallest_group_size` samples. Deriving that number from the
#'   design (the smallest number of samples in any group) instead of a fixed
#'   3 means a 2-vs-2 run no longer demands three expressing samples and a
#'   6-vs-6 run no longer keeps genes expressed in only half a group.
#'
#' @param sample_info Data frame with a `condition` column.
#' @param floor Minimum returned value. Default `2`.
#' @return An integer >= `floor`.
#' @keywords internal
.smallest_group_size <- function(sample_info, floor = 2L) {
  n <- table(as.character(sample_info$condition))
  if (length(n) == 0) return(as.integer(floor))
  max(as.integer(floor), as.integer(min(n)))
}


#' Versions of R and the packages that determine the results
#'
#' @description One record used by the run JSON, the RDS bundle and the
#'   report's session section, so every artifact of a run states the same
#'   versions. Packages that are not installed are reported as `NA`.
#'
#' @return A named list of version strings.
#' @keywords internal
.software_versions <- function() {
  v <- function(pkg) tryCatch(as.character(utils::packageVersion(pkg)),
                              error = function(e) NA_character_)
  list(
    R               = paste(R.version$major, R.version$minor, sep = "."),
    platform        = R.version$platform,
    bisrDE          = v("bisrDE"),
    DESeq2          = v("DESeq2"),
    apeglm          = v("apeglm"),
    edgeR           = v("edgeR"),
    clusterProfiler = v("clusterProfiler"),
    ReactomePA      = v("ReactomePA"),
    msigdbr         = v("msigdbr"),
    ComplexHeatmap  = v("ComplexHeatmap"),
    quarto_r        = v("quarto")
  )
}


#' Where the analysis is running, for the Methods text
#'
#' @description Resolved at analysis time and stored in the run options, so a
#'   report re-rendered elsewhere still names the machine that produced the
#'   numbers. Precedence: the `BISR_PLATFORM_NAME` environment variable (set it
#'   in your job script, e.g. "VCU's High Performance Research Computing
#'   cluster"), else a Slurm description when `SLURM_JOB_ID` is set, else
#'   `Sys.info()` plus the R version.
#'
#' @return A single string.
#' @keywords internal
.run_platform <- function() {
  env <- Sys.getenv("BISR_PLATFORM_NAME", "")
  if (nzchar(env)) return(env)
  si <- Sys.info()
  rv <- paste(R.version$major, R.version$minor, sep = ".")
  job <- Sys.getenv("SLURM_JOB_ID", "")
  if (nzchar(job)) {
    return(sprintf("a Slurm cluster node (%s, job %s), R %s", si[["nodename"]], job, rv))
  }
  sprintf("%s %s (%s), R %s", si[["sysname"]], si[["release"]], si[["machine"]], rv)
}
