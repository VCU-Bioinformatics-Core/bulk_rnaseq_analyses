# Small string + file helpers shared by the pipeline orchestrator.
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.5 spec.

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
