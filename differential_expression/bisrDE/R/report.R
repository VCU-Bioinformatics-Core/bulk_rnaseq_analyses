# Quarto report renderer.
# Phase 5.6 architectural change: parent's `report_generator.R` glue-built a
# fresh Rmd file for every run. Package version uses a real qmd template
# under `inst/qmd/` with native `params:` and child sections under
# `inst/qmd/_sections/`. Eliminates the glue-concatenation bug class.

#' Render the bulk RNA-seq DE report
#'
#' @description Render the package's Quarto template (`inst/qmd/report.qmd`)
#'   to a self-contained HTML file. Inputs come from the RDS bundle saved
#'   by [run_pipeline()].
#'
#' @param rds_path Path to the analysis RDS produced by [run_pipeline()].
#'   Stores `list(results, comparisons, out_dirs, pca_plot, pca_plotly,
#'   pca_3d, annotation)`.
#' @param output_dir Directory where the rendered HTML will be written.
#'   Created if it does not exist.
#' @param brs_ticket Optional BRS ticket identifier. Default `""`. When
#'   non-empty, rendered as a `subtitle:` line under the report title.
#' @param analyst Analyst name displayed in the Overview section.
#'   Default `"Mikail Bala"`.
#' @param report_prefix Prefix for the output HTML filename. Default
#'   `"rnaseq_analysis"`. Combined with a timestamp at render time.
#' @return Absolute path to the rendered HTML file.
#' @details Renders by:
#'   \enumerate{
#'     \item Locating the qmd template via `system.file("qmd/report.qmd",
#'           package = "bisrDE")`.
#'     \item Copying the template + the `_sections/` child folder into a
#'           tempdir.
#'     \item Substituting `\{\{SUBTITLE_LINE\}\}` in the YAML preamble with
#'           either `subtitle: "<brs_ticket>"` (if non-empty) or an empty
#'           string (no stray YAML field).
#'     \item Calling `quarto::quarto_render` with `execute_params` =
#'           `list(rds_path, analyst, brs_ticket, genome_assembly,
#'           sections_dir)`.
#'     \item Moving the rendered HTML to `output_dir` with a timestamped
#'           filename.
#'   }
#'   `genome_assembly` is derived from the RDS's `annotation` slot
#'   (`"human"` -> `"GRCh38 human primary assembly"`,
#'   `"mouse"` -> `"GRCm39 mouse primary assembly"`).
#'
#' @importFrom quarto quarto_render
#' @export
generate_report <- function(rds_path,
                            output_dir,
                            brs_ticket    = "",
                            analyst       = "Mikail Bala",
                            report_prefix = "rnaseq_analysis") {
  if (!file.exists(rds_path)) {
    stop("generate_report: RDS file does not exist: ", rds_path)
  }
  if (is.null(brs_ticket) || is.na(brs_ticket)) brs_ticket <- ""
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  rds_path   <- normalizePath(rds_path,   mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  # Read annotation from the RDS so we can render the genome-assembly line
  # in the manuscript Methods text correctly.
  rds_data   <- readRDS(rds_path)
  annotation <- rds_data[[7]]
  genome_assembly <- switch(
    tolower(as.character(annotation)),
    human = "GRCh38 human primary assembly",
    mouse = "GRCm39 mouse primary assembly",
    paste0(annotation, " primary assembly")
  )

  template_path <- system.file("qmd/report.qmd", package = "bisrDE")
  if (!nzchar(template_path) || !file.exists(template_path)) {
    stop("generate_report: report.qmd template not found in bisrDE package; ",
         "did you run devtools::load_all() or install the package?")
  }
  sections_src <- system.file("qmd/_sections", package = "bisrDE")
  if (!nzchar(sections_src) || !dir.exists(sections_src)) {
    stop("generate_report: qmd/_sections/ child folder not found in bisrDE package")
  }

  # Stage template + child sections in a fresh tempdir; rendering produces
  # auxiliary files (cache, support files) that we don't want polluting
  # `output_dir`. With `embed-resources: true` the final HTML is
  # self-contained, so we just move that one file at the end.
  td <- tempfile("bisrDE_qmd_render_")
  dir.create(td, recursive = TRUE)

  # Substitute the conditional subtitle line in the YAML preamble.
  template_lines <- readLines(template_path)
  subtitle_line <- if (nzchar(brs_ticket)) {
    paste0('subtitle: "', brs_ticket, '"')
  } else {
    ""
  }
  template_lines <- gsub("\\{\\{SUBTITLE_LINE\\}\\}", subtitle_line,
                         template_lines, fixed = FALSE)

  qmd_dest <- file.path(td, "report.qmd")
  writeLines(template_lines, qmd_dest)

  sections_dest <- file.path(td, "_sections")
  dir.create(sections_dest, recursive = TRUE, showWarnings = FALSE)
  for (f in list.files(sections_src, full.names = TRUE)) {
    file.copy(f, sections_dest, overwrite = TRUE)
  }

  # Trailing phase: this runs after run_pipeline() has emitted its "done"
  # event, so a front-end shows the bar at 100% while the report renders.
  .emit_event("phase", name = "rendering report")
  cli::cli_alert_info("Rendering Quarto report (template: {.path {template_path}})")

  quarto::quarto_render(
    input          = qmd_dest,
    output_format  = "html",
    execute_params = list(
      rds_path        = rds_path,
      analyst         = analyst,
      brs_ticket      = brs_ticket,
      genome_assembly = genome_assembly,
      sections_dir    = normalizePath(sections_dest, mustWork = TRUE)
    ),
    quiet = FALSE
  )

  rendered_html <- file.path(td, "report.html")
  if (!file.exists(rendered_html)) {
    stop("generate_report: Quarto did not produce report.html in ", td)
  }

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_file <- file.path(output_dir,
                           paste0(report_prefix, "_", ts, ".html"))
  file.copy(rendered_html, output_file, overwrite = TRUE)

  cli::cli_alert_success("Report rendered: {.path {output_file}}")
  invisible(output_file)
}
