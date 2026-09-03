# cli session-log helpers (from Phase 2's CLI overhaul).
# Source-of-truth during Phase 5.2-5.7 is still the parent's de.R; this
# package copy is a snapshot per the Phase 5.5 spec.

#' Start a tee'd session log for the current pipeline run
#'
#' @description Open a `sink(split = TRUE)` on the output stream so every
#'   `cli::*`, `print()`, and `cat()` call from this point on is mirrored
#'   to BOTH the terminal AND a timestamped session log file.
#'
#' @param log_dir Directory where the session log will be written. Created
#'   if it does not exist.
#' @return A list with `path` (file path of the log) and `con` (open file
#'   connection). The caller MUST close the connection at end-of-run via
#'   [stop_session_log()].
#' @details `cli` routes its output to `stderr()` whenever an output or
#'   message sink is active (cli NEWS #153) — which is why this log used to
#'   come out empty (the whole pipeline narrates through `cli::*`). We call
#'   `cli::start_app(output = "stdout")` to pin cli to stdout so the
#'   `split = TRUE` output sink both shows it live on the terminal AND writes
#'   it to the log. `message()` / `warning()` / `stop()` are still not split.
#'
#' @export
start_session_log <- function(log_dir) {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  path <- file.path(log_dir, paste0(ts, "_session.log"))
  con <- file(path, open = "wt")
  # Pin cli to stdout BEFORE opening the sink (see @details) so the split
  # sink captures cli's narration instead of it escaping to stderr.
  # `.auto_close = FALSE` is essential: start_app() otherwise stops the app
  # when THIS function's frame exits, i.e. before the pipeline logs anything.
  # stop_session_log() pops it explicitly instead.
  cli::start_app(output = "stdout", .auto_close = FALSE)
  sink(con, type = "output", split = TRUE)
  list(path = path, con = con)
}


#' Close a session log started by [start_session_log()]
#'
#' @param log Object returned by [start_session_log()].
#' @return Invisibly `NULL`.
#' @details Defensive: only unwinds the output sink stack if any sinks are
#'   still active (`sink.number() > 0`), and only closes the connection
#'   if it is still open.
#'
#' @export
stop_session_log <- function(log) {
  if (sink.number() > 0) sink(NULL, type = "output")
  cli::stop_app()
  if (!is.null(log$con) && isOpen(log$con)) close(log$con)

  if (!is.null(log$path) && file.exists(log$path)) {
    # Optional colour-preserved HTML copy, rendered from the RAW log while the
    # ANSI escapes are still present (aha needs them). Produces a linkable page
    # a PI can open with warnings still highlighted — something the fixed-size
    # `freeze` screenshot cannot do for a 5000-line log. No-op without aha.
    .render_log_html(log$path)

    # The split sink captured cli's raw ANSI + progress-redraw bytes; strip them
    # so the on-disk log is clean, plain-text and greppable.
    lines <- gsub("\r", "", cli::ansi_strip(readLines(log$path, warn = FALSE)))
    writeLines(lines, log$path)
  }
  invisible(NULL)
}


#' Render a session log to colour-preserved HTML with `aha`
#'
#' @description Convert the raw (still ANSI-coloured) session log into a
#'   standalone HTML page next to it, so warnings/errors keep their colour and
#'   the whole run is linkable from the report. Entirely optional: a no-op when
#'   the `aha` binary is not on `PATH`, and any failure is swallowed — a
#'   cosmetic artifact must never break a finished run.
#'
#' @param log_path Path to the raw session log.
#' @return Invisibly the HTML path, or `NULL` if nothing was written.
#' @keywords internal
.render_log_html <- function(log_path) {
  if (!nzchar(Sys.which("aha"))) return(invisible(NULL))
  html <- sub("\\.log$", ".html", log_path)
  ok <- tryCatch({
    status <- system2(
      "aha",
      args   = c("--black", "--title", shQuote("bisrDE session log"),
                 "-f", shQuote(log_path)),
      stdout = html, stderr = FALSE
    )
    identical(as.integer(status), 0L) && file.exists(html)
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (isTRUE(ok)) {
    cli::cli_alert_success("Session log (HTML): {.path {html}}")
    invisible(html)
  } else {
    if (file.exists(html)) unlink(html)
    invisible(NULL)
  }
}


#' Write a structured JSON run summary (relink-style `.report.json`)
#'
#' @description End-of-run provenance record written alongside the session
#'   log: run metadata, input paths, gene/sample counts, per-contrast DE +
#'   enrichment outcomes, timings, and output paths. Modelled on the
#'   relink-fastq `.report.json`. Failures are non-fatal (a warning only).
#'
#' @keywords internal
.write_report_json <- function(path, runid, brs_ticket, annotation, id_type,
                               counts_path, samplesheet_path, outdir,
                               started, finished, n_genes_input,
                               n_genes_filtered, n_samples, n_samplesheet_rows,
                               comparisons, results, rds_path, session_log,
                               padj = 0.05, lfc = 0.58,
                               options = NULL, software = NULL, qc = NULL) {
  iso <- function(t) format(t, "%Y-%m-%dT%H:%M:%S")
  n_sig_sets <- function(g) {
    if (is.null(g)) return(0L)
    sum(!is.na(g@result$p.adjust) & g@result$p.adjust < padj)
  }

  comp_summaries <- lapply(seq_along(results), function(i) {
    r  <- results[[i]]
    nm <- comparisons[[i]]$name
    if (is.null(r) || is.null(r$deseq)) return(list(name = nm, ok = FALSE))
    d   <- r$deseq
    sig <- sum(!is.na(d$padj) & d$padj < padj & abs(d$log2FoldChange) >= lfc)
    o   <- if (is.null(r$de_options)) list() else r$de_options
    list(
      name        = nm,
      ok          = TRUE,
      n_genes     = nrow(d),
      n_sig       = sig,
      n_padj_na   = sum(is.na(d$padj)),
      design      = o$design,
      reference_level = o$reference_level,
      coefficient = o$coef,
      lfc_shrink  = o$lfc_shrink,
      independent_filtering = o$independent_filtering,
      enrichments = list(GO = !is.null(r$gsea),     KEGG     = !is.null(r$kegg),
                         Reactome = !is.null(r$reactome),
                         Hallmark = !is.null(r$hallmark)),
      n_significant_sets = list(GO = n_sig_sets(r$gsea), KEGG = n_sig_sets(r$kegg),
                                Reactome = n_sig_sets(r$reactome),
                                Hallmark = n_sig_sets(r$hallmark))
    )
  })

  summary <- list(
    runid        = runid,
    brs_ticket   = brs_ticket,
    annotation   = annotation,
    id_type      = id_type,
    started      = iso(started),
    finished     = iso(finished),
    duration_sec = round(as.numeric(difftime(finished, started, units = "secs")), 1),
    inputs       = list(counts = counts_path, samplesheet = samplesheet_path,
                        outdir = outdir),
    n_genes_input      = n_genes_input,
    n_genes_filtered   = n_genes_filtered,
    n_samples          = n_samples,
    n_samplesheet_rows = n_samplesheet_rows,
    thresholds   = list(padj = padj, log2fc = lfc,
                        fold_change = round(2^lfc, 3)),
    contrasts    = vapply(comparisons, function(x) x$name, character(1)),
    comparisons  = comp_summaries,
    options      = options,
    software     = software,
    qc           = if (is.null(qc)) NULL else qc[setdiff(names(qc), c("libsize_M", "detected"))],
    outputs      = list(rds = rds_path, session_log = session_log)
  )

  tryCatch({
    jsonlite::write_json(summary, path, auto_unbox = TRUE, pretty = TRUE,
                         null = "null", na = "null")
    cli::cli_alert_success("Run summary: {.path {path}}")
  }, error = function(e)
    cli::cli_alert_warning("Could not write run summary JSON: {conditionMessage(e)}"))

  invisible(path)
}
