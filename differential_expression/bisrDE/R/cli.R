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
  # The split sink captured cli's raw ANSI + progress-redraw bytes; strip them
  # so the on-disk log is clean, plain-text and greppable.
  if (!is.null(log$path) && file.exists(log$path)) {
    lines <- gsub("\r", "", cli::ansi_strip(readLines(log$path, warn = FALSE)))
    writeLines(lines, log$path)
  }
  invisible(NULL)
}
