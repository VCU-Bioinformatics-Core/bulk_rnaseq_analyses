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
#' @details R's `message()` / `warning()` / `stop()` channel is NOT split
#'   (R's `sink` does not support `split = TRUE` for message-type sinks);
#'   those still print to terminal but are not captured by this log. For
#'   the analytical pipeline this is fine: `cli::*` writes to stdout.
#'
#' @export
start_session_log <- function(log_dir) {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  path <- file.path(log_dir, paste0(ts, "_session.log"))
  con <- file(path, open = "wt")
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
  if (!is.null(log$con) && isOpen(log$con)) close(log$con)
  invisible(NULL)
}
