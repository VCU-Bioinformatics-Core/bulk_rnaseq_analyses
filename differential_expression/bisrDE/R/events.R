# Structured progress events.
#
# When the environment variable BISR_EVENTS_FILE is set, the pipeline appends
# one JSON object per line (NDJSON) describing what it is doing. The Go TUI
# creates that file, exports the variable to the child, and tails it to drive a
# live bubbletea UI — the same "state -> re-render" model Claude Code gets from
# React/Ink, rather than trying to scrape pretty terminal output.
#
# Entirely opt-in: with the variable unset every function here is a no-op, so
# the bash launcher and plain `Rscript de.R` behave exactly as before.
#
# Event shapes (field `t` is the type):
#   {"t":"start","total":30,"comparisons":3,"runid":"..."}
#   {"t":"phase","name":"Sample Exploration QC plots"}
#   {"t":"tick","current":7,"total":30,"i":1,"n":3,
#    "comparison":"A_vs_B","step":"volcano"}
#   {"t":"done","ok":true}

#' Path of the progress-event stream, or `""` when disabled
#' @keywords internal
.events_path <- function() Sys.getenv("BISR_EVENTS_FILE", "")


#' Is structured progress-event emission enabled?
#'
#' @description TRUE when `BISR_EVENTS_FILE` names a writable path. Used both
#'   to gate [.emit_event()] and to suppress the pipeline's own `cli` progress
#'   bar, so exactly one component renders progress at a time.
#' @keywords internal
.events_enabled <- function() nzchar(.events_path())


#' Append one progress event to the NDJSON stream
#'
#' @description No-op unless [.events_enabled()]. Every failure is swallowed:
#'   a telemetry side-channel must never take down an analysis run. The file is
#'   opened in append mode per event (a few dozen per run) so each line is
#'   flushed immediately and visible to a tailing reader.
#'
#' @param t Event type: `"start"`, `"phase"`, `"tick"`, `"done"` or `"error"`.
#' @param ... Additional scalar fields to include in the JSON object.
#' @return Invisibly `NULL`.
#' @keywords internal
.emit_event <- function(t, ...) {
  path <- .events_path()
  if (!nzchar(path)) return(invisible(NULL))

  tryCatch({
    payload <- c(list(t = t), list(...))
    line <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null",
                             digits = NA)
    con <- file(path, open = "at")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    writeLines(as.character(line), con)
  }, error = function(e) NULL, warning = function(w) NULL)

  invisible(NULL)
}
