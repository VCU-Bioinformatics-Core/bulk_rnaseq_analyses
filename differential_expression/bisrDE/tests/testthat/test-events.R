# Structured progress events (bisrDE/R/events.R). The Go TUI tails this NDJSON
# stream to render a live UI, so the contract matters: valid JSON, one object
# per line, and a complete no-op when BISR_EVENTS_FILE is unset.

with_events_file <- function(code) {
  path <- tempfile(fileext = ".ndjson")
  old <- Sys.getenv("BISR_EVENTS_FILE", unset = NA)
  Sys.setenv(BISR_EVENTS_FILE = path)
  on.exit({
    if (is.na(old)) Sys.unsetenv("BISR_EVENTS_FILE") else Sys.setenv(BISR_EVENTS_FILE = old)
  }, add = TRUE)
  force(code)
  path
}

test_that("events are disabled unless BISR_EVENTS_FILE is set", {
  old <- Sys.getenv("BISR_EVENTS_FILE", unset = NA)
  Sys.unsetenv("BISR_EVENTS_FILE")
  on.exit(if (!is.na(old)) Sys.setenv(BISR_EVENTS_FILE = old), add = TRUE)

  expect_false(bisrDE:::.events_enabled())
  # Must be a silent no-op, not an error.
  expect_silent(bisrDE:::.emit_event("tick", current = 1, total = 10))
})

test_that("emitted events are one valid JSON object per line", {
  skip_if_not_installed("jsonlite")
  path <- with_events_file({
    bisrDE:::.emit_event("start", total = 30L, comparisons = 3L, runid = "r1")
    bisrDE:::.emit_event("tick", current = 1L, total = 30L, step = "DESeq2",
                         comparison = "A_vs_B", i = 1L, n = 3L)
    bisrDE:::.emit_event("done", ok = TRUE)
  })

  lines <- readLines(path, warn = FALSE)
  expect_length(lines, 3)

  parsed <- lapply(lines, jsonlite::fromJSON)
  expect_equal(vapply(parsed, function(x) x$t, character(1)),
               c("start", "tick", "done"))

  # Scalars must be unboxed (auto_unbox), not wrapped in 1-element arrays —
  # the Go structs decode scalars, so a boxed value would fail to parse.
  expect_equal(parsed[[1]]$total, 30L)
  expect_equal(parsed[[2]]$step, "DESeq2")
  expect_equal(parsed[[2]]$comparison, "A_vs_B")
  expect_true(isTRUE(parsed[[3]]$ok))
})

test_that("events append rather than truncate", {
  path <- with_events_file({
    for (i in 1:5) bisrDE:::.emit_event("tick", current = i, total = 5L)
  })
  expect_length(readLines(path, warn = FALSE), 5)
})

test_that("emitting never errors even when the path is unwritable", {
  old <- Sys.getenv("BISR_EVENTS_FILE", unset = NA)
  Sys.setenv(BISR_EVENTS_FILE = file.path(tempdir(), "no_such_dir", "e.ndjson"))
  on.exit({
    if (is.na(old)) Sys.unsetenv("BISR_EVENTS_FILE") else Sys.setenv(BISR_EVENTS_FILE = old)
  }, add = TRUE)

  # A telemetry side-channel must never take down an analysis run.
  expect_silent(bisrDE:::.emit_event("tick", current = 1, total = 10))
})
