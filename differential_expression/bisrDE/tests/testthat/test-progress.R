# Progress-bar granularity. cli never redraws on its own — a bar only repaints
# inside cli_progress_update() — so run_analysis() reports sub-steps via a
# `tick` callback and run_pipeline() sizes the bar as
# n_comparisons * .STEPS_PER_COMPARISON. If the two drift apart the bar either
# stalls short of 100% or overshoots, so guard the invariant.

test_that(".STEPS_PER_COMPARISON matches the number of .tick() calls", {
  src <- paste(deparse(body(bisrDE:::run_analysis)), collapse = "\n")
  hits <- gregexpr('\\.tick\\("', src)[[1]]
  n_ticks <- if (hits[1] == -1L) 0L else length(hits)

  expect_equal(n_ticks, bisrDE:::.STEPS_PER_COMPARISON)
})

test_that("run_analysis exposes an optional tick callback defaulting to NULL", {
  fm <- formals(bisrDE:::run_analysis)
  expect_true("tick" %in% names(fm))
  expect_null(eval(fm$tick))
})

test_that(".STEPS_PER_COMPARISON is a positive whole number", {
  n <- bisrDE:::.STEPS_PER_COMPARISON
  expect_true(is.numeric(n) && length(n) == 1L)
  expect_gt(n, 0)
  expect_equal(n, as.integer(n))
})
