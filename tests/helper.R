# helper.R -- assertions shared by every tests/test_*.R.
#
# Sourced, not run. Each test file keeps its own `ok` / `fail` counters in its
# own process, prints a one-line summary and exits non-zero on failure, which is
# what tests/run_all.R aggregates.
#
# The assertions are deliberately few. What makes these tests worth having is
# not the framework but the ORACLES: most of them compare the code against a
# quantity known independently -- a planted duplicate, an algebraic identity, a
# variance that was requested -- rather than against a value the code produced
# earlier. A test that only pins today's output freezes bugs as well as
# behaviour.

ok <- 0L
fail <- 0L

#' TRUE or the test fails.
check <- function(cond, msg) {
  if (isTRUE(cond)) {
    ok <<- ok + 1L
  } else {
    fail <<- fail + 1L
    cat("  FAIL:", msg, "\n")
  }
  invisible(isTRUE(cond))
}

#' Equal to within a tolerance, with the discrepancy printed when it is not.
#'
#' Default tolerance is loose enough for a Monte Carlo quantity and tight enough
#' to catch a wrong formula; pass `tol` explicitly for an algebraic identity,
#' where 1e-10 is the right order.
check_near <- function(actual, expected, tol = 1e-6, msg = "") {
  d <- max(abs(actual - expected))
  if (is.finite(d) && d <= tol) {
    ok <<- ok + 1L
  } else {
    fail <<- fail + 1L
    cat("  FAIL:", msg, "-- max|diff| =", format(d, digits = 3),
        "exceeds", tol, "\n")
  }
  invisible(NULL)
}

#' The expression MUST raise an error.
#'
#' Used for the guards that exist to stop a silent mistake: a sparsity that
#' overruns its own request, a pool larger than its candidate set. A guard that
#' has stopped firing is worse than no guard, because the caller trusts it.
check_error <- function(expr, msg) {
  got <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  check(got, paste0(msg, " (expected an error, got none)"))
}

#' Source project files, quietly, in dependency order.
load_code <- function(...) {
  files <- c(...)
  for (f in files) {
    suppressMessages(invisible(utils::capture.output(
      source(here::here("code", f)))))
  }
}

#' Print the summary and set the exit status. Call once, last.
finish <- function(label) {
  cat(sprintf("\n%s: %d passed, %d failed\n", label, ok, fail))
  if (fail > 0) quit(status = 1)
}
