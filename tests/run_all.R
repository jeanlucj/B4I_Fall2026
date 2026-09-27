# run_all.R -- run every tests/test_*.R in its own process and aggregate.
#
# A fresh process per file is deliberate. The test files source project code and
# set seeds, so running them in one session would let one file's loaded
# functions, RNG state or options decide another file's result -- and a suite
# whose outcome depends on file order is worse than no suite.
#
# Usage:
#   Rscript tests/run_all.R          fast tier, seconds
#   Rscript tests/run_all.R --all    everything, including the MCMC fits
#
# The slow tier fits real chains and takes minutes. It is excluded by default so
# the fast tier stays cheap enough to run after every change, which is the only
# way a suite actually gets run.

here::i_am("tests/run_all.R")

args     <- commandArgs(trailingOnly = TRUE)
run_slow <- "--all" %in% args

slow_files <- c("test_fits.R")          # minutes: BGLR Multitrait and MegaLMM

files <- sort(list.files(here::here("tests"), pattern = "^test_.*[.]R$"))
if (!run_slow) files <- setdiff(files, slow_files)

if (length(files) == 0) stop("no test files found in tests/", call. = FALSE)

cat("Running", length(files), "test file(s)",
    if (run_slow) "(including the slow tier)" else "(fast tier only)", "\n\n")

results <- vapply(files, function(f) {
  cat("---", f, "---\n")
  t0  <- Sys.time()
  out <- system2("Rscript", shQuote(here::here("tests", f)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  secs   <- as.numeric(difftime(Sys.time(), t0, units = "s"))

  # Echo the summary line and anything that failed; the rest of a passing
  # file's chatter (package load messages, MCMC progress) is noise.
  keep <- grep("passed,|FAIL:|SKIP:|Error", out, value = TRUE)
  cat(paste0("  ", keep, collapse = "\n"), "\n")
  cat(sprintf("  (%.1fs)\n\n", secs))

  if (is.null(status)) 0L else as.integer(status)
}, integer(1))

failed <- names(results)[results != 0L]

cat(strrep("=", 60), "\n")
if (length(failed) == 0) {
  cat("all", length(files), "test file(s) passed\n")
  if (!run_slow) {
    cat("the slow tier was skipped:", paste(slow_files, collapse = ", "),
        "-- run with --all\n")
  }
} else {
  cat(length(failed), "of", length(files), "test file(s) FAILED:\n")
  cat(paste0("  ", failed, collapse = "\n"), "\n")
  quit(status = 1)
}
