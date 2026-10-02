# ============================================================
# REFRESH THE NUMBERS IN VALIDATION_DESIGN.md
#
#   Rscript code/validation_report.R                  # latest vintage
#   Rscript code/validation_report.R --vintage 2026-10-02
#   Rscript code/validation_report.R --check          # report drift, write nothing
#
# VALIDATION_DESIGN.md is mostly argument, and an argument should not change
# because a number moved. But it also carries tables that come straight out of
# output/validation/<vintage>/, and those go stale the moment a trial is added
# or the QC screen keeps a different set. On 2026-10-02 the document still
# reported lambda from four folds (oat 1.28) while the vintage on disk had eight
# (oat 0.83) -- a difference large enough to change which species counts as
# adequately powered.
#
# So the data-dependent blocks are DELIMITED and regenerated; everything else is
# left exactly as written:
#
#   <!-- BEGIN GENERATED: lambda -->  ...  <!-- END GENERATED: lambda -->
#   <!-- BEGIN GENERATED: power -->   ...  <!-- END GENERATED: power -->
#   <!-- BEGIN GENERATED: vintage --> ...  <!-- END GENERATED: vintage -->
#
# A block whose markers are missing is reported and skipped rather than guessed
# at, because silently appending a table to a document is worse than leaving it
# out.
#
# Run from validate_refresh.R with --refresh-doc, or on its own.
# ============================================================

library(tidyverse)

here::i_am("code/validation_report.R")

source(here::here("code", "validation_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
check_only <- "--check" %in% args

out_root <- here::here("output", "validation")
doc_file <- here::here("VALIDATION_DESIGN.md")

vintages <- sort(basename(list.dirs(out_root, recursive = FALSE)))
if (length(vintages) == 0) stop("no vintages under ", out_root, call. = FALSE)
vintage <- arg_value("--vintage", utils::tail(vintages, 1))
vin_dir <- file.path(out_root, vintage)
if (!dir.exists(vin_dir)) stop("no such vintage: ", vin_dir, call. = FALSE)

message("reading ", vin_dir)

read_if <- function(name) {
  f <- file.path(vin_dir, name)
  if (file.exists(f)) readr::read_csv(f, show_col_types = FALSE) else NULL
}

cv    <- read_if("crossval_summary.csv")
grid  <- read_if("power_grid.csv")
vint  <- read_if("vintage.csv")
qc    <- if (file.exists(here::here("output", "trial_qc.csv")))
           readr::read_csv(here::here("output", "trial_qc.csv"),
                           show_col_types = FALSE) else NULL

fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)

# ------------------------------------------------------------
# The blocks
# ------------------------------------------------------------

block_vintage <- function() {
  if (is.null(vint)) return(NULL)
  lines <- c(
    sprintf("*Numbers below are from vintage **%s**.*", vintage), "",
    "| species | trials | plots | accessions | eligible | reliability of As |",
    "|---|---|---|---|---|---|",
    sprintf("| %s | %d | %d | %d | %d | %s |",
            vint$species, vint$n_trials, vint$n_plots,
            vint$n_accessions, vint$n_eligible, fmt(vint$reliability_As)))
  if (!is.null(qc)) {
    kept <- qc$studyName[qc$keep]
    drop <- qc$studyName[!qc$keep]
    lines <- c(lines, "",
      sprintf("Trials kept by the QC screen (%d): %s", length(kept),
              paste(kept, collapse = ", ")))
    if (length(drop)) {
      lines <- c(lines,
        sprintf("Dropped (%d): %s", length(drop), paste(drop, collapse = ", ")))
    }
  }
  lines
}

block_lambda <- function() {
  if (is.null(cv)) return(NULL)
  x <- dplyr::filter(cv, set == "informative folds only")
  if (nrow(x) == 0) x <- cv
  c("| | λ | across-fold spread (interaction_frac) | folds | accession-level *r* |",
    "|---|---|---|---|---|",
    sprintf("| %s | **%s** | %s | %d | %s |",
            x$species, fmt(x$lambda_mean), fmt(x$interaction_frac),
            x$n_folds, fmt(x$r_accession_mean)))
}

block_power <- function() {
  if (is.null(grid)) return(NULL)
  # The pool size is a CONFIGURED choice per species -- oat 20, pea 30 by
  # default -- not a formula, so the row to report is the one at that n rather
  # than the nearest to some geometry. Getting this wrong quietly reported pea
  # at oat's pool size.
  n_pool  <- validation_setting("n_per_pool")
  n_loc   <- validation_setting("n_locations")
  g <- grid |>
    dplyr::filter(n == unname(n_pool[species])) |>
    dplyr::arrange(species, P)
  if (nrow(g) == 0) {
    message("  power grid has no rows at the configured pool sizes -- skipped")
    return(NULL)
  }
  # A grid written before lambda and the interaction were measured still sweeps
  # them, and collapsing that to one row per budget would pick an arbitrary
  # lambda and present it as the estimate. Say so instead of guessing.
  swept <- dplyr::n_distinct(g$lambda) > 1 ||
           dplyr::n_distinct(g$interaction_frac) > 1 ||
           ("sided" %in% names(g) && dplyr::n_distinct(g$sided) > 1)
  if (swept) {
    message("  power_grid.csv still sweeps lambda/interaction/sided -- it ",
            "predates the measured values. Re-run code/validate_power.R.")
    return(c(
      "> **Stale.** `power_grid.csv` in this vintage was written before λ and the",
      "> pool × location interaction were measured, so it still sweeps them and",
      "> cannot be collapsed to a single statement. Re-run",
      "> `Rscript code/validate_power.R` and then `code/validation_report.R`."))
  }
  c("One-sided α = 0.05, at the configured pool size, with λ and the",
    "pool × location interaction as measured (see §6).",
    "",
    "| | plots/loc | n per pool | total plots | λ | interaction | power |",
    "|---|---|---|---|---|---|---|",
    sprintf("| **%s** | %d | %d | %d | %s | %s | %s |",
            g$species, as.integer(g$P / n_loc), g$n, g$P,
            fmt(g$lambda), fmt(g$interaction_frac), fmt(g$power)))
}

blocks <- list(vintage = block_vintage(), lambda = block_lambda(),
               power = block_power())

# ------------------------------------------------------------
# Splice
# ------------------------------------------------------------

doc <- readLines(doc_file, warn = FALSE)

splice <- function(doc, key, body) {
  b <- sprintf("<!-- BEGIN GENERATED: %s -->", key)
  e <- sprintf("<!-- END GENERATED: %s -->", key)
  i <- which(doc == b); j <- which(doc == e)
  if (length(i) != 1 || length(j) != 1 || j <= i) {
    message("  no markers for '", key, "' in ", basename(doc_file), " -- skipped")
    return(list(doc = doc, changed = FALSE))
  }
  old <- doc[(i + 1):(j - 1)]
  changed <- !identical(trimws(old), trimws(body))
  list(doc = c(doc[seq_len(i)], body, doc[j:length(doc)]), changed = changed)
}

changed_any <- FALSE
for (k in names(blocks)) {
  if (is.null(blocks[[k]])) {
    message("  no data for '", k, "' in this vintage -- skipped")
    next
  }
  r <- splice(doc, k, blocks[[k]])
  doc <- r$doc
  if (r$changed) {
    changed_any <- TRUE
    message("  '", k, "' block: numbers differ from the document")
  } else {
    message("  '", k, "' block: unchanged")
  }
}

if (check_only) {
  message(if (changed_any)
            "\nDRIFT: VALIDATION_DESIGN.md does not match this vintage. Re-run without --check."
          else "\nVALIDATION_DESIGN.md is up to date with this vintage.")
  quit(save = "no", status = if (changed_any) 1L else 0L)
}

writeLines(doc, doc_file)
message("\nwrote ", basename(doc_file),
        if (changed_any) " (numbers updated)" else " (no change)")
