# ============================================================
# REFRESH THE VALIDATION-TRIAL DESIGN FROM NEW T3 DATA
#
#   Rscript code/validate_refresh.R                 # the whole chain
#   Rscript code/validate_refresh.R --from assemble # resume part-way
#   Rscript code/validate_refresh.R --dry-run       # print the plan, run nothing
#   Rscript code/validate_refresh.R --only validate # just the validation half
#
# Run this when new trials land on T3. It drives the existing scripts in order
# and stops at the first failure, so the thing to read is the step that failed
# rather than this file.
#
# WHAT IT DOES NOT TOUCH: the simulation's parameters. SIM_VAR_SHARES,
# SIM_PR_AS_COR and SIM_RESID_COR in code/sim_config.R are FROZEN at the values
# the six-trial fit gave on 2026-09-21, deliberately, so that simulation results
# stay comparable across data vintages. New variance components will appear in
# output/BGLR_variance_components.csv and will be used by the validation power
# calculation -- which is what should happen -- but nothing here edits
# sim_config.R. If you ever do want to refresh them, that is a separate,
# deliberate act: sim_observed_parameters() prints the current estimates and
# changing the constants invalidates every cached simulation result.
#
# NEEDS T3 CREDENTIALS. Steps 1-2 talk to T3/Oat over BrAPI and read .Renviron.
# Steps 3 onward are offline.
# ============================================================

library(tidyverse)

here::i_am("code/validate_refresh.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
dry_run <- "--dry-run" %in% args
from    <- arg_value("--from", NULL)
only    <- arg_value("--only", NULL)

# The chain, in dependency order. `needs_t3` marks the steps that will fail
# without credentials; `slow` is a warning, not a limit.
STEPS <- tibble::tribble(
  ~stage,      ~script,                             ~needs_t3, ~slow, ~why,
  "discover",  "find_trials_with_B4I_accessions.R",  TRUE,   TRUE,
    "finds the new trials and downloads their observations; without this the new plots do not exist locally",
  "curate",    "curate_oat_accessions.R",            TRUE,   FALSE,
    "new trials may bring accessions that are duplicates of existing lines",
  "curate",    "curate_pea_accessions.R",            TRUE,   FALSE,
    "the same on the pea side",
  "curate",    "curation_report.R",                  FALSE,  FALSE,
    "rewrites CURATION.md so the collapses are documented at this vintage",
  "assemble",  "assemble_B4I_phenotypes.R",          FALSE,  FALSE,
    "rebuilds the plot table the models read",
  "fit",       "BGLR_multi_trait_model.R",           FALSE,  TRUE,
    "the variance components and per-accession effects everything below depends on",
  "validate",  "validate_pool_selection.R",          FALSE,  FALSE,
    "rebuilds the As+/As- pools and diffs them against the previous vintage",
  "validate",  "validate_power.R",                   FALSE,  TRUE,
    "recomputes power from the NEW variance components",
  "validate",  "validate_design.R",                  FALSE,  FALSE,
    "writes the field book and checks its balance",
  "validate",  "validate_crossval.R",                FALSE,  TRUE,
    "re-estimates the attenuation lambda by leave-one-trial-out"
)

STAGE_ORDER <- c("discover", "curate", "assemble", "fit", "validate")

steps <- STEPS
if (!is.null(only)) {
  if (!only %in% STAGE_ORDER) {
    stop("--only must be one of: ", paste(STAGE_ORDER, collapse = ", "),
         call. = FALSE)
  }
  steps <- dplyr::filter(steps, stage == only)
}
if (!is.null(from)) {
  if (!from %in% STAGE_ORDER) {
    stop("--from must be one of: ", paste(STAGE_ORDER, collapse = ", "),
         call. = FALSE)
  }
  keep <- STAGE_ORDER[match(from, STAGE_ORDER):length(STAGE_ORDER)]
  steps <- dplyr::filter(steps, stage %in% keep)
}
if (nrow(steps) == 0) stop("nothing to run", call. = FALSE)

cat("\n", strrep("=", 72), "\nValidation-trial refresh: ", nrow(steps),
    " step(s)\n", strrep("=", 72), "\n", sep = "")
for (i in seq_len(nrow(steps))) {
  cat(sprintf("%2d. [%-9s] %-38s %s%s\n", i, steps$stage[i], steps$script[i],
              if (steps$needs_t3[i]) "T3 " else "   ",
              if (steps$slow[i]) "slow" else ""))
}
cat("\nSimulation parameters in code/sim_config.R are NOT touched by any of",
    "this.\n")

if (dry_run) {
  cat("\n--dry-run: nothing was run. Each step's purpose:\n\n")
  for (i in seq_len(nrow(steps))) {
    cat("  ", steps$script[i], "\n     ", steps$why[i], "\n", sep = "")
  }
  quit(save = "no")
}

# A missing .Renviron is worth catching before a twenty-minute download, not
# after it.
if (any(steps$needs_t3)) {
  env_file <- here::here(".Renviron")
  if (file.exists(env_file)) readRenviron(env_file)
  if (Sys.getenv("T3_USERNAME") == "" || Sys.getenv("T3_PASSWORD") == "") {
    stop("steps ", paste(steps$script[steps$needs_t3], collapse = ", "),
         " need T3 credentials, and T3_USERNAME / T3_PASSWORD are not set. ",
         "Put them in .Renviron at the project root, or run with ",
         "--from assemble if the download is already done.", call. = FALSE)
  }
}

log_dir <- here::here("output", "refresh_logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

run_step <- function(script) {
  path <- here::here("code", script)
  if (!file.exists(path)) stop("missing script: ", path, call. = FALSE)
  log <- file.path(log_dir, paste0(stamp, "_", tools::file_path_sans_ext(script), ".log"))
  cat("\n", strrep("-", 72), "\n", script, "\n", strrep("-", 72), "\n", sep = "")
  t0 <- Sys.time()
  status <- system2("Rscript", shQuote(path), stdout = log, stderr = log)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  tail_n <- function(n) {
    x <- readLines(log, warn = FALSE)
    utils::tail(x, n)
  }
  if (status != 0) {
    cat(paste0("  ", tail_n(25), collapse = "\n"), "\n")
    stop(script, " failed after ", round(mins, 1), " min (exit ", status,
         "). Full log: ", log, call. = FALSE)
  }
  cat(paste0("  ", tail_n(6), collapse = "\n"), "\n")
  cat(sprintf("  ok (%.1f min) -> %s\n", mins, basename(log)))
  mins
}

timings <- purrr::map_dbl(steps$script, run_step)

cat("\n", strrep("=", 72), "\nDone in ", round(sum(timings), 1), " min\n",
    strrep("=", 72), "\n", sep = "")
print(tibble::tibble(step = steps$script, minutes = round(timings, 1)))

# ------------------------------------------------------------
# What changed, since that is the point of refreshing
# ------------------------------------------------------------

vint_dir <- here::here("output", "validation")
if (dir.exists(vint_dir)) {
  vints <- sort(list.dirs(vint_dir, recursive = FALSE, full.names = FALSE))
  cat("\nvalidation vintages on disk:", paste(vints, collapse = ", "), "\n")
  if (length(vints) >= 2) {
    cat("\nThe pool churn between the last two vintages is in\n  ",
        file.path(vint_dir, utils::tail(vints, 1), "pool_diff.csv"), "\n",
        "A pool retaining well under 70% of its members across one data",
        " vintage is\na warning about the EFFECT ESTIMATES, not about the",
        " design.\n", sep = "")
  }
}

cat("\nRead next:\n")
cat("  VALIDATION_DESIGN.md                     the design and what the numbers mean\n")
cat("  output/validation/<today>/pool_summary.csv   who is in which pool\n")
cat("  output/validation/<today>/power_grid.csv     power at the new variance components\n")
cat("  output/validation/<today>/field_book.csv     the trial layout\n")
