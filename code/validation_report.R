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

# PREFER A PLAIN DATE. Archived vintages carry a suffix -- e.g.
# `2026-10-03_5trial-whitelist` -- and a suffix sorts AFTER the bare date, so
# taking the lexicographic last silently reported from the archive: every new
# table came back "no data for ... in this vintage -- skipped" while the old
# ones read "unchanged", which looks like a clean document rather than the
# wrong one. Fall back to the full list only if no bare date exists.
dated <- grep("^\\d{4}-\\d{2}-\\d{2}$", vintages, value = TRUE)
default_vintage <- utils::tail(if (length(dated)) dated else vintages, 1)
vintage <- arg_value("--vintage", default_vintage)
vin_dir <- file.path(out_root, vintage)
if (!dir.exists(vin_dir)) stop("no such vintage: ", vin_dir, call. = FALSE)

message("reading ", vin_dir)

read_if <- function(name) {
  f <- file.path(vin_dir, name)
  if (file.exists(f)) readr::read_csv(f, show_col_types = FALSE) else NULL
}

cv    <- read_if("crossval_summary.csv")
psum  <- read_if("pool_summary.csv")
grid  <- read_if("power_grid.csv")
vint  <- read_if("vintage.csv")
scal  <- read_if("crossval_lambda_scales.csv")
byyr  <- read_if("crossval_lambda_year.csv")
estp  <- read_if("estimand_power.csv")
t1    <- read_if("power_type1.csv")
qc    <- if (file.exists(here::here("output", "trial_qc.csv")))
           readr::read_csv(here::here("output", "trial_qc.csv"),
                           show_col_types = FALSE) else NULL

fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)

#' Is a table still sweeping something that is now measured?
#'
#' Counted WITHIN species, because lambda and interaction_frac are per-species
#' measurements and two species legitimately give two values -- a naive
#' across-the-board count calls a correct table "swept".
#'
#' Shared by both power blocks on purpose. The detector used to live inside
#' block_power() and checked three named columns; a table that gained a fourth
#' swept column would have been published stale without complaint. One helper
#' means a column added to one block cannot be forgotten in the other.
.is_swept <- function(g, cols, by = "species") {
  cols <- intersect(cols, names(g))
  if (length(cols) == 0) return(FALSE)
  per <- g |>
    dplyr::group_by(dplyr::across(dplyr::all_of(intersect(by, names(g))))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(cols), dplyr::n_distinct),
                     .groups = "drop")
  any(vapply(per[cols], \(x) any(x > 1), logical(1)))
}

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

# §2's contrast table was hand-maintained and went stale twice -- it still read
# 14.3 / 16.5 when the vintage on disk produced 15.24 / 18.01 -- so it is
# generated now. pool_summary.csv is the source.
block_pools <- function() {
  if (is.null(psum)) return(NULL)
  idx <- if ("pool_index" %in% names(psum)) unique(psum$pool_index)[1] else NA
  c("| species | n per pool | candidates | θ | ΔAs (g/m²) | ΔPr | mean Pr, As+ | mean Pr, As− |",
    "|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %d | %d | %s | **%s** | %s | %s | %s |",
            psum$species, psum$n_per_pool, psum$n_candidates,
            fmt(psum$theta, 3), fmt(psum$dAs, 2), fmt(psum$dPr, 3),
            fmt(psum$mean_Pr_plus, 2), fmt(psum$mean_Pr_minus, 2)),
    "",
    sprintf("*Selection rule: `%s`.*", if (is.na(idx)) "unstamped" else idx))
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
  # n_locations is swept deliberately, so the headline table pins it at the
  # configured value and the sweep is reported by the `estimands` block.
  if ("n_locations" %in% names(g)) {
    g <- dplyr::filter(g, n_locations == n_loc)
  }
  if (nrow(g) == 0) {
    message("  power grid has no rows at the configured pool sizes -- skipped")
    return(NULL)
  }
  # A grid written before lambda and the interaction were measured still sweeps
  # them, and collapsing that to one row per budget would pick an arbitrary
  # lambda and present it as the estimate. Say so instead of guessing.
  # Count distinct values WITHIN a species: lambda and interaction_frac are
  # per-species measurements, so two species legitimately give two values and
  # the naive across-the-board count calls a correct grid "swept".
  swept <- .is_swept(g, c("lambda", "interaction_frac", "sided",
                          "interaction_mode", "df_mode"))
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

#' The pre-registered analysis, printed from PREREG_MODEL.
#'
#' Generated rather than hand-written so section 4 cannot drift from the object
#' every power number is computed through -- which is exactly what happened
#' when section 4 specified a mixed model while both power routes used a
#' two-stage t-test.
block_analysis <- function() {
  m <- PREREG_MODEL
  out <- c(
    sprintf("*Analysis `%s`, fixed %s. Generated from `PREREG_MODEL` in",
            m$version, m$dated),
    "`code/validation_functions.R`, which is what every power number is",
    "computed through.*", "",
    "One random-effects structure, three instantiations:", "", "```",
    m$base, "```", "",
    "| estimand | role | response | `<FIXED>` | `<AsxE>` |",
    "|---|---|---|---|---|")
  for (nm in names(m$estimands)) {
    e <- m$estimands[[nm]]
    out <- c(out, sprintf("| `%s` | %s | %s | `%s` | `%s` |",
                          nm, e$role, e$response, e$fixed, e$AsxE))
  }
  c(out, "", paste0("**Multiplicity.** ", m$multiplicity), "",
    "**Secondary:**",
    paste0("- ", m$secondary), "",
    "**Considered and not done:**",
    paste0("- ", m$not_done))
}

#' Power for the three estimands -- the headline table.
block_estimands <- function() {
  if (is.null(estp)) return(NULL)
  if (.is_swept(estp, c("prereg", "plots_per_loc"), by = character(0))) {
    return(c("> **Stale.** `estimand_power.csv` mixes analysis versions or plot",
             "> budgets. Re-run `Rscript code/validate_power.R`."))
  }
  g <- estp |>
    dplyr::filter(lambda_source == "global") |>
    dplyr::arrange(estimand, target, n_loc)
  if (nrow(g) == 0) return(NULL)

  lo <- estp |>
    dplyr::filter(lambda_source == "lower95") |>
    dplyr::select(estimand, target, n_loc, power_lo = power)

  g <- dplyr::left_join(g, lo, by = c("estimand", "target", "n_loc"))

  c(sprintf("At **%d plots per location**, one-sided α = 0.05, λ as measured.",
            g$plots_per_loc[1]),
    "Locations are swept because the pool × location term is the largest of the",
    "three variance components and carries only `n_loc − 1` degrees of freedom.",
    "",
    "| estimand | role | target | locations | effect | SE | df | power | power at λ's lower 95% |",
    "|---|---|---|---|---|---|---|---|---|",
    sprintf("| `%s` | %s | %s | %d | %s | %s | %s | **%s** | %s |",
            g$estimand, g$role, g$target, g$n_loc, fmt(g$effect),
            fmt(g$SE, 3), fmt(g$df, 1), fmt(g$power), fmt(g$power_lo)),
    "",
    sprintf("The pool × location term is %s%% of the variance.",
            fmt(mean(g$pct_var_interaction), 0)))
}

#' Lambda on the two scales, and by year.
block_lambda_scales <- function() {
  if (is.null(scal)) return(NULL)
  out <- c(
    "λ is estimated twice: in absolute g/m², which is what the trial is sized",
    "in, and with the held-out trial's response divided by its own SD.",
    "`interaction_frac` is `sd(λ)/|mean(λ)|`, so it is scale-free and the two",
    "rows per species are directly comparable — which is what settles whether",
    "the across-fold spread is interaction or just the trials differing in",
    "spread. `interaction_frac_corrected` additionally removes fold-level",
    "estimation noise, using the accession-clustered standard errors.",
    "",
    "| species | scale | folds | λ | λ in g/m² | SE of λ | interaction_frac | corrected |",
    "|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %d | %s | %s | %s | %s | %s |",
            scal$species, scal$scale, scal$n_folds, fmt(scal$lambda_mean, 3),
            fmt(scal$lambda_gm2_mean, 2), fmt(scal$lambda_se_of_mean, 3),
            fmt(scal$interaction_frac, 3),
            fmt(scal$interaction_frac_corrected, 3)))

  if (!is.null(byyr)) {
    b <- dplyr::filter(byyr, scale == "raw")
    out <- c(out, "",
      "By year, as a **diagnostic** — the global estimate over all folds is what",
      "the power table uses, and nothing in the chain re-sizes on a subset of folds.",
      "",
      "| species | year | folds | λ | SD |",
      "|---|---|---|---|---|",
      sprintf("| %s | %d | %d | %s | %s |", b$species, b$year, b$n_folds,
              fmt(b$lambda_mean, 3), fmt(b$lambda_sd, 3)))
  }
  out
}

#' The eligibility rule, stated once.
block_eligibility <- function() {
  if (is.null(vint) || !"eligibility_rule" %in% names(vint)) return(NULL)
  c("| species | rule | reliability threshold | eligible | of |",
    "|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %d | %d |",
            vint$species, vint$eligibility_rule, fmt(vint$rel_min, 3),
            vint$n_eligible, vint$n_accessions),
    "",
    paste("Candidates must clear that reliability **and** have met at least two",
          "distinct partners in the data the model was fitted to. The partner",
          "floor is not redundant: reliability measures posterior precision,",
          "not identifiability, and an accession grown with a single partner",
          "has its producer and associate effects perfectly aliased while still",
          "scoring well by borrowing from genotyped relatives."))
}

blocks <- list(vintage = block_vintage(), eligibility = block_eligibility(),
               pools = block_pools(), analysis = block_analysis(),
               estimands = block_estimands(), lambda = block_lambda(),
               lambda_scales = block_lambda_scales(), power = block_power())

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
