# ============================================================
# FILTER THE SIMULATION RESULTS TO THE DESIGN THAT WAS RUN
#
#   Rscript code/sim_filter_results.R                        # design-driven
#   Rscript code/sim_filter_results.R --infer-design         # derive it from the reps
#   Rscript code/sim_filter_results.R --results a.csv --design b.csv
#   Rscript code/sim_filter_results.R --in-place             # overwrite, after a backup
#
# WHY THIS EXISTS. sim_run.R rebuilds simulation_results.csv by globbing
# output/simulation/*.rds, so every cache file still on disk joins the table --
# including files written by an earlier grid under different levels or a
# different generator. That is S10 in EVALUATION_SIMULATION.md's catalogue.
#
# THREE THINGS MAKE THIS LESS TRIVIAL THAN IT LOOKS.
#
# 1. The join is two keys, not one. The BGLR half does not depend on the
#    MegaLMM levers and is fitted once per data scenario, so its rows carry
#    K = eigen_variance = fixed_main_effect = NA. Matching everything on the
#    full design key would drop every additive, dge_ige, row_mean and
#    both_means row -- the entire comparator side of the experiment.
#
# 2. `rep` must NOT be part of the key. simulation_design.csv is rewritten by
#    every sim_run.R invocation, including `--combine`, and it records whatever
#    --reps THAT invocation was given. So a `--combine` run without --reps after
#    a --reps 5 grid leaves a rep-1 design describing 5-rep results, and keying
#    on rep would delete four fifths of the data. A design authorises SETTINGS;
#    replication is a separate axis.
#
# 3. The fraction is not portable. sim_design() seeds optFederov, but its
#    search still depends on the RNG stream, the AlgDesign version and the
#    platform, so regenerating the design locally does NOT reproduce the subset
#    a cluster chose. The design that describes a set of results is the one
#    written beside them, not one recomputed afterwards. When that file has
#    been lost or overwritten, --infer-design recovers the setting list from
#    the replicates themselves.
#
# WHAT THIS CANNOT DO. The results CSV carries no seed column, while the cache
# filename does. Two runs that assign different seeds to the same (scenario,
# rep) -- which is what happens when --reps is added after the fact, since rep
# is the inner index of the seed numbering -- both land in the CSV as rep 1 and
# are indistinguishable here. The script counts them and refuses to guess.
# ============================================================

library(tidyverse)

here::i_am("code/sim_filter_results.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
out_dir      <- here::here("output")
results_path <- arg_value("--results", file.path(out_dir, "simulation_results.csv"))
design_path  <- arg_value("--design",  file.path(out_dir, "simulation_design.csv"))
in_place     <- "--in-place" %in% args
infer        <- "--infer-design" %in% args

if (!file.exists(results_path)) stop("not found: ", results_path, call. = FALSE)
results <- readr::read_csv(results_path, show_col_types = FALSE)
message("results: ", nrow(results), " rows from ", basename(results_path))

SETTING <- c("scenario", "K", "eigen_variance", "fixed_main_effect")
r6 <- function(x) round(as.numeric(x), 6)
results <- dplyr::mutate(results, eigen_variance = r6(eigen_variance))
is_mm   <- stringr::str_starts(results$model, "megalmm")

# ------------------------------------------------------------
# STAGE 1: is this row from the CURRENT generator at all?
#
# The most reliable authority, and the only one always available, is
# sim_config.R -- versioned in the repo, and the definition of what the current
# grid is. Three independent signatures, because any one of them could in
# principle be coincidence and all three together cannot:
#
#   * the scenario NAME format. The bivariate rewrite added the GxE axis, so a
#     current name ends in _g<gxe_cor*100>; a pre-rewrite name does not.
#   * the LEVELS. Sparsity moved from 5/15/45% to 1.6/4.8/16/48% and
#     eigen_variance from 0.2/0.5/0.8 to 0.25/0.75.
#   * the MODEL names. The old single-trait scorer emitted `oat_mean` and
#     `oat_plus_pea` baselines, which the bivariate scorer does not produce.
#
# This is what a design file cannot do for you here: an earlier grid run with
# the same --reps passes any "present in every replicate" test, and an earlier
# grid's own design file is long overwritten.
# ------------------------------------------------------------

source(here::here("code", "dge_ige_functions.R"))
source(here::here("code", "sim_config.R"))

CURRENT_MODELS <- c("additive", "dge_ige", "row_mean", "both_means",
                    "megalmm", "megalmm_U")

current <- results |>
  dplyr::mutate(
    .name_ok    = stringr::str_detect(scenario, "_g\\d+$"),
    .sparsity_ok = r6(sparsity) %in% r6(SIM_LEVELS$sparsity),
    .nacc_ok     = n_acc %in% SIM_LEVELS$n_acc,
    .ev_ok       = is.na(eigen_variance) |
                     eigen_variance %in% r6(SIM_MEGALMM_LEVELS$eigen_variance),
    .K_ok        = is.na(K) | K %in% SIM_MEGALMM_LEVELS$K,
    .model_ok    = model %in% CURRENT_MODELS)

stale <- dplyr::filter(current, !(.name_ok & .sparsity_ok & .nacc_ok &
                                    .ev_ok & .K_ok & .model_ok))
current <- dplyr::filter(current, .name_ok, .sparsity_ok, .nacc_ok,
                         .ev_ok, .K_ok, .model_ok) |>
  dplyr::select(-dplyr::starts_with("."))

cat("\n=== stage 1: is the row from the current generator? ===\n")
cat(nrow(current), "current,", nrow(stale), "from an earlier grid\n")
if (nrow(stale) > 0) {
  cat("\nwhat marked the earlier-grid rows as stale:\n")
  print(tibble::tibble(
    signature = c("old scenario-name format (no _g suffix)",
                  "sparsity not a current level",
                  "n_acc not a current level",
                  "eigen_variance not a current level",
                  "K not a current level",
                  "model name the current scorer never emits"),
    rows = c(sum(!stale$.name_ok), sum(!stale$.sparsity_ok),
             sum(!stale$.nacc_ok), sum(!stale$.ev_ok),
             sum(!stale$.K_ok), sum(!stale$.model_ok))))
  cat("\ntheir levels:\n")
  print(dplyr::count(stale, n_acc, sparsity, name = "rows"), n = 20)
  cat("\ntheir models:\n")
  print(dplyr::count(stale, model, name = "rows"), n = 20)
}

is_mm <- stringr::str_starts(current$model, "megalmm")

# ------------------------------------------------------------
# STAGE 2: within the current grid, was this MegaLMM setting in the fraction?
#
# This needs the design, and the design is the part that goes missing. Two
# sources, in order of trustworthiness:
#
#   --design PATH     the simulation_design.csv written BESIDE these results.
#                     Authoritative when it is the right one, and the script
#                     checks that it is before trusting it.
#   --infer-design    the setting list the replicates unanimously share. Sound
#                     once stage 1 has removed the earlier grid, because a
#                     coherent run gives every replicate the same settings.
#
# Nothing here keys on `rep`. simulation_design.csv is rewritten by every
# sim_run.R invocation including --combine, recording whatever --reps THAT one
# was given, so a bare --combine after a --reps 5 grid leaves a rep-1 design
# describing 5-rep results. A design authorises settings; replication is a
# separate axis.
# ------------------------------------------------------------

if (infer) {
  seen <- current[is_mm, ] |>
    dplyr::distinct(rep, dplyr::across(dplyr::all_of(SETTING)))
  n_reps <- dplyr::n_distinct(seen$rep)
  votes <- dplyr::count(seen, dplyr::across(dplyr::all_of(SETTING)),
                        name = "reps_seen")
  design  <- dplyr::filter(votes, reps_seen == n_reps) |>
    dplyr::select(dplyr::all_of(SETTING))
  partial <- dplyr::filter(votes, reps_seen < n_reps)
  cat("\n=== stage 2: inferred design ===\n")
  cat(nrow(design), "setting(s) present in all", n_reps, "replicate(s);",
      nrow(partial), "in only some (treated as leftovers)\n")
} else {
  if (!file.exists(design_path)) stop("not found: ", design_path, call. = FALSE)
  d <- readr::read_csv(design_path, show_col_types = FALSE)
  design <- d |>
    dplyr::mutate(eigen_variance = r6(eigen_variance)) |>
    dplyr::distinct(dplyr::across(dplyr::all_of(SETTING)))
  cat("\n=== stage 2: design from", basename(design_path), "===\n")
  cat(nrow(d), "runs ->", nrow(design), "distinct setting(s)\n")
  if (dplyr::n_distinct(d$rep) < dplyr::n_distinct(current$rep)) {
    cat("NOTE: the design records", dplyr::n_distinct(d$rep),
        "replicate(s) but the results carry", dplyr::n_distinct(current$rep),
        "--\n  rep is ignored in the match, as it must be.\n")
  }
  got <- current[is_mm, ] |>
    dplyr::distinct(dplyr::across(dplyr::all_of(SETTING)))
  matched <- nrow(dplyr::semi_join(got, design, by = SETTING))
  if (nrow(got) > 0 && matched < 0.5 * nrow(got)) {
    warning("this design matches only ", matched, " of ", nrow(got),
            " MegaLMM settings in the results, so it is almost certainly from ",
            "a DIFFERENT fraction -- optFederov's search is not reproducible ",
            "across machines or AlgDesign versions. Fetch the design written ",
            "beside these results, or use --infer-design.",
            call. = FALSE, immediate. = TRUE)
  }
}

scenarios <- unique(design$scenario)

# ------------------------------------------------------------
# Filter
# ------------------------------------------------------------

kept <- dplyr::bind_rows(
  # BGLR half and baselines: authorised by the scenario alone
  current[!is_mm, ] |> dplyr::filter(scenario %in% scenarios),
  # MegaLMM: scenario AND the three levers
  current[is_mm, ] |> dplyr::semi_join(design, by = SETTING)
) |>
  dplyr::arrange(scenario, rep, model)

cat("\n=== kept ===\n")
cat(nrow(kept), "of", nrow(results), "rows,",
    dplyr::n_distinct(kept$scenario), "scenarios,",
    dplyr::n_distinct(kept$rep), "replicate(s)\n")
print(dplyr::count(kept, model, name = "rows"))
cat("\nrows per replicate:\n"); print(dplyr::count(kept, rep, name = "rows"))

# ------------------------------------------------------------
# What went, and why -- the two reasons mean different things
# ------------------------------------------------------------

gone <- dplyr::anti_join(current, kept, by = names(current))
cat("\n=== dropped ===\n")
cat(nrow(stale), "row(s) at stage 1 (earlier grid)\n")
cat(nrow(gone), "row(s) at stage 2 (current grid, setting not in the fraction)\n")
cat(nrow(stale) + nrow(gone), "of", nrow(results), "total\n")

# ------------------------------------------------------------
# The ambiguity this cannot resolve
# ------------------------------------------------------------

dup <- kept |>
  dplyr::count(model, scenario, rep, K, eigen_variance, fixed_main_effect,
               name = "copies") |>
  dplyr::filter(copies > 1)
if (nrow(dup) > 0) {
  cat("\n=== WARNING:", nrow(dup),
      "(model, scenario, rep, setting) combination(s) appear more than once ===\n")
  cat("Two cache files under DIFFERENT seeds both call themselves this",
      "replicate.\nThat is what adding --reps to an existing run does: rep is",
      "the inner index\nof the seed numbering, so rep 1 of --reps 5 is not rep",
      "1 of --reps 1.\n")
  cat("The CSV carries no seed, so they cannot be told apart here. Either\n",
      " - drop the affected replicate, or\n",
      " - clear output/simulation/ on the cluster and re-combine.\n")
  print(dplyr::count(dup, rep, name = "affected_combinations"))
}

# Done LAST, after the diagnostics: the anti_join above keys on every column
# of `current`, so dropping columns before it would make that join fail.
# Columns that belonged to the OLD scorer survive as all-NA, because the
# combined table is a bind_rows union over cache files with different schemas.
# Dropping them is part of the same clean-up: an all-NA r_* column invites
# someone to wonder why their metric is empty.
empty <- names(kept)[purrr::map_lgl(kept, ~ all(is.na(.x)))]
if (length(empty) > 0) {
  cat("\n=== dropping", length(empty), "all-NA column(s) left by the old scorer ===\n")
  cat(paste(empty, collapse = ", "), "\n")
  kept <- dplyr::select(kept, -dplyr::all_of(empty))
}

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

if (in_place) {
  backup <- paste0(tools::file_path_sans_ext(results_path), "_unfiltered.csv")
  file.copy(results_path, backup, overwrite = TRUE)
  readr::write_csv(kept, results_path)
  message("\nwrote ", results_path, " (original kept at ", basename(backup), ")")
} else {
  out <- paste0(tools::file_path_sans_ext(results_path), "_filtered.csv")
  readr::write_csv(kept, out)
  message("\nwrote ", out)
}
