# ============================================================
# RUN THE INTERACTION-FOCUSED SIMULATION
#
#   Rscript code/sim_int_run.R --check                 # sanity check, first
#   Rscript code/sim_int_run.R --pilot                 # four cheap runs
#   Rscript code/sim_int_run.R                         # 120 scenarios x 3 reps
#   Rscript code/sim_int_run.R --reps 5                # additive: reuses 1-3
#   Rscript code/sim_int_run.R --filter "sparsity >= 0.09"
#   Rscript code/sim_int_run.R --task 3 --ntasks 20    # one slice, job array
#   Rscript code/sim_int_run.R --combine               # rebuild the CSV
#
# The question: WHERE are the two boundaries at which MegaLMM starts to beat the
# Kronecker kernel at recovering the oat x pea interaction? The main design
# (code/sim_run.R) established that both exist -- a density threshold between
# 4.8% and 16% observed, and a rank threshold between 1 and 5 -- but could not
# locate either, because its MegaLMM sweep is a D-optimal fraction and so cannot
# support questions conditioned on a particular configuration.
#
# This design is a FULL FACTORIAL over the five factors that matter, with the
# pinning switch crossed rather than fixed. See code/sim_int_config.R for what
# was changed and why.
#
# Everything downstream is shared with the main simulation: the same generator,
# the same fitters, the same scorer, the same surfaces. Only the design and the
# output names differ -- so the two are comparable, and nothing is reimplemented.
#
# Outputs: output/simulation_int/<scenario>_rep<k>_s<seed>_<half>.rds
#          output/simulation_int_results.csv
#          output/simulation_int_design.csv
# ============================================================

library(tidyverse)

here::i_am("code/sim_int_run.R")

source(here::here("code", "dge_ige_functions.R"))
source(here::here("code", "sim_config.R"))        # SIM_SCORE_SET, SIM_BGLR_*
source(here::here("code", "sim_int_config.R"))
source(here::here("code", "sim_generate.R"))
source(here::here("code", "sim_fit.R"))
source(here::here("code", "megalmm_setup.R"))

out_dir   <- here::here("output")
cache_dir <- here::here("output", "simulation_int")
run_dir   <- here::here("output", "simulation_int_runs")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
has_flag <- function(flag) flag %in% args

n_reps       <- as.integer(arg_value("--reps", as.character(SIM_INT_REPS)))
check        <- has_flag("--check")
pilot        <- has_flag("--pilot")
combine_only <- has_flag("--combine")
refresh      <- has_flag("--refresh")
task         <- as.integer(arg_value("--task", NA))
ntasks       <- as.integer(arg_value("--ntasks", NA))
filter_expr  <- arg_value("--filter", NULL)

# ------------------------------------------------------------
# Sanity check
#
# Same purpose as the main design's: a negative result is only worth reporting
# if the pipeline can produce a positive one. Here the positive control is the
# one this design exists to measure -- a dense, rank-1 scenario, where MegaLMM
# should beat the Kronecker kernel on the interaction. If it does not, the
# comparison is mis-wired and every boundary this design reports is an artefact.
# ------------------------------------------------------------

if (check) {
  message("dense rank-1 scenario: MegaLMM should BEAT dge_ige on r_int")
  grms <- sim_grms(120L, seed = 120L)
  set.seed(1L)
  sim <- simulate_experiment(grms$G_oat, grms$G_pea, sparsity = 0.48,
                             n_factors = 1, interaction_pct = 0.20, n_envs = 1)
  sp <- prepare_scenario(sim)

  d <- fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, TRUE,
                   nIter = 3000L, burnIn = 600L, seed = 1L)
  mm <- fit_megalmm_both(sp$train, sim$G_oat, sim$G_pea,
                         runID = file.path(tempdir(), "int_check"),
                         K = SIM_INT_MEGALMM$K,
                         eigen_variance = SIM_INT_MEGALMM$eigen_variance,
                         fixed_main_effect = FALSE)
  sd_ <- score_predictions(d,  sim, sp$idx, "dge_ige", train = sp$train)
  sm  <- score_predictions(mm, sim, sp$idx, "megalmm", train = sp$train)
  print(dplyr::bind_rows(sd_, sm) |>
          dplyr::select(model, r_int_oat, r_int_pea, r_gma_oat, r_gma_pea))

  ok_int <- sm$r_int_oat > sd_$r_int_oat
  ok_gma <- sd_$r_gma_oat > 0.7
  if (!ok_int) {
    stop("MegaLMM did NOT beat dge_ige on r_int in the dense rank-1 case ",
         sprintf("(%.3f vs %.3f). The comparison is mis-wired.",
                 sm$r_int_oat, sd_$r_int_oat), call. = FALSE)
  }
  if (!ok_gma) {
    stop(sprintf("dge_ige recovered GMA at only %.3f in a dense scenario.",
                 sd_$r_gma_oat), call. = FALSE)
  }
  message(sprintf("OK: r_int megalmm %.3f > dge_ige %.3f; dge_ige r_gma %.3f",
                  sm$r_int_oat, sd_$r_int_oat, sd_$r_gma_oat))
  quit(save = "no")
}

# ------------------------------------------------------------
# The design
# ------------------------------------------------------------

scenarios <- sim_int_grid(n_reps = n_reps)
runs      <- sim_int_runs(scenarios)

message(nrow(scenarios), " data scenario(s), ", nrow(runs),
        " MegaLMM run(s) -- FULL factorial, ", n_reps, " replicate(s)")

if (pilot) {
  scenarios <- scenarios |>
    dplyr::filter(rep == 1, n_acc == 200, sparsity == 0.09, n_envs == 1) |>
    dplyr::slice_head(n = 2)
  runs <- dplyr::semi_join(runs, scenarios, by = c("scenario", "rep"))
}

if (!is.null(filter_expr)) {
  scenarios <- dplyr::filter(scenarios, !!rlang::parse_expr(filter_expr))
  runs <- dplyr::semi_join(runs, scenarios, by = c("scenario", "rep"))
}

if (nrow(scenarios) == 0) { message("nothing to run"); quit(save = "no") }

dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

# Only a run writes the design. `--combine` fits nothing, so it has no business
# recording one -- the mistake that left a one-replicate design beside
# five-replicate results in the main simulation.
if (!combine_only) {
  readr::write_csv(runs, file.path(out_dir, "simulation_int_design.csv"))
}

if (!is.na(task) && !is.na(ntasks)) {
  keep <- which((seq_len(nrow(scenarios)) - 1L) %% ntasks == (task - 1L))
  scenarios <- scenarios[keep, ]
  runs <- dplyr::semi_join(runs, scenarios, by = c("scenario", "rep"))
  message("task ", task, " of ", ntasks, ": ", nrow(scenarios), " scenario(s)")
  if (nrow(scenarios) == 0) quit(save = "no")
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

# The scoring scheme is part of the cache key.
#
# Files written before 30 September 2026 hold scores computed on a HELD-OUT 20%
# of the observed cells, with the models fitted at 0.8x the labelled sparsity.
# They are not comparable with anything written since, and a globbed CSV would
# mix the two silently -- which is the class of bug this project has already been
# bitten by twice. Marking the scheme in the filename means old files simply miss
# and are recomputed, and the combine glob below will not pick them up.
SIM_CACHE_SCHEME <- "v2"

cache_path <- function(scenario, rep, seed, suffix) {
  file.path(cache_dir,
            paste0(scenario, "_rep", rep, "_s", seed, "_",
                   SIM_CACHE_SCHEME, "_", suffix, ".rds"))
}

run_one <- function(sc_row, sc_runs) {
  design_cols <- tibble::tibble(
    scenario = sc_row$scenario, rep = sc_row$rep, seed = sc_row$seed,
    n_acc = sc_row$n_acc, sparsity = sc_row$sparsity,
    n_factors = sc_row$n_factors, interaction_pct = sc_row$interaction_pct,
    n_envs = sc_row$n_envs, gxe_cor = sc_row$gxe_cor)

  bglr_file <- cache_path(sc_row$scenario, sc_row$rep, sc_row$seed, "bglr")
  mm_file <- function(K, eigen_variance, fixed_main_effect) {
    cache_path(sc_row$scenario, sc_row$rep, sc_row$seed,
               sprintf("mm_K%d_ev%02d_fx%d", K, round(eigen_variance * 100),
                       as.integer(fixed_main_effect)))
  }
  wanted <- purrr::pmap_chr(
    dplyr::select(sc_runs, K, eigen_variance, fixed_main_effect), mm_file)

  if (!refresh && file.exists(bglr_file) && all(file.exists(wanted))) {
    message("  cached: ", sc_row$scenario, " rep ", sc_row$rep)
    return(invisible(NULL))
  }

  message("\n### ", sc_row$scenario, " rep ", sc_row$rep, " ###")

  grms <- sim_grms(sc_row$n_acc, seed = sc_row$n_acc)
  set.seed(sc_row$seed)
  sim <- simulate_experiment(
    G_oat = grms$G_oat, G_pea = grms$G_pea,
    sparsity = sc_row$sparsity, n_factors = sc_row$n_factors,
    interaction_pct = sc_row$interaction_pct, n_envs = sc_row$n_envs,
    gxe_cor = sc_row$gxe_cor)

  if (refresh || !file.exists(bglr_file)) {
    # `additive` is deliberately not fitted -- see SIM_INT_BGLR_MODELS
    b <- run_scenario_bglr(sim, SIM_SCORE_SET, seed = sc_row$seed,
                           models = SIM_INT_BGLR_MODELS)
    saveRDS(dplyr::bind_cols(b, design_cols[rep(1, nrow(b)), ]), bglr_file)
  }

  purrr::pwalk(dplyr::select(sc_runs, K, eigen_variance, fixed_main_effect),
               function(K, eigen_variance, fixed_main_effect) {
    f <- mm_file(K, eigen_variance, fixed_main_effect)
    if (!refresh && file.exists(f)) return(invisible(NULL))
    rd <- file.path(run_dir, paste0(sc_row$scenario, "_", sc_row$rep, "_",
                                    K, "_", round(eigen_variance * 100), "_",
                                    as.integer(fixed_main_effect)))
    dir.create(rd, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(rd, recursive = TRUE), add = TRUE)
    out <- run_scenario_megalmm(sim, K = K, eigen_variance = eigen_variance,
                                score_set = SIM_SCORE_SET,
                                seed = sc_row$seed, run_dir = rd,
                                fixed_main_effect = fixed_main_effect)$scores
    saveRDS(dplyr::bind_cols(out, design_cols[rep(1, nrow(out)), ]), f)
  })
  invisible(NULL)
}

if (!combine_only) {
  invisible(purrr::pmap(scenarios, function(...) {
    sc_row <- tibble::tibble(...)
    sc_runs <- dplyr::filter(runs, scenario == sc_row$scenario,
                             rep == sc_row$rep)
    if (nrow(sc_runs) == 0) return(NULL)
    run_one(sc_row, sc_runs)
  }))
}

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

cache_files <- list.files(cache_dir,
                          pattern = paste0("_", SIM_CACHE_SCHEME,
                                            "_(bglr|mm_K[0-9]+_ev[0-9]+_fx[01])\\.rds$"),
                          full.names = TRUE)
if (length(cache_files) == 0) { message("no cached results yet"); quit(save = "no") }

results <- purrr::map(cache_files, function(f) {
  x <- readRDS(f)
  if (!"seed" %in% names(x) || all(is.na(x$seed))) {
    x$seed <- as.integer(stringr::str_match(basename(f), "_s(\\d+)_")[, 2])
  }
  x
}) |> purrr::list_rbind()

readr::write_csv(results, file.path(out_dir, "simulation_int_results.csv"))
message("\n", nrow(results), " rows -> output/simulation_int_results.csv")

# ------------------------------------------------------------
# The two boundaries, straight out of the cache
# ------------------------------------------------------------

fz <- function(r) atanh(pmax(pmin(r, 1 - 1e-8), -1 + 1e-8))

gap <- results |>
  dplyr::filter(model %in% c("dge_ige", "megalmm")) |>
  dplyr::select(model, scenario, rep, n_acc, sparsity, n_factors,
                interaction_pct, n_envs, fixed_main_effect,
                r_int_oat, r_int_pea, r_gma_oat, r_gma_pea) |>
  tidyr::pivot_longer(dplyr::starts_with("r_"), names_to = "response",
                      values_to = "r") |>
  dplyr::filter(!is.na(r)) |>
  dplyr::mutate(z = fz(r))

bg <- gap |> dplyr::filter(model == "dge_ige") |>
  dplyr::select(scenario, rep, response, z_bglr = z)
mm <- gap |> dplyr::filter(model == "megalmm") |>
  dplyr::select(-model) |> dplyr::rename(z_mm = z)

paired <- dplyr::inner_join(mm, bg, by = c("scenario", "rep", "response")) |>
  dplyr::mutate(d = z_mm - z_bglr)

if (nrow(paired) > 0) {
  cat("\n=== MegaLMM minus dge_ige on r_int_oat, by density and rank ===\n")
  cat("(positive favours MegaLMM; `wins` is the share of scenarios)\n\n")
  print(paired |>
          dplyr::filter(response == "r_int_oat") |>
          dplyr::group_by(sparsity, n_factors, fixed_main_effect) |>
          dplyr::summarise(n = dplyr::n(), mean_d = round(mean(d), 3),
                           wins = paste0(round(100 * mean(d > 0)), "%"),
                           .groups = "drop") |>
          tidyr::pivot_wider(names_from = n_factors,
                             values_from = c(mean_d, wins, n)),
        n = 30, width = 200)

  cat("\n=== does pinning help? by response and density ===\n")
  print(paired |>
          dplyr::group_by(response, sparsity, fixed_main_effect) |>
          dplyr::summarise(mean_r = round(mean(tanh(z_mm)), 3), .groups = "drop") |>
          tidyr::pivot_wider(names_from = fixed_main_effect,
                             values_from = mean_r,
                             names_prefix = "pin_") |>
          dplyr::mutate(gain_from_pinning = round(pin_TRUE - pin_FALSE, 3)),
        n = 40)

  readr::write_csv(paired, file.path(out_dir, "simulation_int_paired.csv"))
}
