# ============================================================
# RUN THE SIMULATION GRID
#
# Compares the MegaLMM factor framework with the DGE-IGE framework across the
# factorial in code/sim_config.R. Each scenario is cached, so the grid can be
# run in pieces, interrupted, and resumed.
#
#   Rscript code/sim_run.R                  # the whole grid, 1 replicate
#   Rscript code/sim_run.R --reps 5         # 5 replicates
#   Rscript code/sim_run.R --pilot          # 4 fast scenarios, to check it runs
#   Rscript code/sim_run.R --filter "n_acc == 200 & n_envs == 1"
#   Rscript code/sim_run.R --extended       # sparsity out to 25% and 50%
#   Rscript code/sim_run.R --check          # sanity check, see below
#   Rscript code/sim_run.R --refresh        # ignore the cache
#   Rscript code/sim_run.R --task 3 --ntasks 20   # slice 3 of 20, for a job array
#   Rscript code/sim_run.R --trace          # also record a convergence trace
#   Rscript code/sim_run.R --combine        # just rebuild the CSV from the cache
#
# Outputs: output/simulation/<scenario>_rep<k>.rds   one per scenario
#          output/simulation_results.csv             everything, combined
#          output/simulation_summary.png             the headline comparison
# ============================================================

library(tidyverse)

here::i_am("code/sim_run.R")

source(here::here("code", "sim_config.R"))
source(here::here("code", "sim_generate.R"))
source(here::here("code", "sim_fit.R"))
source(here::here("code", "megalmm_setup.R"))

out_dir   <- here::here("output")
cache_dir <- here::here("output", "simulation")
run_dir   <- here::here("output", "simulation_runs")

# ------------------------------------------------------------
# Arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
has_flag <- function(flag) flag %in% args

n_reps  <- as.integer(arg_value("--reps", "1"))
check   <- has_flag("--check")
extended <- has_flag("--extended")
trace   <- has_flag("--trace")
combine_only <- has_flag("--combine")
task    <- as.integer(arg_value("--task", NA))
ntasks  <- as.integer(arg_value("--ntasks", NA))
refresh <- has_flag("--refresh")
filter_expr <- arg_value("--filter", NULL)
pilot   <- has_flag("--pilot")

# ------------------------------------------------------------
# The grid
# ------------------------------------------------------------

# ------------------------------------------------------------
# Sanity check
#
# A negative result is only worth reporting if the pipeline can produce a
# positive one. This fits a DENSE, strongly structured scenario in which
# MegaLMM should clearly recover the simulated signal, and stops if it does
# not. Run it before trusting any sweep: the failure mode this guards
# against -- a wiring mistake that makes every MegaLMM number approximately
# zero -- looks exactly like "the method does not work here".
# ------------------------------------------------------------

if (check) {
  message("sanity check: dense matrix, 100 x 100 at 50% observed")
  grms <- sim_grms(100, seed = 100)
  set.seed(7)
  sim <- simulate_experiment(grms$G_oat, grms$G_pea,
                             sparsity = 0.50, n_factors = 1,
                             interaction_pct = 0.20, n_envs = 1)
  scores <- run_scenario(sim, run_dir = file.path(run_dir, "check"), seed = 7)
  print(as.data.frame(scores), row.names = FALSE, digits = 3)

  mm <- scores$r_total[scores$model == "megalmm"]
  ad <- scores$r_total[scores$model == "additive"]
  if (is.na(mm) || mm < 0.5) {
    stop("MegaLMM reached r = ", round(mm, 3), " on a dense, strongly ",
         "structured matrix where it should exceed 0.5. The pipeline is ",
         "broken, not the method.", call. = FALSE)
  }
  message("OK: MegaLMM r = ", round(mm, 3), ", additive r = ", round(ad, 3))
  quit(save = "no")
}

grid <- sim_grid(levels = if (extended) SIM_LEVELS_EXTENDED else SIM_LEVELS,
                 n_reps = n_reps)

if (pilot) {
  # The cheapest corner of each axis that still varies the thing being
  # tested: smallest panel, one environment, both interaction ranks that
  # matter, at the sparsity the real experiment sits at.
  grid <- grid |>
    dplyr::filter(n_acc == 200, n_envs == 1, sparsity == 0.03,
                  interaction_pct %in% c(0, 0.20)) |>
    dplyr::slice_head(n = 4)
}

if (!is.null(filter_expr)) {
  grid <- dplyr::filter(grid, !!rlang::parse_expr(filter_expr))
}

message(nrow(grid), " scenario(s) to run")
if (nrow(grid) == 0) quit(save = "no")

dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

# MegaLMM settings are crossed with the data design but cached apart: they
# change nothing about the simulated experiment, so refitting the BGLR half
# for each of them would waste most of the compute.
mm_grid <- tidyr::expand_grid(!!!SIM_MEGALMM_LEVELS)

# A job array hands each task a slice of the data scenarios. Slicing by
# scenario rather than by row keeps a scenario's two halves in one task, so a
# simulation is never generated twice.
if (!is.na(task) && !is.na(ntasks)) {
  keep <- which((seq_len(nrow(grid)) - 1L) %% ntasks == (task - 1L))
  grid <- grid[keep, ]
  message("task ", task, " of ", ntasks, ": ", nrow(grid), " scenario(s)")
  if (nrow(grid) == 0) quit(save = "no")
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

cache_path <- function(scenario, rep, suffix) {
  file.path(cache_dir, paste0(scenario, "_rep", rep, "_", suffix, ".rds"))
}

run_one <- function(scenario, seed, n_acc, sparsity, n_factors,
                    interaction_pct, n_envs, rep) {

  design <- tibble::tibble(scenario = scenario, rep = rep, n_acc = n_acc,
                           sparsity = sparsity, n_factors = n_factors,
                           interaction_pct = interaction_pct, n_envs = n_envs)

  bglr_file <- cache_path(scenario, rep, "bglr")
  mm_files <- purrr::pmap_chr(mm_grid, \(K, eigen_variance)
    cache_path(scenario, rep, sprintf("mm_K%d_ev%02d", K, round(eigen_variance * 100))))

  if (!refresh && file.exists(bglr_file) && all(file.exists(mm_files))) {
    message("  cached: ", scenario, " rep ", rep)
    return(dplyr::bind_rows(readRDS(bglr_file),
                            purrr::map(mm_files, readRDS) |> purrr::list_rbind()))
  }

  message("\n### ", scenario, " rep ", rep, " ###")

  grms <- sim_grms(n_acc, seed = n_acc)   # same panel for a given size
  set.seed(seed)
  sim <- simulate_experiment(
    G_oat = grms$G_oat, G_pea = grms$G_pea,
    sparsity = sparsity, n_factors = n_factors,
    interaction_pct = interaction_pct, n_envs = n_envs
  )
  message("  ", nrow(sim$obs), " observations over ", n_acc, " x ", n_acc,
          " cells (", round(nrow(sim$obs) / n_acc, 1), " per pea column)")

  # --- BGLR half ---
  if (refresh || !file.exists(bglr_file)) {
    bglr <- run_scenario_bglr(sim, seed = seed) |>
      dplyr::bind_cols(design[rep(1, 4), ])
    saveRDS(bglr, bglr_file)
  } else {
    bglr <- readRDS(bglr_file)
  }

  # --- MegaLMM half, once per setting ---
  mm <- purrr::pmap(mm_grid, function(K, eigen_variance) {
    f <- cache_path(scenario, rep, sprintf("mm_K%d_ev%02d", K, round(eigen_variance * 100)))
    if (!refresh && file.exists(f)) return(readRDS(f))

    message("  megalmm K = ", K, ", eigenvectors to ", eigen_variance)
    d <- file.path(run_dir, paste0(scenario, "_rep", rep, "_K", K))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)

    res <- run_scenario_megalmm(
      sim, K = K, eigen_variance = eigen_variance, seed = seed,
      run_dir = d, n_chunks = if (trace) SIM_TRACE_CHUNKS else 1L
    )
    unlink(d, recursive = TRUE)

    out <- dplyr::bind_cols(res$scores, design[rep(1, nrow(res$scores)), ])
    if (!is.null(res$trace)) {
      saveRDS(dplyr::bind_cols(res$trace, design[rep(1, nrow(res$trace)), ],
                               tibble::tibble(K = K, eigen_variance = eigen_variance)),
              sub("\\.rds$", "_trace.rds", f))
    }
    saveRDS(out, f)
    out
  }) |> purrr::list_rbind()

  dplyr::bind_rows(bglr, mm)
}

if (!combine_only) {
  invisible(grid |>
    dplyr::select(scenario, seed, n_acc, sparsity, n_factors, interaction_pct,
                  n_envs, rep) |>
    purrr::pmap(run_one))
}

# Always rebuild the combined table from the cache rather than from this run:
# with a job array, no single process sees every scenario.
results <- list.files(cache_dir, pattern = "_(bglr|mm_K[0-9]+_ev[0-9]+)\\.rds$",
                      full.names = TRUE) |>
  purrr::map(readRDS) |>
  purrr::list_rbind()

if (nrow(results) == 0) {
  message("no cached results yet")
  quit(save = "no")
}

readr::write_csv(results, file.path(out_dir, "simulation_results.csv"))

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

# Best MegaLMM setting per scenario, so the headline tables compare like with
# like rather than averaging over settings that are being swept
best_mm <- results |>
  dplyr::filter(model == "megalmm") |>
  dplyr::group_by(scenario, rep) |>
  dplyr::slice_max(r_interaction, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

results_best <- dplyr::bind_rows(dplyr::filter(results, model != "megalmm"), best_mm)

cat("\n=== MegaLMM settings, averaged over the design ===\n")
results |>
  dplyr::filter(model == "megalmm") |>
  dplyr::group_by(K, eigen_variance) |>
  dplyr::summarise(r_total = mean(r_total),
                   r_interaction = mean(r_interaction, na.rm = TRUE),
                   seconds = mean(seconds, na.rm = TRUE), .groups = "drop") |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\n=== Accuracy for the true total genetic value ===\n")
results_best |>
  dplyr::group_by(n_factors, interaction_pct, sparsity, n_acc, n_envs, model) |>
  dplyr::summarise(r = mean(r_total), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = r) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\n=== Accuracy for the interaction alone ===\n")
results_best |>
  dplyr::filter(n_factors > 0) |>
  dplyr::group_by(n_factors, interaction_pct, sparsity, n_acc, n_envs, model) |>
  dplyr::summarise(r = mean(r_interaction), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = r) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

# --- where does each framework win? ---
head_to_head <- results_best |>
  dplyr::filter(model %in% c("dge_ige", "megalmm")) |>
  dplyr::select(scenario, rep, n_acc, sparsity, n_factors, interaction_pct,
                n_envs, model, r_total, r_interaction) |>
  tidyr::pivot_wider(names_from = model,
                     values_from = c(r_total, r_interaction)) |>
  dplyr::mutate(
    total_gain = r_total_megalmm - r_total_dge_ige,
    int_gain   = r_interaction_megalmm - r_interaction_dge_ige
  )

cat("\n=== MegaLMM minus DGE-IGE (positive = MegaLMM ahead) ===\n")
head_to_head |>
  dplyr::group_by(n_factors, sparsity, n_acc) |>
  dplyr::summarise(total_gain = mean(total_gain),
                   int_gain = mean(int_gain, na.rm = TRUE), .groups = "drop") |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

p <- results_best |>
  dplyr::filter(n_factors > 0) |>
  dplyr::mutate(
    sparsity = factor(paste0(sparsity * 100, "% observed")),
    panel = factor(paste0(n_acc, " x ", n_acc)),
    factors = factor(paste0(n_factors, " factor", ifelse(n_factors > 1, "s", "")))
  ) |>
  ggplot2::ggplot(ggplot2::aes(sparsity, r_interaction, colour = model,
                               group = model)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  ggplot2::stat_summary(fun = mean, geom = "line") +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 2) +
  ggplot2::facet_grid(panel ~ factors) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
  ggplot2::labs(
    title = "Recovering the oat x pea interaction",
    subtitle = "correlation with the simulated interaction at held-out cells",
    x = NULL, y = "correlation", colour = NULL
  )

ggplot2::ggsave(file.path(out_dir, "simulation_summary.png"), p,
                width = 10, height = 6, dpi = 150)

message("\nwrote output/simulation_results.csv and output/simulation_summary.png")
