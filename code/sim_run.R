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

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

run_one <- function(scenario, seed, n_acc, sparsity, n_factors,
                    interaction_pct, n_envs, rep) {

  cache_file <- file.path(cache_dir, paste0(scenario, "_rep", rep, ".rds"))
  if (file.exists(cache_file) && !refresh) {
    message("  cached: ", scenario, " rep ", rep)
    return(readRDS(cache_file))
  }

  message("\n### ", scenario, " rep ", rep, " ###")

  grms <- sim_grms(n_acc, seed = n_acc)   # same panel for a given size

  set.seed(seed)
  sim <- simulate_experiment(
    G_oat = grms$G_oat, G_pea = grms$G_pea,
    sparsity = sparsity, n_factors = n_factors,
    interaction_pct = interaction_pct, n_envs = n_envs
  )
  message("  ", nrow(sim$obs), " observations over ",
          n_acc, " x ", n_acc, " cells")

  scenario_run_dir <- file.path(run_dir, paste0(scenario, "_rep", rep))
  dir.create(scenario_run_dir, showWarnings = FALSE, recursive = TRUE)

  scores <- run_scenario(sim, run_dir = scenario_run_dir, seed = seed)

  result <- scores |>
    dplyr::mutate(scenario = scenario, rep = rep, n_acc = n_acc,
                  sparsity = sparsity, n_factors = n_factors,
                  interaction_pct = interaction_pct, n_envs = n_envs,
                  .before = 1)

  saveRDS(result, cache_file)
  unlink(scenario_run_dir, recursive = TRUE)   # MegaLMM run state, regenerable
  result
}

results <- grid |>
  dplyr::select(scenario, seed, n_acc, sparsity, n_factors, interaction_pct,
                n_envs, rep) |>
  purrr::pmap(run_one) |>
  purrr::list_rbind()

readr::write_csv(results, file.path(out_dir, "simulation_results.csv"))

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n=== Accuracy for the true total genetic value ===\n")
results |>
  dplyr::group_by(n_factors, interaction_pct, sparsity, n_acc, n_envs, model) |>
  dplyr::summarise(r = mean(r_total), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = r) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\n=== Accuracy for the interaction alone ===\n")
results |>
  dplyr::filter(n_factors > 0) |>
  dplyr::group_by(n_factors, interaction_pct, sparsity, n_acc, n_envs, model) |>
  dplyr::summarise(r = mean(r_interaction), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = r) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

# --- where does each framework win? ---
head_to_head <- results |>
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

p <- results |>
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
