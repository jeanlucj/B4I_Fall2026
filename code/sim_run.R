# ============================================================
# RUN THE SIMULATION DESIGN
#
# Compares the MegaLMM factor framework with the bivariate DGE-IGE framework
# across the design in code/sim_config.R. Every fit is cached, so the design can
# be run in pieces, interrupted, and resumed.
#
#   Rscript code/sim_run.R                  # the design, 1 replicate
#   Rscript code/sim_run.R --reps 5         # 5 replicates
#   Rscript code/sim_run.R --full           # every candidate, not the fraction
#   Rscript code/sim_run.R --runs 200       # a bigger D-optimal fraction
#   Rscript code/sim_run.R --pilot          # a handful of fast runs, to check it works
#   Rscript code/sim_run.R --filter "n_acc == 200 & environment == 'one'"
#   Rscript code/sim_run.R --check          # sanity check, see below
#   Rscript code/sim_run.R --refresh        # ignore the cache
#   Rscript code/sim_run.R --task 3 --ntasks 20   # slice 3 of 20, for a job array
#   Rscript code/sim_run.R --combine        # just rebuild the CSV from the cache
#
# TWO HALVES, CACHED APART. The DGE-IGE half depends only on the simulated
# data, so it is fitted once per data scenario. The MegaLMM half additionally
# depends on K, eigen_variance and fixed_main_effect, and is now fitted in both
# orientations, so it is where the compute goes -- and it is the half the
# D-optimal fraction cuts. By default 120 data scenarios and 150 MegaLMM runs,
# against 960 runs for the full crossing.
#
# HOW TO READ THE RESULTS. Because most scenarios carry only some MegaLMM
# settings, the answer is not a table of cell means. Fit the design model --
# main effects and two-way interactions in the composite factors -- to the
# outcome columns. simulation_design_effects.csv does that as a first pass.
#
# Outputs: output/simulation/<scenario>_rep<k>_s<seed>_<half>.rds
#          output/simulation_results.csv            everything, combined
#          output/simulation_design.csv             the design that was run
#          output/simulation_design_effects.csv     the design model, fitted
#          output/simulation_summary.png            the headline comparison
# ============================================================

library(tidyverse)

here::i_am("code/sim_run.R")

source(here::here("code", "dge_ige_functions.R"))   # read_grm, grm_factor
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
full    <- has_flag("--full")
n_runs  <- as.integer(arg_value("--runs", as.character(SIM_DESIGN_RUNS)))
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
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  scores <- dplyr::bind_rows(
    run_scenario_bglr(sim, seed = 7),
    run_scenario_megalmm(sim, K = SIM_MEGALMM_K,
                         eigen_variance = SIM_EIGEN_VARIANCE, seed = 7,
                         run_dir = file.path(run_dir, "check"))$scores
  )
  key <- c("model", "r_total_oat", "r_gma_oat", "r_int_oat",
           "r_total_pea", "r_oat_prod", "r_oat_assoc",
           "r_pea_prod", "r_pea_assoc")
  print(as.data.frame(scores[, intersect(key, names(scores))]),
        row.names = FALSE, digits = 3)

  pick <- function(model, col) scores[[col]][scores$model == model][1]

  # Three things must hold on a dense, strongly structured matrix, and each
  # catches a different kind of wiring mistake.
  fails <- character(0)

  mm <- pick("megalmm", "r_total_oat"); ad <- pick("additive", "r_total_oat")
  if (is.na(mm) || mm < 0.5) {
    fails <- c(fails, sprintf(
      "MegaLMM reached r_total_oat = %.3f where it should exceed 0.5", mm))
  }

  # BOTH traits must work. A single-trait bug would leave the pea side empty
  # or at chance while the oat side looked fine.
  mm_pea <- pick("megalmm", "r_total_pea")
  if (is.na(mm_pea) || mm_pea < 0.5) {
    fails <- c(fails, sprintf(
      "MegaLMM reached r_total_pea = %.3f; the pea trait is not being fitted",
      mm_pea))
  }

  # All four effects must be recovered. The two ASSOCIATE effects are the ones
  # that come from the other orientation, so a transpose error or a failure to
  # run the second fit shows up here and nowhere else.
  for (col in c("r_oat_prod", "r_oat_assoc", "r_pea_prod", "r_pea_assoc")) {
    v <- pick("dge_ige", col)
    if (is.na(v) || v < 0.5) {
      fails <- c(fails, sprintf("dge_ige %s = %.3f, should exceed 0.5", col, v))
    }
  }

  if (length(fails) > 0) {
    stop("the pipeline is broken, not the method:\n  ",
         paste(fails, collapse = "\n  "), call. = FALSE)
  }
  message("OK: MegaLMM r_total oat ", round(mm, 3), " / pea ",
          round(mm_pea, 3), ", additive oat ", round(ad, 3),
          "; all four effects recovered")
  quit(save = "no")
}

design <- sim_design(
  levels = if (extended) SIM_LEVELS_EXTENDED else SIM_LEVELS,
  n_runs = n_runs, n_reps = n_reps, full = full)

runs      <- design$runs
scenarios <- design$scenarios

message(nrow(scenarios), " data scenario(s), ", nrow(runs), " MegaLMM run(s)",
        if (!full) sprintf(" -- D-optimal fraction, model rank %d of %d terms",
                           design$model_rank, design$model_terms) else
                  " -- full crossing")
if (!full && design$model_rank < design$model_terms) {
  warning("the chosen design cannot estimate every main effect and two-way ",
          "interaction (rank ", design$model_rank, " of ", design$model_terms,
          " terms). Raise --runs.", call. = FALSE)
}

if (pilot) {
  # The cheapest corner that still varies what is being tested: smallest panel,
  # one environment, the sparsest matrix, and both ends of the interaction axis.
  runs <- runs |>
    dplyr::filter(n_acc == 200, environment == "one", sparsity == 0.016,
                  interaction %in% c("none", "f1_i20")) |>
    dplyr::slice_head(n = 4)
  scenarios <- dplyr::semi_join(scenarios, runs, by = c("scenario", "rep"))
}

if (!is.null(filter_expr)) {
  runs <- dplyr::filter(runs, !!rlang::parse_expr(filter_expr))
  scenarios <- dplyr::semi_join(scenarios, runs, by = c("scenario", "rep"))
}

if (nrow(runs) == 0) { message("nothing to run"); quit(save = "no") }

dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(runs, file.path(out_dir, "simulation_design.csv"))

# A job array hands each task a slice of the DATA SCENARIOS, not of the runs,
# so a scenario's simulated experiment is generated once and both halves stay
# in the same task.
if (!is.na(task) && !is.na(ntasks)) {
  keep <- which((seq_len(nrow(scenarios)) - 1L) %% ntasks == (task - 1L))
  scenarios <- scenarios[keep, ]
  runs <- dplyr::semi_join(runs, scenarios, by = c("scenario", "rep"))
  message("task ", task, " of ", ntasks, ": ", nrow(scenarios),
          " scenario(s), ", nrow(runs), " run(s)")
  if (nrow(scenarios) == 0) quit(save = "no")
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

# The seed is part of the key, not just of the contents. `seed` is
# SIM_BASE_SEED + the row number over the WHOLE grid, so adding a level to any
# axis renumbers every row and changes every seed -- while the scenario name,
# which encodes only the design, stays the same. Without the seed here a
# changed grid silently reuses results generated under a seed it would never
# have chosen. With it, those files simply miss and are recomputed.
cache_path <- function(scenario, rep, seed, suffix) {
  file.path(cache_dir,
            paste0(scenario, "_rep", rep, "_s", seed, "_", suffix, ".rds"))
}

#' One data scenario: generate it, fit the DGE-IGE half, and fit whichever
#' MegaLMM settings the design asked for on this scenario.
run_scenario <- function(sc_row, sc_runs) {
  design_cols <- tibble::tibble(
    scenario = sc_row$scenario, rep = sc_row$rep, n_acc = sc_row$n_acc,
    sparsity = sc_row$sparsity, interaction = sc_row$interaction,
    environment = sc_row$environment, n_factors = sc_row$n_factors,
    interaction_pct = sc_row$interaction_pct, n_envs = sc_row$n_envs,
    gxe_cor = sc_row$gxe_cor
  )
  seed <- sc_row$seed
  bglr_file <- cache_path(sc_row$scenario, sc_row$rep, seed, "bglr")
  # argument names must match the columns pmap() hands in
  mm_file <- function(K, eigen_variance, fixed_main_effect) {
    cache_path(sc_row$scenario, sc_row$rep, seed,
               sprintf("mm_K%d_ev%02d_fx%d", K, round(eigen_variance * 100),
                       as.integer(fixed_main_effect)))
  }
  wanted <- purrr::pmap_chr(
    dplyr::select(sc_runs, K, eigen_variance, fixed_main_effect), mm_file)

  if (!refresh && file.exists(bglr_file) && all(file.exists(wanted))) {
    message("  cached: ", sc_row$scenario, " rep ", sc_row$rep)
    return(dplyr::bind_rows(readRDS(bglr_file),
                            purrr::map(wanted, readRDS) |> purrr::list_rbind()))
  }

  message("\n### ", sc_row$scenario, " rep ", sc_row$rep, " ###")

  grms <- sim_grms(sc_row$n_acc, seed = sc_row$n_acc)  # one panel per size
  set.seed(seed)
  sim <- simulate_experiment(
    G_oat = grms$G_oat, G_pea = grms$G_pea,
    sparsity = sc_row$sparsity, n_factors = sc_row$n_factors,
    interaction_pct = sc_row$interaction_pct, n_envs = sc_row$n_envs,
    gxe_cor = sc_row$gxe_cor
  )
  message("  ", nrow(sim$obs), " observations over ", sc_row$n_acc, " x ",
          sc_row$n_acc, " cells (",
          round(nrow(sim$obs) / sc_row$n_acc, 1), " per column)")

  # --- DGE-IGE half: once per scenario, both traits, Multitrait ---
  if (refresh || !file.exists(bglr_file)) {
    bglr <- run_scenario_bglr(sim, seed = seed) |>
      dplyr::bind_cols(design_cols[rep(1, 4), ])
    saveRDS(bglr, bglr_file)
  } else {
    bglr <- readRDS(bglr_file)
  }

  # --- MegaLMM half: only the settings this scenario was assigned ---
  mm <- purrr::pmap(
    dplyr::select(sc_runs, K, eigen_variance, fixed_main_effect),
    function(K, eigen_variance, fixed_main_effect) {
      f <- mm_file(K, eigen_variance, fixed_main_effect)
      if (!refresh && file.exists(f)) return(readRDS(f))

      message("  megalmm K = ", K, ", eigenvectors to ", eigen_variance,
              ", fixed main effect = ", fixed_main_effect,
              " (both orientations)")
      d <- file.path(run_dir, paste0(sc_row$scenario, "_rep", sc_row$rep,
                                     "_K", K, "_fx", as.integer(fixed_main_effect)))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)

      res <- run_scenario_megalmm(
        sim, K = K, eigen_variance = eigen_variance, seed = seed,
        run_dir = d, fixed_main_effect = fixed_main_effect)
      unlink(d, recursive = TRUE)

      out <- dplyr::bind_cols(res$scores, design_cols[rep(1, nrow(res$scores)), ])
      saveRDS(out, f)
      out
    }) |> purrr::list_rbind()

  dplyr::bind_rows(bglr, mm)
}

if (!combine_only) {
  invisible(purrr::pmap(scenarios, function(...) {
    sc_row <- tibble::tibble(...)
    sc_runs <- dplyr::filter(runs, scenario == sc_row$scenario,
                             rep == sc_row$rep)
    if (nrow(sc_runs) == 0) return(NULL)
    run_scenario(sc_row, sc_runs)
  }))
}

# Always rebuild the combined table from the cache rather than from this run:
# with a job array, no single process sees every scenario.
results <- list.files(cache_dir, pattern = "_(bglr|mm_K[0-9]+_ev[0-9]+_fx[01])\\.rds$",
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
#
# The design is fractional, so most scenarios carry only some MegaLMM settings
# and a table of cell means would be comparing unlike with unlike. The summary
# is therefore a MODEL fitted to the outcomes -- main effects and two-way
# interactions in the composite factors -- plus marginal means for reading.
#
# This also removes a bias the old summary had: it took the best of twelve
# MegaLMM settings per scenario, judged on the same held-out cells it then
# reported, and compared that maximum against a DGE-IGE fitted once. There is
# no per-scenario maximum to take any more.
# ------------------------------------------------------------

outcomes <- c("r_total_oat", "r_gma_oat", "r_int_oat",
              "r_total_pea", "r_gma_pea", "r_int_pea",
              "r_oat_prod", "r_oat_assoc", "r_pea_prod", "r_pea_assoc",
              "r_oat_gma", "r_pea_gma")
outcomes <- intersect(outcomes, names(results))

cat("\n=== Marginal means by model ===\n")
results |>
  dplyr::group_by(model) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(outcomes),
                                 \(x) mean(x, na.rm = TRUE)),
                   n = dplyr::n(), .groups = "drop") |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\n=== Where each framework wins, by sparsity and interaction rank ===\n")
results |>
  dplyr::filter(model %in% c("additive", "dge_ige", "megalmm")) |>
  dplyr::group_by(sparsity, interaction, model) |>
  dplyr::summarise(r_total = mean(r_total_oat, na.rm = TRUE),
                   r_int = mean(r_int_oat, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = c(r_total, r_int)) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\n=== Producer effects are kinship-shrunk, associate effects are not ===\n")
cat("    MegaLMM carries the row species' producer effect in latent structure\n")
cat("    and the column species' associate effect in an unshrunk intercept, so\n")
cat("    the associate columns should fall away faster as the matrix thins.\n\n")
results |>
  dplyr::filter(model %in% c("dge_ige", "megalmm")) |>
  dplyr::group_by(model, sparsity) |>
  dplyr::summarise(dplyr::across(dplyr::any_of(c("r_oat_prod", "r_pea_prod",
                                                 "r_oat_assoc", "r_pea_assoc")),
                                 \(x) mean(x, na.rm = TRUE)),
                   .groups = "drop") |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

# --- the design model ---
#
# Fitted only to the MegaLMM rows, because the design's MegaLMM factors (K,
# eigen_variance, fixed_main_effect) do not apply to the other models.
design_vars <- c("n_acc", "sparsity", "interaction", "environment",
                 "K", "eigen_variance", "fixed_main_effect")

mm_rows <- dplyr::filter(results, model == "megalmm")

effects <- if (nrow(mm_rows) > length(design_vars) * 3) {
  purrr::map(outcomes, \(y) {
    d <- mm_rows |>
      dplyr::select(dplyr::all_of(c(y, design_vars))) |>
      tidyr::drop_na() |>
      dplyr::mutate(dplyr::across(dplyr::all_of(design_vars), factor))
    if (nrow(d) < 20 || dplyr::n_distinct(d[[y]]) < 5) return(NULL)
    fit <- stats::lm(stats::reformulate(paste0("(", paste(design_vars, collapse = " + "), ")^2"),
                                        response = y), data = d)
    a <- stats::anova(fit)
    tibble::tibble(outcome = y, term = rownames(a), df = a$Df,
                   mean_sq = a$`Mean Sq`, F = a$`F value`, p = a$`Pr(>F)`) |>
      dplyr::filter(term != "Residuals") |>
      dplyr::arrange(dplyr::desc(F))
  }) |> purrr::compact() |> purrr::list_rbind()
} else NULL

if (!is.null(effects) && nrow(effects) > 0) {
  readr::write_csv(effects, file.path(out_dir, "simulation_design_effects.csv"))
  cat("\n=== Design model: the largest effects on each outcome ===\n")
  effects |>
    dplyr::group_by(outcome) |>
    dplyr::slice_head(n = 3) |>
    dplyr::ungroup() |>
    as.data.frame() |>
    print(row.names = FALSE, digits = 3)
} else {
  message("\ntoo few MegaLMM runs cached to fit the design model yet")
}

p <- results |>
  dplyr::filter(model %in% c("additive", "dge_ige", "megalmm"),
                interaction != "none") |>
  dplyr::mutate(sparsity = factor(paste0(sparsity * 100, "%")),
                panel = factor(paste0(n_acc, " x ", n_acc))) |>
  tidyr::pivot_longer(c(r_total_oat, r_int_oat, r_oat_assoc),
                      names_to = "metric", values_to = "r") |>
  ggplot2::ggplot(ggplot2::aes(sparsity, r, colour = model, group = model)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  ggplot2::stat_summary(fun = mean, geom = "line", na.rm = TRUE) +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 2, na.rm = TRUE) +
  ggplot2::facet_grid(metric ~ panel, scales = "free_y") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
  ggplot2::labs(
    title = "Total genetic value, interaction, and the oat associate effect",
    subtitle = paste0("averaged over the design; ",
                      "r_oat_assoc is the effect MegaLMM holds in an ",
                      "unshrunk intercept"),
    x = NULL, y = "correlation with the truth", colour = NULL
  )

ggplot2::ggsave(file.path(out_dir, "simulation_summary.png"), p,
                width = 10, height = 7, dpi = 150)

message("\nwrote output/simulation_results.csv, simulation_design.csv and ",
        "simulation_summary.png")
