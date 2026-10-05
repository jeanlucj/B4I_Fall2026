# ============================================================
# THE SIMULATION'S VARIANCE COMPONENTS, TABULATED
#
#   Rscript code/sim_variance_components.R
#
# One row per trait x interaction level x component: how much variance each
# effect gets in the simulated data, as a share of the trait's total
# PHENOTYPIC variance and of its GENETIC variance (producer + associate +
# interaction; the residual has no genetic share).
#
# The parameters are those in code/sim_config.R (SIM_VAR_SHARES,
# SIM_PR_AS_COR, SIM_RESID_COR) and the interaction levels in
# code/sim_int_config.R. They come from the 2026-09-21 fit and are STALE against
# output/BGLR_variance_components.csv -- see docs/B4I_followups.md.
#
# The simulation's scale is unit total variance per trait, so `variance` equals
# `frac_phenotypic` by construction; both are kept so the file reads either way.
# The numbers are the budget sim_variance_budget() hands to simulate_experiment(),
# which rescales each draw to hit its target variance exactly.
#
# Output: output/simulation_variance_components.csv
# ============================================================

library(tidyverse)

here::i_am("code/sim_variance_components.R")

source(here::here("code", "dge_ige_functions.R"))
source(here::here("code", "sim_config.R"))
source(here::here("code", "sim_int_config.R"))
source(here::here("code", "sim_generate.R"))

budget_rows <- function(interaction_pct) {
  V <- sim_variance_budget(interaction_pct)
  tibble::tribble(
    ~trait, ~component, ~variance,
    "oat_yield", "producer (oat)",    V[["oat_prod"]],
    "oat_yield", "associate (pea)",   V[["pea_assoc"]],
    "oat_yield", "interaction",       interaction_pct,
    "oat_yield", "residual",          V[["e_oat"]],
    "pea_yield", "producer (pea)",    V[["pea_prod"]],
    "pea_yield", "associate (oat)",   V[["oat_assoc"]],
    "pea_yield", "interaction",       interaction_pct,
    "pea_yield", "residual",          V[["e_pea"]]
  ) |>
    dplyr::mutate(interaction_pct = interaction_pct, .before = 1)
}

tab <- SIM_INT_LEVELS$interaction_pct |>
  purrr::map(budget_rows) |>
  purrr::list_rbind() |>
  dplyr::mutate(
    is_genetic = component != "residual",
    .by = c(interaction_pct, trait),
    total = sum(variance),
    genetic_total = sum(variance[is_genetic]),
    frac_phenotypic = variance / total,
    frac_genetic = dplyr::if_else(is_genetic, variance / genetic_total,
                                  NA_real_)
  ) |>
  dplyr::select(interaction_pct, trait, component, variance,
                frac_phenotypic, frac_genetic, genetic_total) |>
  dplyr::mutate(dplyr::across(c(variance, frac_phenotypic, frac_genetic,
                                genetic_total), \(x) round(x, 4)))

readr::write_csv(tab, here::here("output", "simulation_variance_components.csv"))
print(tab, n = Inf)
