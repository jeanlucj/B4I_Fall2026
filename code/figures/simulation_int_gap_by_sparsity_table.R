# ============================================================
# Reshape the n_factors facet of the simulation gap-by-sparsity summary into a
# wide table: sparsity in rows; n_factors in the main columns; r_megalmm and
# r_dge_ige nested within each n_factors.
#
#   Rscript code/figures/simulation_int_gap_by_sparsity_table.R
# ============================================================
library(tidyverse)

here::i_am("code/figures/simulation_int_gap_by_sparsity_table.R")

tab <- here::here("output", "simulation_int",
                  "simulation_int_gap_by_sparsity_r_int_oat.csv") |>
  read_csv(show_col_types = FALSE) |>
  filter(facet == "n_factors") |>
  mutate(n_factors = as.integer(level)) |>
  select(sparsity, n_factors, r_megalmm, r_dge_ige) |>
  arrange(n_factors) |>
  pivot_wider(
    names_from = n_factors,
    values_from = c(r_megalmm, r_dge_ige),
    names_glue = "n_factors_{n_factors}_{.value}"
  ) |>
  # Group by n_factors, megalmm then dge_ige within each
  select(sparsity, starts_with("n_factors_1_"), starts_with("n_factors_3_"),
         starts_with("n_factors_5_")) |>
  arrange(sparsity)

write_csv(tab, here::here("output", "figures",
                          "simulation_int_gap_by_sparsity_table.csv"))
