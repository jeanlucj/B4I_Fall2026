# ============================================================
# FLOOR FOR THE FIRST-FACTOR CORRELATIONS in docs/interaction_decomposition.md
# section 4.4: what a RANDOM direction in the centred rank-30 oat basis scores
# against a simulated true subspace of 1, 3 or 5 factors at n_acc = 400.
#
#   Rscript code/first_factor_null.R
#
# The floor depends on how many true factors there are: the multiple correlation
# with a k-dimensional span rises with k even for a direction that knows
# nothing, so a k = 5 result cannot be read against the k = 1 floor.
# Output: output/simulation_decomp_first_factor_null.csv
# ============================================================
library(tidyverse)
here::i_am("code/first_factor_null.R")
for (f in c("dge_ige_functions.R","sim_config.R","sim_int_config.R","sim_generate.R","interaction_decomp.R")) source(here::here("code", f))
g <- sim_grms(400, seed = 400)
A <- grm_basis(g$G_oat, rank = 30)
Ac <- sweep(A, 2, colMeans(A), "-")
res <- map_dfr(c(1, 3, 5), \(k) map_dfr(1:60, \(i) {
  set.seed(1000 * k + i)
  sim <- simulate_experiment(g$G_oat, g$G_pea, sparsity = 0.016, n_factors = k,
                             interaction_pct = 0.1, n_envs = 1, gxe_cor = 1)
  v <- Ac %*% rnorm(ncol(Ac))
  ff <- first_factor_cors(v, sim$truth$U_oat)
  tibble(k = k, multiple = ff$multiple, mean_cor = mean(ff$cors), max_cor = max(ff$cors))
}))
null_summary <- res |> group_by(k) |> summarise(across(c(multiple, mean_cor, max_cor), list(med = median, p95 = \(x) quantile(x, .95))))
write_csv(null_summary, here::here("output", "simulation_decomp_first_factor_null.csv"))
print(null_summary, width = 200)
