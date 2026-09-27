# test_generator.R
#
# The bivariate simulator in code/sim_generate.R. Every scenario label is a
# claim -- "1.6% observed", "rank 1", "20% interaction variance", "genetic
# correlation 0.6 across environments" -- and if the generator does not deliver
# what the label says then every result read off the design is mislabelled. So
# these tests compare the realised truth against what was REQUESTED, which the
# generator itself never checks.
#
# Run: Rscript tests/test_generator.R

library(tidyverse)
here::i_am("tests/test_generator.R")
source(here::here("tests", "helper.R"))
load_code("sim_config.R", "dge_ige_functions.R", "sim_generate.R")

set.seed(2)

# A small synthetic panel, so the tests need neither data/ nor minutes.
make_panel <- function(n, seed = 1) {
  withr::with_seed(seed, {
    mk <- function(tag) {
      X <- matrix(rnorm(n * n), n, n)
      G <- tcrossprod(X) / n
      diag(G) <- diag(G) + 0.1               # keep it comfortably positive
      dimnames(G) <- list(paste0(tag, seq_len(n)), paste0(tag, seq_len(n)))
      G
    }
    list(G_oat = mk("o"), G_pea = mk("p"))
  })
}
panel <- make_panel(60)

# ------------------------------------------------------------
# 1. draw_effect and draw_effect_pair hit their targets EXACTLY
#
# Both rescale empirically, so a single replicate has the variance and
# covariance the design asked for rather than only having them in expectation.
# That is a design choice worth pinning: it removes a source of between-replicate
# noise, and if it silently stopped holding the scenario labels would drift.
# ------------------------------------------------------------

L <- grm_factor(panel$G_oat)
v <- draw_effect(L, 0.37)
check_near(stats::var(v), 0.37, tol = 1e-10, "draw_effect hits its target variance")
check_near(mean(v), 0, tol = 1e-10, "draw_effect is centred")
check(all(draw_effect(L, 0) == 0), "a zero target gives exactly zero")

Sig <- sigma_from_cor(0.4, 0.9, -0.3)
check_near(Sig, matrix(c(0.4, -0.3 * sqrt(0.36), -0.3 * sqrt(0.36), 0.9), 2, 2),
           tol = 1e-12, "sigma_from_cor builds the right matrix")

pr <- draw_effect_pair(L, Sig)
check(identical(dim(pr), c(nrow(L), 2L)), "draw_effect_pair returns n x 2")
check_near(stats::cov(pr), Sig, tol = 1e-9,
           "draw_effect_pair hits the requested covariance exactly")
check_near(stats::cor(pr[, 1], pr[, 2]), -0.3, tol = 1e-9,
           "and therefore the requested correlation")

# a zero-variance column must not blow up chol()
check_near(stats::var(draw_effect_pair(L, sigma_from_cor(0.5, 0, 0))[, 2]), 0,
           tol = 1e-12, "a zero-variance member of the pair comes back as zero")

# ------------------------------------------------------------
# 2. sample_combinations: the sparsity label is exact
#
# The guarantee of SIM_MIN_PER_ACC observations per accession can overrun a
# request that sits on the floor, which would mislabel the scenario. Asserted in
# the function; asserted here too, at every level actually used.
# ------------------------------------------------------------

# At the REAL panel sizes, because the floor scales with n: a 60-panel cannot
# support 1.6% at all (its floor is 5%). sample_combinations() needs no
# relationship matrix, so testing the real configuration costs nothing.
for (s in SIM_LEVELS$sparsity) {
  for (n in SIM_LEVELS$n_acc) {
    cl <- sample_combinations(n, n, s)
    check(nrow(cl) == round(s * n * n),
          sprintf("sparsity %.3f at n=%d returns exactly the requested count", s, n))
    check(min(table(cl$oat)) >= SIM_MIN_PER_ACC,
          sprintf("every oat meets the minimum (s=%.3f, n=%d)", s, n))
    check(min(table(cl$pea)) >= SIM_MIN_PER_ACC,
          sprintf("every pea meets the minimum (s=%.3f, n=%d)", s, n))
    check(!any(duplicated(cl$cell)), "no cell is observed twice")
  }
}

# below the floor it must refuse, not thin silently
check_error(sample_combinations(200, 200, 0.001),
            "an impossible sparsity is refused")
# ON the floor the forced minimum overruns the request, so the count assertion
# must fire -- this is the guard that stops a scenario mislabelling itself
check_error(sample_combinations(200, 200, SIM_MIN_PER_ACC * 200 / 200^2),
            "a sparsity exactly on the floor is refused rather than mislabelled")
# and 1.6% must clear the floor with room to spare, which is why it was chosen
check(0.016 > SIM_MIN_PER_ACC * 200 / 200^2,
      "the lowest swept sparsity clears the 200-panel floor")

# ------------------------------------------------------------
# 3. simulate_experiment: the variance budget and the four effects
# ------------------------------------------------------------

sim <- simulate_experiment(panel$G_oat, panel$G_pea, sparsity = 0.30,
                           n_factors = 1, interaction_pct = 0.20, n_envs = 1)
tr <- sim$truth

for (trait in c("oat", "pea")) {
  check_near(sum(tr$V[[trait]]), 1, tol = 1e-10,
             sprintf("the %s-yield budget sums to 1", trait))
}

check_near(stats::var(tr$oat_prod),  tr$V$oat[["producer"]],  tol = 1e-9,
           "var(oat producer) is the oat budget's producer share")
check_near(stats::var(tr$pea_assoc), tr$V$oat[["associate"]], tol = 1e-9,
           "var(pea associate) is the oat budget's associate share")
check_near(stats::var(tr$pea_prod),  tr$V$pea[["producer"]],  tol = 1e-9,
           "var(pea producer) is the pea budget's producer share")
check_near(stats::var(tr$oat_assoc), tr$V$pea[["associate"]], tol = 1e-9,
           "var(oat associate) is the pea budget's associate share")

check_near(stats::cor(tr$oat_prod, tr$oat_assoc), SIM_PR_AS_COR[["oat"]],
           tol = 1e-9, "the oat producer-associate correlation is as requested")
check_near(stats::cor(tr$pea_prod, tr$pea_assoc), SIM_PR_AS_COR[["pea"]],
           tol = 1e-9, "the pea producer-associate correlation is as requested")

check(qr(tr$I_oat)$rank == 1, "the oat interaction has the requested rank")
check(qr(tr$I_pea)$rank == 1, "the pea interaction has the requested rank")
check(abs(stats::cor(as.vector(tr$I_oat), as.vector(tr$I_pea))) < 0.15,
      "the two interaction surfaces are independent, as chosen")
check_near(stats::var(as.vector(tr$I_oat)), 0.20, tol = 1e-9,
           "the interaction variance is the requested share")

check(all(c("y_oat", "y_pea") %in% names(sim$obs)),
      "both yields are generated")
check(nrow(sim$obs) == round(0.30 * 60 * 60), "the right number of plots")

# rank 0 must mean NO interaction at all, not a small one
sim0 <- simulate_experiment(panel$G_oat, panel$G_pea, 0.30, 0, 0, 1)
check(all(sim0$truth$I_oat == 0) && all(sim0$truth$I_pea == 0),
      "n_factors = 0 gives an exactly zero interaction")
check_near(sum(sim0$truth$V$oat), 1, tol = 1e-10,
           "and the budget still sums to 1 without it")

# rank 5 must actually be rank 5
sim5 <- simulate_experiment(panel$G_oat, panel$G_pea, 0.30, 5, 0.20, 1)
check(qr(sim5$truth$I_oat)$rank == 5, "n_factors = 5 gives rank 5")

# ------------------------------------------------------------
# 4. Genotype x environment
#
# gxe_cor = 1 must be a no-op, and the realised across-environment correlation
# must match what was asked for. Without the second check the axis could be
# present in name only.
# ------------------------------------------------------------

sim_stable <- simulate_experiment(panel$G_oat, panel$G_pea, 0.30, 1, 0.20, 10,
                                  gxe_cor = 1)
check(is.null(sim_stable$truth$gxe$oat_prod),
      "gxe_cor = 1 draws no environment-specific effects at all")
check_near(sim_stable$truth$V$oat[["producer_env"]], 0, tol = 1e-12,
           "and records zero environment-specific variance")
check(sim_stable$truth$realised_gxe_cor == 1,
      "the realised correlation is exactly 1")

sim_gxe <- simulate_experiment(panel$G_oat, panel$G_pea, 0.30, 1, 0.20, 10,
                               gxe_cor = 0.6)
check_near(sim_gxe$truth$realised_gxe_cor, 0.6, tol = 0.1,
           "gxe_cor = 0.6 realises a correlation near 0.6")
check(!is.null(sim_gxe$truth$gxe$oat_prod),
      "and does draw environment-specific effects")
check_near(sum(sim_gxe$truth$V$oat), 1, tol = 1e-10,
           "splitting off GxE leaves the budget summing to 1")
check_near(sim_gxe$truth$V$oat[["producer"]] +
             sim_gxe$truth$V$oat[["producer_env"]],
           sim_stable$truth$V$oat[["producer"]], tol = 1e-10,
           "and the TOTAL producer variance is unchanged by rho")

# with one environment there is nothing to be specific to, so it is ignored
sim_1env <- simulate_experiment(panel$G_oat, panel$G_pea, 0.30, 1, 0.20, 1,
                                gxe_cor = 0.6)
check(sim_1env$settings$gxe_cor == 1,
      "gxe_cor is forced to 1 when there is only one environment")

check_error(simulate_experiment(panel$G_oat, panel$G_pea, 0.30, 1, 0.20, 10,
                                gxe_cor = 1.5),
            "a gxe_cor outside [0,1] is refused")

finish("generator tests")
