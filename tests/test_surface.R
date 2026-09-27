# test_surface.R
#
# The surface decomposition that every metric in the simulation rests on, and the
# masking that decides which cells are scored.
#
# Every model returns two oat x pea surfaces and each reported effect is a margin
# of one of them, so if the decomposition is wrong then r_gma, r_interaction and
# all four effect-recovery numbers are wrong together, in a way no single
# comparison would reveal. These are algebraic identities, so the tolerances are
# machine precision, not "close enough".
#
# Run: Rscript tests/test_surface.R

library(tidyverse)
here::i_am("tests/test_surface.R")
source(here::here("tests", "helper.R"))
load_code("sim_config.R", "dge_ige_functions.R", "sim_generate.R", "sim_fit.R")

set.seed(3)

# Non-square on purpose: a transposed row/column mean would conform on a square
# matrix and pass.
M <- outer(rnorm(5), rnorm(8), "+") + matrix(rnorm(40), 5, 8)

# ------------------------------------------------------------
# 1. The decomposition is exact and complementary
# ------------------------------------------------------------

check_near(additive_part(M) + interaction_part(M), M, tol = 1e-12,
           "M == additive_part(M) + interaction_part(M)")
check_near(rowMeans(interaction_part(M)), rep(0, nrow(M)), tol = 1e-12,
           "interaction_part has zero row means")
check_near(colMeans(interaction_part(M)), rep(0, ncol(M)), tol = 1e-12,
           "interaction_part has zero column means")
check_near(rowMeans(additive_part(M)), rowMeans(M), tol = 1e-12,
           "additive_part keeps the row margin -- which is why it IS the effect")
check_near(colMeans(additive_part(M)), colMeans(M), tol = 1e-12,
           "additive_part keeps the column margin")

# an additive surface must residualise to EXACTLY zero: that is what makes
# r_interaction come back NA for the additive model rather than a small number
A <- outer(rnorm(5), rnorm(8), "+")
check_near(interaction_part(A), matrix(0, 5, 8), tol = 1e-12,
           "a purely additive surface has no interaction part")
check(stats::sd(as.vector(interaction_part(A))) < 1e-12,
      "and therefore zero variance, so a correlation against it is undefined")

# idempotence: decomposing a part returns it unchanged
check_near(additive_part(additive_part(M)), additive_part(M), tol = 1e-12,
           "additive_part is idempotent")
check_near(interaction_part(interaction_part(M)), interaction_part(M), tol = 1e-12,
           "interaction_part is idempotent")
# and the two parts are orthogonal in the usual inner product
check_near(sum(additive_part(M) * interaction_part(M)), 0, tol = 1e-9,
           "the two parts are orthogonal")

# ------------------------------------------------------------
# 2. Which margin is which effect
#
# The naming convention is the thing most easily got backwards, and getting it
# backwards would swap producer and associate silently. Built from known effects
# so the answer is decided in advance.
# ------------------------------------------------------------

oat_prod <- rnorm(5); pea_assoc <- rnorm(8)
surface_oat <- outer(oat_prod, pea_assoc, "+")
check_near(rowMeans(surface_oat) - mean(rowMeans(surface_oat)),
           oat_prod - mean(oat_prod), tol = 1e-12,
           "the row margin of the oat-yield surface IS the oat producer effect")
check_near(colMeans(surface_oat) - mean(colMeans(surface_oat)),
           pea_assoc - mean(pea_assoc), tol = 1e-12,
           "the column margin of the oat-yield surface IS the pea associate effect")

# ------------------------------------------------------------
# 3. Masking: determinism, floors, disjointness
#
# Both halves of a scenario call split_observations() with the same seed, and
# they are fitted in different processes. If it ever stopped being deterministic
# every head-to-head comparison would silently be between models scored on
# different data.
# ------------------------------------------------------------

make_panel <- function(n, tag) {
  X <- matrix(rnorm(n * n), n, n)
  G <- tcrossprod(X) / n; diag(G) <- diag(G) + 0.1
  dimnames(G) <- list(paste0(tag, seq_len(n)), paste0(tag, seq_len(n)))
  G
}
sim <- simulate_experiment(make_panel(50, "o"), make_panel(50, "p"),
                           sparsity = 0.30, n_factors = 1,
                           interaction_pct = 0.20, n_envs = 1)

a <- split_observations(sim, 0.2, seed = 7)
b <- split_observations(sim, 0.2, seed = 7)
check(identical(a$held$cell, b$held$cell),
      "the same seed gives the same split -- both halves see one split")
check(identical(a$train$cell, b$train$cell), "and the same training set")

c2 <- split_observations(sim, 0.2, seed = 8)
check(!identical(a$held$cell, c2$held$cell), "a different seed gives a different split")

check(length(intersect(a$train$cell, a$held$cell)) == 0,
      "train and held are disjoint")
check(nrow(a$train) + nrow(a$held) == nrow(sim$obs),
      "and together they are everything")
check(min(table(a$train$oat)) >= SIM_FLOOR_OBS,
      "every oat keeps at least the floor in training")
check(min(table(a$train$pea)) >= SIM_FLOOR_OBS,
      "every pea keeps at least the floor in training")
check(abs(nrow(a$held) / nrow(sim$obs) - 0.2) < 0.02,
      "the held-out share is close to the fraction asked for")

check(all(c("y_oat_std", "y_pea_std") %in% names(a$train)),
      "both traits are standardised")

# standardising within environment must leave each environment centred
sim10 <- simulate_experiment(make_panel(50, "o"), make_panel(50, "p"),
                             0.30, 1, 0.20, 10)
st <- standardize_within_env(sim10$obs)
per_env <- st |>
  dplyr::group_by(env) |>
  dplyr::summarise(m_oat = mean(y_oat_std), m_pea = mean(y_pea_std),
                   s_oat = stats::sd(y_oat_std), .groups = "drop")
check_near(per_env$m_oat, rep(0, nrow(per_env)), tol = 1e-10,
           "each environment's oat yield is centred")
check_near(per_env$m_pea, rep(0, nrow(per_env)), tol = 1e-10,
           "each environment's pea yield is centred")
check_near(per_env$s_oat, rep(1, nrow(per_env)), tol = 1e-8,
           "and scaled to unit SD")

# ------------------------------------------------------------
# 4. score_predictions: shape, and that a perfect predictor scores 1
#
# Feeding it the truth is the cheapest possible end-to-end oracle: if the
# plumbing from surfaces to metrics is right, every correlation must be 1.
# ------------------------------------------------------------

tr <- sim$truth
perfect <- list(
  surface = list(
    oat = outer(tr$oat_prod, tr$pea_assoc, "+") + tr$I_oat,
    pea = outer(tr$oat_assoc, tr$pea_prod, "+") + tr$I_pea),
  interaction = NULL)

sc <- score_predictions(perfect, sim, a$held, "perfect")
for (col in c("r_total_oat", "r_gma_oat", "r_int_oat",
              "r_total_pea", "r_gma_pea", "r_int_pea",
              "r_oat_prod", "r_oat_assoc", "r_pea_prod", "r_pea_assoc",
              "r_oat_gma", "r_pea_gma")) {
  check(col %in% names(sc), paste("score_predictions reports", col))
  check_near(sc[[col]], 1, tol = 1e-6,
             paste("a perfect predictor scores 1 on", col))
}

# and a pure-noise predictor must score near zero, not near one
noise <- list(surface = list(oat = matrix(rnorm(50 * 50), 50, 50,
                                         dimnames = dimnames(tr$I_oat)),
                             pea = matrix(rnorm(50 * 50), 50, 50,
                                          dimnames = dimnames(tr$I_pea))),
              interaction = NULL)
sn <- score_predictions(noise, sim, a$held, "noise")
check(abs(sn$r_total_oat) < 0.2, "a noise predictor scores near zero on r_total_oat")
check(abs(sn$r_oat_prod) < 0.3, "and near zero on r_oat_prod")

finish("surface tests")
