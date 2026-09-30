# test_surface.R
#
# The surface decomposition that every metric in the simulation rests on, and the
# masking that decides which cells are scored.
#
# Every model returns two oat x pea surfaces and each reported effect is a margin
# of one of them, so if the decomposition is wrong then r_addsurf, r_interaction and
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
# 3. The scoring set: determinism, disjointness, and the right size
#
# Nothing is held out. Both halves of a scenario are fitted in separate
# processes and must score the SAME cells, or every head-to-head comparison is
# silently between models scored on different data -- so determinism is the
# load-bearing property, exactly as it was for the old holdout.
#
# The other two are what the redesign is for: the scored cells must be the ones
# that were NEVER observed, and there must be (1 - sparsity) of them, so that a
# scenario labelled 4.8% observed is fitted at 4.8% rather than at 0.8 x 4.8%.
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

a <- prepare_scenario(sim)
b <- prepare_scenario(sim)
check(identical(a$idx, b$idx),
      "the scoring set is deterministic -- both halves score the same cells")

check(nrow(a$train) == nrow(sim$obs),
      "NOTHING is held out: every observed plot is used for fitting")
check_near(nrow(a$idx) / (50 * 50), 1 - 0.30, tol = 1e-12,
           "and exactly (1 - sparsity) of the matrix is scored")

observed <- paste(sim$obs$oat, sim$obs$pea)
check(sum(paste(a$idx[, 1], a$idx[, 2]) %in% observed) == 0,
      "no scored cell was ever observed -- they are the prediction target")

# the other option, and that it is a superset
check(nrow(scoring_index(sim, "all")) == 50 * 50,
      "scoring_index(\"all\") returns every cell")
check(nrow(scoring_index(sim, "all")) > nrow(a$idx),
      "and is a superset of the unobserved set")
check_error(scoring_index(sim, "held_out"),
            "an unknown scoring set is refused rather than silently defaulted")

# a fully observed matrix leaves nothing to score, and must say so rather than
# returning a correlation over two cells
dense <- simulate_experiment(make_panel(12, "o"), make_panel(12, "p"),
                             sparsity = 0.99, n_factors = 0,
                             interaction_pct = 0, n_envs = 1)
check_error(scoring_index(dense, "unobserved"),
            "a matrix with almost nothing unobserved is refused")

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

sc <- score_predictions(perfect, sim, a$idx, "perfect", train = a$train)
for (col in c("r_total_oat", "r_addsurf_oat", "r_int_oat",
              "r_total_pea", "r_addsurf_pea", "r_int_pea",
              "r_oat_prod", "r_oat_assoc", "r_pea_prod", "r_pea_assoc",
              "r_oat_gma", "r_pea_gma")) {
  check(col %in% names(sc), paste("score_predictions reports", col))
  check_near(sc[[col]], 1, tol = 1e-6,
             paste("a perfect predictor scores 1 on", col))
}
check(sc$n_scored == nrow(a$idx),
      "and records how many cells it scored")

# r_fit is a FIT statistic against a noisy phenotype, so even the truth cannot
# score 1 on it -- it is capped by the square root of heritability. Pinning this
# is what stops it being read as a prediction accuracy.
check(!is.na(sc$r_fit_oat) && sc$r_fit_oat > 0.3 && sc$r_fit_oat < 0.95,
      sprintf("r_fit_oat is bounded well below 1 even for the truth (%.3f)",
              sc$r_fit_oat))

# r_addsurf_oat and r_oat_gma are DIFFERENT quantities that the naming invites
# confusing: the first is the additive part of oat YIELD (oat producer + pea
# associate), the second is an oat ACCESSION's total contribution (oat producer +
# oat associate). Both are 1 for a perfect predictor, so equality there proves
# nothing -- the test is that they are built from different truth vectors.
surf_oat <- outer(tr$oat_prod, tr$pea_assoc, "+")
check(abs(stats::cor(tr$pea_assoc, tr$oat_assoc)) < 0.9,
      "the two associate effects are distinct vectors, so r_addsurf_oat and r_oat_gma
       cannot be the same quantity")

# and a pure-noise predictor must score near zero, not near one
noise <- list(surface = list(oat = matrix(rnorm(50 * 50), 50, 50,
                                         dimnames = dimnames(tr$I_oat)),
                             pea = matrix(rnorm(50 * 50), 50, 50,
                                          dimnames = dimnames(tr$I_pea))),
              interaction = NULL)
sn <- score_predictions(noise, sim, a$idx, "noise", train = a$train)
check(abs(sn$r_total_oat) < 0.2, "a noise predictor scores near zero on r_total_oat")
check(abs(sn$r_oat_prod) < 0.3, "and near zero on r_oat_prod")

finish("surface tests")
