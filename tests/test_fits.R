# test_fits.R -- SLOW tier, excluded from tests/run_all.R unless --all is given.
#
# Everything above this file tests the plumbing: identities, labels, balance.
# This one tests the claim the whole simulation exists to make -- that a model
# given kinship and a partner structure recovers the four effects better than
# not modelling them at all -- by fitting real chains to a small dense scenario
# where the truth is known.
#
# Three things make it a test rather than a demonstration:
#
#   * the ORACLE is the simulated truth, not a stored number, so it does not
#     freeze whatever the code happens to do today;
#   * it includes a NULL: a scenario with no genetic signal must score near
#     zero. A recovery test alone passes on code that leaks the answer;
#   * the thresholds are deliberately loose. They are set to catch a broken
#     model -- a swapped effect, a lost kinship term, a transposed surface --
#     not to pin a performance level. A tight threshold on a short MCMC chain
#     is a test that fails for the wrong reason.
#
# Chains are far shorter than production (SIM_BGLR_NITER), so the numbers here
# are not the simulation's numbers. Minutes, not seconds.
#
# Run: Rscript tests/test_fits.R

library(tidyverse)
here::i_am("tests/test_fits.R")
source(here::here("tests", "helper.R"))
load_code("sim_config.R", "dge_ige_functions.R", "megalmm_setup.R",
          "sim_generate.R", "sim_fit.R")

# Seeded at the top: these thresholds are loose but not infinitely loose, and an
# unseeded suite that fails one run in twenty teaches people to ignore it.
set.seed(20260927L)

# Short chains, and a panel small enough that Multitrait finishes.
N_ITER <- 3000L
BURN   <- 600L
N_ACC  <- 40L

make_panel <- function(n, tag, seed) {
  withr::with_seed(seed, {
    X <- matrix(stats::rnorm(n * n), n, n)
    G <- tcrossprod(X) / n
    diag(G) <- diag(G) + 0.1
    dimnames(G) <- list(paste0(tag, seq_len(n)), paste0(tag, seq_len(n)))
    G
  })
}
G_oat <- make_panel(N_ACC, "o", 11)
G_pea <- make_panel(N_ACC, "p", 12)

# Dense on purpose: 48% observed is the easiest cell in the design, so a
# failure here is a failure of the model and not of identifiability.
sim <- simulate_experiment(G_oat, G_pea, sparsity = 0.48, n_factors = 1,
                           interaction_pct = 0.20, n_envs = 1)
sp  <- prepare_scenario(sim)

check(nrow(sp$train) == nrow(sim$obs),
      "nothing is held out: the models are fitted at the labelled sparsity")
check_near(nrow(sp$idx) / (N_ACC * N_ACC), 1 - 0.48, tol = 1e-12,
           "and scored on exactly the unobserved share of the matrix")

cat("fitting the additive model...\n")
add <- fit_dge_ige(sp$train, G_oat, G_pea, with_interaction = FALSE,
                   nIter = N_ITER, burnIn = BURN, seed = 1L)
cat("fitting the DGE-IGE model with the interaction term...\n")
dge <- fit_dge_ige(sp$train, G_oat, G_pea, with_interaction = TRUE,
                   kron_rank = 12L, nIter = N_ITER, burnIn = BURN, seed = 1L)

base <- baseline_predictions(sp$train, N_ACC, N_ACC,
                            rownames(G_oat), rownames(G_pea))
# `both_means` is the baseline to beat. `row_mean` has no column effect at all,
# so its column-margin metrics are NA by construction rather than poor -- which
# is correct, but makes it the wrong comparator.

s_add  <- score_predictions(add, sim, sp$idx, "additive", train = sp$train)
s_dge  <- score_predictions(dge, sim, sp$idx, "dge_ige", train = sp$train)
s_mean <- score_predictions(base$both_means, sim, sp$idx, "both_means", train = sp$train)
s_row  <- score_predictions(base$row_mean, sim, sp$idx, "row_mean", train = sp$train)
check(is.na(s_row$r_pea_assoc) && is.na(s_row$r_pea_prod),
      "the row-mean baseline has no column effect, so its column metrics are NA")
check(!is.na(s_mean$r_pea_assoc) && !is.na(s_mean$r_pea_prod),
      "and the both-means baseline does have them")

# ------------------------------------------------------------
# 1. Shape: both surfaces, right dimensions, right names
#
# A transposed surface is the failure mode that would still produce plausible
# numbers, because the panels here are square. Checked on the names, not the
# dimensions, so a square panel cannot hide it.
# ------------------------------------------------------------

for (nm in c("oat", "pea")) {
  check(identical(dim(add$surface[[nm]]), c(N_ACC, N_ACC)),
        paste("the", nm, "surface is n_oat x n_pea"))
  check(identical(rownames(add$surface[[nm]]), rownames(G_oat)),
        paste("the", nm, "surface has OATS on the rows"))
  check(identical(colnames(add$surface[[nm]]), rownames(G_pea)),
        paste("and PEAS on the columns"))
}
check(!is.null(dge$interaction), "the interaction term returns a surface")
check(is.null(add$interaction),
      "and the additive model returns none, rather than a zero one")

# ------------------------------------------------------------
# 2. Recovery: all four effects, better than the cell mean
#
# The four are asserted separately. A model that recovered producers and lost
# associates would still look fine on r_total, and the associate effects are
# the ones the project is about.
# ------------------------------------------------------------

for (col in c("r_oat_prod", "r_oat_assoc", "r_pea_prod", "r_pea_assoc")) {
  check(s_dge[[col]] > 0.3,
        sprintf("DGE-IGE recovers %s (got %.2f)", col, s_dge[[col]]))
  # At 48% observed the cell-mean baseline is already very good at the MAIN
  # effects -- measured at r = 0.89 on r_pea_prod against the model's 0.89 -- so
  # the assertion is "not worse", not "better". Where the model has to win is
  # the sparse cells and the interaction, tested elsewhere.
  check(s_dge[[col]] >= s_mean[[col]] - 0.05,
        sprintf("and is no worse than the cell-mean baseline on %s (%.2f vs %.2f)",
                col, s_dge[[col]], s_mean[[col]]))
}
for (col in c("r_total_oat", "r_total_pea", "r_addsurf_oat", "r_addsurf_pea")) {
  check(s_dge[[col]] > 0.3, sprintf("and recovers %s (got %.2f)", col, s_dge[[col]]))
}

# the interaction term must buy something on the interaction, and nothing much
# on the additive part -- that is what makes it an interaction term
check(!is.na(s_dge$r_int_oat) && s_dge$r_int_oat > 0,
      sprintf("the interaction term correlates with the true interaction (%.2f)",
              s_dge$r_int_oat))
check(is.na(s_add$r_int_oat) || abs(s_add$r_int_oat) < 1e-8,
      "the additive model has no interaction to score")
check(abs(s_dge$r_addsurf_oat - s_add$r_addsurf_oat) < 0.25,
      "adding the interaction does not move the additive part much")

# ------------------------------------------------------------
# 3. The variance components point at the simulated values
#
# Loose: a short chain on 40 accessions will not pin a variance ratio. What is
# asserted is the ORDERING and the SIGN, which is what a swapped or mis-scaled
# term would break.
# ------------------------------------------------------------

vc <- add$varcomp
check(all(c("var_Pr", "var_As") %in% names(vc)),
      "the fit reports producer and associate variances")
g_oat_row <- dplyr::filter(vc, term == "G_oat")
g_pea_row <- dplyr::filter(vc, term == "G_pea")
check(nrow(g_oat_row) == 1 && nrow(g_pea_row) == 1,
      "one row per species")
check(g_oat_row$var_Pr > 0 && g_pea_row$var_Pr > 0,
      "both producer variances are positive")

# the simulated producer-associate correlations are both negative, and are not
# symmetric between species -- the fit should not invent a positive one
check(g_oat_row$cor_PrAs < 0.4 && g_pea_row$cor_PrAs < 0.4,
      sprintf("neither fitted cor_PrAs is strongly positive (%.2f, %.2f)",
              g_oat_row$cor_PrAs, g_pea_row$cor_PrAs))

# ------------------------------------------------------------
# 4. The NULL: the fit must not correlate with a truth it never saw
#
# The load-bearing test of the pair. A recovery test alone passes on code that
# leaks the answer -- if `score_predictions` reached into `sim$truth` on both
# sides, every correlation above would be 1 and nothing would complain.
#
# Note what the null is NOT. Simulating a scenario with near-zero genetic
# variance does not work: correlation is scale-free, so with 48% of the matrix
# observed and G known, the fit recovers the DIRECTION of a 0.1%-variance effect
# at r = 0.5 while shrinking its size to nothing. Measured, not assumed --
# var_shares of 0.001 gave r_pea_prod = 0.56. A small variance is not a null.
#
# The null that does work is an independent truth: fit on one scenario, score
# against another drawn with a different seed. Same code path, same scorer, and
# the only way to score above noise is a leak.
# ------------------------------------------------------------

# One independent truth is not enough to be a null. With 40 accessions and a G
# whose effective df is ~22, two INDEPENDENT draws correlate at |r| up to 0.55 by
# chance (measured: median 0.13, max 0.55 over 40 pairs). So the null is averaged
# over several draws, where it concentrates on zero and a leak cannot hide.
null_cols <- c("r_oat_prod", "r_oat_assoc", "r_pea_prod", "r_pea_assoc",
               "r_addsurf_oat", "r_addsurf_pea")
null_runs <- purrr::map(seq_len(8), \(i) {
  other <- simulate_experiment(G_oat, G_pea, sparsity = 0.48, n_factors = 1,
                               interaction_pct = 0.20, n_envs = 1)
  # the same fit and the same scored cells; only the TRUTH differs
  score_predictions(dge, other, sp$idx, "independent_truth", train = sp$train)
}) |> purrr::list_rbind()

for (col in null_cols) {
  m <- mean(abs(null_runs[[col]]))
  check(m < 0.3,
        sprintf("independent truths average near zero on %s (mean |r| = %.2f)",
                col, m))
  check(m < abs(s_dge[[col]]) / 2,
        sprintf("and far below the real truth on %s (%.2f vs %.2f)",
                col, m, s_dge[[col]]))
}

# a permutation of the real truth, which controls for its marginal distribution
perm <- sim
perm$truth$oat_prod  <- sim$truth$oat_prod[sample(N_ACC)]
perm$truth$pea_assoc <- sim$truth$pea_assoc[sample(N_ACC)]
s_perm <- score_predictions(dge, perm, sp$idx, "permuted", train = sp$train)
check(abs(s_perm$r_oat_prod) < s_dge$r_oat_prod,
      "a permuted truth scores worse than the real one on r_oat_prod")
check(abs(s_perm$r_pea_assoc) < s_dge$r_pea_assoc,
      "and on r_pea_assoc")

# and the fit really was estimated from data: a scenario with almost no genetic
# variance must produce almost no FITTED genetic variance, even though its
# correlations stay moderate for the reason above
cat("fitting a near-null scenario for its variance components...\n")
flat <- c(producer = 0.001, associate = 0.001, interaction = 0,
          residual = 0.998)
null_sim <- simulate_experiment(
  G_oat, G_pea, sparsity = 0.48, n_factors = 0, interaction_pct = 0, n_envs = 1,
  var_shares = list(oat = flat, pea = flat))
null_sp  <- prepare_scenario(null_sim)
null_fit <- fit_dge_ige(null_sp$train, G_oat, G_pea, with_interaction = FALSE,
                        nIter = N_ITER, burnIn = BURN, seed = 2L)
gen <- function(v) dplyr::filter(v, term %in% c("G_oat", "G_pea"))
nv <- gen(null_fit$varcomp); rv <- gen(vc)
check(all(nv$var_Pr < rv$var_Pr) && all(nv$var_As < rv$var_As),
      sprintf("a near-null scenario fits less genetic variance than a real one (%.3f vs %.3f on var_Pr)",
              max(nv$var_Pr), max(rv$var_Pr)))
check(max(nv$var_Pr, nv$var_As) < 0.15,
      sprintf("and little of it in absolute terms (max %.3f on a unit-variance scale)",
              max(nv$var_Pr, nv$var_As)))

# ------------------------------------------------------------
# 5. MegaLMM in both orientations, and the GMA assembly
#
# The assembly is the part of Part 4 with no algebraic check available: each
# species' GMA comes from a margin of one orientation plus a margin of the
# other, so it can only be verified against the truth by fitting both.
# ------------------------------------------------------------

if (!requireNamespace("MegaLMM", quietly = TRUE)) {
  cat("  SKIP: MegaLMM is not installed\n")
} else {
  cat("fitting MegaLMM in both orientations...\n")
  run_dir <- file.path(tempdir(), "test_mm")
  mm <- fit_megalmm_both(sp$train, G_oat, G_pea, runID = run_dir,
                         K = 5L, eigen_variance = 0.75,
                         fixed_main_effect = FALSE)

  check(!is.null(mm$surface$oat) && !is.null(mm$surface$pea),
        "both orientations return a surface")
  check(identical(rownames(mm$surface$oat), rownames(G_oat)) &&
          identical(colnames(mm$surface$oat), rownames(G_pea)),
        "and both are returned oat-rows x pea-cols, as the scorer assumes")
  check(identical(dim(mm$surface$pea), dim(mm$surface$oat)),
        "so the two surfaces are conformable")

  s_mm <- score_predictions(mm, sim, sp$idx, "megalmm", train = sp$train)
  for (col in c("r_oat_prod", "r_oat_assoc", "r_pea_prod", "r_pea_assoc")) {
    check(s_mm[[col]] > 0.2,
          sprintf("MegaLMM recovers %s (got %.2f)", col, s_mm[[col]]))
  }

  # The GMA is producer + associate, and each comes from a different
  # orientation. If the assembly is right it must beat the single-orientation
  # margin it is built from; if the orientations were mixed up it will not.
  # Both surfaces are oat-rows x pea-cols, so BOTH oat effects are row margins:
  # the producer off the oat-yield surface, the associate off the pea-yield one.
  # colMeans(surface$pea) is the PEA's producer effect, not the oat's associate.
  oat_gma_true <- sim$truth$oat_prod + sim$truth$oat_assoc
  from_oat_side <- rowMeans(mm$surface$oat)
  assembled     <- rowMeans(mm$surface$oat) + rowMeans(mm$surface$pea)
  check(stats::cor(assembled, oat_gma_true) >
          stats::cor(from_oat_side, oat_gma_true),
        sprintf("assembling across orientations beats one orientation (%.2f vs %.2f)",
                stats::cor(assembled, oat_gma_true),
                stats::cor(from_oat_side, oat_gma_true)))
  check_near(s_mm$r_oat_gma, stats::cor(assembled, oat_gma_true), tol = 1e-6,
             "and r_oat_gma is exactly that assembled correlation")
}

finish("fit tests")
