# ============================================================
# SIMULATION CONFIGURATION
#
# The factorial design for comparing the MegaLMM factor framework with the
# DGE-IGE (producer-associate) framework, and the parameters that ground it
# in the B4I data.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

# ------------------------------------------------------------
# The factorial
# ------------------------------------------------------------

SIM_LEVELS <- list(
  # Panel size: oat accessions x pea accessions (square panels)
  n_acc = c(200, 400),

  # Share of the n_oat x n_pea combinations actually observed, as 1.6 / 4.8 /
  # 16 / 48% -- two pairs of the same digits a decade apart.
  #
  # The earlier 5 / 15 / 45% was chosen to straddle the crossover where the
  # factor model starts to win, and it does. The low end was added afterwards
  # for a different reason: the fixed-loading main-effect factor
  # (SIM_MEGALMM_LEVELS below) is meant to help most exactly where the matrix
  # is thinnest, because that is where reconstructing a main effect from free
  # factors fails. Testing it only at densities where it is not needed answers
  # the wrong question.
  #
  # 1.6% is the floor worth using. SIM_MIN_PER_ACC = 3 makes 1.50% the lowest
  # feasible density at n_acc = 200 (0.75% at 400), and hugging it makes
  # sample_combinations() overrun the request -- so 1.6% keeps enough slack
  # that every level returns the exact number of observations asked for.
  # B4I itself sits at 1.1% before trimming and 3.2% after.
  sparsity = c(0.016, 0.048, 0.16, 0.48),

  # Rank of the true oat x pea interaction. 0 = no interaction at all;
  # 1 = a single axis; 5 = five equally important axes. This is the axis the
  # whole comparison turns on (see BACKGROUND in SIMULATION.md).
  n_factors = c(0, 1, 5),

  # Interaction variance as a share of the total. Meaningless when
  # n_factors = 0, and those cells are dropped from the grid below.
  interaction_pct = c(0.10, 0.20),

  # One physical environment, or ten each holding a tenth of the combinations
  n_envs = c(1, 10),

  # Genetic correlation of an accession's producer and associate effects
  # between two environments -- ordinary GxE, distinct from the oat x pea
  # interaction. 1 makes every effect perfectly stable, which is what the
  # generator did before this axis existed; 0.6 puts 40% of each effect's
  # variance into environment-specific deviations.
  #
  # It matters because it is the attenuation the validation trial turns on:
  # an effect estimated in one set of environments is only worth having if it
  # survives into another. code/validate_crossval.R measures it on the real
  # trials and found the across-fold spread of the calibration slope to be
  # 0.57-0.69 of its mean, which is roughly where 0.6 sits.
  #
  # Unidentifiable with a single environment, so those cells are collapsed to
  # 1 in the grid below, the same way interaction_pct is at n_factors = 0.
  gxe_cor = c(1.0, 0.6)
)

# MegaLMM settings, swept as part of the design rather than crossed with it.
#
# They change nothing about the simulated experiment, so the expensive DGE-IGE
# half must not be refitted for each of them -- the two halves are cached
# separately. Since MegaLMM is now fitted in BOTH orientations, the sweep is
# also where most of the compute goes, which is why sim_design() takes a
# D-optimal fraction of it rather than running every combination.
#
# K: MegaLMM shrinks surplus factors through the ARD prior, so erring high is
# meant to be cheap. Five is about the most one could hope to detect in a real
# intercrop experiment and ten is past the point of usefulness, so the sweep
# brackets the practical range rather than exploring beyond it.
#
# eigen_variance: how much of the column species' genetic variance the
# covariates offered to Lambda should span. Irrelevant covariates are likewise
# meant to be shrunk away, so this tests that claim as much as it tunes
# anything. Two levels, well apart, is enough to see whether it matters.
#
# fixed_main_effect: give the row species a first factor with loadings pinned
# at 1, so its main effect has a term of its own instead of having to be
# reconstructed from latent factors. See docs/MegaLMM_sparsity_challenge.md for
# why that is the obvious thing to try, and setup_megalmm_state() for the three
# settings it drags along with it.
SIM_MEGALMM_LEVELS <- list(
  K = c(5, 10),
  eigen_variance = c(0.25, 0.75),
  fixed_main_effect = c(FALSE, TRUE)
)

#' The design grid, with the cells that would be duplicates collapsed.
#'
#' Two axes are meaningless at one end of another axis, and their levels would
#' otherwise run the same simulation twice:
#'
#'   * `interaction_pct` has nothing to scale when `n_factors = 0`;
#'   * `gxe_cor` has nothing to vary across when `n_envs = 1`.
#'
#' Both are pinned to a single value there and the duplicates dropped.
#'
#' Sparsity is written into the scenario name in tenths of a percent, so 1.6%
#' and 2% do not collide as they would at whole percents.
sim_grid <- function(levels = SIM_LEVELS, n_reps = 1L) {
  grid <- tidyr::expand_grid(!!!levels)

  grid <- grid |>
    dplyr::mutate(
      interaction_pct = dplyr::if_else(n_factors == 0, 0, interaction_pct),
      gxe_cor         = dplyr::if_else(n_envs == 1, 1, gxe_cor)
    ) |>
    dplyr::distinct()

  tidyr::expand_grid(grid, rep = seq_len(n_reps)) |>
    dplyr::mutate(
      scenario = sprintf("n%d_sp%03d_f%d_i%02d_e%02d_g%03d",
                         n_acc, round(sparsity * 1000), n_factors,
                         round(interaction_pct * 100), n_envs,
                         round(gxe_cor * 100)),
      seed = SIM_BASE_SEED + dplyr::row_number(),
      .before = 1
    )
}

SIM_BASE_SEED <- 20260921L

# ------------------------------------------------------------
# The fractional design
#
# Crossing the data-generating grid with the MegaLMM sweep and running every
# cell is 960 combinations, and with MegaLMM now fitted in BOTH orientations
# that is 1,920 fits per replicate. Most of it is redundant: what we want from
# the sweep is the main effect of each lever and the two-way interactions
# between them. Three-way interactions are not interpretable anyway.
#
# TWO AXES HAVE TO BE RECODED FIRST. `interaction_pct` only exists when
# n_factors > 0 and `gxe_cor` only when n_envs > 1, so as separate factors they
# are NESTED rather than crossed, and a main-effects-plus-two-way model is
# singular on the full grid -- not just on a fraction of it. Folding each
# nested pair into one composite factor makes the design properly crossed:
#
#   interaction : none, f1_i10, f1_i20, f5_i10, f5_i20        (5 levels)
#   environment : one, ten_stable, ten_gxe                    (3 levels)
#
# With n_acc (2), sparsity (4), K (2), eigen_variance (2) and
# fixed_main_effect (2) that is 960 candidates and 82 parameters for main
# effects plus all two-way interactions, and the full factorial has rank 82.
# A D-optimal subset of 150 runs keeps full rank at D-efficiency ~0.6 and
# touches all 120 data scenarios, so the DGE-IGE half loses no coverage.
#
# HOW THE RESULTS ARE THEN READ. Not as a table of cell means -- most cells
# have only some MegaLMM settings. Fit the design model to the outcomes and
# report main effects and two-way interactions. That also removes the
# best-of-twelve-settings selection bias the old summary had, because there is
# no longer a per-scenario maximum to take.
# ------------------------------------------------------------

# How the composite levels map back to what the generator needs
SIM_INTERACTION_MAP <- tibble::tribble(
  ~interaction, ~n_factors, ~interaction_pct,
  "none",       0,          0,
  "f1_i10",     1,          0.10,
  "f1_i20",     1,          0.20,
  "f5_i10",     5,          0.10,
  "f5_i20",     5,          0.20
)

SIM_ENVIRONMENT_MAP <- tibble::tribble(
  ~environment,  ~n_envs, ~gxe_cor,
  "one",         1,       1.0,
  "ten_stable",  10,      1.0,
  "ten_gxe",     10,      0.6
)

#' The full candidate set, in both codings.
#'
#' One row per (data scenario x MegaLMM setting), with the composite factors
#' the design model uses alongside the underlying levels the generator needs.
sim_candidates <- function(levels = SIM_LEVELS,
                           mm_levels = SIM_MEGALMM_LEVELS) {
  tidyr::expand_grid(
    n_acc = levels$n_acc,
    sparsity = levels$sparsity,
    interaction = SIM_INTERACTION_MAP$interaction,
    environment = SIM_ENVIRONMENT_MAP$environment,
    K = mm_levels$K,
    eigen_variance = mm_levels$eigen_variance,
    fixed_main_effect = mm_levels$fixed_main_effect
  ) |>
    dplyr::left_join(SIM_INTERACTION_MAP, by = "interaction") |>
    dplyr::left_join(SIM_ENVIRONMENT_MAP, by = "environment")
}

#' Everything a run needs: the data scenarios, and the runs chosen from them.
#'
#' @param n_runs Size of the D-optimal subset. 150 is the default because it is
#'   where the model first becomes full rank with room to spare AND every data
#'   scenario is still touched.
#' @param full TRUE ignores the fraction and returns every candidate, for when
#'   the compute is available and the cell means are wanted directly.
#' @return list(scenarios, runs, efficiency, model_rank, model_terms)
sim_design <- function(levels = SIM_LEVELS, mm_levels = SIM_MEGALMM_LEVELS,
                       n_runs = SIM_DESIGN_RUNS, n_reps = 1L,
                       seed = SIM_DESIGN_SEED, full = FALSE) {

  cand <- sim_candidates(levels, mm_levels)
  design_vars <- c("n_acc", "sparsity", "interaction", "environment",
                   "K", "eigen_variance", "fixed_main_effect")

  chosen <- if (full) {
    cand
  } else {
    if (!requireNamespace("AlgDesign", quietly = TRUE)) {
      stop("the fractional design needs the AlgDesign package; install it or ",
           "call sim_design(full = TRUE)", call. = FALSE)
    }
    d <- as.data.frame(lapply(cand[design_vars], factor))
    withr::with_seed(seed, {
      opt <- AlgDesign::optFederov(~ .^2, data = d, nTrials = n_runs,
                                   criterion = "D", nRepeats = 20,
                                   maxIteration = 200)
    })
    cand[opt$rows, ]
  }

  X <- stats::model.matrix(~ .^2,
                           data = as.data.frame(lapply(chosen[design_vars], factor)))

  # the data scenarios the chosen runs need generated, numbered and seeded once
  scenarios <- chosen |>
    dplyr::distinct(n_acc, sparsity, interaction, environment,
                    n_factors, interaction_pct, n_envs, gxe_cor) |>
    dplyr::arrange(n_acc, sparsity, interaction, environment) |>
    dplyr::mutate(
      scenario = sprintf("n%d_sp%03d_f%d_i%02d_e%02d_g%03d",
                         n_acc, round(sparsity * 1000), n_factors,
                         round(interaction_pct * 100), n_envs,
                         round(gxe_cor * 100)),
      .before = 1
    )

  scenarios <- tidyr::expand_grid(scenarios, rep = seq_len(n_reps)) |>
    dplyr::mutate(seed = SIM_BASE_SEED + dplyr::row_number())

  runs <- chosen |>
    dplyr::left_join(dplyr::select(scenarios, scenario, rep, seed,
                                   n_acc, sparsity, interaction, environment),
                     by = c("n_acc", "sparsity", "interaction", "environment"),
                     relationship = "many-to-many")

  list(scenarios = scenarios, runs = runs,
       efficiency = if (full) NA_real_ else
         attr(chosen, "Dea") %||% NA_real_,
       model_rank = qr(X)$rank, model_terms = ncol(X))
}

# Runs in the fractional design, and the seed the search used. 150 is not a
# round number chosen for tidiness: below about 120 the model loses rank, and
# above 200 the gain in D-efficiency stops paying for the compute.
SIM_DESIGN_RUNS <- 150L
SIM_DESIGN_SEED <- 20260925L


# The same design with the sparsity axis resolved more finely, for locating
# the crossover precisely. Stays clear of the 1.50% floor at n_acc = 200.
#   Rscript code/sim_run.R --extended
SIM_LEVELS_EXTENDED <- modifyList(
  SIM_LEVELS,
  list(sparsity = c(0.016, 0.032, 0.048, 0.08, 0.16, 0.32, 0.48, 0.60)))

# ------------------------------------------------------------
# Parameters taken from the B4I data
#
# Recompute with sim_observed_parameters() if the data changes; these are the
# values it returns for the six B4I intercrop trials as of 2026-09-21.
# ------------------------------------------------------------

# Each TRAIT's variance split between the producer effect of the species whose
# yield it is, the associate effect of the other species, and residual -- from
# the fitted bivariate model. The two traits are deliberately NOT symmetric,
# because the data is not:
#
#   oat yield : oat producer 395.1 / pea associate 299.3 / residual 1181.5
#   pea yield : pea producer 185.9 / oat associate 105.3 / residual  572.9
#
# Read "producer" as "of the species this trait belongs to" and "associate" as
# "of the other species" throughout.
SIM_VAR_SHARES <- list(
  oat = c(producer = 0.211, associate = 0.160, residual = 0.630),
  pea = c(producer = 0.215, associate = 0.121, residual = 0.663)
)

# Within-species correlation between a species' producer effect (on its own
# yield) and its associate effect (on its partner's). This is the off-diagonal
# of Sigma_oat and Sigma_pea -- the covariance that motivates fitting the two
# traits jointly, and the reason the simulation needs both.
#
# These are the GENETIC correlations from the fit (-0.065, -0.234), not the
# correlations between the BLUPs (-0.405, -0.441). The BLUP correlation is
# inflated because the two effects are estimated with correlated errors: on any
# plot both contribute and they trade off.
SIM_PR_AS_COR <- c(oat = -0.065, pea = -0.234)

# Residual correlation between the two yields on the same plot: competition,
# once trial and block are removed.
SIM_RESID_COR <- -0.122

# Between-environment heterogeneity. Environment scale factors are drawn
# log-normal, so this is the SD of log(within-trial SD) across trials.
#
#   0.669  all six B4I trials
#   0.167  excluding B4I_2025_AL
#
# AL was a near-total crop failure (mean 8.3 g/m2 against 127-383 elsewhere)
# and drives the larger figure on its own. The default is the honest "as
# observed" value; set SIM_ENV_LOG_SD to 0.167 to ask what happens in a season
# where nothing fails.
SIM_ENV_LOG_SD <- 0.669

# SD of log(trial mean), for the environment means. 1.355 as observed.
SIM_ENV_MEAN_LOG_SD <- 1.355

# Grand mean of oat yield, g/m2
SIM_GRAND_MEAN <- 167.3

#' Recompute the parameters above from the real phenotype table.
#'
#' Returns the numbers, and prints the figures with and without the failed
#' trial so the choice of SIM_ENV_LOG_SD stays a visible decision.
sim_observed_parameters <- function(
    pheno_file = here::here("output", "B4I_intercrop_pheno.rds"),
    varcomp_file = here::here("output", "BGLR_variance_components.csv")) {

  pheno <- readRDS(pheno_file) |> dplyr::filter(!is.na(oat_yield))

  per_trial <- pheno |>
    dplyr::group_by(studyName) |>
    dplyr::summarise(n = dplyr::n(), mean = mean(oat_yield),
                     sd = stats::sd(oat_yield), .groups = "drop")

  spread <- function(d) c(env_log_sd = stats::sd(log(d$sd)),
                          env_mean_log_sd = stats::sd(log(d$mean)))

  v <- readr::read_csv(varcomp_file, show_col_types = FALSE)
  g <- dplyr::filter(v, component == "genetic")
  r <- dplyr::filter(v, component == "residual")

  # Which component belongs to which trait. For oat yield the producer is the
  # OAT's effect on its own yield (G_oat's var_Pr) and the associate is the
  # PEA's effect on the oat (G_pea's var_As); for pea yield it is the mirror.
  comp <- function(term, field) mean(g[[field]][g$term == term])
  budget <- list(
    oat = c(producer = comp("G_oat", "var_Pr"),
            associate = comp("G_pea", "var_As"),
            residual  = mean(r$var_oat, na.rm = TRUE)),
    pea = c(producer = comp("G_pea", "var_Pr"),
            associate = comp("G_oat", "var_As"),
            residual  = mean(r$var_pea, na.rm = TRUE))
  )

  list(
    per_trial   = per_trial,
    grand_mean  = mean(pheno$oat_yield),
    all_trials  = spread(per_trial),
    no_failure  = spread(dplyr::filter(per_trial, sd > 0.5 * stats::median(per_trial$sd))),
    # absolute variances, and the same as shares summing to 1 within each trait
    var_absolute = budget,
    var_shares  = purrr::map(budget, \(b) b / sum(b)),
    # the off-diagonals: within-species producer-associate, and residual
    pr_as_cor = c(oat = mean(g$cor_PrAs[g$term == "G_oat"]),
                  pea = mean(g$cor_PrAs[g$term == "G_pea"])),
    resid_cor = mean(r$cor_pea_oat, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# Fitting
# ------------------------------------------------------------

# Fraction of observed cells held out for scoring
SIM_CV_FRACTION <- 0.20

# Minimum observations an oat or pea must retain after masking
SIM_FLOOR_OBS <- 2

# Minimum observations every accession gets by design. Below 3 a pea
# environment can end up with a single oat, which has no estimable residual
# variance and takes MegaLMM's sampler to NaN. Matches the trimming the real
# analysis applies after the fact.
SIM_MIN_PER_ACC <- 3L

# MegaLMM
SIM_MEGALMM_K          <- 10   # default when not swept
SIM_MEGALMM_BURN_ROUND <- 5
SIM_MEGALMM_BURN_ITER  <- 60
SIM_MEGALMM_SAMPLE     <- 250
SIM_MEGALMM_FIT_X_FROM <- 4

# Share of pea genetic variance covered by the eigenvectors offered to
# MegaLMM as environmental covariates
SIM_EIGEN_VARIANCE <- 0.80   # default when not swept

# Split the MegaLMM sampling into this many chunks and score after each, to
# see whether accuracy is still climbing when the chain stops. Answers the
# chain-length question without a sweep: one run reports its own trace.
SIM_TRACE_CHUNKS <- 5L

# BGLR
SIM_BGLR_NITER  <- 6000
SIM_BGLR_BURNIN <- 1000

# The DGE-IGE interaction term is a Kronecker kernel over observed
# combinations, which is n_obs x n_obs. Rather than eigen-decomposing that at
# 19,200 observations, the term is built from the leading eigenvectors of
# G_oat and G_pea: their Kronecker product spans the same space and the
# leading block captures most of it. This is the rank used per species, so
# the interaction basis has SIM_KRON_RANK^2 columns.
#
# Set to NA to build the exact kernel instead, which is only tractable at the
# smaller and sparser cells; sim_fit.R falls back to the approximation with a
# warning above SIM_KRON_EXACT_MAX observations.
SIM_KRON_RANK      <- 30L
SIM_KRON_EXACT_MAX <- 3000L
