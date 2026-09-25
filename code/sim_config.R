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

# MegaLMM settings swept separately from the data-generating design, because
# they change nothing about the simulated experiment and so must not cause the
# expensive BGLR fits to be repeated. sim_run.R crosses these with SIM_LEVELS
# and caches the two halves apart.
#
# K: MegaLMM shrinks surplus factors through the ARD prior, so erring high is
# meant to be cheap. Five is about the most one could hope to detect in a real
# intercrop experiment and ten is past the point of usefulness, so the sweep
# brackets the practical range rather than exploring beyond it.
#
# eigen_variance: how much pea genetic variance the covariates offered to
# Lambda should span. Irrelevant covariates are likewise meant to be shrunk
# away, so this tests that claim as much as it tunes anything.
# fixed_main_effect: give MegaLMM a first factor with loadings pinned at 1,
# so the oat main effect has a term of its own instead of having to be
# reconstructed from latent factors. See docs/MegaLMM_sparsity_challenge.md
# for why that is the obvious thing to try, and setup_megalmm_state() for the
# three settings it drags along with it.
#
# Crossing all three doubles the MegaLMM half to 12 settings and the whole
# grid to 720 fits. The BGLR half is cached apart and is not refitted, so the
# added cost is 6 MegaLMM fits per scenario at roughly 10-40 s each.
SIM_MEGALMM_LEVELS <- list(
  K = c(5, 10),
  eigen_variance = c(0.20, 0.50, 0.80),
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

# Variance of oat yield split between producer, associate and residual, from
# the fitted bivariate model: 395.1 / 299.3 / 1181.5.
SIM_VAR_SHARES <- c(producer = 0.211, associate = 0.160, residual = 0.630)

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
  g <- v |> dplyr::filter(component == "genetic")
  vPr <- mean(g$var_Pr[g$term == "G_oat"])
  vAs <- mean(g$var_As[g$term == "G_pea"])
  vE  <- mean(v$var_oat[v$component == "residual"], na.rm = TRUE)
  shares <- c(producer = vPr, associate = vAs, residual = vE)
  shares <- shares / sum(shares)

  list(
    per_trial   = per_trial,
    grand_mean  = mean(pheno$oat_yield),
    all_trials  = spread(per_trial),
    no_failure  = spread(dplyr::filter(per_trial, sd > 0.5 * stats::median(per_trial$sd))),
    var_shares  = shares
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
