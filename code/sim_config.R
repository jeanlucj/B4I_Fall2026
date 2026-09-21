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

  # Share of the n_oat x n_pea combinations actually observed. The B4I
  # experiment sits at 1.1% before trimming and 3.2% after, so 3% is roughly
  # where we are and 6/12% are the densities a redesign could buy.
  # 5 / 15 / 45%. The first pass used 3 / 6 / 12%, which is where the B4I
  # experiment sits, but a pilot showed all three are deep in the region
  # where the factor model has nothing to work with: at 20% observed MegaLMM
  # loses badly to an additive model and at 50% it wins. Three levels that
  # straddle that crossover say more than three that agree it has not
  # happened yet. 5% still brackets the real experiment.
  sparsity = c(0.05, 0.15, 0.45),

  # Rank of the true oat x pea interaction. 0 = no interaction at all;
  # 1 = a single axis; 5 = five equally important axes. This is the axis the
  # whole comparison turns on (see BACKGROUND in SIMULATION.md).
  n_factors = c(0, 1, 5),

  # Interaction variance as a share of the total. Meaningless when
  # n_factors = 0, and those cells are dropped from the grid below.
  interaction_pct = c(0.10, 0.20),

  # One physical environment, or ten each holding a tenth of the combinations
  n_envs = c(1, 10)
)

#' The design grid, with the redundant no-interaction cells collapsed.
#'
#' With `n_factors = 0` there is no interaction, so `interaction_pct` has
#' nothing to scale and the two levels would be the same simulation run twice.
#' Those cells are kept once, at 0.
sim_grid <- function(levels = SIM_LEVELS, n_reps = 1L) {
  grid <- tidyr::expand_grid(!!!levels)

  grid <- grid |>
    dplyr::mutate(
      interaction_pct = dplyr::if_else(n_factors == 0, 0, interaction_pct)
    ) |>
    dplyr::distinct()

  tidyr::expand_grid(grid, rep = seq_len(n_reps)) |>
    dplyr::mutate(
      scenario = sprintf("n%d_sp%02d_f%d_i%02d_e%02d",
                         n_acc, round(sparsity * 100), n_factors,
                         round(interaction_pct * 100), n_envs),
      seed = SIM_BASE_SEED + dplyr::row_number(),
      .before = 1
    )
}

SIM_BASE_SEED <- 20260921L

# The same design with the sparsity axis resolved more finely, including the
# 3% the real experiment sits at, for locating the crossover precisely.
#   Rscript code/sim_run.R --extended
SIM_LEVELS_EXTENDED <- modifyList(
  SIM_LEVELS, list(sparsity = c(0.03, 0.05, 0.10, 0.15, 0.25, 0.45, 0.60)))

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
SIM_MEGALMM_K          <- 10
SIM_MEGALMM_BURN_ROUND <- 5
SIM_MEGALMM_BURN_ITER  <- 60
SIM_MEGALMM_SAMPLE     <- 250
SIM_MEGALMM_FIT_X_FROM <- 4

# Share of pea genetic variance covered by the eigenvectors offered to
# MegaLMM as environmental covariates
SIM_EIGEN_VARIANCE <- 0.80

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
