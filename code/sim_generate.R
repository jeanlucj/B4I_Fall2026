# ============================================================
# SIMULATING AN OAT x PEA INTERCROP EXPERIMENT
#
# Generates oat yield for a set of observed oat x pea combinations under a
# known truth, so that the MegaLMM and DGE-IGE frameworks can be scored
# against something rather than against each other.
#
# The model, for oat accession i grown with pea accession j in environment k:
#
#   y_ijk = mu_k + s_k * ( Pr_i + As_j + I_ij + e_ijk )
#
#   Pr ~ N(0, G_oat * V_Pr)         producer: what the oat does for itself
#   As ~ N(0, G_pea * V_As)         associate: what the pea does to the oat
#   I_ij = sum_f u_if * lam_jf      interaction, rank n_factors
#     with u_.f ~ N(0, G_oat), lam_.f ~ N(0, G_pea)
#   e ~ N(0, V_e)
#
# The interaction is the point of the exercise. Drawing the oat scores from
# G_oat and the pea loadings from G_pea makes a single factor's covariance
# exactly G_oat[i,i'] * G_pea[j,j'] -- so as n_factors grows the interaction
# converges on the full Kronecker structure the DGE-IGE model assumes, while
# at n_factors = 1 it is as low-rank as MegaLMM could hope for. The design
# therefore slides smoothly between the two frameworks' home ground instead
# of favouring one by construction.
#
# Environments carry two separable things.
#
#   * Mean and variance heterogeneity: mu_k and s_k. s_k scales the whole
#     signal, so on its own it changes the spread and nothing else -- no
#     accession changes rank.
#   * Genotype x environment, controlled by `gxe_cor`. Each producer and
#     associate effect splits into a stable share and an environment-specific
#     one, so an accession's effect correlates `gxe_cor` with itself in
#     another environment. At gxe_cor = 1 there is none, which is what this
#     generator did before the argument existed.
#
# Keeping them apart matters: the first is a nuisance the analysis removes by
# standardising within environment, the second is the attenuation that decides
# whether an effect estimated here means anything there. The leave-one-trial-
# out work in code/validate_crossval.R measures the second on the real data.
#
# grm_factor() and read_grm() come from code/dge_ige_functions.R, which
# code/sim_run.R sources: the simulation draws effects with the same
# decomposition the DGE-IGE model fits with.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

#' Draw a vector with covariance G, scaled to variance `target_var`.
#'
#' Scaling is empirical -- the realised vector is rescaled to the target --
#' so a single replicate has the variance the design asked for rather than
#' only having it in expectation.
draw_effect <- function(L, target_var) {
  v <- as.vector(L %*% stats::rnorm(ncol(L)))
  if (target_var <= 0) return(rep(0, length(v)))
  v <- v - mean(v)
  v * sqrt(target_var / stats::var(v))
}

#' Which oat x pea combinations were observed.
#'
#' A simple random sample of cells leaves some accessions with one
#' observation or none, and a pea environment holding one oat has no
#' estimable residual variance: MegaLMM's ARD sampler answers that with NaN
#' rather than an error. Every accession is therefore guaranteed
#' `min_per_acc` observations before the remainder is filled at random. That
#' is also the more realistic design -- nobody grows an entry once -- and it
#' keeps `sparsity` meaning what it says, since the panel is never trimmed
#' afterwards.
sample_combinations <- function(n_oat, n_pea, sparsity,
                                min_per_acc = SIM_MIN_PER_ACC) {
  n_cells <- n_oat * n_pea
  n_obs <- round(sparsity * n_cells)
  needed <- min_per_acc * max(n_oat, n_pea)
  if (n_obs < needed) {
    stop("sparsity ", sparsity, " gives ", n_obs, " observations, fewer than ",
         "the ", needed, " needed to see every accession ", min_per_acc,
         " times", call. = FALSE)
  }

  cell_of <- function(oat, pea) (pea - 1L) * n_oat + oat

  # Rounds of a random oat-to-pea matching: each round gives every oat and
  # every pea one observation
  cells <- integer(0)
  for (r in seq_len(min_per_acc)) {
    pea <- if (n_pea >= n_oat) sample(n_pea, n_oat) else sample(rep_len(seq_len(n_pea), n_oat))
    cells <- union(cells, cell_of(seq_len(n_oat), pea))
    if (n_pea > n_oat) {
      oat <- sample(rep_len(seq_len(n_oat), n_pea))
      cells <- union(cells, cell_of(oat, seq_len(n_pea)))
    }
  }

  # Rounds can collide, so top up anyone still short
  top_up <- function(cells) {
    oat_of <- ((cells - 1L) %% n_oat) + 1L
    pea_of <- ((cells - 1L) %/% n_oat) + 1L
    oat_n <- tabulate(oat_of, n_oat)
    pea_n <- tabulate(pea_of, n_pea)
    short_oat <- which(oat_n < min_per_acc)
    short_pea <- which(pea_n < min_per_acc)
    for (i in short_oat) {
      add <- setdiff(cell_of(i, sample(n_pea)), cells)
      cells <- c(cells, utils::head(add, min_per_acc - oat_n[i]))
    }
    for (j in short_pea) {
      add <- setdiff(cell_of(sample(n_oat), j), cells)
      cells <- c(cells, utils::head(add, min_per_acc - pea_n[j]))
    }
    cells
  }
  cells <- top_up(cells)

  remaining <- setdiff(seq_len(n_cells), cells)
  extra <- sample(remaining, max(0, n_obs - length(cells)))
  cells <- sort(c(cells, extra))

  # The guaranteed minimum can overrun the request when sparsity sits on the
  # floor: `extra` is then empty and there is nothing to trim against, so the
  # achieved sparsity quietly exceeds the one in the scenario's name. Every
  # level in SIM_LEVELS clears the floor with slack, so this should never fire
  # -- which is exactly why it is worth asserting rather than recording.
  if (length(cells) != n_obs) {
    stop("sparsity ", sparsity, " asked for ", n_obs, " observations but the ",
         "minimum of ", min_per_acc, " per accession forces ", length(cells),
         ". The scenario would be mislabelled; raise sparsity or lower ",
         "min_per_acc.", call. = FALSE)
  }

  tibble::tibble(
    cell = cells,
    oat  = ((cells - 1L) %% n_oat) + 1L,
    pea  = ((cells - 1L) %/% n_oat) + 1L
  )
}

#' Simulate one experiment.
#'
#' @param G_oat,G_pea Relationship matrices, already subset to the panel size.
#' @param sparsity Share of the n_oat x n_pea cells observed.
#' @param n_factors Rank of the interaction; 0 for none.
#' @param interaction_pct Interaction variance as a share of the total.
#' @param n_envs Physical environments; combinations are split evenly.
#' @param gxe_cor Genetic correlation of an accession's producer and associate
#'   effects between two environments. 1 makes every effect perfectly stable,
#'   which is what this function did before the argument existed; 0.6 puts 40%
#'   of each effect's variance into environment-specific deviations. Ignored
#'   when `n_envs == 1`, where it is unidentifiable.
#' @param var_shares Producer / associate / residual shares of the remainder.
#' @param env_log_sd SD of log environment scale factor.
#' @param env_mean_log_sd SD of log environment mean.
#' @return list(obs, truth, G_oat, G_pea, settings)
simulate_experiment <- function(G_oat, G_pea, sparsity, n_factors,
                                interaction_pct, n_envs, gxe_cor = 1,
                                var_shares = SIM_VAR_SHARES,
                                env_log_sd = SIM_ENV_LOG_SD,
                                env_mean_log_sd = SIM_ENV_MEAN_LOG_SD,
                                grand_mean = SIM_GRAND_MEAN) {
  if (gxe_cor < 0 || gxe_cor > 1) {
    stop("gxe_cor must be a correlation in [0, 1], not ", gxe_cor, call. = FALSE)
  }
  # With one environment there is nothing for an effect to be specific TO
  rho <- if (n_envs > 1) gxe_cor else 1
  n_oat <- nrow(G_oat)
  n_pea <- nrow(G_pea)

  L_oat <- grm_factor(G_oat)
  L_pea <- grm_factor(G_pea)

  # Variance budget: the interaction takes its share, the rest is split in
  # the proportions the real data shows
  shares <- var_shares / sum(var_shares)
  V_I  <- interaction_pct
  rest <- 1 - V_I
  V_Pr <- rest * shares[["producer"]]
  V_As <- rest * shares[["associate"]]
  V_e  <- rest * shares[["residual"]]

  # Each effect splits into a stable share rho and an environment-specific
  # share 1 - rho, so the TOTAL producer and associate variances are V_Pr and
  # V_As whatever rho is, and the budget still sums to 1. At rho = 1 these two
  # calls are exactly the old ones.
  Pr <- draw_effect(L_oat, V_Pr * rho)
  As <- draw_effect(L_pea, V_As * rho)

  # Interaction: n_factors axes of equal importance
  if (n_factors > 0 && V_I > 0) {
    U   <- vapply(seq_len(n_factors), \(f) draw_effect(L_oat, 1), numeric(n_oat))
    Lam <- vapply(seq_len(n_factors), \(f) draw_effect(L_pea, 1), numeric(n_pea))
    I_mat <- (U %*% t(Lam)) / sqrt(n_factors)
    I_mat <- I_mat * sqrt(V_I / stats::var(as.vector(I_mat)))
  } else {
    U <- Lam <- NULL
    I_mat <- matrix(0, n_oat, n_pea)
  }

  cells <- sample_combinations(n_oat, n_pea, sparsity)

  # Environments: each takes an equal slice of the combinations at random
  env <- if (n_envs > 1) {
    sample(rep_len(seq_len(n_envs), nrow(cells)))
  } else {
    rep(1L, nrow(cells))
  }
  env_scale <- if (n_envs > 1) exp(stats::rnorm(n_envs, 0, env_log_sd)) else 1
  env_mean  <- if (n_envs > 1) {
    grand_mean * exp(stats::rnorm(n_envs, 0, env_mean_log_sd))
  } else grand_mean

  genetic <- Pr[cells$oat] + As[cells$pea] +
    I_mat[cbind(cells$oat, cells$pea)]
  resid <- stats::rnorm(nrow(cells), 0, sqrt(V_e))

  # --- genotype x environment ---
  #
  # Drawn LAST, and only when it is asked for, so that rho = 1 consumes no
  # random numbers and reproduces the pre-GxE generator bit for bit. Each
  # environment gets its own draw with variance V * (1 - rho), so an
  # accession's effect correlates rho with itself in another environment while
  # its total variance stays V.
  Pr_env <- As_env <- NULL
  gxe <- rep(0, nrow(cells))
  if (rho < 1) {
    Pr_env <- vapply(seq_len(n_envs), \(k) draw_effect(L_oat, V_Pr * (1 - rho)),
                     numeric(n_oat))
    As_env <- vapply(seq_len(n_envs), \(k) draw_effect(L_pea, V_As * (1 - rho)),
                     numeric(n_pea))
    gxe <- Pr_env[cbind(cells$oat, env)] + As_env[cbind(cells$pea, env)]
    genetic <- genetic + gxe
  }

  # Realised, rather than requested: the mean correlation between an
  # accession's effect in one environment and in another. Recorded so the axis
  # can be checked instead of trusted (level S4).
  realised_gxe_cor <- if (rho < 1 && n_envs > 1) {
    per_env <- Pr + Pr_env            # n_oat x n_envs, oat effect by environment
    mean(stats::cor(per_env)[upper.tri(stats::cor(per_env))])
  } else if (n_envs > 1) 1 else NA_real_

  obs <- cells |>
    dplyr::mutate(
      env       = env,
      oat_name  = rownames(G_oat)[oat],
      pea_name  = rownames(G_pea)[pea],
      producer  = Pr[oat],
      associate = As[pea],
      interaction = I_mat[cbind(oat, pea)],
      gxe       = gxe,
      genetic   = genetic,
      # The scale factor multiplies signal and noise alike, so it is variance
      # heterogeneity only. Any genotype x environment interaction comes from
      # the gxe term above, not from this.
      y = env_mean[env] + env_scale[env] * (genetic + resid)
    )

  list(
    obs = obs,
    # producer and associate are the STABLE parts -- what another environment
    # could predict, and therefore what the models are scored against. The
    # environment-specific parts are by construction unpredictable.
    truth = list(producer = Pr, associate = As, interaction = I_mat,
                 producer_env = Pr_env, associate_env = As_env,
                 U = U, Lambda = Lam,
                 env_scale = env_scale, env_mean = env_mean,
                 gxe_cor = rho, realised_gxe_cor = realised_gxe_cor,
                 V = c(producer = V_Pr * rho, associate = V_As * rho,
                       producer_env = V_Pr * (1 - rho),
                       associate_env = V_As * (1 - rho),
                       interaction = V_I, residual = V_e)),
    G_oat = G_oat, G_pea = G_pea,
    settings = list(n_oat = n_oat, n_pea = n_pea, sparsity = sparsity,
                    n_factors = n_factors, interaction_pct = interaction_pct,
                    n_envs = n_envs, gxe_cor = rho)
  )
}

#' Real GRMs, subset to a panel size, as the backbone of the simulation.
#'
#' Using the actual B4I relationship matrices rather than simulated ones
#' keeps the relatedness structure -- including its unevenness -- realistic.
#' Accessions are taken at random so a panel is not systematically the most
#' or least related part of the collection.
sim_grms <- function(n_acc,
                     oat_file = here::here("data", "GRM_Avena.rds"),
                     pea_file = here::here("data", "GRM_Pisum.rds"),
                     seed = 1L) {
  G_oat <- read_grm(oat_file)
  G_pea <- read_grm(pea_file)

  if (nrow(G_oat) < n_acc || nrow(G_pea) < n_acc) {
    stop("panel of ", n_acc, " requested but the GRMs hold ", nrow(G_oat),
         " oat and ", nrow(G_pea), " pea accessions", call. = FALSE)
  }

  withr::with_seed(seed, {
    oat_keep <- sort(sample(rownames(G_oat), n_acc))
    pea_keep <- sort(sample(rownames(G_pea), n_acc))
  })

  list(G_oat = G_oat[oat_keep, oat_keep, drop = FALSE],
       G_pea = G_pea[pea_keep, pea_keep, drop = FALSE])
}
