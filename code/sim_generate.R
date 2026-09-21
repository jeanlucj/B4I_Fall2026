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
# Environments carry mean and variance heterogeneity only: s_k scales the
# whole signal, so genetic correlations between environments stay 1 and the
# only GxE is in the spread. That is what the B4I trials show once the
# failed site is taken into account, and it keeps the environment axis from
# quietly introducing a second kind of interaction.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

#' Square-root factor of a relationship matrix, dropping null directions.
grm_factor <- function(G, tol = 1e-8) {
  e <- eigen(G, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  L <- sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), "*")
  rownames(L) <- rownames(G)
  L
}

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
#' @param var_shares Producer / associate / residual shares of the remainder.
#' @param env_log_sd SD of log environment scale factor.
#' @param env_mean_log_sd SD of log environment mean.
#' @return list(obs, truth, G_oat, G_pea, settings)
simulate_experiment <- function(G_oat, G_pea, sparsity, n_factors,
                                interaction_pct, n_envs,
                                var_shares = SIM_VAR_SHARES,
                                env_log_sd = SIM_ENV_LOG_SD,
                                env_mean_log_sd = SIM_ENV_MEAN_LOG_SD,
                                grand_mean = SIM_GRAND_MEAN) {
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

  Pr <- draw_effect(L_oat, V_Pr)
  As <- draw_effect(L_pea, V_As)

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

  obs <- cells |>
    dplyr::mutate(
      env       = env,
      oat_name  = rownames(G_oat)[oat],
      pea_name  = rownames(G_pea)[pea],
      producer  = Pr[oat],
      associate = As[pea],
      interaction = I_mat[cbind(oat, pea)],
      genetic   = genetic,
      # The scale factor multiplies signal and noise alike: variance
      # heterogeneity, with genetic correlations across environments still 1
      y = env_mean[env] + env_scale[env] * (genetic + resid)
    )

  list(
    obs = obs,
    truth = list(producer = Pr, associate = As, interaction = I_mat,
                 U = U, Lambda = Lam,
                 env_scale = env_scale, env_mean = env_mean,
                 V = c(producer = V_Pr, associate = V_As,
                       interaction = V_I, residual = V_e)),
    G_oat = G_oat, G_pea = G_pea,
    settings = list(n_oat = n_oat, n_pea = n_pea, sparsity = sparsity,
                    n_factors = n_factors, interaction_pct = interaction_pct,
                    n_envs = n_envs)
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
  read_grm <- function(p) {
    g <- readRDS(p)
    if (is.list(g) && !is.null(g$G)) g$G else g
  }
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
