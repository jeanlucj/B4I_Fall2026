# ============================================================
# SIMULATING AN OAT x PEA INTERCROP EXPERIMENT
#
# Generates BOTH yields for a set of observed oat x pea combinations under a
# known truth, so that the MegaLMM and DGE-IGE frameworks can be scored against
# something rather than against each other.
#
# The model, for oat i with pea j in environment k:
#
#   y_oat = mu_oat_k + s_oat_k ( oatPr_i + peaAs_j + I_oat_ij + gxe + e_oat )
#   y_pea = mu_pea_k + s_pea_k ( peaPr_j + oatAs_i + I_pea_ij + gxe + e_pea )
#
#   (oatPr, oatAs) ~ N(0, Sigma_oat (x) G_oat)
#   (peaPr, peaAs) ~ N(0, Sigma_pea (x) G_pea)
#   (e_oat, e_pea) ~ N(0, R)
#   I_oat, I_pea     two INDEPENDENT surfaces, each of rank n_factors
#
# Four genetic effects, not two. A species' producer effect (on its own yield)
# and its associate effect (on its partner's) are correlated properties of the
# same genotype, and that covariance -- Sigma's off-diagonal -- is the whole
# reason the real analysis fits the two yields jointly. An earlier version of
# this file simulated oat yield alone, which left the oat's associate effect and
# the pea's producer effect out of existence and forced the DGE-IGE comparator
# to be univariate.
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

#' Draw a CORRELATED PAIR of effects, both with covariance G, with exactly the
#' requested 2 x 2 covariance between them.
#'
#' This is what makes the simulation bivariate. A species' producer effect (on
#' its own yield) and its associate effect (on its partner's) are two
#' correlated properties of the same genotype -- `Sigma (x) G` in the
#' producer-associate model -- and that covariance is the reason the real
#' analysis fits the two yields jointly rather than separately.
#'
#' The realised pair is whitened and then re-coloured, so its empirical
#' covariance IS `Sigma` rather than merely having it in expectation. Same
#' reasoning as draw_effect()'s empirical rescaling, extended to two columns.
#'
#' @param L Square-root factor of the relationship matrix, from grm_factor().
#' @param Sigma 2 x 2 target covariance, (producer, associate).
#' @return An n x 2 matrix, columns (producer, associate).
draw_effect_pair <- function(L, Sigma) {
  n <- nrow(L)
  if (all(diag(Sigma) <= 0)) return(matrix(0, n, 2))

  Z <- L %*% matrix(stats::rnorm(ncol(L) * 2), ncol(L), 2)
  Z <- sweep(Z, 2, colMeans(Z), "-")

  # whiten: empirical covariance exactly I, so the re-colouring below is exact
  S <- stats::cov(Z)
  Zw <- Z %*% solve(chol(S))

  # a zero-variance column would make chol(Sigma) fail; handle it by scaling
  # the columns separately and only then introducing the correlation
  if (any(diag(Sigma) <= 0)) {
    out <- sweep(Zw, 2, sqrt(pmax(diag(Sigma), 0)), "*")
    return(out)
  }
  Zw %*% chol(Sigma)
}

#' Build a 2 x 2 covariance from two variances and a correlation.
sigma_from_cor <- function(var1, var2, cor12) {
  cv <- cor12 * sqrt(var1 * var2)
  matrix(c(var1, cv, cv, var2), 2, 2)
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

#' Simulate one bivariate experiment.
#'
#' Both yields, and therefore all four genetic effects. For oat i with pea j in
#' environment k:
#'
#'   y_oat = mu_oat_k + s_oat_k ( oatPr_i + peaAs_j + I_oat_ij
#'                                + oatPr_ik + peaAs_jk + e_oat )
#'   y_pea = mu_pea_k + s_pea_k ( peaPr_j + oatAs_i + I_pea_ij
#'                                + peaPr_jk + oatAs_ik + e_pea )
#'
#' with (oatPr, oatAs) ~ N(0, Sigma_oat (x) G_oat) and likewise for pea, so a
#' species' effect on its own yield and its effect on its partner's are
#' correlated -- the covariance the producer-associate model exists to
#' estimate, and the reason the fitted comparator must be Multitrait.
#'
#' Earlier versions simulated oat yield alone, which meant the oat's associate
#' effect and the pea's producer effect did not exist and the DGE-IGE
#' comparator had to be univariate.
#'
#' @param G_oat,G_pea Relationship matrices, already subset to the panel size.
#' @param sparsity Share of the n_oat x n_pea cells observed.
#' @param n_factors Rank of each interaction surface; 0 for none.
#' @param interaction_pct Interaction variance as a share of each trait's total.
#' @param n_envs Physical environments; combinations are split evenly.
#' @param gxe_cor Genetic correlation of an effect between two environments.
#'   1 makes every effect perfectly stable. Ignored when `n_envs == 1`, where
#'   it is unidentifiable.
#' @param var_shares Per-trait producer / associate / residual shares. Read
#'   "producer" as "of the species whose yield this is".
#' @param pr_as_cor Within-species producer-associate correlation, per species.
#' @param resid_cor Residual correlation between the two yields on a plot.
#' @return list(obs, truth, G_oat, G_pea, settings)
simulate_experiment <- function(G_oat, G_pea, sparsity, n_factors,
                                interaction_pct, n_envs, gxe_cor = 1,
                                var_shares = SIM_VAR_SHARES,
                                pr_as_cor = SIM_PR_AS_COR,
                                resid_cor = SIM_RESID_COR,
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

  # ---- variance budget, one per trait ----
  #
  # Each trait's interaction takes its share off the top and the rest is split
  # in the proportions the real data shows. A species' producer variance comes
  # from its OWN trait's budget and its associate variance from the other
  # trait's, because that is where each effect acts.
  V_I <- interaction_pct
  rest <- 1 - V_I
  sh <- lapply(var_shares, \(x) x / sum(x))

  V <- c(
    oat_prod  = rest * sh$oat[["producer"]],   # oat on oat yield
    pea_assoc = rest * sh$oat[["associate"]],  # pea on oat yield
    e_oat     = rest * sh$oat[["residual"]],
    pea_prod  = rest * sh$pea[["producer"]],   # pea on pea yield
    oat_assoc = rest * sh$pea[["associate"]],  # oat on pea yield
    e_pea     = rest * sh$pea[["residual"]]
  )

  # ---- the four genetic effects, as two correlated pairs ----
  oat_pair <- draw_effect_pair(
    L_oat, sigma_from_cor(V[["oat_prod"]] * rho, V[["oat_assoc"]] * rho,
                          pr_as_cor[["oat"]]))
  pea_pair <- draw_effect_pair(
    L_pea, sigma_from_cor(V[["pea_prod"]] * rho, V[["pea_assoc"]] * rho,
                          pr_as_cor[["pea"]]))

  oat_prod <- oat_pair[, 1]; oat_assoc <- oat_pair[, 2]
  pea_prod <- pea_pair[, 1]; pea_assoc <- pea_pair[, 2]

  # ---- two interaction surfaces, independent ----
  #
  # Independent by choice: nothing says a pairing that suits the oat also suits
  # the pea. Fitting the real model with the specific-combination term puts the
  # correlation at about +0.20, but with 1,869 of 2,059 combinations in a single
  # plot that is not separable from plot quality -- see SIMULATION.md.
  draw_interaction <- function() {
    if (!(n_factors > 0 && V_I > 0)) {
      return(list(I = matrix(0, n_oat, n_pea), U = NULL, Lambda = NULL))
    }
    U   <- vapply(seq_len(n_factors), \(f) draw_effect(L_oat, 1), numeric(n_oat))
    Lam <- vapply(seq_len(n_factors), \(f) draw_effect(L_pea, 1), numeric(n_pea))
    I <- (U %*% t(Lam)) / sqrt(n_factors)
    I <- I * sqrt(V_I / stats::var(as.vector(I)))
    list(I = I, U = U, Lambda = Lam)
  }
  int_oat <- draw_interaction()
  int_pea <- draw_interaction()

  # ---- which cells, and in which environment ----
  cells <- sample_combinations(n_oat, n_pea, sparsity)

  env <- if (n_envs > 1) {
    sample(rep_len(seq_len(n_envs), nrow(cells)))
  } else {
    rep(1L, nrow(cells))
  }
  # each trait gets its own environment means and scales: a site that is good
  # for oats is not necessarily good for peas
  env_scale_oat <- if (n_envs > 1) exp(stats::rnorm(n_envs, 0, env_log_sd)) else 1
  env_scale_pea <- if (n_envs > 1) exp(stats::rnorm(n_envs, 0, env_log_sd)) else 1
  env_mean_oat  <- if (n_envs > 1) {
    grand_mean * exp(stats::rnorm(n_envs, 0, env_mean_log_sd))
  } else grand_mean
  env_mean_pea  <- if (n_envs > 1) {
    grand_mean * exp(stats::rnorm(n_envs, 0, env_mean_log_sd))
  } else grand_mean

  # ---- genotype x environment on all four effects ----
  gxe_env <- list(oat_prod = NULL, oat_assoc = NULL,
                  pea_prod = NULL, pea_assoc = NULL)
  gxe_oat <- gxe_pea <- rep(0, nrow(cells))
  if (rho < 1) {
    per_env <- function(L, target) {
      vapply(seq_len(n_envs), \(k) draw_effect(L, target), numeric(nrow(L)))
    }
    gxe_env$oat_prod  <- per_env(L_oat, V[["oat_prod"]]  * (1 - rho))
    gxe_env$oat_assoc <- per_env(L_oat, V[["oat_assoc"]] * (1 - rho))
    gxe_env$pea_prod  <- per_env(L_pea, V[["pea_prod"]]  * (1 - rho))
    gxe_env$pea_assoc <- per_env(L_pea, V[["pea_assoc"]] * (1 - rho))

    gxe_oat <- gxe_env$oat_prod[cbind(cells$oat, env)] +
               gxe_env$pea_assoc[cbind(cells$pea, env)]
    gxe_pea <- gxe_env$pea_prod[cbind(cells$pea, env)] +
               gxe_env$oat_assoc[cbind(cells$oat, env)]
  }

  realised_gxe_cor <- if (rho < 1 && n_envs > 1) {
    pe <- oat_prod + gxe_env$oat_prod
    mean(stats::cor(pe)[upper.tri(stats::cor(pe))])
  } else if (n_envs > 1) 1 else NA_real_

  # ---- correlated residuals: the two yields on one plot ----
  R <- sigma_from_cor(V[["e_oat"]], V[["e_pea"]], resid_cor)
  E <- matrix(stats::rnorm(nrow(cells) * 2), nrow(cells), 2) %*% chol(R)

  genetic_oat <- oat_prod[cells$oat] + pea_assoc[cells$pea] +
    int_oat$I[cbind(cells$oat, cells$pea)] + gxe_oat
  genetic_pea <- pea_prod[cells$pea] + oat_assoc[cells$oat] +
    int_pea$I[cbind(cells$oat, cells$pea)] + gxe_pea

  obs <- cells |>
    dplyr::mutate(
      env      = env,
      oat_name = rownames(G_oat)[oat],
      pea_name = rownames(G_pea)[pea],
      # the scale factors multiply signal and noise alike, so they are variance
      # heterogeneity only; any GxE comes from the gxe terms above
      y_oat = env_mean_oat[env] + env_scale_oat[env] * (genetic_oat + E[, 1]),
      y_pea = env_mean_pea[env] + env_scale_pea[env] * (genetic_pea + E[, 2]),
      genetic_oat = genetic_oat, genetic_pea = genetic_pea
    )

  list(
    obs = obs,
    # The four effects here are the STABLE parts -- what another environment
    # could predict, and therefore what the models are scored against. The
    # environment-specific parts are by construction unpredictable.
    truth = list(
      oat_prod = oat_prod, oat_assoc = oat_assoc,
      pea_prod = pea_prod, pea_assoc = pea_assoc,
      I_oat = int_oat$I, I_pea = int_pea$I,
      U_oat = int_oat$U, Lambda_oat = int_oat$Lambda,
      U_pea = int_pea$U, Lambda_pea = int_pea$Lambda,
      gxe = gxe_env,
      env_scale = list(oat = env_scale_oat, pea = env_scale_pea),
      env_mean  = list(oat = env_mean_oat,  pea = env_mean_pea),
      gxe_cor = rho, realised_gxe_cor = realised_gxe_cor,
      # per trait, each summing to 1
      V = list(
        oat = c(producer = V[["oat_prod"]] * rho,
                associate = V[["pea_assoc"]] * rho,
                producer_env = V[["oat_prod"]] * (1 - rho),
                associate_env = V[["pea_assoc"]] * (1 - rho),
                interaction = V_I, residual = V[["e_oat"]]),
        pea = c(producer = V[["pea_prod"]] * rho,
                associate = V[["oat_assoc"]] * rho,
                producer_env = V[["pea_prod"]] * (1 - rho),
                associate_env = V[["oat_assoc"]] * (1 - rho),
                interaction = V_I, residual = V[["e_pea"]])
      ),
      Sigma = list(
        oat = sigma_from_cor(V[["oat_prod"]], V[["oat_assoc"]], pr_as_cor[["oat"]]),
        pea = sigma_from_cor(V[["pea_prod"]], V[["pea_assoc"]], pr_as_cor[["pea"]]),
        R = R
      )
    ),
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
