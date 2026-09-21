# ============================================================
# FITTING THE TWO FRAMEWORKS TO A SIMULATED EXPERIMENT
#
# Three models, all asked the same question: predict oat yield in oat x pea
# combinations they were not shown.
#
#   additive   DGE-IGE without an interaction term: producer + associate.
#              The model actually used on the B4I data, and the floor the
#              others must clear to have found any interaction at all.
#   dge_ige    additive plus the specific-combination term, with covariance
#              G_oat (x) G_pea. The proposal's Eqn 2 in full.
#   megalmm    the factor model on the oat x pea matrix, with pea GRM
#              eigenvectors as environmental covariates.
#
# Both frameworks see identical data, identical masking and identical
# preprocessing, and are scored on the same held-out cells against the
# simulated truth rather than against each other.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

# ------------------------------------------------------------
# Preprocessing and masking
# ------------------------------------------------------------

#' Centre and scale within environment.
#'
#' Environments differ in mean and in spread by construction, and neither is
#' of interest. Standardising within environment is what the real analysis
#' does and removes both, at the cost of estimating a mean and an SD from
#' however few observations an environment holds -- which is exactly the
#' penalty the n_envs axis is there to expose.
standardize_within_env <- function(obs) {
  obs |>
    dplyr::group_by(env) |>
    dplyr::mutate(
      y_std = if (dplyr::n() > 2 && stats::sd(y) > 0) {
        (y - mean(y)) / stats::sd(y)
      } else {
        y - mean(y)
      }
    ) |>
    dplyr::ungroup()
}

#' Hold out a share of the observations, keeping a floor in every oat and pea.
#'
#' The floor is not tidiness: an oat or pea left with nothing takes MegaLMM's
#' ARD sampler to NaN deep inside sample_Lambda_prec_ARD rather than to a
#' clean error. Same rule as the real analysis in code/megalmm_oat_pea.R.
mask_observations <- function(obs, fraction, floor_obs = SIM_FLOOR_OBS) {
  oat_left <- table(obs$oat)
  pea_left <- table(obs$pea)
  order_try <- sample(nrow(obs))
  target <- round(fraction * nrow(obs))

  held <- logical(nrow(obs))
  n_held <- 0L
  for (k in order_try) {
    if (n_held >= target) break
    i <- as.character(obs$oat[k]); j <- as.character(obs$pea[k])
    if (oat_left[[i]] > floor_obs && pea_left[[j]] > floor_obs) {
      held[k] <- TRUE
      n_held <- n_held + 1L
      oat_left[[i]] <- oat_left[[i]] - 1L
      pea_left[[j]] <- pea_left[[j]] - 1L
    }
  }
  dplyr::mutate(obs, held = held)
}

# ------------------------------------------------------------
# Design matrices
# ------------------------------------------------------------

grm_basis <- function(G, rank = NULL, variance = NULL) {
  e <- eigen(G, symmetric = TRUE)
  vals <- pmax(e$values, 0)
  k <- if (!is.null(rank)) {
    min(rank, sum(vals > 1e-8 * max(vals)))
  } else {
    which(cumsum(vals) / sum(vals) >= variance)[1]
  }
  V <- sweep(e$vectors[, seq_len(k), drop = FALSE], 2, sqrt(vals[seq_len(k)]), "*")
  rownames(V) <- rownames(G)
  V
}

#' Row-wise Kronecker basis for the specific-combination term.
#'
#' The exact term needs the n_obs x n_obs kernel G_oat[i,i'] * G_pea[j,j'],
#' which at 19,200 observations is a 2.9 GB matrix to eigen-decompose. The
#' same space is spanned by products of the two species' own eigenvectors, so
#' the leading `rank` of each are combined instead: column (a-1)*rank + b of
#' the result is A[,a] * B[,b], and a fitted coefficient vector reshapes to
#' a rank x rank matrix with the full interaction surface A %*% Beta %*% t(B).
#'
#' That reshaping is what makes prediction for every cell cheap; building the
#' basis for all n_oat * n_pea cells directly would not be.
kron_basis <- function(A, B, oat_idx, pea_idx) {
  k <- ncol(A)
  A[oat_idx, rep(seq_len(k), each = k), drop = FALSE] *
    B[pea_idx, rep(seq_len(k), times = k), drop = FALSE]
}

# ------------------------------------------------------------
# The models
# ------------------------------------------------------------

#' DGE-IGE, with or without the specific-combination term.
#'
#' @return list(total, interaction) full n_oat x n_pea prediction surfaces.
fit_dge_ige <- function(train, G_oat, G_pea, with_interaction,
                        kron_rank = SIM_KRON_RANK,
                        nIter = SIM_BGLR_NITER, burnIn = SIM_BGLR_BURNIN) {
  n_oat <- nrow(G_oat); n_pea <- nrow(G_pea)

  L_oat <- grm_basis(G_oat, variance = 0.999)
  L_pea <- grm_basis(G_pea, variance = 0.999)

  ETA <- list(
    oat = list(X = L_oat[train$oat, , drop = FALSE], model = "BRR"),
    pea = list(X = L_pea[train$pea, , drop = FALSE], model = "BRR")
  )

  A <- B <- NULL
  if (with_interaction) {
    A <- grm_basis(G_oat, rank = kron_rank)
    B <- grm_basis(G_pea, rank = kron_rank)
    ETA$int <- list(X = kron_basis(A, B, train$oat, train$pea), model = "BRR")
  }

  fit <- BGLR::BGLR(y = train$y_std, ETA = ETA, nIter = nIter, burnIn = burnIn,
                    verbose = FALSE, saveAt = file.path(tempdir(), "sim_"))

  pr <- as.vector(L_oat %*% fit$ETA$oat$b)
  as_ <- as.vector(L_pea %*% fit$ETA$pea$b)
  total <- outer(pr, as_, "+")

  interaction <- matrix(0, n_oat, n_pea)
  if (with_interaction) {
    Beta <- matrix(fit$ETA$int$b, nrow = ncol(A), ncol = ncol(B), byrow = TRUE)
    interaction <- A %*% Beta %*% t(B)
    total <- total + interaction
  }

  dimnames(total) <- dimnames(interaction) <- list(rownames(G_oat), rownames(G_pea))
  list(total = total, interaction = interaction)
}

#' MegaLMM on the oat x pea matrix.
fit_megalmm <- function(train, G_oat, G_pea, runID,
                        K = SIM_MEGALMM_K,
                        eigen_variance = SIM_EIGEN_VARIANCE) {
  n_oat <- nrow(G_oat); n_pea <- nrow(G_pea)

  Y <- matrix(NA_real_, n_oat, n_pea,
              dimnames = list(rownames(G_oat), rownames(G_pea)))
  Y[cbind(train$oat, train$pea)] <- train$y_std

  # Columns and rows with nothing left cannot be sampled; they keep their
  # place in the returned surface but are dropped from the fit
  keep_col <- colSums(!is.na(Y)) > 0
  keep_row <- rowSums(!is.na(Y)) > 0
  Y_fit <- Y[keep_row, keep_col, drop = FALSE]

  eig <- grm_basis(G_pea[keep_col, keep_col, drop = FALSE],
                   variance = eigen_variance)
  X_Env <- cbind(intercept = 1, eig)
  X_Env_groups <- c(1, rep(2, ncol(eig)))

  accNames <- tibble::tibble(germplasmName = rownames(Y_fit))

  state <- setup_megalmm_state(
    accNames = accNames, wideData = Y_fit,
    kinMat = G_oat[keep_row, keep_row, drop = FALSE],
    envCov = list(X_Env = X_Env, X_Env_groups = X_Env_groups,
                  newEnv = character(0)),
    runID = runID, K = K
  )
  state <- run_megalmm(state, burn_rounds = SIM_MEGALMM_BURN_ROUND,
                       burn_iter = SIM_MEGALMM_BURN_ITER,
                       sample_iter = SIM_MEGALMM_SAMPLE,
                       fit_X_from = SIM_MEGALMM_FIT_X_FROM)
  post <- megalmm_posterior(state, accessions = accNames$germplasmName)

  total <- matrix(0, n_oat, n_pea,
                  dimnames = list(rownames(G_oat), rownames(G_pea)))
  total[keep_row, keep_col] <- post$U
  list(total = total, interaction = NULL)
}

# ------------------------------------------------------------
# Scoring
# ------------------------------------------------------------

#' Strip oat and pea main effects from a full prediction surface.
#'
#' Neither framework hands back an "interaction" on the same terms: the
#' DGE-IGE model has an explicit term, MegaLMM's factors carry main effects
#' and interaction together. Removing row and column means from the whole
#' surface puts both on the same footing, and does the same to the truth.
#' An additive model residualises to exactly zero, which is the correct
#' answer for it.
interaction_part <- function(M) {
  M - rowMeans(M) - rep(colMeans(M), each = nrow(M)) + mean(M)
}

score_predictions <- function(pred, sim, held, label) {
  idx <- cbind(held$oat, held$pea)

  truth_total <- sim$truth$producer[held$oat] +
    sim$truth$associate[held$pea] +
    sim$truth$interaction[idx]

  truth_int <- interaction_part(sim$truth$interaction)[idx]
  pred_int  <- interaction_part(pred$total)[idx]

  safe_cor <- function(a, b) {
    if (stats::sd(a) < 1e-10 || stats::sd(b) < 1e-10) return(NA_real_)
    stats::cor(a, b)
  }

  tibble::tibble(
    model = label,
    r_total       = safe_cor(pred$total[idx], truth_total),
    r_interaction = safe_cor(pred_int, truth_int),
    r_observed    = safe_cor(pred$total[idx], held$y_std)
  )
}

#' Margin-only baselines, as in the real analysis.
baseline_predictions <- function(train, n_oat, n_pea) {
  grand <- mean(train$y_std)
  oat_mean <- tapply(train$y_std, train$oat, mean)
  pea_mean <- tapply(train$y_std, train$pea, mean)

  om <- rep(grand, n_oat); om[as.integer(names(oat_mean))] <- oat_mean
  pm <- rep(grand, n_pea); pm[as.integer(names(pea_mean))] <- pea_mean

  list(
    oat_mean     = list(total = matrix(om, n_oat, n_pea)),
    oat_plus_pea = list(total = outer(om, pm, "+") - grand)
  )
}

# ------------------------------------------------------------
# One scenario, end to end
# ------------------------------------------------------------

run_scenario <- function(sim, cv_fraction = SIM_CV_FRACTION,
                         run_dir = tempdir(), seed = 1L,
                         models = c("additive", "dge_ige", "megalmm")) {
  set.seed(seed)

  obs <- sim$obs |> standardize_within_env() |> mask_observations(cv_fraction)
  train <- dplyr::filter(obs, !held)
  held  <- dplyr::filter(obs, held)

  n_oat <- nrow(sim$G_oat); n_pea <- nrow(sim$G_pea)

  preds <- list()
  timing <- list()

  if ("additive" %in% models) {
    t0 <- Sys.time()
    preds$additive <- fit_dge_ige(train, sim$G_oat, sim$G_pea, FALSE)
    timing$additive <- as.numeric(difftime(Sys.time(), t0, units = "s"))
  }
  if ("dge_ige" %in% models) {
    t0 <- Sys.time()
    preds$dge_ige <- fit_dge_ige(train, sim$G_oat, sim$G_pea, TRUE)
    timing$dge_ige <- as.numeric(difftime(Sys.time(), t0, units = "s"))
  }
  if ("megalmm" %in% models) {
    t0 <- Sys.time()
    preds$megalmm <- fit_megalmm(train, sim$G_oat, sim$G_pea,
                                 runID = file.path(run_dir, "megalmm"))
    timing$megalmm <- as.numeric(difftime(Sys.time(), t0, units = "s"))
  }

  preds <- c(preds, baseline_predictions(train, n_oat, n_pea))

  scores <- purrr::imap(preds, \(p, nm) score_predictions(p, sim, held, nm)) |>
    purrr::list_rbind()

  scores |>
    dplyr::mutate(
      seconds = unlist(timing)[model] |> unname(),
      n_train = nrow(train), n_held = nrow(held),
      .after = model
    )
}
