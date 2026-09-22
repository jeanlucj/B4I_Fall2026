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
#'
#' Returns `Eta_mean`, the predicted phenotype, NOT `U`. `U` is a genetic
#' value and excludes the per-column intercept, which is where MegaLMM keeps
#' the pea main effect -- so scoring `U` against a truth that contains the pea
#' effect asks it to predict a component it structurally cannot hold. `U` is
#' returned alongside for the comparisons where that is the quantity wanted.
#'
#' @param n_chunks Split the sampling into this many pieces and score after
#'   each, giving a convergence trace. 1 disables it.
#' @param score_fn Called as score_fn(Eta) after each chunk.
fit_megalmm <- function(train, G_oat, G_pea, runID,
                        K = SIM_MEGALMM_K,
                        eigen_variance = SIM_EIGEN_VARIANCE,
                        n_chunks = 1L, score_fn = NULL,
                        fixed_main_effect = FALSE) {
  n_oat <- nrow(G_oat); n_pea <- nrow(G_pea)

  Y <- matrix(NA_real_, n_oat, n_pea,
              dimnames = list(rownames(G_oat), rownames(G_pea)))
  Y[cbind(train$oat, train$pea)] <- train$y_std

  keep_col <- colSums(!is.na(Y)) > 0
  keep_row <- rowSums(!is.na(Y)) > 0
  Y_fit <- Y[keep_row, keep_col, drop = FALSE]

  # With a fixed factor MegaLMM's own per-column scaling is switched off, so Y
  # is put on a sensible scale here instead -- once, globally, not per column.
  if (fixed_main_effect) {
    Y_fit <- Y_fit / stats::sd(Y_fit, na.rm = TRUE)
  }

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
    runID = runID, K = K, fixed_main_effect = fixed_main_effect
  )

  place <- function(M) {
    out <- matrix(0, n_oat, n_pea,
                  dimnames = list(rownames(G_oat), rownames(G_pea)))
    out[keep_row, keep_col] <- M
    out
  }

  trace <- NULL
  if (n_chunks > 1L && !is.null(score_fn)) {
    state <- run_megalmm(state, burn_rounds = SIM_MEGALMM_BURN_ROUND,
                         burn_iter = SIM_MEGALMM_BURN_ITER,
                         sample_iter = 0, fit_X_from = SIM_MEGALMM_FIT_X_FROM)
    per <- ceiling(SIM_MEGALMM_SAMPLE / n_chunks)
    trace <- purrr::map(seq_len(n_chunks), \(k) {
      state <<- MegaLMM::sample_MegaLMM(state, per, verbose = FALSE)
      state <<- MegaLMM::save_posterior_chunk(state)
      Eta <- MegaLMM::get_posterior_mean(
        MegaLMM::load_posterior_param(state, "Eta_mean"))
      tibble::tibble(iterations = k * per, score_fn(place(Eta)))
    }) |> purrr::list_rbind()
  } else {
    state <- run_megalmm(state, burn_rounds = SIM_MEGALMM_BURN_ROUND,
                         burn_iter = SIM_MEGALMM_BURN_ITER,
                         sample_iter = SIM_MEGALMM_SAMPLE,
                         fit_X_from = SIM_MEGALMM_FIT_X_FROM)
  }

  post <- megalmm_posterior(state, accessions = accNames$germplasmName)

  # The fixed factor's loadings should be constant across columns. Checked
  # rather than assumed: the feature is lightly exercised upstream, and a
  # silent failure here would look like the unconstrained model.
  fixed_ok <- NA
  if (fixed_main_effect && !is.null(post$Lambda)) {
    fixed_ok <- stats::sd(post$Lambda[1, ]) < 1e-8
    if (!fixed_ok) {
      warning("the first factor's loadings are not constant (sd = ",
              signif(stats::sd(post$Lambda[1, ]), 3),
              "); Lambda_fixed did not hold", call. = FALSE)
    }
  }

  list(total = place(post$Eta_mean), U = place(post$U),
       interaction = NULL, trace = trace,
       Lambda = post$Lambda, U_F = post$U_F, fixed_ok = fixed_ok)
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
# One scenario, in two halves
#
# The BGLR models depend only on the simulated data; MegaLMM additionally
# depends on K and on how many pea eigenvectors it is offered. Sweeping those
# two would refit the expensive BGLR half six times over for no reason, so the
# halves are run and cached separately and joined afterwards.
# ------------------------------------------------------------

#' The train/held split, derived from the seed so both halves see the same one.
split_observations <- function(sim, cv_fraction = SIM_CV_FRACTION, seed = 1L) {
  set.seed(seed)
  obs <- sim$obs |> standardize_within_env() |> mask_observations(cv_fraction)
  list(train = dplyr::filter(obs, !held), held = dplyr::filter(obs, held))
}

run_scenario_bglr <- function(sim, cv_fraction = SIM_CV_FRACTION, seed = 1L) {
  sp <- split_observations(sim, cv_fraction, seed)
  n_oat <- nrow(sim$G_oat); n_pea <- nrow(sim$G_pea)

  t0 <- Sys.time()
  add <- fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, FALSE)
  t_add <- as.numeric(difftime(Sys.time(), t0, units = "s"))

  t0 <- Sys.time()
  dge <- fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, TRUE)
  t_dge <- as.numeric(difftime(Sys.time(), t0, units = "s"))

  preds <- c(list(additive = add, dge_ige = dge),
             baseline_predictions(sp$train, n_oat, n_pea))

  purrr::imap(preds, \(p, nm) score_predictions(p, sim, sp$held, nm)) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      seconds = c(t_add, t_dge, NA_real_, NA_real_),
      n_train = nrow(sp$train), n_held = nrow(sp$held), .after = model
    )
}

run_scenario_megalmm <- function(sim, K, eigen_variance,
                                 cv_fraction = SIM_CV_FRACTION, seed = 1L,
                                 run_dir = tempdir(), n_chunks = 1L,
                                 fixed_main_effect = FALSE) {
  sp <- split_observations(sim, cv_fraction, seed)

  score_fn <- if (n_chunks > 1L) {
    function(Eta) {
      score_predictions(list(total = Eta), sim, sp$held, "megalmm") |>
        dplyr::select(-model)
    }
  } else NULL

  t0 <- Sys.time()
  mm <- fit_megalmm(sp$train, sim$G_oat, sim$G_pea,
                    runID = file.path(run_dir, "megalmm"),
                    K = K, eigen_variance = eigen_variance,
                    n_chunks = n_chunks, score_fn = score_fn,
                    fixed_main_effect = fixed_main_effect)
  t_mm <- as.numeric(difftime(Sys.time(), t0, units = "s"))

  scores <- dplyr::bind_rows(
    score_predictions(mm, sim, sp$held, "megalmm"),
    score_predictions(list(total = mm$U), sim, sp$held, "megalmm_U")
  ) |>
    dplyr::mutate(seconds = c(t_mm, NA_real_),
                  n_train = nrow(sp$train), n_held = nrow(sp$held),
                  K = K, eigen_variance = eigen_variance,
                  fixed_main_effect = fixed_main_effect,
                  .after = model)

  # How well does the model recover the oat main effect? This is the quantity
  # the fixed factor exists to rescue, and the row average is the bar it has
  # to clear, so both are recorded alongside the held-out accuracy.
  row_mean <- tapply(sp$train$y_std, sp$train$oat, mean)
  oats <- as.integer(names(row_mean))
  main_effect <- tibble::tibble(
    r_mainfactor_truePr = if (!is.null(mm$U_F)) {
      stats::cor(mm$U_F[, 1], sim$truth$producer[seq_len(nrow(mm$U_F))])
    } else NA_real_,
    r_rowmean_truePr = stats::cor(row_mean, sim$truth$producer[oats]),
    fixed_ok = mm$fixed_ok
  )

  list(scores = dplyr::bind_cols(scores, main_effect[rep(1, nrow(scores)), ]),
       trace = mm$trace)
}
