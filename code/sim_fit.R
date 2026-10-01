# ============================================================
# FITTING THE TWO FRAMEWORKS TO A SIMULATED EXPERIMENT
#
# Three models, all asked the same question about BOTH yields: predict oat and
# pea yield in oat x pea combinations they were not shown.
#
#   additive   the bivariate producer-associate model without an interaction
#              term. The model actually used on the B4I data, and the floor the
#              others must clear to have found any interaction at all.
#   dge_ige    the same plus the specific-combination term, covariance
#              G_oat (x) G_pea, built at low rank. The proposal's Eqn 2 in full.
#   megalmm    the factor model, fitted TWICE -- once with oats as rows and pea
#              yield's partner as columns, once transposed.
#
# Both DGE-IGE variants are fitted by fit_producer_associate() in
# code/dge_ige_functions.R, the same function the production analysis and the
# leave-one-trial-out cross-validation use. The comparator is therefore the
# model the project actually runs, not a univariate stand-in for it.
#
# WHY MEGALMM IS RUN TWICE. MegaLMM has a per-column intercept and no per-row
# one, so in any one orientation the ROW species' producer effect is carried by
# kinship-shrunk latent structure while the COLUMN species' associate effect
# sits in an unshrunk fixed intercept. One orientation therefore estimates one
# species' producer effect well and the other species' associate effect badly.
# Running both and assembling gives every effect from the orientation that
# treats it best:
#
#   oat GMA = rowMeans(surface_oat) + rowMeans(surface_pea)
#   pea GMA = colMeans(surface_pea) + colMeans(surface_oat)
#
# EVERY MODEL RETURNS THE SAME THING: two full oat x pea surfaces, one per
# trait, both in oat-rows x pea-cols layout. That is what makes the metrics
# comparable -- each effect is a margin of a surface, taken the same way from
# every framework, rather than each model's own idea of what it estimated.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

# ------------------------------------------------------------
# Preprocessing and masking
# ------------------------------------------------------------

#' Centre and scale each trait within environment.
#'
#' Environments differ in mean and in spread by construction, and neither is of
#' interest. Standardising within environment is what the real analysis does and
#' removes both, at the cost of estimating a mean and an SD from however few
#' observations an environment holds -- which is exactly the penalty the n_envs
#' axis is there to expose. Each trait is standardised separately, because a
#' site that is good for oats need not be good for peas.
standardize_within_env <- function(obs) {
  scale_one <- function(y, n) {
    if (n > 2 && stats::sd(y) > 0) (y - mean(y)) / stats::sd(y) else y - mean(y)
  }
  obs |>
    dplyr::group_by(env) |>
    dplyr::mutate(
      y_oat_std = scale_one(y_oat, dplyr::n()),
      y_pea_std = scale_one(y_pea, dplyr::n())
    ) |>
    dplyr::ungroup()
}

#' Which cells to score, as a two-column index matrix.
#'
#' NO OBSERVATIONS ARE HELD OUT. In a simulation the truth is known for every
#' one of the n_oat x n_pea cells, so the natural test set is the cells that were
#' never observed at all -- which is what the model is actually for, and which
#' the design already provides (1 - sparsity) of. Holding out a share of the
#' OBSERVED cells instead had three costs:
#'
#'   * it made the nominal sparsity a lie. Holding out 20% meant the models were
#'     fitted at 0.8x the labelled density -- a scenario labelled 4.8% observed
#'     trained at 3.84% -- with nothing in the output saying so.
#'   * it threw away most of the available evaluation. At 10% observed there are
#'     72 held-out cells against 3,240 never-observed ones, so the per-cell
#'     metrics were estimated on 2% of the matrix while a 45x larger and more
#'     relevant set went unused. Measured: r_int came out 0.017 on the held-out
#'     cells and 0.144 on the never-observed ones, the former being mostly noise.
#'   * it spent 20% of the data to answer a question the simulation does not need
#'     answered. Guarding against overfitting is what a holdout is for, and here
#'     the criterion is a noise-free truth that was never the fitting target.
#'
#' @param set "unobserved" scores the cells with no plot -- the prediction
#'   target, and the default. "all" scores every cell, which additionally
#'   includes the cells the model was fitted to and so is a blend of prediction
#'   and fit; it differs materially only where the matrix is dense.
scoring_index <- function(sim, set = SIM_SCORE_SET) {
  n_oat <- nrow(sim$G_oat); n_pea <- nrow(sim$G_pea)
  observed <- matrix(FALSE, n_oat, n_pea)
  observed[cbind(sim$obs$oat, sim$obs$pea)] <- TRUE
  keep <- switch(set,
    unobserved = !observed,
    all        = matrix(TRUE, n_oat, n_pea),
    stop("SIM_SCORE_SET must be \"unobserved\" or \"all\", not \"", set, "\"",
         call. = FALSE))
  if (sum(keep) < 3) {
    stop("only ", sum(keep), " cell(s) to score; a correlation needs at least 3",
         call. = FALSE)
  }
  which(keep, arr.ind = TRUE, useNames = FALSE)
}

# ------------------------------------------------------------
# The models
#
# grm_basis() and kron_basis(), which build the design matrices, live in
# code/dge_ige_functions.R beside the other GRM machinery -- the DGE-IGE model
# uses them too.
# ------------------------------------------------------------

#' The bivariate producer-associate model, with or without the
#' specific-combination term.
#'
#' A thin wrapper on fit_producer_associate() in code/dge_ige_functions.R: it
#' puts the simulated observations into the plot-table shape that function
#' expects, and turns its named effect vectors back into the two surfaces every
#' model in this file returns.
#'
#' @return list(surface = list(oat, pea), interaction = list(oat, pea) or NULL)
fit_dge_ige <- function(train, G_oat, G_pea, with_interaction,
                        kron_rank = SIM_KRON_RANK,
                        nIter = SIM_BGLR_NITER, burnIn = SIM_BGLR_BURNIN,
                        seed = 1L, save_effects = FALSE, saveAt = NULL,
                        storage_mode = "double") {
  oat_names <- rownames(G_oat); pea_names <- rownames(G_pea)

  dat <- tibble::tibble(
    oatAcc   = oat_names[train$oat],
    peaAcc   = pea_names[train$pea],
    oatYield = train$y_oat_std,
    peaYield = train$y_pea_std,
    trialF       = droplevels(factor(train$env)),
    # one block per environment: the simulation has no within-environment
    # blocking, and a constant block factor would be aliased with the trial
    blockNumberF = droplevels(factor(train$env))
  )

  f <- fit_producer_associate(
    dat, G_oat, G_pea, seed = seed, nIter = nIter, burnIn = burnIn,
    fit_mix_term = with_interaction,
    kron_rank = if (with_interaction) kron_rank else NA,
    saveAt = saveAt %||% file.path(tempdir(), "sim_dge_"),
    save_effects = save_effects, storage_mode = storage_mode
  )

  # Pad back to the full panel: an accession with no training observation gets
  # zero rather than being dropped, so the surfaces are always n_oat x n_pea.
  pad <- function(v, names_in, all_names) {
    out <- stats::setNames(numeric(length(all_names)), all_names)
    out[names_in] <- v
    out
  }
  oat_prod  <- pad(f$effects$oat$PrEff, f$effects$oat$accession, oat_names)
  oat_assoc <- pad(f$effects$oat$AsEff, f$effects$oat$accession, oat_names)
  pea_prod  <- pad(f$effects$pea$PrEff, f$effects$pea$accession, pea_names)
  pea_assoc <- pad(f$effects$pea$AsEff, f$effects$pea$accession, pea_names)

  # oat yield = oat's producer + pea's associate; pea yield is the mirror.
  # Both surfaces are oat-rows x pea-cols.
  surface_oat <- outer(oat_prod, pea_assoc, "+")
  surface_pea <- outer(oat_assoc, pea_prod, "+")

  interaction <- NULL
  if (with_interaction && !is.null(f$interaction)) {
    pad_mat <- function(M) {
      out <- matrix(0, length(oat_names), length(pea_names),
                    dimnames = list(oat_names, pea_names))
      out[rownames(M), colnames(M)] <- M
      out
    }
    interaction <- list(oat = pad_mat(f$interaction$oat),
                        pea = pad_mat(f$interaction$pea))
    surface_oat <- surface_oat + interaction$oat
    surface_pea <- surface_pea + interaction$pea
  }

  # The Kronecker bases, padded the same way. Bilinearity survives padding --
  # pad_mat(A Beta B') == pad_rows(A) Beta pad_rows(B)' -- so the decomposition
  # can work from these and still describe the PADDED surface, which is the one
  # score_predictions() scores. Decomposing the unpadded surface would centre
  # over a different panel and report shares of a different quantity.
  kron <- NULL
  if (with_interaction && !is.null(f$kron)) {
    pad_rows <- function(M, all_names) {
      out <- matrix(0, length(all_names), ncol(M),
                    dimnames = list(all_names, NULL))
      out[rownames(M), ] <- M
      out
    }
    kron <- f$kron
    kron$A <- pad_rows(f$kron$A, oat_names)
    kron$B <- pad_rows(f$kron$B, pea_names)
  }

  list(surface = list(oat = surface_oat, pea = surface_pea),
       interaction = interaction, varcomp = f$varcomp,
       kron = kron, mcmc = f$mcmc,
       beta = if (!is.null(f$kron)) f$fit$ETA$G_mix$beta else NULL)
}

#' MegaLMM on one orientation of the matrix.
#'
#' Orientation-agnostic: `row_idx`, `col_idx` and `y` say which species is on
#' the rows, which is on the columns, and which yield fills the cells. Called
#' twice per scenario by fit_megalmm_both().
#'
#' Returns `Eta_mean`, the predicted phenotype, NOT `U`. `U` is a genetic value
#' and excludes the per-column intercept, which is where MegaLMM keeps the
#' column species' main effect -- so scoring `U` against a truth that contains
#' that effect asks it to predict a component it structurally cannot hold.
#'
#' @param G_row,G_col Relationship matrices for the row and column species.
#'   Only the row species gets kinship inside MegaLMM; the column species is
#'   represented by its GRM eigenvectors as covariates on the Lambda prior.
fit_megalmm_one <- function(row_idx, col_idx, y, G_row, G_col, runID,
                            K = SIM_MEGALMM_K,
                            eigen_variance = SIM_EIGEN_VARIANCE,
                            fixed_main_effect = FALSE) {
  n_row <- nrow(G_row); n_col <- nrow(G_col)

  Y <- matrix(NA_real_, n_row, n_col,
              dimnames = list(rownames(G_row), rownames(G_col)))
  Y[cbind(row_idx, col_idx)] <- y

  keep_col <- colSums(!is.na(Y)) > 0
  keep_row <- rowSums(!is.na(Y)) > 0
  Y_fit <- Y[keep_row, keep_col, drop = FALSE]

  # With a fixed factor MegaLMM's own per-column scaling is switched off, so Y
  # is put on a sensible scale here instead -- once, globally, not per column.
  if (fixed_main_effect) {
    Y_fit <- Y_fit / stats::sd(Y_fit, na.rm = TRUE)
  }

  eig <- grm_basis(G_col[keep_col, keep_col, drop = FALSE],
                   variance = eigen_variance)
  X_Env <- cbind(intercept = 1, eig)
  X_Env_groups <- c(1, rep(2, ncol(eig)))

  accNames <- tibble::tibble(germplasmName = rownames(Y_fit))

  state <- setup_megalmm_state(
    accNames = accNames, wideData = Y_fit,
    kinMat = G_row[keep_row, keep_row, drop = FALSE],
    envCov = list(X_Env = X_Env, X_Env_groups = X_Env_groups,
                  newEnv = character(0)),
    runID = runID, K = K, fixed_main_effect = fixed_main_effect
  )
  state <- run_megalmm(state, burn_rounds = SIM_MEGALMM_BURN_ROUND,
                       burn_iter = SIM_MEGALMM_BURN_ITER,
                       sample_iter = SIM_MEGALMM_SAMPLE,
                       fit_X_from = SIM_MEGALMM_FIT_X_FROM)
  post <- megalmm_posterior(state, accessions = accNames$germplasmName)

  # Pad dropped rows and columns back with zero so the surface is always the
  # full panel. With SIM_MIN_PER_ACC = 3 and nothing held out, nothing should
  # be dropped; `n_dropped` is reported so that assumption is visible.
  place <- function(M) {
    out <- matrix(0, n_row, n_col,
                  dimnames = list(rownames(G_row), rownames(G_col)))
    out[keep_row, keep_col] <- M
    out
  }

  # The pinned factor's loadings should be constant across columns. Checked
  # rather than assumed: a silent failure looks like the unconstrained model.
  fixed_ok <- NA
  if (fixed_main_effect && !is.null(post$Lambda)) {
    fixed_ok <- stats::sd(post$Lambda[1, ]) < 1e-8
    if (!fixed_ok) {
      warning("the first factor's loadings are not constant (sd = ",
              signif(stats::sd(post$Lambda[1, ]), 3),
              "); Lambda_fixed did not hold", call. = FALSE)
    }
  }

  list(surface = place(post$Eta_mean), U = place(post$U),
       Lambda = post$Lambda, U_F = post$U_F, fixed_ok = fixed_ok,
       n_dropped = sum(!keep_row) + sum(!keep_col))
}

#' MegaLMM in both orientations, returning the two surfaces in one layout.
#'
#' The pea-side fit has peas on the rows, so its surface is transposed on the
#' way out: both come back oat-rows x pea-cols, like every other model here.
fit_megalmm_both <- function(train, G_oat, G_pea, runID,
                             K = SIM_MEGALMM_K,
                             eigen_variance = SIM_EIGEN_VARIANCE,
                             fixed_main_effect = FALSE) {
  oat_side <- fit_megalmm_one(
    row_idx = train$oat, col_idx = train$pea, y = train$y_oat_std,
    G_row = G_oat, G_col = G_pea, runID = file.path(runID, "oat_side"),
    K = K, eigen_variance = eigen_variance,
    fixed_main_effect = fixed_main_effect)

  pea_side <- fit_megalmm_one(
    row_idx = train$pea, col_idx = train$oat, y = train$y_pea_std,
    G_row = G_pea, G_col = G_oat, runID = file.path(runID, "pea_side"),
    K = K, eigen_variance = eigen_variance,
    fixed_main_effect = fixed_main_effect)

  list(
    surface = list(oat = oat_side$surface, pea = t(pea_side$surface)),
    interaction = NULL,
    oat_side = oat_side, pea_side = pea_side,
    fixed_ok = c(oat = oat_side$fixed_ok, pea = pea_side$fixed_ok),
    n_dropped = oat_side$n_dropped + pea_side$n_dropped
  )
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

#' The additive half of the same split: M = additive_part(M) + interaction_part(M).
#'
#' Every model here hands back a full oat x pea surface, so the additive part
#' of that surface is the model's estimate of "producer + associate", whatever
#' machinery produced it. Nothing MegaLMM-specific is needed: its main effects
#' are spread across all the factors and the column intercepts, and averaging
#' over a margin collects them regardless.
#'
#' In particular `rowMeans(M)` is the oat main effect of the surface -- an
#' oat's average value over the pea environments -- which is what makes it a
#' better estimate of the true producer effect than the first factor's scores.
#' Expanding MegaLMM's own decomposition shows why:
#'
#'   rowMeans(U_F %*% Lambda + U_R) = U_F %*% rowMeans(Lambda) + rowMeans(U_R)
#'
#' every factor contributes, weighted by its mean loading, and nothing forces
#' the free factors' mean loadings to zero.
additive_part <- function(M) {
  outer(rowMeans(M), colMeans(M), "+") - mean(M)
}

score_predictions <- function(pred, sim, idx, label, train = NULL) {
  t <- sim$truth

  safe_cor <- function(a, b) {
    if (length(a) < 3) return(NA_real_)
    if (stats::sd(a) < 1e-10 || stats::sd(b) < 1e-10) return(NA_real_)
    stats::cor(a, b)
  }

  # ---- per-trait metrics, at the scored cells ----
  #
  # The truth is the STABLE effects plus the interaction. Environment-specific
  # deviations are unpredictable by construction and correctly count against
  # every model.
  per_trait <- function(surface, prod_vec, assoc_vec, I_true) {
    truth_gma   <- prod_vec[idx[, 1]] + assoc_vec[idx[, 2]]
    truth_total <- truth_gma + I_true[idx]
    c(
      total       = safe_cor(surface[idx], truth_total),
      gma         = safe_cor(additive_part(surface)[idx], truth_gma),
      interaction = safe_cor(interaction_part(surface)[idx],
                             interaction_part(I_true)[idx])
    )
  }

  oat_trait <- per_trait(pred$surface$oat, t$oat_prod, t$pea_assoc, t$I_oat)
  pea_trait <- per_trait(pred$surface$pea, t$oat_assoc, t$pea_prod, t$I_pea)

  # ---- fit against the observed phenotype ----
  #
  # NOT a cross-validation number. Nothing is held out, so this is measured on
  # the cells the model was fitted to and is a goodness of FIT, which is why it
  # is named r_fit_ and not r_obs_. It is kept because a model whose surface does
  # not track its own training data has gone wrong in a way the truth-based
  # metrics can hide.
  fit_r <- c(oat = NA_real_, pea = NA_real_)
  if (!is.null(train) && nrow(train) > 2) {
    oidx <- cbind(train$oat, train$pea)
    fit_r["oat"] <- safe_cor(pred$surface$oat[oidx], train$y_oat_std)
    fit_r["pea"] <- safe_cor(pred$surface$pea[oidx], train$y_pea_std)
  }

  # ---- effect recovery, from the margins of the surfaces ----
  #
  # Which margin gives which effect follows from the layout. Both surfaces are
  # oat-rows x pea-cols, so:
  #   oat-yield surface: row margin = the OAT's producer, col margin = the
  #                      PEA's associate
  #   pea-yield surface: row margin = the OAT's associate, col margin = the
  #                      PEA's producer
  # These are full-panel margins, one value per accession, and always were --
  # an accession effect has no held-out version.
  oat_prod_hat  <- rowMeans(pred$surface$oat)
  pea_assoc_hat <- colMeans(pred$surface$oat)
  oat_assoc_hat <- rowMeans(pred$surface$pea)
  pea_prod_hat  <- colMeans(pred$surface$pea)

  tibble::tibble(
    model = label,
    n_scored = nrow(idx),
    # oat yield
    r_total_oat = oat_trait[["total"]], r_addsurf_oat = oat_trait[["gma"]],
    r_int_oat = oat_trait[["interaction"]], r_fit_oat = fit_r[["oat"]],
    # pea yield
    r_total_pea = pea_trait[["total"]], r_addsurf_pea = pea_trait[["gma"]],
    r_int_pea = pea_trait[["interaction"]], r_fit_pea = fit_r[["pea"]],
    # the four effects
    r_oat_prod  = safe_cor(oat_prod_hat,  t$oat_prod),
    r_oat_assoc = safe_cor(oat_assoc_hat, t$oat_assoc),
    r_pea_prod  = safe_cor(pea_prod_hat,  t$pea_prod),
    r_pea_assoc = safe_cor(pea_assoc_hat, t$pea_assoc),
    # general mixing ability per SPECIES: producer + associate, each taken from
    # whichever surface carries it. Not the same quantity as r_addsurf_oat above,
    # which is the additive part of oat YIELD -- see SIMULATION_GLOSSARY.md.
    r_oat_gma = safe_cor(oat_prod_hat + oat_assoc_hat,
                         t$oat_prod + t$oat_assoc),
    r_pea_gma = safe_cor(pea_prod_hat + pea_assoc_hat,
                         t$pea_prod + t$pea_assoc)
  )
}

#' Margin-only baselines, as in the real analysis.
#'
#' The floor the models have to clear: each accession's own training mean, with
#' no borrowing of any kind. Returned in the same two-surface shape as the
#' models so it goes through the same scorer.
baseline_predictions <- function(train, n_oat, n_pea, oat_names, pea_names) {
  margins <- function(y) {
    grand <- mean(y)
    om <- rep(grand, n_oat); pm <- rep(grand, n_pea)
    o <- tapply(y, train$oat, mean); p <- tapply(y, train$pea, mean)
    om[as.integer(names(o))] <- o
    pm[as.integer(names(p))] <- p
    list(om = om, pm = pm, grand = grand)
  }
  mo <- margins(train$y_oat_std)
  mp <- margins(train$y_pea_std)

  named <- function(M) { dimnames(M) <- list(oat_names, pea_names); M }

  list(
    # the row species' own mean only
    row_mean = list(surface = list(
      oat = named(matrix(mo$om, n_oat, n_pea)),
      pea = named(matrix(mp$om, n_oat, n_pea))), interaction = NULL),
    # both margins
    both_means = list(surface = list(
      oat = named(outer(mo$om, mo$pm, "+") - mo$grand),
      pea = named(outer(mp$om, mp$pm, "+") - mp$grand)), interaction = NULL)
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

#' Everything a scenario's fits need: the training table and the cells to score.
#'
#' All observed plots are used for fitting -- see scoring_index() for why nothing
#' is held out. Deterministic, so the two halves of a scenario, run in separate
#' processes, always score the same cells.
prepare_scenario <- function(sim, score_set = SIM_SCORE_SET) {
  list(train = standardize_within_env(sim$obs),
       idx   = scoring_index(sim, score_set))
}

run_scenario_bglr <- function(sim, score_set = SIM_SCORE_SET, seed = 1L,
                              models = c("additive", "dge_ige")) {
  sp <- prepare_scenario(sim, score_set)
  oat_names <- rownames(sim$G_oat); pea_names <- rownames(sim$G_pea)

  timed <- function(expr) {
    t0 <- Sys.time()
    v <- force(expr)
    list(v = v, secs = as.numeric(difftime(Sys.time(), t0, units = "s")))
  }

  fits <- list()
  secs <- numeric(0)
  if ("additive" %in% models) {
    a <- timed(fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, FALSE, seed = seed))
    fits$additive <- a$v; secs <- c(secs, a$secs)
  }
  if ("dge_ige" %in% models) {
    d <- timed(fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, TRUE, seed = seed))
    fits$dge_ige <- d$v; secs <- c(secs, d$secs)
  }

  # The margin-only baselines are free -- no fitting, just training means -- so
  # they are always included as the floor every model has to clear.
  base <- baseline_predictions(sp$train, nrow(sim$G_oat), nrow(sim$G_pea),
                               oat_names, pea_names)
  preds <- c(fits, base)

  purrr::imap(preds, \(pp, nm)
              score_predictions(pp, sim, sp$idx, nm, train = sp$train)) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      seconds = c(secs, rep(NA_real_, length(base))),
      n_train = nrow(sp$train), .after = model
    )
}

run_scenario_megalmm <- function(sim, K, eigen_variance,
                                 score_set = SIM_SCORE_SET, seed = 1L,
                                 run_dir = tempdir(),
                                 fixed_main_effect = FALSE) {
  sp <- prepare_scenario(sim, score_set)

  t0 <- Sys.time()
  mm <- fit_megalmm_both(sp$train, sim$G_oat, sim$G_pea,
                         runID = run_dir, K = K,
                         eigen_variance = eigen_variance,
                         fixed_main_effect = fixed_main_effect)
  t_mm <- as.numeric(difftime(Sys.time(), t0, units = "s"))

  # `U` excludes the per-column intercept, which is where the column species'
  # main effect lives, so it is scored alongside rather than instead.
  mm_U <- list(surface = list(oat = mm$oat_side$U, pea = t(mm$pea_side$U)),
               interaction = NULL)

  scores <- dplyr::bind_rows(
    score_predictions(mm,   sim, sp$idx, "megalmm",   train = sp$train),
    score_predictions(mm_U, sim, sp$idx, "megalmm_U", train = sp$train)
  ) |>
    dplyr::mutate(seconds = c(t_mm, NA_real_),
                  n_train = nrow(sp$train),
                  K = K, eigen_variance = eigen_variance,
                  fixed_main_effect = fixed_main_effect,
                  .after = model)

  # DIAGNOSTICS of the fixed factor, not estimates of a main effect -- for that
  # see r_oat_prod and r_pea_prod in score_predictions(), which take margins of
  # the whole surface and so collect the main effect wherever the model put it.
  #
  # The pinned factor's scores are compared with the truth for the species on
  # the ROWS of that orientation, since that is whose main effect it carries.
  mainfactor <- function(side, G_row, truth_vec) {
    if (is.null(side$U_F)) return(NA_real_)
    kept <- match(rownames(side$U_F), rownames(G_row))
    if (anyNA(kept)) return(NA_real_)
    r <- stats::cor(side$U_F[, 1], truth_vec[kept])
    # a FREE factor's sign is arbitrary, so only its magnitude means anything
    if (fixed_main_effect) r else abs(r)
  }

  row_mean_r <- function(y, idx, truth_vec) {
    m <- tapply(y, idx, mean)
    stats::cor(m, truth_vec[as.integer(names(m))])
  }

  diagnostics <- tibble::tibble(
    r_mainfactor_oat = mainfactor(mm$oat_side, sim$G_oat, sim$truth$oat_prod),
    r_mainfactor_pea = mainfactor(mm$pea_side, sim$G_pea, sim$truth$pea_prod),
    r_rowmean_oat = row_mean_r(sp$train$y_oat_std, sp$train$oat,
                               sim$truth$oat_prod),
    r_rowmean_pea = row_mean_r(sp$train$y_pea_std, sp$train$pea,
                               sim$truth$pea_prod),
    fixed_ok_oat = mm$fixed_ok[["oat"]], fixed_ok_pea = mm$fixed_ok[["pea"]],
    n_dropped = mm$n_dropped
  )

  list(scores = dplyr::bind_cols(scores, diagnostics[rep(1, nrow(scores)), ]),
       trace = NULL)
}
