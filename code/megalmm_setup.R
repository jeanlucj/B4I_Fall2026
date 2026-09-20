# ============================================================
# MegaLMM MODEL CONSTRUCTION FOR THE OAT x PEA MATRIX
#
# Adapted from BarleyWinterSurvival/code/08_megalmm.R (which itself came from
# DSFAS_Recommending/code/setupMegaLMMstate.R), keeping the fixes made there:
#
#   * MegaLMM_control() and MegaLMM_priors() are called through do.call() with
#     a list of VALUES. They read their arguments with match.call() + eval(),
#     so a variable passed by name is looked up in the wrong frame and the call
#     dies with "object 'K' not found".
#   * The missing-data map is chosen with a bounds check, defaulting to the
#     last (least-grouped) candidate. The list's length depends on the data, so
#     a hard-coded index silently means something different on different runs.
#   * U rows come back labelled "<level>::<term>"; they are undecorated and
#     then CHECKED against the names that went in, because a silent mismatch
#     would attach every prediction to the wrong accession.
#
# What is new here is the environmental-covariate half, from
# DSFAS_Recommending/code/WSSW_CV0.R: the pea covariates enter as `X` in the
# Lambda prior so that a pea environment's factor loadings can be predicted
# from what kind of pea it is. With ~5 oats per pea environment that is not a
# refinement, it is the only reason the model has anything to work with.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

#' Build a MegaLMM_state for an oat x pea matrix.
#'
#' @param accNames One-column tibble of germplasmName in the ROW ORDER of
#'   `wideData`; MegaLMM evaluates the formula in this data frame.
#' @param wideData Oat accession x pea environment matrix (NA = never grown).
#' @param kinMat Oat relationship matrix; NULL gives the identity.
#' @param envCov NULL, or a list with `X_Env` (all environments x covariates,
#'   rownames = environment names), `X_Env_groups` (one group label per
#'   covariate column) and `newEnv` (environments to hold out entirely; may be
#'   `character(0)`, in which case the covariates still inform Lambda but no
#'   new-environment prediction is set up).
#' @param runID Folder for MegaLMM's run state (regeneratable; gitignored).
#' @param K Number of factors.
#' @param whichNAmap Which candidate missing-data map to use; NULL = the last.
setup_megalmm_state <- function(accNames, wideData, kinMat = NULL,
                                envCov = NULL, runID = "megalmm_run",
                                K = 15, whichNAmap = NULL, verbose = FALSE) {
  stopifnot(nrow(accNames) == nrow(wideData))

  if (is.null(kinMat)) {
    kinMat <- diag(nrow(accNames))
    rownames(kinMat) <- colnames(kinMat) <- accNames$germplasmName
  }

  # Split the covariates into training and held-out environments, and drop the
  # held-out columns from Y. With no held-out environments the covariates are
  # still attached; they just have nothing to extrapolate to.
  predict_new_env <- FALSE
  if (!is.null(envCov)) {
    newEnv <- envCov$newEnv %||% character(0)
    oldEnv <- setdiff(colnames(wideData), newEnv)
    envCov$X_Env_Train <- envCov$X_Env[oldEnv, , drop = FALSE]
    envCov$X_Env_Test  <- envCov$X_Env[newEnv, , drop = FALSE]
    wideData <- wideData[, oldEnv, drop = FALSE]
    predict_new_env <- length(newEnv) > 0
  }

  run_parameters <- do.call(MegaLMM::MegaLMM_control, list(
    h2_divisions = 20,
    burn = 0,          # burn-in is done manually in run_megalmm()
    thin = 2,
    K = K
  ))

  MegaLMM_state <- MegaLMM::setup_model_MegaLMM(
    Y        = wideData,
    formula  = ~ (1 | germplasmName),
    data     = accNames,
    relmat   = list(germplasmName = kinMat),
    run_parameters = run_parameters,
    run_ID   = runID
  )

  Lambda_prior <- list(
    sampler   = MegaLMM::sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2, rate = 1),
    delta_2   = list(shape = 3, rate = 1),
    delta_iterations_factor = 100
  )

  if (!is.null(envCov)) {
    Lambda_prior <- c(Lambda_prior, list(
      X       = envCov$X_Env_Train,
      X_group = envCov$X_Env_groups,
      # Start with Lambda free of the covariates and turn this on partway
      # through burn-in, so the factors settle before X is asked to explain them
      fit_X   = FALSE,
      Lambda_beta_var_shape = 3,
      Lambda_beta_var_rate  = 1
    ))
  }

  priors <- do.call(MegaLMM::MegaLMM_priors, list(
    tot_Y_var = list(V = 0.5, nu = 5),
    tot_F_var = list(V = 18 / 20, nu = 20),
    h2_priors_resids_fun  = function(h2s, n) 1,
    h2_priors_factors_fun = function(h2s, n) 1,
    Lambda_prior = Lambda_prior
  ))
  MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)

  if (predict_new_env) {
    MegaLMM_state$priors$X_Env_Test <- envCov$X_Env_Test
  }

  maps <- MegaLMM::make_Missing_data_map(
    MegaLMM_state,
    max_NA_groups = ncol(wideData) + 1,
    verbose = verbose
  )
  n_maps <- length(maps$Missing_data_map_list)

  # The list runs from most-grouped (index 1) to least. The default here is
  # the MOST grouped, the opposite of the choice that suits a dense matrix.
  # With about five oats per pea environment there is no per-column missingness
  # pattern worth preserving, and the least-grouped map does not merely run
  # slowly: it makes a group per column, one of which comes back empty, and
  # initialize_MegaLMM() dies in `ncol(ZL) < nrow(ZL[x, ]) * 0.9` with
  # "argument is of length zero". Step towards more grouping if the requested
  # map fails, so a denser matrix can still ask for a finer one.
  idx <- if (is.null(whichNAmap)) 1L else whichNAmap
  if (idx < 1 || idx > n_maps) {
    stop("whichNAmap = ", idx, " but make_Missing_data_map() returned ",
         n_maps, " candidate map(s).", call. = FALSE)
  }

  initialized <- NULL
  for (try_idx in seq(idx, 1L)) {
    candidate <- MegaLMM::set_Missing_data_map(
      MegaLMM_state, maps$Missing_data_map_list[[try_idx]])
    initialized <- tryCatch({
      candidate <- MegaLMM::initialize_variables_MegaLMM(candidate)
      MegaLMM::initialize_MegaLMM(candidate, verbose = verbose)
    }, error = function(e) {
      message("  missing-data map ", try_idx, " failed to initialize (",
              conditionMessage(e), "); trying a more grouped map")
      NULL
    })
    if (!is.null(initialized)) {
      message("  missing-data map ", try_idx, " of ", n_maps)
      break
    }
  }
  if (is.null(initialized)) {
    stop("no missing-data map initialized; the matrix may be too sparse",
         call. = FALSE)
  }
  MegaLMM_state <- initialized

  sample_params <- c("Lambda", "F_h2", "resid_h2", "tot_Eta_prec", "B1")
  if (!is.null(envCov)) {
    sample_params <- c(sample_params,
                       "U_F", "F", "B2_F", "Lambda_beta", "Lambda_beta_var")
  }
  MegaLMM_state$Posterior$posteriorSample_params <- sample_params
  MegaLMM_state$Posterior$posteriorMean_params   <- "Eta_mean"

  posterior_functions <- list(
    # The factor part alone, and the factor part plus the environment-specific
    # genetic residual. The second is the CV2 target: an oat's value in a pea
    # environment it was never grown in, from the peas it WAS grown with.
    U_CV2noU_R = "U_F %*% Lambda",
    U_CV2      = "U_F %*% Lambda + U_R",
    G  = "t(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])",
    R  = "t(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])",
    h2 = "(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])"
  )
  if (predict_new_env) {
    posterior_functions <- c(posterior_functions, list(
      U_CV0   = "U_F %*% Lambda_beta %*% t(X_Env_Test)",
      Eta_CV0 = "F %*% Lambda_beta %*% t(X_Env_Test)"
    ))
  }
  MegaLMM_state$Posterior$posteriorFunctions <- posterior_functions

  MegaLMM::clear_Posterior(MegaLMM_state)
}

#' Burn in, then sample.
#'
#' Burn-in is manual and in rounds. Each round before `fit_X_from` re-orders
#' the factors largest-to-smallest; from that round on, the covariates are
#' switched into the Lambda prior instead. Re-ordering restarts the chain, so
#' samples collected before it are discarded.
run_megalmm <- function(MegaLMM_state, burn_rounds = 8, burn_iter = 100,
                        sample_iter = 400, fit_X_from = 7, verbose = FALSE) {
  uses_X <- !is.null(MegaLMM_state$priors$Lambda_prior$X)

  message("burn-in: ", burn_rounds, " rounds x ", burn_iter, " iterations",
          if (uses_X) paste0(" (covariates switched on at round ", fit_X_from, ")") else "")

  for (i in seq_len(burn_rounds)) {
    if (uses_X && i >= fit_X_from) {
      MegaLMM_state$priors$Lambda_prior$fit_X <- TRUE
    } else {
      MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state,
                                                drop_cor_threshold = 0.6)
    }
    MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
    MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, burn_iter,
                                             verbose = verbose)
    message("  round ", i, " of ", burn_rounds, " done")
  }

  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  message("sampling: ", sample_iter, " iterations")
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, sample_iter,
                                           verbose = verbose)
  MegaLMM::save_posterior_chunk(MegaLMM_state)
}

#' Posterior means worth keeping.
megalmm_posterior <- function(MegaLMM_state, accessions = NULL) {
  .pm <- function(p) MegaLMM::get_posterior_mean(
    MegaLMM::load_posterior_param(MegaLMM_state, p))

  have <- names(MegaLMM_state$Posterior$posteriorFunctions)
  G <- .pm("G")

  out <- list(
    Lambda   = .pm("Lambda"),
    U        = .pm("U_CV2"),
    U_noU_R  = .pm("U_CV2noU_R"),
    Eta_mean = MegaLMM::load_posterior_param(MegaLMM_state, "Eta_mean"),
    G = G, R = .pm("R"), h2 = .pm("h2"),
    G_cor = stats::cov2cor(G)
  )
  if ("U_CV0" %in% have) out$U_CV0 <- .pm("U_CV0")

  out$U       <- .undecorate_rows(out$U, accessions)
  out$U_noU_R <- .undecorate_rows(out$U_noU_R, accessions)
  out
}

# MegaLMM labels U rows "<level>::<term>". Strip the suffix, then verify
# against the names that went in rather than trusting the order.
.undecorate_rows <- function(m, accessions) {
  if (is.null(m)) return(m)
  rownames(m) <- stringr::str_remove(rownames(m), "::.*$")
  if (!is.null(accessions)) {
    if (!identical(rownames(m), as.character(accessions))) {
      stop("MegaLMM returned U rows that do not match the accessions passed ",
           "in. Predictions would be attached to the wrong accessions.",
           call. = FALSE)
    }
  }
  m
}
