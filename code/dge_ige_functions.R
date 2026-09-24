# ============================================================
# SHARED DGE-IGE MODEL MACHINERY
#
# The pieces of the bivariate producer-associate model that do not depend on
# which experiment is being fitted: reading and collapsing relationship
# matrices, factoring them so a BRR term puts the coefficients on the
# accession scale, building incidence matrices, and pulling the covariance
# components and per-accession effects back out of a fit.
#
# Used by code/BGLR_multi_trait_model.R for the B4I experiment and by
# code/pilot/ for the earlier pilot, so the two fit the same model rather
# than two models that resemble each other.
#
# Why a BRR on Z L rather than an RKHS on Z G Z' is argued in BACKGROUND.md;
# the short version is that an RKHS term's coefficients live in the
# eigenvector basis of the plot kernel and carry no accession identity, while
# L %*% beta is the accession effect, named.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

# Accessions that share an analysis name are one genotype, so their GRM
# rows and columns are averaged into a single row and column.  Averaging
# the relationships of identical lines is the same as taking the GRM of
# their averaged marker profiles.
collapse_grm <- function(G, name_file) {
  if (!file.exists(name_file)) return(G)

  name_map <- readr::read_csv(name_file, show_col_types = FALSE)
  if (!all(c("germplasmName", "analysis_name") %in% names(name_map))) return(G)

  new_name <- rownames(G)
  hit <- match(new_name, name_map$germplasmName)
  new_name[!is.na(hit)] <- name_map$analysis_name[hit[!is.na(hit)]]

  if (identical(new_name, rownames(G))) return(G)

  f <- factor(new_name, levels = unique(new_name))
  A <- stats::model.matrix(~ 0 + f)
  colnames(A) <- levels(f)
  A <- sweep(A, 2, colSums(A), "/")

  Gc <- t(A) %*% G %*% A
  dimnames(Gc) <- list(colnames(A), colnames(A))

  message("  collapsed ", nrow(G), " -> ", nrow(Gc), " lines using ",
          basename(name_file))
  Gc
}

read_grm <- function(path) {
  if (!file.exists(path)) {
    stop("GRM not found:\n  ", path,
         "\nRun code/create_GRMs_T3.R first.", call. = FALSE)
  }
  g <- readRDS(path)
  if (is.list(g) && !is.null(g$G)) g$G else g
}

grm_factor <- function(G, tol = 1e-8) {
  e <- eigen(G, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  if (sum(keep) < ncol(G)) {
    message("  ", ncol(G) - sum(keep), " of ", ncol(G),
            " dimensions dropped as numerically null")
  }
  L <- sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), "*")
  rownames(L) <- rownames(G)
  L
}

# Incidence matrices, with factor levels forced to the GRM row order
incidence <- function(x, levels) {
  Z <- stats::model.matrix(~ 0 + factor(x, levels = levels))
  colnames(Z) <- levels
  Z
}

#' @param trait_names Column names of the response, in order. Passed
#'   explicitly rather than read from a global `Y`, because the pilot fits
#'   the same model to a different experiment.
#' @param roles Named vector mapping "Pr" and "As" to trait names for this
#'   term: for the oat kernel the effect on oat yield is the producer effect
#'   and the effect on pea yield the associate effect; for pea it is reversed.
covariance_components <- function(fit, term, trait_names, roles) {
  Omega <- fit$ETA[[term]]$Cov$Omega
  dimnames(Omega) <- list(trait_names, trait_names)
  tibble::tibble(
    term        = term,
    var_Pr      = Omega[roles[["Pr"]], roles[["Pr"]]],
    var_As      = Omega[roles[["As"]], roles[["As"]]],
    cov_PrAs    = Omega[roles[["Pr"]], roles[["As"]]],
    cor_PrAs    = cov_PrAs / sqrt(var_Pr * var_As)
  )
}

accession_effects <- function(fit, term, L, trait_names, roles) {
  g <- L %*% fit$ETA[[term]]$beta
  dimnames(g) <- list(rownames(L), trait_names)

  tibble::as_tibble(g, rownames = "accession") |>
    dplyr::transmute(
      accession,
      PrEff = .data[[roles[["Pr"]]]],
      AsEff = .data[[roles[["As"]]]],
      GMA   = PrEff + AsEff   # GMA = Pr + As, docs/B4I_Proposal_Models.docx
    )
}

# ------------------------------------------------------------
# Fitting the bivariate model
#
# Factored out of code/BGLR_multi_trait_model.R so that the cross-validation
# in code/validate_crossval.R refits the SAME model on each training fold
# rather than a second model that resembles it.  A cross-validation whose
# training fits differ from the production fit measures the wrong thing.
#
# Trait order is fixed here as (peaYield, oatYield) and the roles are assigned
# from it: for the oat kernel the effect on oat yield is the producer effect
# and the effect on pea yield the associate effect; for pea it is reversed.
# That ordering is the single highest-consequence convention in the project --
# swapping it turns every producer into an associate with no error -- so it
# lives in one place.
# ------------------------------------------------------------

DGE_IGE_TRAITS <- c("peaYield", "oatYield")

DGE_IGE_ROLES <- list(
  G_oat = c(Pr = "oatYield", As = "peaYield"),
  G_pea = c(Pr = "peaYield", As = "oatYield"),
  G_mix = c(Pr = "oatYield", As = "peaYield")
)

#' Fit the bivariate producer-associate model to a plot table.
#'
#' @param dat Plot table with oatAcc, peaAcc, oatYield, peaYield, trialF,
#'   blockNumberF.  Factors should already be droplevels()'d.
#' @param G_oat,G_pea Relationship matrices, already collapsed onto the
#'   analysis names.  Subset to the accessions in `dat` here.
#' @param fit_mix_term The specific-combination term.  Off by default: it is
#'   estimable only when combinations are replicated, and building its kernel
#'   costs an eigendecomposition the size of the number of combinations.
#' @return list(fit, L_oat, L_pea, effects, varcomp, resid_cov, accessions)
fit_producer_associate <- function(dat, G_oat, G_pea, seed = 12567,
                                   nIter = 20000, burnIn = 3000, thin = 10,
                                   fit_mix_term = FALSE, saveAt = NULL,
                                   verbose = FALSE) {

  oatAccs <- sort(unique(dat$oatAcc))
  peaAccs <- sort(unique(dat$peaAcc))

  missing_oat <- setdiff(oatAccs, rownames(G_oat))
  missing_pea <- setdiff(peaAccs, rownames(G_pea))
  if (length(missing_oat) || length(missing_pea)) {
    stop("accessions phenotyped but absent from the GRMs -- oat: ",
         length(missing_oat), ", pea: ", length(missing_pea), call. = FALSE)
  }

  Go <- G_oat[oatAccs, oatAccs, drop = FALSE]
  Gp <- G_pea[peaAccs, peaAccs, drop = FALSE]

  Y <- as.matrix(dat[, DGE_IGE_TRAITS])
  colnames(Y) <- DGE_IGE_TRAITS
  stopifnot(!anyNA(Y))

  # Full trial dummies and NO separate intercept: a full dummy set plus an
  # intercept is rank-deficient, and BGLR samples along the ridge rather than
  # erroring, so nothing means anything and the chain mixes badly.
  incTrials <- stats::model.matrix(~ 0 + trialF, dat)
  colnames(incTrials) <- levels(dat$trialF)

  # Blocks are informative only where a trial has more than one. A trial with
  # a single block gives a factor constant within it, perfectly aliased with
  # its own fixed trial effect.
  blocks_per_trial <- dat |>
    dplyr::distinct(trialF, blockNumberF) |>
    dplyr::count(trialF, name = "n_blocks")

  informative <- dat |>
    dplyr::left_join(blocks_per_trial, by = "trialF") |>
    dplyr::filter(n_blocks > 1) |>
    dplyr::pull(blockNumberF) |>
    unique() |>
    as.character()

  incBlocks <- stats::model.matrix(~ 0 + blockNumberF, dat)
  colnames(incBlocks) <- levels(dat$blockNumberF)
  incBlocks <- incBlocks[, colnames(incBlocks) %in% informative, drop = FALSE]

  Z_oat <- incidence(dat$oatAcc, rownames(Go))
  Z_pea <- incidence(dat$peaAcc, rownames(Gp))
  stopifnot(nrow(Z_oat) == nrow(Y), nrow(Z_pea) == nrow(Y))

  L_oat <- grm_factor(Go)
  L_pea <- grm_factor(Gp)

  ETA <- list(
    trial = list(X = incTrials,        model = "FIXED"),
    G_pea = list(X = Z_pea %*% L_pea,  model = "BRR"),
    G_oat = list(X = Z_oat %*% L_oat,  model = "BRR")
  )
  if (ncol(incBlocks) > 0) {
    ETA <- append(ETA, list(block = list(X = incBlocks, model = "BRR")),
                  after = 1)
  }

  L_mix <- NULL
  if (fit_mix_term) {
    mixIDs <- sort(unique(paste(dat$oatAcc, dat$peaAcc, sep = "::")))
    parts <- stringr::str_split_fixed(mixIDs, stringr::fixed("::"), 2)
    G_mix <- Go[parts[, 1], parts[, 1], drop = FALSE] *
             Gp[parts[, 2], parts[, 2], drop = FALSE]
    dimnames(G_mix) <- list(mixIDs, mixIDs)
    L_mix <- grm_factor(G_mix)
    Z_mix <- incidence(paste(dat$oatAcc, dat$peaAcc, sep = "::"), mixIDs)
    ETA$G_mix <- list(X = Z_mix %*% L_mix, model = "BRR")
  }

  set.seed(seed)
  fit <- BGLR::Multitrait(
    y = Y, ETA = ETA, intercept = FALSE,
    resCov = list(df0 = 4, S0 = NULL, type = "UN"),
    nIter = nIter, burnIn = burnIn, thin = thin,
    saveAt = saveAt %||% file.path(tempdir(), "dge_ige_"),
    verbose = verbose
  )

  terms <- intersect(names(DGE_IGE_ROLES), names(ETA))
  varcomp <- purrr::map(terms, \(tm)
    covariance_components(fit, tm, DGE_IGE_TRAITS, DGE_IGE_ROLES[[tm]])) |>
    purrr::list_rbind()

  R <- fit$resCov$R
  dimnames(R) <- list(DGE_IGE_TRAITS, DGE_IGE_TRAITS)

  list(
    fit = fit, L_oat = L_oat, L_pea = L_pea, L_mix = L_mix,
    effects = list(
      oat = accession_effects(fit, "G_oat", L_oat, DGE_IGE_TRAITS,
                              DGE_IGE_ROLES$G_oat),
      pea = accession_effects(fit, "G_pea", L_pea, DGE_IGE_TRAITS,
                              DGE_IGE_ROLES$G_pea)
    ),
    varcomp = varcomp,
    resid_cov = tibble::tibble(
      var_pea = R["peaYield", "peaYield"], var_oat = R["oatYield", "oatYield"],
      cov_pea_oat = R["peaYield", "oatYield"],
      cor_pea_oat = R["peaYield", "oatYield"] /
                    sqrt(R["peaYield", "peaYield"] * R["oatYield", "oatYield"])
    ),
    accessions = list(oat = oatAccs, pea = peaAccs),
    n_plots = nrow(dat), n_blocks_fitted = ncol(incBlocks)
  )
}

rank_summary <- function(effects, id) {
  effects |>
    dplyr::group_by(seed) |>
    dplyr::mutate(rank_GMA = dplyr::min_rank(dplyr::desc(GMA))) |>
    dplyr::group_by(.data[[id]]) |>
    dplyr::summarise(
      mean_rank = mean(rank_GMA),
      sd_rank   = stats::sd(rank_GMA),
      mean_GMA  = mean(GMA),
      sd_GMA    = stats::sd(GMA),
      mean_Pr   = mean(PrEff),
      mean_As   = mean(AsEff),
      .groups   = "drop"
    ) |>
    dplyr::arrange(mean_rank)
}
