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

# A full set of dummy columns for a factor, no intercept.
#
# model.matrix(~ 0 + f) refuses a single-level factor -- "contrasts can be
# applied only to factors with 2 or more levels" -- even though the answer is
# an unambiguous column of ones. That case never arises in the six-trial
# production fit, but it does for a simulated experiment with one environment,
# and it is the caller's problem either way.
dummy_matrix <- function(f) {
  f <- droplevels(as.factor(f))
  if (nlevels(f) < 2) {
    Z <- matrix(1, length(f), 1, dimnames = list(NULL, levels(f)))
    return(Z)
  }
  Z <- stats::model.matrix(~ 0 + f)
  colnames(Z) <- levels(f)
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
# Low-rank bases for design matrices
#
# Moved here from code/sim_fit.R: the simulation and the DGE-IGE model both
# need them, and the specific-combination term cannot be built at scale
# without them. The index arithmetic in kron_basis() and the reshape that
# undoes it are derived in docs/specific-combination_kronecker.md -- read that
# before touching either, because getting it wrong fails silently.
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
#'   estimable only when combinations are replicated.
#' @param kron_rank How to build that term. `NA` uses the EXACT kernel
#'   `G_oat[i,i'] * G_pea[j,j']` over observed combinations, which costs an
#'   eigendecomposition the size of the number of combinations -- fine at the
#'   2,059 of the real experiment, impossible at the 19,200 of a dense
#'   simulation. An integer instead builds the low-rank row-wise Kronecker
#'   basis from the leading `kron_rank` eigenvectors of each species, giving a
#'   `kron_rank^2`-column design matrix whose fitted coefficients reshape to a
#'   full interaction surface. See docs/specific-combination_kronecker.md.
#' @return list(fit, L_oat, L_pea, effects, varcomp, resid_cov, accessions,
#'   interaction) -- `interaction` is a per-trait pair of full oat x pea
#'   surfaces when the mix term was fitted at low rank, else NULL.
fit_producer_associate <- function(dat, G_oat, G_pea, seed = 12567,
                                   nIter = 20000, burnIn = 3000, thin = 10,
                                   fit_mix_term = FALSE, kron_rank = NA,
                                   saveAt = NULL, verbose = FALSE) {

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
  incTrials <- dummy_matrix(dat$trialF)

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

  incBlocks <- dummy_matrix(dat$blockNumberF)
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
  kron_A <- kron_B <- NULL
  oat_idx <- match(dat$oatAcc, oatAccs)
  pea_idx <- match(dat$peaAcc, peaAccs)

  if (fit_mix_term) {
    if (is.na(kron_rank)) {
      # exact kernel over observed combinations
      mixIDs <- sort(unique(paste(dat$oatAcc, dat$peaAcc, sep = "::")))
      parts <- stringr::str_split_fixed(mixIDs, stringr::fixed("::"), 2)
      G_mix <- Go[parts[, 1], parts[, 1], drop = FALSE] *
               Gp[parts[, 2], parts[, 2], drop = FALSE]
      dimnames(G_mix) <- list(mixIDs, mixIDs)
      L_mix <- grm_factor(G_mix)
      Z_mix <- incidence(paste(dat$oatAcc, dat$peaAcc, sep = "::"), mixIDs)
      ETA$G_mix <- list(X = Z_mix %*% L_mix, model = "BRR")
    } else {
      # low-rank row-wise Kronecker basis: one column per (oat direction x pea
      # direction) pair, so the term is kron_rank^2 columns wide however many
      # combinations there are
      kron_A <- grm_basis(Go, rank = kron_rank)
      kron_B <- grm_basis(Gp, rank = kron_rank)
      ETA$G_mix <- list(X = kron_basis(kron_A, kron_B, oat_idx, pea_idx),
                        model = "BRR")
    }
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
    # Full oat x pea interaction surfaces, one per trait, when the mix term was
    # built at low rank. The coefficients come back as kron_rank^2 x 2 (one
    # column per trait); each column reshapes row-major to kron_rank x
    # kron_rank and expands to A %*% Beta %*% t(B). byrow = TRUE is not
    # optional -- see docs/specific-combination_kronecker.md.
    interaction = if (fit_mix_term && !is.na(kron_rank)) {
      b <- fit$ETA$G_mix$beta
      surf <- function(trait) {
        Beta <- matrix(b[, match(trait, DGE_IGE_TRAITS)],
                       nrow = ncol(kron_A), ncol = ncol(kron_B), byrow = TRUE)
        out <- kron_A %*% Beta %*% t(kron_B)
        dimnames(out) <- list(oatAccs, peaAccs)
        out
      }
      list(oat = surf("oatYield"), pea = surf("peaYield"))
    } else NULL,
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
