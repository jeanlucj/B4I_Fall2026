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
