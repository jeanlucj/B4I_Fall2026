# ============================================================
# BUILD THE INPUTS FOR THE OAT x PEA MegaLMM ANALYSIS
#
# Reframes the intercrop experiment from the oat's point of view: every pea
# accession is an "environment" in which oats were grown.  That gives an
# oat-accession x pea-environment matrix of oat grain yield, which is what
# MegaLMM factors.  The matrix is very sparse -- most oat x pea pairs were
# never grown together -- and whether the factor model can carry that
# sparsity is the question the analysis exists to answer.
#
# Three ways of filling a cell are supported, chosen with `fill`:
#
#   raw           the oat yield as measured
#   centered      oat yield minus its trial mean
#   standardized  centered, then divided by the trial's standard deviation
#
# Cells with more than one plot take the mean of the plots.
#
# Pea "environments" get covariates, which is what lets MegaLMM borrow
# strength across peas that were never paired with the same oats:
#   * pea grain yield BLUE across the whole experiment
#   * pea phenology BLUE (flowering by default; see `phenology_trait`)
#   * eigenvectors of the pea GRM, enough to reach `eigen_variance`
#
# Outputs: output/megalmm_inputs_<fill>.rds
# ============================================================

library(tidyverse)

here::i_am("code/megalmm_build_inputs.R")

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

pheno_file <- here::here("output", "B4I_intercrop_pheno.rds")
obs_file   <- here::here("output", "B4I_observations.rds")
cache_dir  <- here::here("output", "trial_cache")
out_dir    <- here::here("output")

grm_files <- c(
  oat = here::here("output", "GRM_Avena.rds"),
  pea = here::here("output", "GRM_Pisum.rds")
)

analysis_name_files <- c(
  oat = here::here("output", "oat_analysis_names.csv"),
  pea = here::here("output", "pea_analysis_names.csv")
)

# Build one input set per filling rule, so they can be compared on equal terms
fills <- c("raw", "centered", "standardized")

# A pea environment with one or two oats in it carries no estimable residual
# variance, and MegaLMM's ARD sampler produces NaN rather than failing
# cleanly on such a column. Trim iteratively until every environment has at
# least `min_obs_per_env` oats and every oat at least `min_obs_per_oat`
# environments.
min_obs_per_env <- 3
min_obs_per_oat <- 2

# Pea phenology covariate.  Flowering and maturity BLUEs correlate 0.68, so
# using both adds little; flowering is the default because it was scored in
# all six trials while maturity was not scored in AL.
phenology_trait <- "flower"   # "flower" or "maturity"

# Take as many pea-GRM eigenvectors as are needed to reach this share of the
# pea genetic variance.  The pea spectrum is flat (PC1 is only 11% of the
# total), so 80% costs about 99 eigenvectors.
eigen_variance <- 0.80

trait_names <- c(
  oat_yield    = "Grain yield - g/m2|CO_350:0000260",
  pea_yield    = "Pea Grain Yield - g/m2|CO_xxx:0003008",
  flower       = "Pea Flowering Date - 10% - Julian Day|CO_xxx:0003014",
  maturity     = "Pea Maturity Date - Julian Day|CO_xxx:0003015"
)

# ------------------------------------------------------------
# The oat x pea matrix
# ------------------------------------------------------------

#' Fill a cell value from the plot-level yields.
#'
#' Centering and standardizing are done WITHIN TRIAL, before anything is
#' averaged, because the trial is the unit that carries the management and
#' season effects.  Standardizing additionally puts every trial on the same
#' scale, which matters here because the trials differ several-fold in mean
#' yield and MegaLMM gives each column one residual variance.
scale_within_trial <- function(pheno, fill) {
  out <- pheno |>
    dplyr::filter(!is.na(oat_yield)) |>
    dplyr::group_by(studyName)

  out <- switch(
    fill,
    raw          = dplyr::mutate(out, cell_value = oat_yield),
    centered     = dplyr::mutate(out, cell_value = oat_yield - mean(oat_yield)),
    standardized = dplyr::mutate(
      out,
      cell_value = (oat_yield - mean(oat_yield)) / stats::sd(oat_yield)
    ),
    stop("fill must be one of raw, centered, standardized", call. = FALSE)
  )

  dplyr::ungroup(out)
}

#' Drop thin rows and columns until every one meets its minimum.
#'
#' Iterative because dropping a sparse column can take a row below its own
#' minimum and vice versa.
trim_matrix <- function(Y, min_env, min_oat) {
  repeat {
    keep_col <- colSums(!is.na(Y)) >= min_env
    Y <- Y[, keep_col, drop = FALSE]
    keep_row <- rowSums(!is.na(Y)) >= min_oat
    Y <- Y[keep_row, , drop = FALSE]
    if (all(keep_col) && all(keep_row)) break
  }
  Y
}

#' Oat accessions (rows) x pea accessions (columns), NA where never grown.
build_matrix <- function(pheno, fill) {
  cells <- scale_within_trial(pheno, fill) |>
    dplyr::group_by(germplasmName, intercropGermplasmName) |>
    dplyr::summarise(
      value   = mean(cell_value),
      n_plots = dplyr::n(),
      .groups = "drop"
    )

  wide <- cells |>
    dplyr::select(germplasmName, intercropGermplasmName, value) |>
    tidyr::pivot_wider(
      names_from  = intercropGermplasmName,
      values_from = value
    ) |>
    dplyr::arrange(germplasmName)

  mat <- as.matrix(dplyr::select(wide, -germplasmName))
  rownames(mat) <- wide$germplasmName
  mat <- mat[, sort(colnames(mat)), drop = FALSE]

  list(Y = mat, cells = cells)
}

# ------------------------------------------------------------
# Pea covariates
# ------------------------------------------------------------

#' BLUE per pea accession, adjusted for trial.
#'
#' Trial and accession both enter as fixed effects and the accession value is
#' averaged over trials, so an accession seen only in a high-yielding trial is
#' not credited with that trial's mean.  With one trial there is nothing to
#' adjust and the accession mean is returned.
pea_blue <- function(plot_data, value_col) {
  x <- plot_data |>
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(peaAcc)) |>
    dplyr::mutate(peaAcc = factor(peaAcc), studyName = factor(studyName))

  if (nlevels(x$studyName) < 2) {
    return(
      x |>
        dplyr::group_by(peaAcc) |>
        dplyr::summarise(blue = mean(.data[[value_col]]), .groups = "drop")
    )
  }

  fit <- stats::lm(
    stats::reformulate(c("studyName", "peaAcc"), response = value_col),
    data = x
  )

  grid <- tidyr::expand_grid(
    peaAcc    = levels(x$peaAcc),
    studyName = levels(x$studyName)
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), factor))

  grid$fitted <- stats::predict(fit, grid)

  grid |>
    dplyr::group_by(peaAcc) |>
    dplyr::summarise(blue = mean(fitted), .groups = "drop")
}

#' Pea trait BLUEs, on the analysis names so they join to the matrix columns.
pea_trait_blues <- function(obs_file, cache_dir, trait_names, name_file) {
  units <- list.files(cache_dir, "^units_", full.names = TRUE) |>
    purrr::map(readRDS) |>
    purrr::list_rbind() |>
    dplyr::select(observationUnitDbId, studyName, peaAcc = intercropGermplasmName)

  if (file.exists(name_file)) {
    nm <- readr::read_csv(name_file, show_col_types = FALSE)
    units <- units |>
      dplyr::left_join(dplyr::select(nm, germplasmName, analysis_name),
                       by = c("peaAcc" = "germplasmName")) |>
      dplyr::mutate(peaAcc = dplyr::coalesce(analysis_name, peaAcc)) |>
      dplyr::select(-analysis_name)
  }

  obs <- readRDS(obs_file) |>
    dplyr::filter(observationVariableName %in% trait_names,
                  !is.na(value), value != "") |>
    dplyr::mutate(
      trait = names(trait_names)[match(observationVariableName, trait_names)],
      value = suppressWarnings(as.numeric(value))
    ) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::inner_join(units, by = "observationUnitDbId")

  c("pea_yield", "flower", "maturity") |>
    purrr::map(\(tr) {
      d <- dplyr::filter(obs, trait == tr)
      if (nrow(d) == 0) return(NULL)
      pea_blue(dplyr::rename(d, !!tr := value), tr) |>
        dplyr::rename(!!tr := blue)
    }) |>
    purrr::compact() |>
    purrr::reduce(dplyr::full_join, by = "peaAcc")
}

#' Eigenvectors of the pea GRM covering `variance` of the genetic variance.
grm_eigenvectors <- function(G, variance) {
  e <- eigen(G, symmetric = TRUE)
  vals <- pmax(e$values, 0)
  k <- which(cumsum(vals) / sum(vals) >= variance)[1]

  V <- e$vectors[, seq_len(k), drop = FALSE]
  rownames(V) <- rownames(G)
  colnames(V) <- paste0("PC", seq_len(k))

  message("  ", k, " eigenvectors reach ", round(100 * variance), "% of the ",
          "pea genetic variance (PC1 alone is ",
          round(100 * vals[1] / sum(vals), 1), "%)")
  V
}

# Collapse a GRM onto the analysis names, as the model script does
collapse_grm <- function(G, name_file) {
  if (!file.exists(name_file)) return(G)
  nm <- readr::read_csv(name_file, show_col_types = FALSE)
  new_name <- rownames(G)
  hit <- match(new_name, nm$germplasmName)
  new_name[!is.na(hit)] <- nm$analysis_name[hit[!is.na(hit)]]
  if (identical(new_name, rownames(G))) return(G)
  f <- factor(new_name, levels = unique(new_name))
  A <- stats::model.matrix(~ 0 + f)
  colnames(A) <- levels(f)
  A <- sweep(A, 2, colSums(A), "/")
  Gc <- t(A) %*% G %*% A
  dimnames(Gc) <- list(colnames(A), colnames(A))
  Gc
}

read_grm <- function(path) {
  g <- readRDS(path)
  if (is.list(g) && !is.null(g$G)) g$G else g
}

# ------------------------------------------------------------
# Driver
# ------------------------------------------------------------

pheno <- readRDS(pheno_file)

G_oat_full <- collapse_grm(read_grm(grm_files[["oat"]]), analysis_name_files[["oat"]])
G_pea_full <- collapse_grm(read_grm(grm_files[["pea"]]), analysis_name_files[["pea"]])

blues <- pea_trait_blues(obs_file, cache_dir, trait_names,
                         analysis_name_files[["pea"]])

if (all(c("flower", "maturity") %in% names(blues))) {
  r <- stats::cor(blues$flower, blues$maturity, use = "complete.obs")
  message("pea flowering and maturity BLUEs correlate ", round(r, 3),
          "; using ", phenology_trait)
}

build_one <- function(fill) {
  message("\n=== ", fill, " ===")

  Y <- build_matrix(pheno, fill)$Y
  message("matrix: ", nrow(Y), " oat x ", ncol(Y), " pea, ",
          sum(!is.na(Y)), " filled cells (",
          round(100 * mean(!is.na(Y)), 2), "%)")

  # Keep only lines that have both a phenotype and a genotype
  Y <- Y[intersect(rownames(Y), rownames(G_oat_full)),
         intersect(colnames(Y), rownames(G_pea_full)), drop = FALSE]

  before <- dim(Y)
  Y <- trim_matrix(Y, min_obs_per_env, min_obs_per_oat)
  message("after trimming to >= ", min_obs_per_env, " oats per pea and >= ",
          min_obs_per_oat, " peas per oat: ", nrow(Y), " x ", ncol(Y),
          " (was ", before[1], " x ", before[2], "), ",
          sum(!is.na(Y)), " cells (", round(100 * mean(!is.na(Y)), 2), "%)")

  G_oat <- G_oat_full[rownames(Y), rownames(Y), drop = FALSE]
  G_pea <- G_pea_full[colnames(Y), colnames(Y), drop = FALSE]

  eig <- grm_eigenvectors(G_pea, eigen_variance)

  cov_tbl <- tibble::tibble(peaAcc = colnames(Y)) |>
    dplyr::left_join(
      dplyr::select(blues, peaAcc, dplyr::all_of(c("pea_yield", phenology_trait))),
      by = "peaAcc"
    )

  n_imputed <- sum(!stats::complete.cases(cov_tbl))
  cov_tbl <- cov_tbl |>
    dplyr::mutate(dplyr::across(
      -peaAcc,
      \(x) tidyr::replace_na(x, mean(x, na.rm = TRUE))
    ))
  if (n_imputed > 0) {
    message("  ", n_imputed, " pea environment(s) had no BLUE and took the mean")
  }

  X_pheno <- cov_tbl |>
    dplyr::select(-peaAcc) |>
    as.matrix() |>
    scale()
  rownames(X_pheno) <- cov_tbl$peaAcc

  X_Env <- cbind(intercept = 1, X_pheno, eig[colnames(Y), , drop = FALSE])

  # Group 1 is the intercept, group 2 the pea phenotypes, group 3 the markers.
  # MegaLMM's ARD prior shrinks each group separately, so 98 eigenvectors
  # cannot swamp the two phenotypic covariates by sheer number.
  X_Env_groups <- c(1, rep(2, ncol(X_pheno)), rep(3, ncol(eig)))

  message("  X_Env: ", nrow(X_Env), " environments x ", ncol(X_Env),
          " covariates (", paste(table(X_Env_groups), collapse = "/"),
          " in groups 1/2/3)")

  stopifnot(identical(rownames(X_Env), colnames(Y)))

  out <- list(
    fill = fill, Y = Y, G_oat = G_oat, G_pea = G_pea,
    X_Env = X_Env, X_Env_groups = X_Env_groups, blues = blues,
    settings = list(phenology_trait = phenology_trait,
                    eigen_variance  = eigen_variance,
                    n_eigenvectors  = ncol(eig),
                    min_obs_per_env = min_obs_per_env,
                    min_obs_per_oat = min_obs_per_oat)
  )

  out_file <- file.path(out_dir, paste0("megalmm_inputs_", fill, ".rds"))
  saveRDS(out, out_file)
  message("  wrote ", basename(out_file))
  invisible(out)
}

invisible(purrr::walk(fills, build_one))
