# ============================================================
# MegaLMM ON THE OAT x PEA MATRIX
#
# Fits the factor model in which every pea accession is an environment for
# the oats, and answers the question the analysis exists for: can a factor
# model carry a matrix this sparse?
#
# The test is cross-validation on the observed cells. A random share of the
# filled cells is masked, the model is fitted to the rest, and the masked
# cells are predicted. Three comparisons make the answer interpretable:
#
#   * the three filling rules (raw / centered / standardized) against each
#     other, which settles empirically which values belong in the matrix;
#   * with and without the pea covariates, which is what tells us whether the
#     extended MegaLMM buys anything here;
#   * against margin-only baselines -- the oat's own mean, the pea's own mean,
#     and both -- because a factor model that cannot beat a row mean has not
#     found genotype-by-genotype structure, whatever its internal fit.
#
# The chosen configuration is then refitted on all the data and saved.
#
# Outputs: output/megalmm_cv_results.csv    the comparison above
#          output/megalmm_fit_<fill>.rds    posterior of the final fit
#          output/megalmm_cv_accuracy.png   the comparison, plotted
# ============================================================

library(tidyverse)

here::i_am("code/megalmm_oat_pea.R")

source(here::here("code", "megalmm_setup.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

out_dir <- here::here("output")
run_dir <- here::here("output", "megalmm_runs")

fills <- c("raw", "centered", "standardized")

# The configuration refitted on everything at the end. Chosen from the
# cross-validation below, not in advance: `centered` gives both the best
# MegaLMM accuracy and by far the best baselines, because subtracting the
# trial mean removes the between-trial differences that otherwise dominate
# every cell, while keeping yields on their natural g/m2 scale.
final_fill <- "centered"

K <- 10                  # factors; the ARD prior shrinks unused ones away

# Cross-validation
cv_fraction <- 0.20      # share of filled cells masked
cv_reps     <- 3
cv_seed     <- 20260920

# Shorter chains for the CV sweep, longer for the final fit
cv_burn_rounds  <- 5
cv_burn_iter    <- 60
cv_sample_iter  <- 250

final_burn_rounds <- 8
final_burn_iter   <- 100
final_sample_iter <- 400

fit_X_from <- 4          # burn-in round at which the covariates switch on

# ------------------------------------------------------------
# Cross-validation
# ------------------------------------------------------------

#' Mask a random share of the filled cells, leaving every row and column
#' with at least `floor_obs` training cells.
#'
#' Masking is over filled cells, not over the whole matrix, so the held-out
#' set is always something we can score against. The floor is not cosmetic:
#' an oat or pea left with nothing takes MegaLMM's ARD sampler to NaN, which
#' surfaces as "missing value where TRUE/FALSE needed" several frames deep in
#' `sample_Lambda_prec_ARD`, not as a clean complaint about empty data. Cells
#' are therefore taken in random order and kept only while both their row and
#' their column can spare one.
mask_cells <- function(Y, fraction, seed, floor_obs = 2) {
  filled <- which(!is.na(Y))
  idx <- arrayInd(filled, dim(Y))

  row_left <- rowSums(!is.na(Y))
  col_left <- colSums(!is.na(Y))

  set.seed(seed)
  order_try <- sample(seq_along(filled))
  target <- round(fraction * length(filled))

  held <- integer(0)
  for (k in order_try) {
    if (length(held) >= target) break
    i <- idx[k, 1]; j <- idx[k, 2]
    if (row_left[i] > floor_obs && col_left[j] > floor_obs) {
      held <- c(held, filled[k])
      row_left[i] <- row_left[i] - 1L
      col_left[j] <- col_left[j] - 1L
    }
  }

  Y_train <- Y
  Y_train[held] <- NA

  list(Y_train = Y_train, held = held, truth = Y[held],
       achieved = length(held) / length(filled))
}

#' Margin-only predictions, as the floor the factor model has to clear.
margin_baselines <- function(Y_train, held) {
  row_mean <- rowMeans(Y_train, na.rm = TRUE)
  col_mean <- colMeans(Y_train, na.rm = TRUE)
  grand    <- mean(Y_train, na.rm = TRUE)

  idx <- arrayInd(held, dim(Y_train))
  r <- row_mean[idx[, 1]]
  c <- col_mean[idx[, 2]]

  # An oat or pea whose training cells all vanished falls back to the grand mean
  r[is.na(r)] <- grand
  c[is.na(c)] <- grand

  list(
    oat_mean      = r,
    pea_mean      = c,
    oat_plus_pea  = r + c - grand
  )
}

#' One fit: mask, fit, predict the masked cells.
#'
#' @return a tibble of accuracy for the model and each baseline.
cv_once <- function(inputs, fraction, seed, use_covariates, label,
                    burn_rounds, burn_iter, sample_iter) {
  Y <- inputs$Y
  split <- mask_cells(Y, fraction, seed)

  message("  held out ", length(split$held), " of ", sum(!is.na(Y)),
          " cells (", round(100 * split$achieved, 1), "% of the ",
          round(100 * fraction), "% asked for)")
  stopifnot(all(rowSums(!is.na(split$Y_train)) > 0),
            all(colSums(!is.na(split$Y_train)) > 0))

  accNames <- tibble::tibble(germplasmName = rownames(Y))
  envCov <- if (use_covariates) {
    list(X_Env = inputs$X_Env, X_Env_groups = inputs$X_Env_groups,
         newEnv = character(0))
  } else NULL

  runID <- file.path(run_dir, paste0("cv_", label, "_", seed))

  state <- setup_megalmm_state(
    accNames = accNames, wideData = split$Y_train, kinMat = inputs$G_oat,
    envCov = envCov, runID = runID, K = K
  )
  state <- run_megalmm(state, burn_rounds = burn_rounds, burn_iter = burn_iter,
                       sample_iter = sample_iter, fit_X_from = fit_X_from)

  post <- megalmm_posterior(state, accessions = accNames$germplasmName)

  stopifnot(identical(dim(post$U), dim(Y)))

  # Score all three of MegaLMM's candidate targets rather than picking one.
  # Eta_mean is the predicted phenotype and is the like-for-like comparison
  # with the baselines; the two U matrices are genetic values, which omit the
  # per-environment intercept and so are at a disadvantage across columns.
  # Reporting all three removes the doubt that a weak result is just the wrong
  # target.
  megalmm_preds <- list(
    MegaLMM_Eta      = post$Eta_mean[split$held],
    MegaLMM_U        = post$U[split$held],
    MegaLMM_U_noU_R  = post$U_noU_R[split$held]
  )

  base <- margin_baselines(split$Y_train, split$held)

  scores <- c(megalmm_preds, base) |>
    purrr::imap(\(p, nm) tibble::tibble(
      predictor = nm,
      r = suppressWarnings(stats::cor(p, split$truth, use = "complete.obs")),
      rmse = sqrt(mean((p - split$truth)^2, na.rm = TRUE))
    )) |>
    purrr::list_rbind()

  scores |>
    dplyr::mutate(fill = inputs$fill, covariates = use_covariates,
                  seed = seed, n_held = length(split$held), .before = 1)
}

# ------------------------------------------------------------
# Driver
# ------------------------------------------------------------

dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- fills |>
  rlang::set_names() |>
  purrr::map(\(f) readRDS(file.path(out_dir, paste0("megalmm_inputs_", f, ".rds"))))

Y0 <- inputs[[1]]$Y
message("matrix: ", nrow(Y0), " oat accessions x ", ncol(Y0),
        " pea environments, ", sum(!is.na(Y0)), " filled cells (",
        round(100 * mean(!is.na(Y0)), 2), "%)")
message("covariates: ", ncol(inputs[[1]]$X_Env), " per pea environment")

grid <- tidyr::expand_grid(
  fill = fills,
  use_covariates = c(TRUE, FALSE),
  rep = seq_len(cv_reps)
) |>
  dplyr::mutate(seed = cv_seed + rep)

message("\ncross-validation: ", nrow(grid), " fits (", cv_fraction * 100,
        "% of cells held out)")

cv <- purrr::pmap(grid, function(fill, use_covariates, rep, seed) {
  message("\n--- ", fill, " | covariates = ", use_covariates, " | rep ", rep, " ---")
  cv_once(
    inputs = inputs[[fill]], fraction = cv_fraction, seed = seed,
    use_covariates = use_covariates,
    label = paste0(fill, "_", ifelse(use_covariates, "cov", "nocov")),
    burn_rounds = cv_burn_rounds, burn_iter = cv_burn_iter,
    sample_iter = cv_sample_iter
  )
}) |>
  purrr::list_rbind()

readr::write_csv(cv, file.path(out_dir, "megalmm_cv_results.csv"))

summary_tbl <- cv |>
  dplyr::group_by(fill, covariates, predictor) |>
  dplyr::summarise(
    mean_r  = mean(r), sd_r = stats::sd(r),
    mean_rmse = mean(rmse), .groups = "drop"
  ) |>
  dplyr::arrange(fill, covariates, dplyr::desc(mean_r))

cat("\n=== Cross-validated accuracy (correlation with held-out cells) ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE, digits = 3)

# --- does the factor model beat the margins, and do covariates help? ---
model_only <- summary_tbl |>
  dplyr::filter(stringr::str_starts(predictor, "MegaLMM")) |>
  dplyr::group_by(fill, covariates) |>
  dplyr::slice_max(mean_r, n = 1) |>
  dplyr::ungroup() |>
  dplyr::rename(best_megalmm = predictor)
best_base  <- summary_tbl |>
  dplyr::filter(!stringr::str_starts(predictor, "MegaLMM")) |>
  dplyr::group_by(fill, covariates) |>
  dplyr::slice_max(mean_r, n = 1) |>
  dplyr::ungroup()

comparison <- model_only |>
  dplyr::select(fill, covariates, best_megalmm, megalmm_r = mean_r) |>
  dplyr::left_join(
    dplyr::select(best_base, fill, covariates,
                  best_baseline = predictor, baseline_r = mean_r),
    by = c("fill", "covariates")
  ) |>
  dplyr::mutate(gain = megalmm_r - baseline_r)

cat("\n=== MegaLMM against the best margin-only baseline ===\n")
print(as.data.frame(comparison), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# Final fit on everything
# ------------------------------------------------------------

message("\n=== final fit: ", final_fill, ", all cells ===")

inp <- inputs[[final_fill]]
accNames <- tibble::tibble(germplasmName = rownames(inp$Y))

state <- setup_megalmm_state(
  accNames = accNames, wideData = inp$Y, kinMat = inp$G_oat,
  envCov = list(X_Env = inp$X_Env, X_Env_groups = inp$X_Env_groups,
                newEnv = character(0)),
  runID = file.path(run_dir, paste0("final_", final_fill)), K = K
)
state <- run_megalmm(state, burn_rounds = final_burn_rounds,
                     burn_iter = final_burn_iter,
                     sample_iter = final_sample_iter, fit_X_from = fit_X_from)

post <- megalmm_posterior(state, accessions = accNames$germplasmName)

cat("\n=== Final fit ===\n")
cat("Lambda: ", nrow(post$Lambda), " factors x ", ncol(post$Lambda),
    " pea environments\n", sep = "")

# How much of each factor is actually used: the ARD prior shrinks surplus
# factors towards zero, so this says whether K was generous enough.
factor_size <- rowSums(post$Lambda^2)
cat("factor sizes (sum of squared loadings):\n")
print(round(factor_size, 3))

cat("\nper-environment heritability implied by the factor model:\n")
print(summary(as.vector(post$h2)))

cat("\ngenetic correlations between pea environments:\n")
print(summary(post$G_cor[upper.tri(post$G_cor)]))

saveRDS(
  list(posterior = post, inputs_settings = inp$settings, fill = final_fill,
       K = K, cv = cv),
  file.path(out_dir, paste0("megalmm_fit_", final_fill, ".rds"))
)

# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

p_cv <- cv |>
  dplyr::mutate(
    covariates = ifelse(covariates, "with pea covariates", "no covariates"),
    predictor  = forcats::fct_reorder(predictor, r)
  ) |>
  ggplot2::ggplot(ggplot2::aes(predictor, r, colour = fill)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 3,
                        position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
                        position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::facet_wrap(~ covariates) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Predicting held-out oat x pea cells",
    subtitle = paste0(cv_fraction * 100, "% of filled cells masked, ",
                      cv_reps, " replicates; mean +/- se"),
    x = NULL, y = "correlation with held-out values", colour = "cell value"
  )

ggplot2::ggsave(file.path(out_dir, "megalmm_cv_accuracy.png"), p_cv,
                width = 10, height = 5, dpi = 150)

message("\nwrote:\n  output/megalmm_cv_results.csv",
        "\n  output/megalmm_fit_", final_fill, ".rds",
        "\n  output/megalmm_cv_accuracy.png")
