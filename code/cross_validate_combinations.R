# ============================================================
# CROSS-VALIDATION OF THE BIVARIATE MODEL, OVER NEW COMBINATIONS
#
#   Rscript code/cross_validate_combinations.R --quick      # 2 folds, short chains
#   Rscript code/cross_validate_combinations.R              # 5 folds
#   Rscript code/cross_validate_combinations.R --folds 10 --kron-rank 80
#
# THE USE CASE THIS MEASURES. Oat and pea accessions have been evaluated in
# SOME combinations; we want to predict combinations that have not been grown.
# So the thing held out is a COMBINATION, not a plot and not an accession, and
# the constraint is that both of its components stay in training -- in different
# combinations. An accession the model has never seen is a separate question
# (that is accession-wise CV, item 7 of docs/B4I_followups.md) and a much harder
# one; this is the question the programme actually faces.
#
# WHY THAT CHANGES THE MASKING. Holding out plots at random does not work: 90.8%
# of combinations occur in a single plot, so a random plot holdout mostly leaves
# the combination itself in training through its other plots, and where it does
# not, the component accessions may vanish too. So:
#
#   * the unit is the combination -- ALL plots of a held-out combination go out
#     together, or the model has seen that exact pairing;
#   * a combination is maskable only if both components have at least
#     `min_train_combos` OTHER combinations remaining in training;
#   * folds are built greedily under that constraint, and the number of
#     combinations it refuses to mask is reported rather than hidden.
#
# 1,904 of 2,059 combinations are maskable at min_train_combos = 1. The 155 that
# are not belong to the 99 oat and 84 pea accessions with a single partner, whose
# producer and associate effects are aliased anyway -- so their exclusion is a
# property of the design, not of this script.
#
# WHAT IS COMPARED. The bivariate producer-associate model with and without the
# specific-combination term, plus the two margin-only baselines. The interaction
# uses the LOW-RANK Kronecker basis, not the exact kernel: the exact kernel is
# defined only over OBSERVED combinations and so cannot predict a held-out one
# at all, which is the whole point here. `--kron-rank` sets the rank.
#
# HOW ACCURACY IS SCORED. Predictions are genetic values; the observations carry
# trial and block effects. Both sides are centred within block before
# correlating, which removes the field effects without needing to reconstruct
# fitted coefficients for them. Blocks hold 50 plots at the median, so little
# genetic signal goes with them.
#
# Outputs: output/crossval/<date>/predictions.csv
#          output/crossval/<date>/accuracy.csv        overall, per model and trait
#          output/crossval/<date>/accuracy_by_support.csv
#          output/crossval/<date>/folds.csv
# ============================================================

library(tidyverse)

here::i_am("code/cross_validate_combinations.R")

source(here::here("code", "dge_ige_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
quick      <- "--quick" %in% args
n_folds    <- as.integer(arg_value("--folds", if (quick) "2" else "5"))
kron_rank  <- as.integer(arg_value("--kron-rank", "50"))
min_train  <- as.integer(arg_value("--min-train-combos", "1"))
n_iter     <- as.integer(arg_value("--n-iter", if (quick) "3000" else "12000"))
burn_in    <- as.integer(arg_value("--burn-in", if (quick) "600" else "2000"))
seed       <- as.integer(arg_value("--seed", "20261002"))

out_dir <- here::here("output", "crossval", format(Sys.Date(), "%Y-%m-%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------

grms <- list(
  oat = collapse_grm(read_grm(here::here("data", "GRM_Avena.rds")),
                     here::here("output", "oat_analysis_names.csv")),
  pea = collapse_grm(read_grm(here::here("data", "GRM_Pisum.rds")),
                     here::here("output", "pea_analysis_names.csv"))
)
dat <- b4i_plot_table(grms = grms)

combos <- dat |>
  dplyr::count(oatAcc, peaAcc, mixID, name = "plots")
oat_support <- dplyr::count(combos, oatAcc, name = "oat_combos")
pea_support <- dplyr::count(combos, peaAcc, name = "pea_combos")
combos <- combos |>
  dplyr::left_join(oat_support, by = "oatAcc") |>
  dplyr::left_join(pea_support, by = "peaAcc")

message(nrow(dat), " plots | ", nrow(combos), " combinations | ",
        dplyr::n_distinct(dat$oatAcc), " oat, ",
        dplyr::n_distinct(dat$peaAcc), " pea")
message(sprintf("%.1f%% of combinations occur in a single plot",
                100 * mean(combos$plots == 1)))

# ------------------------------------------------------------
# Folds
#
# Greedy under the support constraint: walk the combinations in random order and
# assign one to the current fold only if masking it would leave both components
# with at least `min_train` other combinations in training. A combination that
# can never be masked is reported, not silently dropped.
# ------------------------------------------------------------

assign_folds <- function(combos, n_folds, min_train, seed) {
  set.seed(seed)
  n <- nrow(combos)
  fold <- rep(NA_integer_, n)
  # how many combinations each accession has left available to mask
  for (k in seq_len(n_folds)) {
    oat_left <- table(combos$oatAcc)
    pea_left <- table(combos$peaAcc)
    # within a fold, masking is simultaneous, so count against the fold only
    oat_used <- stats::setNames(integer(length(oat_left)), names(oat_left))
    pea_used <- stats::setNames(integer(length(pea_left)), names(pea_left))
    target <- ceiling(sum(is.na(fold)) / (n_folds - k + 1))
    n_in <- 0L
    for (i in sample(which(is.na(fold)))) {
      if (n_in >= target) break
      o <- combos$oatAcc[i]; p <- combos$peaAcc[i]
      if (oat_left[[o]] - oat_used[[o]] - 1L >= min_train &&
          pea_left[[p]] - pea_used[[p]] - 1L >= min_train) {
        fold[i] <- k
        oat_used[[o]] <- oat_used[[o]] + 1L
        pea_used[[p]] <- pea_used[[p]] + 1L
        n_in <- n_in + 1L
      }
    }
  }
  fold
}

combos$fold <- assign_folds(combos, n_folds, min_train, seed)
n_never <- sum(is.na(combos$fold))
message(n_folds, " folds | ", sum(!is.na(combos$fold)), " combinations maskable | ",
        n_never, " never masked (single-partner components)")
print(dplyr::count(combos, fold, name = "combinations"))
readr::write_csv(combos, file.path(out_dir, "folds.csv"))

# ------------------------------------------------------------
# One fold
# ------------------------------------------------------------

MODELS <- c("additive", "dge_ige")

run_fold <- function(k) {
  held_mix <- combos$mixID[!is.na(combos$fold) & combos$fold == k]
  train <- dplyr::filter(dat, !mixID %in% held_mix)
  held  <- dplyr::filter(dat,  mixID %in% held_mix)

  # Blocks and trials present in the held-out set must also be in training, or
  # its field effects are not estimable. Masking combinations cannot normally
  # remove a whole block (14 blocks, 50+ plots each), but it is cheap to check.
  stopifnot(all(as.character(held$trialF) %in% as.character(train$trialF)))

  train <- dplyr::mutate(train, trialF = droplevels(trialF),
                         blockNumberF = droplevels(blockNumberF))

  # support measured on THIS fold's training set -- the quantity accuracy is
  # expected to depend on
  tr_combos <- dplyr::distinct(train, oatAcc, peaAcc)
  support <- held |>
    dplyr::distinct(mixID, oatAcc, peaAcc) |>
    dplyr::left_join(dplyr::count(tr_combos, oatAcc, name = "oat_train_combos"),
                     by = "oatAcc") |>
    dplyr::left_join(dplyr::count(tr_combos, peaAcc, name = "pea_train_combos"),
                     by = "peaAcc") |>
    dplyr::mutate(dplyr::across(dplyr::ends_with("_train_combos"),
                                ~ tidyr::replace_na(.x, 0L)))

  message("\n--- fold ", k, ": ", nrow(train), " train / ", nrow(held),
          " held plots over ", length(held_mix), " combinations ---")

  purrr::map(MODELS, function(mdl) {
    t0 <- Sys.time()
    f <- fit_producer_associate(
      train, grms$oat, grms$pea, seed = seed + k,
      nIter = n_iter, burnIn = burn_in,
      fit_mix_term = (mdl == "dge_ige"),
      kron_rank = if (mdl == "dge_ige") kron_rank else NA,
      saveAt = file.path(tempdir(), paste0("cv_", mdl, "_", k, "_")))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "s"))

    pull <- function(sp, col) {
      e <- f$effects[[sp]]
      stats::setNames(e[[col]], e$accession)
    }
    oat_pr <- pull("oat", "PrEff"); oat_as <- pull("oat", "AsEff")
    pea_pr <- pull("pea", "PrEff"); pea_as <- pull("pea", "AsEff")
    get <- function(v, nm) tidyr::replace_na(unname(v[nm]), 0)

    # oat yield = oat producer + pea associate; pea yield is the mirror
    pred_oat <- get(oat_pr, held$oatAcc) + get(pea_as, held$peaAcc)
    pred_pea <- get(pea_pr, held$peaAcc) + get(oat_as, held$oatAcc)

    if (mdl == "dge_ige" && !is.null(f$interaction)) {
      io <- f$interaction$oat; ip <- f$interaction$pea
      idx <- cbind(match(held$oatAcc, rownames(io)),
                   match(held$peaAcc, colnames(io)))
      ok <- !is.na(idx[, 1]) & !is.na(idx[, 2])
      add_o <- numeric(nrow(held)); add_p <- numeric(nrow(held))
      add_o[ok] <- io[idx[ok, , drop = FALSE]]
      add_p[ok] <- ip[idx[ok, , drop = FALSE]]
      pred_oat <- pred_oat + add_o
      pred_pea <- pred_pea + add_p
    }

    tibble::tibble(
      fold = k, model = mdl, seconds = secs,
      mixID = held$mixID, oatAcc = held$oatAcc, peaAcc = held$peaAcc,
      blockNumberF = as.character(held$blockNumberF),
      obs_oat = held$oatYield, obs_pea = held$peaYield,
      pred_oat = pred_oat, pred_pea = pred_pea
    ) |>
      dplyr::left_join(dplyr::select(support, mixID, oat_train_combos,
                                     pea_train_combos), by = "mixID")
  }) |> purrr::list_rbind() |>
    # the margin-only baselines: each accession's own training mean, no
    # borrowing of any kind. Free, and the floor the models must clear.
    dplyr::bind_rows(
      local({
        om <- train |> dplyr::group_by(oatAcc) |>
          dplyr::summarise(m = mean(oatYield), .groups = "drop")
        pm <- train |> dplyr::group_by(peaAcc) |>
          dplyr::summarise(m = mean(peaYield), .groups = "drop")
        tibble::tibble(
          fold = k, model = "own_mean", seconds = NA_real_,
          mixID = held$mixID, oatAcc = held$oatAcc, peaAcc = held$peaAcc,
          blockNumberF = as.character(held$blockNumberF),
          obs_oat = held$oatYield, obs_pea = held$peaYield,
          pred_oat = tidyr::replace_na(om$m[match(held$oatAcc, om$oatAcc)],
                                       mean(train$oatYield)),
          pred_pea = tidyr::replace_na(pm$m[match(held$peaAcc, pm$peaAcc)],
                                       mean(train$peaYield))
        ) |>
          dplyr::left_join(dplyr::select(support, mixID, oat_train_combos,
                                         pea_train_combos), by = "mixID")
      })
    )
}

preds <- purrr::map(seq_len(n_folds), run_fold) |> purrr::list_rbind()
readr::write_csv(preds, file.path(out_dir, "predictions.csv"))

# ------------------------------------------------------------
# Accuracy
#
# Centred within block on BOTH sides, so trial and block effects drop out
# without their coefficients having to be reconstructed.
# ------------------------------------------------------------

centred <- preds |>
  dplyr::group_by(model, blockNumberF) |>
  dplyr::mutate(dplyr::across(c(obs_oat, obs_pea, pred_oat, pred_pea),
                              ~ .x - mean(.x, na.rm = TRUE))) |>
  dplyr::ungroup()

safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 10 || stats::sd(a[ok]) < 1e-10 || stats::sd(b[ok]) < 1e-10)
    return(NA_real_)
  stats::cor(a[ok], b[ok])
}

accuracy <- centred |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    n_plots = dplyr::n(),
    r_oat = safe_cor(pred_oat, obs_oat),
    r_pea = safe_cor(pred_pea, obs_pea),
    seconds = round(mean(seconds, na.rm = TRUE)), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(r_oat))

cat("\n", strrep("=", 72),
    "\nAccuracy on held-out COMBINATIONS (block-centred)\n",
    strrep("=", 72), "\n", sep = "")
print(accuracy)
readr::write_csv(accuracy, file.path(out_dir, "accuracy.csv"))

# per fold, to see whether any single fold is carrying the answer
cat("\nby fold:\n")
print(centred |> dplyr::group_by(model, fold) |>
        dplyr::summarise(r_oat = round(safe_cor(pred_oat, obs_oat), 3),
                         r_pea = round(safe_cor(pred_pea, obs_pea), 3),
                         .groups = "drop") |>
        tidyr::pivot_wider(names_from = fold,
                           values_from = c(r_oat, r_pea)), width = 200)

# ------------------------------------------------------------
# Accuracy as a function of how much the components were seen
#
# The question behind the question: how many OTHER combinations does an
# accession need before its contribution to a new one can be predicted?
# ------------------------------------------------------------

bin_support <- function(x) {
  cut(x, breaks = c(0, 1, 2, 4, 8, Inf),
      labels = c("1", "2", "3-4", "5-8", "9+"), right = TRUE)
}

by_support <- centred |>
  dplyr::mutate(
    oat_bin = bin_support(oat_train_combos),
    pea_bin = bin_support(pea_train_combos),
    min_bin = bin_support(pmin(oat_train_combos, pea_train_combos))) |>
  dplyr::group_by(model, support = min_bin) |>
  dplyr::summarise(n_plots = dplyr::n(),
                   r_oat = round(safe_cor(pred_oat, obs_oat), 3),
                   r_pea = round(safe_cor(pred_pea, obs_pea), 3),
                   .groups = "drop")

cat("\n", strrep("=", 72),
    "\nAccuracy by training support (the SCARCER of the two components)\n",
    strrep("=", 72), "\n", sep = "")
print(by_support, n = 40)
readr::write_csv(by_support, file.path(out_dir, "accuracy_by_support.csv"))

cat("\nand by the OAT's support alone:\n")
print(centred |>
        dplyr::mutate(support = bin_support(oat_train_combos)) |>
        dplyr::group_by(model, support) |>
        dplyr::summarise(n = dplyr::n(),
                         r_oat = round(safe_cor(pred_oat, obs_oat), 3),
                         .groups = "drop") |>
        tidyr::pivot_wider(names_from = support, values_from = c(n, r_oat)),
      width = 200)

message("\nwrote ", out_dir)
