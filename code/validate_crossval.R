# ============================================================
# LEAVE-ONE-TRIAL-OUT: DOES A PREDICTED ASSOCIATE EFFECT MATERIALISE
# IN AN ENVIRONMENT THE MODEL NEVER SAW?
#
# This is the decision input for the validation trial.  Everything in
# VALIDATION_DESIGN.md is conditional on lambda -- how much of a predicted
# associate contrast actually appears in new plots -- and lambda moves power
# by about 40 points across its plausible range, more than every other choice
# combined.  It can be estimated from trials already in hand, so it should be,
# before seed is ordered.
#
# For each held-out trial: refit the bivariate model on the others using the
# same fit_producer_associate() the production script uses, then ask the
# held-out trial three questions.
#
#   1. CALIBRATION.  Regress the held-out partner yield on the training
#      associate BLUP.  A slope of 1 means the effects are delivered in full;
#      the slope IS lambda.
#   2. PREDICTIVE ABILITY.  Correlate each accession's adjusted mean in the
#      held-out trial with its training BLUP.
#   3. THE REHEARSAL.  Build As+/As- pools from the training fit alone, then
#      measure the contrast those pools actually show in the held-out trial.
#      This is the validation trial in miniature, on data we already have, and
#      it is the most direct estimate of what the real one would see.
#
# NOT every fold answers the same question.  The five 2025 trials share nearly
# all their accessions, so holding one out asks "same lines, new environment"
# -- which is what the validation trial will do.  B4I_2026_IL is almost
# disjoint from the rest, so holding it out asks "new germplasm" instead.
# Both are reported; only the first is the lambda the design needs.
#
# Outputs: output/validation/<vintage>/crossval_folds.csv
#          output/validation/<vintage>/crossval_summary.csv
#          output/validation/<vintage>/crossval_accession.csv
#          output/validation/<vintage>/crossval.png
# ============================================================

library(tidyverse)

here::i_am("code/validate_crossval.R")

source(here::here("code", "dge_ige_functions.R"))
source(here::here("code", "validation_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

out_root <- here::here("output", "validation")
vintage  <- format(Sys.Date(), "%Y-%m-%d")

pheno_file <- here::here("output", "B4I_intercrop_pheno.rds")

grm_files <- c(oat = here::here("data", "GRM_Avena.rds"),
               pea = here::here("data", "GRM_Pisum.rds"))
analysis_name_files <- c(oat = here::here("output", "oat_analysis_names.csv"),
                         pea = here::here("output", "pea_analysis_names.csv"))

# Shorter chains than the production fit: six of them are needed and the
# quantity of interest is a regression slope, not a variance component.
cv_nIter  <- 6000
cv_burnIn <- 1000
cv_seed   <- 12567

# A fold counts as "same lines, new environment" -- the question the
# validation trial asks -- when at least this share of the held-out trial's
# accessions also appear in the training trials.
same_lines_threshold <- 0.5

min_partners <- validation_setting("min_partners")
pr_quantile  <- validation_setting("pr_quantile")

# Pools for the rehearsal are capped by how many eligible accessions the
# held-out trial actually contains.
rehearsal_max_n <- 25L

# ============================================================
# Driver
# ============================================================

out_dir <- file.path(out_root, vintage)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

G <- list(
  oat = collapse_grm(read_grm(grm_files[["oat"]]), analysis_name_files[["oat"]]),
  pea = collapse_grm(read_grm(grm_files[["pea"]]), analysis_name_files[["pea"]])
)

pheno <- readRDS(pheno_file) |>
  dplyr::filter(!is.na(oat_yield), !is.na(pea_yield)) |>
  dplyr::transmute(
    trial    = studyName,
    oatAcc   = as.character(germplasmName),
    peaAcc   = as.character(intercropGermplasmName),
    oatYield = oat_yield, peaYield = pea_yield,
    block    = paste(studyYear, studyName, blockNumber)
  )

trials <- sort(unique(pheno$trial))
message(length(trials), " trials, ", nrow(pheno), " plots")

as_model_frame <- function(d) {
  d |>
    dplyr::mutate(trialF = droplevels(factor(trial)),
                  blockNumberF = droplevels(factor(block)))
}

# ------------------------------------------------------------
# One fold
# ------------------------------------------------------------

run_fold <- function(held) {
  train <- as_model_frame(dplyr::filter(pheno, trial != held))
  test  <- dplyr::filter(pheno, trial == held)

  message("\n--- holding out ", held, " (", nrow(test), " plots, training on ",
          nrow(train), ") ---")

  # BGLR prints its hyperparameter setup to stdout whatever `verbose` says,
  # which buries the fold summaries; keep the return value, drop the chatter.
  invisible(utils::capture.output(
    fit <- fit_producer_associate(train, G$oat, G$pea, seed = cv_seed,
                                  nIter = cv_nIter, burnIn = cv_burnIn)
  ))

  oat_b <- dplyr::select(fit$effects$oat, acc = accession, Pr = PrEff, As = AsEff)
  pea_b <- dplyr::select(fit$effects$pea, acc = accession, Pr = PrEff, As = AsEff)

  # How much of the held-out trial does the training fit actually know about?
  cover_oat <- mean(unique(test$oatAcc) %in% oat_b$acc)
  cover_pea <- mean(unique(test$peaAcc) %in% pea_b$acc)

  dat <- test |>
    dplyr::inner_join(dplyr::rename(oat_b, oat_Pr = Pr, oat_As = As),
                      by = c("oatAcc" = "acc")) |>
    dplyr::inner_join(dplyr::rename(pea_b, pea_Pr = Pr, pea_As = As),
                      by = c("peaAcc" = "acc")) |>
    # the trial mean is not something the model predicts, so remove it
    dplyr::mutate(
      peaYield = peaYield - mean(peaYield),
      oatYield = oatYield - mean(oatYield)
    )

  # --- 1. calibration: slope of realised on predicted ---
  # An accession's ASSOCIATE effect is read on the PARTNER's yield, and the
  # partner's own producer effect is adjusted for rather than left in the
  # residual.
  calib <- function(response, as_col, partner_pr_col, species) {
    if (nrow(dat) < 20) return(NULL)
    m <- stats::lm(stats::reformulate(c(as_col, partner_pr_col), response),
                   data = dat)
    s <- summary(m)$coefficients
    tibble::tibble(
      species = species,
      lambda = s[as_col, "Estimate"], lambda_se = s[as_col, "Std. Error"],
      lambda_p = s[as_col, "Pr(>|t|)"],
      producer_slope = s[partner_pr_col, "Estimate"],
      n_plots = nrow(dat)
    )
  }

  calibration <- dplyr::bind_rows(
    calib("peaYield", "oat_As", "pea_Pr", "oat"),
    calib("oatYield", "pea_As", "oat_Pr", "pea")
  )

  # --- 2. predictive ability at the accession level ---
  # Adjust each plot for the partner's predicted producer effect, then average
  # to the focal accession: the closest analogue to what the validation trial
  # measures.
  acc_level <- function(response, focal, partner_pr_col, as_col, species) {
    d <- dat |>
      dplyr::mutate(adj = .data[[response]] -
                      stats::coef(stats::lm(
                        stats::reformulate(partner_pr_col, response),
                        data = dat))[2] * .data[[partner_pr_col]]) |>
      dplyr::group_by(acc = .data[[focal]]) |>
      dplyr::summarise(realised = mean(adj), n_plots = dplyr::n(),
                       predicted = dplyr::first(.data[[as_col]]),
                       .groups = "drop") |>
      dplyr::mutate(species = species, held_out = held)
    d
  }

  acc_oat <- acc_level("peaYield", "oatAcc", "pea_Pr", "oat_As", "oat")
  acc_pea <- acc_level("oatYield", "peaAcc", "oat_Pr", "pea_As", "pea")

  ability <- dplyr::bind_rows(acc_oat, acc_pea) |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      r_accession = stats::cor(predicted, realised),
      n_accessions = dplyr::n(), .groups = "drop"
    )

  # --- 3. the rehearsal: pools built on training data only ---
  rehearse <- function(species, acc_tbl, blups, focal, response,
                       partner_pr_col) {
    partner_col <- if (focal == "oatAcc") "peaAcc" else "oatAcc"

    # Eligibility is a property of the TRAINING data -- how well an
    # accession's effect could be estimated -- not of the held-out trial.
    # Within any single trial an accession meets only about two partners, so
    # filtering on the test trial would reject almost everyone.
    partners <- train |>
      dplyr::distinct(.data[[focal]], .data[[partner_col]]) |>
      dplyr::count(.data[[focal]], name = "n_partners") |>
      dplyr::rename(acc = 1)

    plots_per <- dplyr::count(test, acc = .data[[focal]], name = "n_plots")

    cand <- blups |>
      dplyr::inner_join(partners, by = "acc") |>
      dplyr::inner_join(plots_per, by = "acc") |>   # must be IN the held-out trial
      dplyr::filter(n_partners >= min_partners)

    n <- min(rehearsal_max_n, floor(nrow(cand) / 3))
    if (n < 5) return(NULL)

    # same constrained rule the real selection uses. PEV is set to 0 because
    # the rehearsal only needs the pools, not their power.
    fake_inp <- list(species = species,
                     accessions = dplyr::mutate(cand, eligible = TRUE),
                     PEV_As = 0)
    bp <- tryCatch(build_pools(fake_inp, n, pr_quantile = pr_quantile),
                   error = function(e) NULL)
    if (is.null(bp)) return(NULL)

    pools <- dplyr::select(bp$pools, acc, pool)
    d <- dat |>
      dplyr::inner_join(pools, by = stats::setNames("acc", focal)) |>
      dplyr::mutate(adj = .data[[response]] -
                      stats::coef(stats::lm(
                        stats::reformulate(partner_pr_col, response),
                        data = dat))[2] * .data[[partner_pr_col]]) |>
      dplyr::group_by(acc = .data[[focal]], pool) |>
      dplyr::summarise(y = mean(adj), .groups = "drop")

    plus <- d$y[d$pool == "As+"]; minus <- d$y[d$pool == "As-"]
    if (length(plus) < 3 || length(minus) < 3) return(NULL)
    tt <- stats::t.test(plus, minus, var.equal = TRUE)

    tibble::tibble(
      species = species, pool_n = n,
      predicted_contrast = bp$dAs,
      realised_contrast = as.numeric(tt$estimate[1] - tt$estimate[2]),
      contrast_se = tt$stderr,
      lambda_pool = as.numeric(tt$estimate[1] - tt$estimate[2]) / bp$dAs,
      p_one_sided = if (as.numeric(tt$estimate[1] - tt$estimate[2]) > 0)
        tt$p.value / 2 else 1 - tt$p.value / 2
    )
  }

  rehearsal <- dplyr::bind_rows(
    rehearse("oat", NULL, oat_b, "oatAcc", "peaYield", "pea_Pr"),
    rehearse("pea", NULL, pea_b, "peaAcc", "oatYield", "oat_Pr")
  )

  # The rehearsal can come back empty when a held-out trial holds too few
  # eligible accessions to build pools from; the fold is still informative for
  # calibration, so carry it with the rehearsal columns missing rather than
  # failing the join.
  if (nrow(rehearsal) == 0) {
    rehearsal <- tibble::tibble(
      species = character(0), pool_n = integer(0),
      predicted_contrast = numeric(0), realised_contrast = numeric(0),
      contrast_se = numeric(0), lambda_pool = numeric(0),
      p_one_sided = numeric(0))
  }

  # A trial cannot show an effect bigger than its own variation. The B4I
  # trials differ roughly 40-fold in mean yield -- AL was a near-total crop
  # failure at 8 g/m2 against 127-383 elsewhere -- so a slope measured on the
  # absolute g/m2 scale is not comparable across them. Record what each fold
  # had to work with, and report the slope on a standardised scale too.
  quality <- tibble::tibble(
    species = c("oat", "pea"),
    # the response for a species' associate effect is its PARTNER's yield
    response_mean = c(mean(test$peaYield), mean(test$oatYield)),
    response_sd   = c(stats::sd(test$peaYield), stats::sd(test$oatYield))
  )

  folds <- calibration |>
    dplyr::left_join(ability, by = "species") |>
    dplyr::left_join(rehearsal, by = "species") |>
    dplyr::left_join(quality, by = "species") |>
    dplyr::mutate(
      held_out = held,
      coverage = dplyr::if_else(species == "oat", cover_oat, cover_pea),
      question = dplyr::if_else(coverage >= same_lines_threshold,
                                "same lines, new environment",
                                "new germplasm"),
      .before = 1
    )

  list(folds = folds, accessions = dplyr::bind_rows(acc_oat, acc_pea))
}

results <- purrr::map(trials, run_fold)

folds <- purrr::map(results, "folds") |> purrr::list_rbind()
acc_all <- purrr::map(results, "accessions") |> purrr::list_rbind()

readr::write_csv(folds, file.path(out_dir, "crossval_folds.csv"))
readr::write_csv(acc_all, file.path(out_dir, "crossval_accession.csv"))

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n\n=== Fold by fold ===\n")
print(folds |>
        dplyr::select(held_out, species, question, coverage, lambda, lambda_se,
                      r_accession, lambda_pool, realised_contrast,
                      predicted_contrast) |>
        as.data.frame(), row.names = FALSE, digits = 3)

same <- dplyr::filter(folds, question == "same lines, new environment")

# A trial whose response barely varies cannot exhibit an effect of any size,
# so its slope is noise rather than evidence about lambda. Flag those rather
# than letting them drag the average.
sd_floor <- 0.4 * stats::median(same$response_sd)
same <- dplyr::mutate(same, informative = response_sd >= sd_floor)

if (any(!same$informative)) {
  message("\n", sum(!same$informative), " fold(s) excluded from the pooled ",
          "lambda for having too little variation to show an effect:\n  ",
          paste(unique(same$held_out[!same$informative]), collapse = ", "))
}

summarise_folds <- function(d, label) {
  d |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      set = label,
      n_folds = dplyr::n(),
      lambda_mean = mean(lambda), lambda_median = stats::median(lambda),
      lambda_sd = stats::sd(lambda),
      # the across-fold spread of the slope is the As x location interaction
      # the power calculation asks for, on the same scale
      interaction_frac = stats::sd(lambda) / abs(mean(lambda)),
      r_accession_mean = mean(r_accession),
      lambda_pool_mean = mean(lambda_pool, na.rm = TRUE),
      lambda_pool_sd = stats::sd(lambda_pool, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::relocate(set, .after = species)
}

summary_tbl <- dplyr::bind_rows(
  summarise_folds(same, "all same-line folds"),
  summarise_folds(dplyr::filter(same, informative), "informative folds only")
)

readr::write_csv(summary_tbl, file.path(out_dir, "crossval_summary.csv"))

cat("\n=== Lambda, over the 'same lines, new environment' folds ===\n")
cat("    This is the number VALIDATION_DESIGN.md section 5 is conditional on.\n\n")
print(as.data.frame(summary_tbl), row.names = FALSE, digits = 3)

newg <- dplyr::filter(folds, question == "new germplasm")
if (nrow(newg) > 0) {
  cat("\n=== For contrast: folds that hold out mostly NEW germplasm ===\n")
  cat("    A different and harder question than the validation trial asks.\n\n")
  print(newg |> dplyr::select(held_out, species, coverage, lambda, r_accession) |>
          as.data.frame(), row.names = FALSE, digits = 3)
}

cat("\n=== What this means for the validation trial ===\n")

inputs <- validation_inputs(min_partners = min_partners)
n_loc  <- validation_setting("n_locations")
p_loc  <- validation_setting("plots_per_location")
n_pool <- validation_setting("n_per_pool")

measured <- purrr::map(unique(summary_tbl$species), \(sp) {
  s <- dplyr::filter(summary_tbl, species == sp, set == "informative folds only")
  inp <- inputs[[sp]]
  bp  <- build_pools(inp, n_pool[[sp]], pr_quantile = pr_quantile)

  # The mean of a handful of fold slopes is itself uncertain; carry that
  # through rather than treating lambda as known.
  se_lambda <- s$lambda_sd / sqrt(s$n_folds)
  lo <- max(0, s$lambda_mean - 1.96 * se_lambda)
  hi <- s$lambda_mean + 1.96 * se_lambda

  pw <- function(lam, int) contrast_power(
    bp$dAs, bp$sigma2_within, n_pool[[sp]], inp$sigma2_e, p_loc * n_loc,
    lambda = lam, sided = 1, interaction_frac = int, n_loc = n_loc)$power

  tibble::tibble(
    species = sp, n = n_pool[[sp]], plots = p_loc * n_loc,
    lambda = s$lambda_mean, lambda_se = se_lambda,
    lambda_lo = lo, lambda_hi = hi,
    interaction_frac = s$interaction_frac,
    power_at_lambda = pw(s$lambda_mean, 0),
    power_with_interaction = pw(s$lambda_mean, s$interaction_frac),
    power_at_lambda_lo = pw(lo, s$interaction_frac)
  )
}) |> purrr::list_rbind()

readr::write_csv(measured, file.path(out_dir, "crossval_power.csv"))

cat("Power at the MEASURED lambda, rather than at an assumed one\n")
cat("(", p_loc, " plots x ", n_loc, " locations, one-sided alpha 0.05):\n\n", sep = "")
print(as.data.frame(measured), row.names = FALSE, digits = 3)

cat("\n")
for (i in seq_len(nrow(measured))) {
  m <- measured[i, ]
  cat("  ", m$species, ": lambda = ", round(m$lambda, 2),
      " (95% CI ", round(m$lambda_lo, 2), " to ", round(m$lambda_hi, 2),
      " over ", "folds)\n", sep = "")
  verdict <- dplyr::case_when(
    m$power_with_interaction >= 0.8 ~
      "adequately powered even allowing for the fold-to-fold spread",
    m$power_with_interaction >= 0.6 ~
      "marginal once the fold-to-fold spread is allowed for",
    TRUE ~ "under-powered once the fold-to-fold spread is allowed for"
  )
  cat("      power ", round(m$power_at_lambda, 2), " ignoring the spread, ",
      round(m$power_with_interaction, 2), " allowing for it -- ", verdict,
      "\n", sep = "")
}

cat("\n  Caveats worth carrying into any decision:\n",
    "   * lambda is estimated from a handful of folds and its own ",
    "confidence interval is wide.\n",
    "   * The across-fold spread is treated here as genuine associate x ",
    "environment\n     interaction, but part of it is estimation noise in ",
    "each fold's slope.\n",
    "   * The model assumes associate effects are constant in absolute g/m2, ",
    "while the\n     trials differ several-fold in mean yield. See the fold ",
    "table's response_sd.\n", sep = "")

# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

p <- acc_all |>
  dplyr::left_join(dplyr::distinct(folds, held_out, species, question),
                   by = c("held_out", "species")) |>
  ggplot2::ggplot(ggplot2::aes(predicted, realised)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  ggplot2::geom_point(ggplot2::aes(colour = question), alpha = 0.5, size = 1) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black",
                       linewidth = 0.6, formula = y ~ x) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red",
                       linetype = 3) +
  ggplot2::facet_grid(species ~ held_out, scales = "free") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(
    title    = "Does a predicted associate effect show up in a trial the model never saw?",
    subtitle = paste0("vintage ", vintage,
                      "; each point an accession. Red dotted = perfect ",
                      "calibration (slope 1); black = fitted slope (lambda)."),
    x = "predicted associate effect (training trials, g/m2)",
    y = "realised partner yield, adjusted (held-out trial, g/m2)",
    colour = NULL
  )

ggplot2::ggsave(file.path(out_dir, "crossval.png"), p,
                width = 13, height = 6, dpi = 150)

message("\nwrote ", out_dir)
