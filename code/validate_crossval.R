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
# Outputs: output/validation/<vintage>/crossval_jackknife.csv  what each trial contributes
#          output/validation/<vintage>/crossval_folds.csv
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

# THE SAME FILTER CHAIN THE PRODUCTION FIT USES. This read the raw phenotype
# file with no trial QC until 2026-10-03, so it ran a fold for every trial
# including the crop failure at AL, and AL sat in the TRAINING set of every
# other fold. Production discarded it. So lambda described a training policy
# the production BLUPs did not share -- SELF_CRITIQUE.md finding A2. With
# b4i_fit_frame() the folds and the fit see the same plots.
pheno <- b4i_fit_frame(pheno_file = pheno_file) |>
  dplyr::transmute(
    trial    = studyName,
    year     = studyYear,
    oatAcc, peaAcc, oatYield, peaYield,
    block    = as.character(blockNumberF)
  )

trials <- sort(unique(pheno$trial))
message(length(trials), " trials, ", nrow(pheno), " plots")
message("  ", paste(trials, collapse = ", "))

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
  #
  # TWO SCALES. `raw` is the slope in g/m2, which is what the power calculation
  # works in. `z` divides the held-out trial's response by its own SD first.
  # The two answer different questions, and section 6 of VALIDATION_DESIGN.md
  # has flagged the difference as unresolved: the model assumes associate
  # effects are constant in absolute g/m2, while the trials differ several-fold
  # in spread. A slope measured in absolute units across environments that
  # differ that much will vary across folds for reasons that are not genotype x
  # environment at all. interaction_frac is sd(lambda)/|mean(lambda)|, which is
  # scale-free, so the two scales' interaction_frac ARE comparable -- and the
  # comparison is the test.
  #
  # TWO LEVELS. The plot-level slope is unbiased but its SE treats plots as
  # independent, when accessions recur across plots and plots sit in blocks.
  # VALIDATION_DESIGN.md section 7b uses lambda_p as a trial-exclusion
  # criterion, so an anticonservative SE there is a decision error, not just a
  # cosmetic one. The accession-level fit is the honest one for inference.
  calib <- function(response, as_col, partner_pr_col, species, scale = "raw") {
    if (nrow(dat) < 20) return(NULL)
    d <- dat
    if (scale == "z") {
      sdy <- stats::sd(d[[response]])
      if (!is.finite(sdy) || sdy <= 0) return(NULL)
      d[[response]] <- d[[response]] / sdy
    }
    m <- stats::lm(stats::reformulate(c(as_col, partner_pr_col), response),
                   data = d)
    co <- summary(m)$coefficients

    # Cluster-robust SE by focal accession: a sandwich with the residuals summed
    # within accession. Written out rather than taken from a package, because
    # the repo has no sandwich dependency and the formula is three lines.
    X  <- stats::model.matrix(m)
    u  <- stats::residuals(m)
    cl <- d[[if (identical(species, "oat")) "oatAcc" else "peaAcc"]]
    XtX_inv <- chol2inv(chol(crossprod(X)))
    meat <- Reduce(`+`, lapply(split(seq_along(u), cl), function(ix) {
      xu <- crossprod(X[ix, , drop = FALSE], u[ix])
      tcrossprod(xu)
    }))
    G  <- length(unique(cl)); k <- ncol(X); nn <- nrow(X)
    adj <- (G / (G - 1)) * ((nn - 1) / (nn - k))
    V  <- adj * XtX_inv %*% meat %*% XtX_inv
    # chol2inv() DROPS DIMNAMES, so diag(V) comes back unnamed and indexing it
    # by the coefficient name silently returned NA -- which propagated into
    # every noise-corrected interaction. Index by position instead, and assert.
    j <- match(as_col, colnames(X))
    stopifnot("the associate predictor is not in the model matrix" = !is.na(j))
    se_cl <- sqrt(V[j, j])
    stopifnot("the clustered variance is not positive" = is.finite(se_cl))
    t_cl  <- co[as_col, "Estimate"] / se_cl

    tibble::tibble(
      species = species, scale = scale,
      lambda = co[as_col, "Estimate"],
      # On the `z` scale the response was divided by the held-out trial's SD,
      # so the slope is in SDs of partner yield per g/m2 of predicted effect
      # and its MAGNITUDE is not comparable with the raw slope. Multiplying by
      # that SD puts it back in g/m2. This is cosmetic -- interaction_frac is
      # sd/|mean| and so is unchanged by any constant rescaling -- but a table
      # showing lambda = 0.02 beside lambda = 0.77 invites the wrong reading.
      lambda_gm2 = co[as_col, "Estimate"] *
        (if (scale == "z") stats::sd(dat[[response]]) else 1),
      # the plot-level SE, retained but NOT to be used for inference
      lambda_se_plot = co[as_col, "Std. Error"],
      lambda_p_plot  = co[as_col, "Pr(>|t|)"],
      # clustered by focal accession: this is the one to quote
      lambda_se = se_cl,
      lambda_p  = 2 * stats::pt(abs(t_cl), df = G - 1, lower.tail = FALSE),
      se_inflation = se_cl / co[as_col, "Std. Error"],
      producer_slope = co[partner_pr_col, "Estimate"],
      n_plots = nrow(d), n_clusters = G
    )
  }

  calibration <- purrr::map(c("raw", "z"), \(sc) dplyr::bind_rows(
    calib("peaYield", "oat_As", "pea_Pr", "oat", sc),
    calib("oatYield", "pea_As", "oat_Pr", "pea", sc)
  )) |> purrr::list_rbind()

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
    # THE REHEARSAL NEEDS THE POOLS, NOT THEIR POWER, so the uncertainty
    # columns are set to zero rather than measured. Measuring them per fold
    # would mean streaming the coefficient draws from eight separate refits to
    # compute a quantity this function never uses. Eligibility is handed over
    # as a precomputed TRUE for the same reason: `cand` has already been
    # filtered on the training data's partner counts just above, which is the
    # right filter here -- a per-fold reliability is not available and the
    # question the rehearsal asks does not need one.
    fake_inp <- list(species = species,
                     accessions = dplyr::mutate(cand,
                       eligible  = TRUE,
                       GMA       = Pr + As,
                       PEV_As_i  = 0, PEV_Pr_i = 0, PEV_GMA_i = 0,
                       rel_As_i  = 1, rel_GMA_i = 1),
                     PEV_As = 0,
                     eligibility = list(rule = "partners", rel_min = 0,
                                        min_partners = min_partners,
                                        partner_floor = min_partners))
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
      year     = dplyr::first(test$year),
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
# response_sd is the held-out trial's own spread in g/m2, so this filter is a
# statement about the RAW scale; it is applied to both scales so the two are
# summarised over the same set of folds and their interaction_frac stays
# comparable. With the trial QC screen now applied upstream, the crop-failure
# trial never reaches here, so this may well exclude nothing -- which is
# itself worth saying, because it means lambda is no longer conditional on a
# fold-inclusion rule chosen by looking at lambda.
sd_floor <- 0.4 * stats::median(same$response_sd)
same <- dplyr::mutate(same, informative = response_sd >= sd_floor)

if (any(!same$informative)) {
  message("\n", sum(!same$informative), " fold(s) excluded from the pooled ",
          "lambda for having too little variation to show an effect:\n  ",
          paste(unique(same$held_out[!same$informative]), collapse = ", "))
}

summarise_folds <- function(d, label) {
  d |>
    dplyr::group_by(species, scale) |>
    dplyr::summarise(
      set = label,
      n_folds = dplyr::n(),
      lambda_mean = mean(lambda), lambda_median = stats::median(lambda),
      lambda_gm2_mean = mean(lambda_gm2),
      lambda_sd = stats::sd(lambda),
      # The mean of a handful of fold slopes is itself uncertain; this is the
      # SE of the GLOBAL estimate, which is what the power table is
      # conditioned on.
      lambda_se_of_mean = stats::sd(lambda) / sqrt(dplyr::n()),
      # The across-fold spread of the slope is the As x location interaction
      # the power calculation asks for. Scale-free, so it is comparable
      # between the raw and standardised scales -- which is the whole point of
      # fitting both.
      interaction_frac = stats::sd(lambda) / abs(mean(lambda)),
      # PART OF THAT SPREAD IS JUST ESTIMATION NOISE IN EACH FOLD'S SLOPE.
      # Under a constant true lambda, E[var(lambda_hat)] = var_true +
      # mean(se^2), so subtracting the mean squared SE leaves the genuine
      # across-environment variance. The SE used is the CLUSTERED one -- the
      # plot-level SE is too small, which would under-correct and leave the
      # interaction looking larger than it is.
      #
      # Report a zero as "not distinguishable from zero", never as "no
      # interaction": the correction can overshoot, and max(0, .) hides that.
      mean_lambda_se2 = mean(lambda_se^2),
      interaction_frac_corrected =
        sqrt(max(0, stats::var(lambda) - mean(lambda_se^2))) / abs(mean(lambda)),
      interaction_floored = stats::var(lambda) <= mean(lambda_se^2),
      r_accession_mean = mean(r_accession),
      lambda_pool_mean = mean(lambda_pool, na.rm = TRUE),
      lambda_pool_sd = stats::sd(lambda_pool, na.rm = TRUE),
      se_inflation_mean = mean(se_inflation),
      .groups = "drop"
    ) |>
    dplyr::relocate(set, .after = scale)
}

summary_tbl <- dplyr::bind_rows(
  summarise_folds(same, "all same-line folds"),
  summarise_folds(dplyr::filter(same, informative), "informative folds only")
)

readr::write_csv(summary_tbl, file.path(out_dir, "crossval_summary.csv"))

# ------------------------------------------------------------
# What each trial contributes: the jackknife
#
# interaction_frac is a property of the SET of folds -- sd(lambda)/|mean| -- so
# no single trial has one, and "does this trial widen the interaction?" cannot
# be read off any column. It is answered by removing that trial's fold and
# recomputing: `narrows_by` is how much the spread falls without it.
#
# This is what makes the decision rule in VALIDATION_DESIGN.md section 7b
# followable. A trial whose own fold has lambda near zero is not predicting
# itself out of sample, and it is also what widens the spread -- the two
# symptoms travel together, because a fold far from the others is both.
# ------------------------------------------------------------

jackknife <- same |>
  dplyr::group_by(species, scale) |>
  dplyr::mutate(
    lambda_all           = mean(lambda),
    interaction_all      = stats::sd(lambda) / abs(mean(lambda)),
    lambda_without       = purrr::map_dbl(dplyr::row_number(),
                                          \(i) mean(lambda[-i])),
    interaction_without  = purrr::map_dbl(dplyr::row_number(),
                             \(i) stats::sd(lambda[-i]) / abs(mean(lambda[-i]))),
    narrows_by           = interaction_all - interaction_without,
    raises_lambda_by     = lambda_without - lambda_all) |>
  dplyr::ungroup() |>
  dplyr::select(held_out, species, scale, year, informative, lambda, lambda_p,
                r_accession, lambda_all, lambda_without, raises_lambda_by,
                interaction_all, interaction_without, narrows_by) |>
  dplyr::arrange(species, scale, dplyr::desc(narrows_by))

readr::write_csv(jackknife, file.path(out_dir, "crossval_jackknife.csv"))

cat("\n=== What each trial contributes (leave-its-fold-out) ===\n")
cat("    narrows_by > 0 means the across-fold spread FALLS without this trial,\n",
    "    i.e. this trial is what widens the As x location interaction.\n\n", sep = "")
print(jackknife |>
        dplyr::select(held_out, species, lambda, lambda_p, r_accession,
                      raises_lambda_by, narrows_by) |>
        as.data.frame(), row.names = FALSE, digits = 3)

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

# ------------------------------------------------------------
# Lambda by year -- A DIAGNOSTIC, NOT AN ESTIMATE
#
# The 5-trial whitelist fit showed a striking split: oat folds held out from
# 2025 gave lambda around 1.4 and those from 2026 around 0.6. That is very
# likely an artifact of the training sets rather than a fact about years. Four
# of that fit's five trials were from 2025, so a held-out 2025 trial was
# predicted by closely related siblings and a held-out 2026 trial was not.
# With four trials from each year and every fold training on seven, the
# asymmetry should largely go.
#
# It is reported because if the split SURVIVES that rebalancing it is a real
# finding about year-to-year transfer and belongs in the write-up. It is NOT
# the number the trial is sized on: nothing in this chain re-sizes on a subset
# of folds, and "predict 2026 from 2025" is not a more valid question than its
# reverse.
#
# The year comes from the plot table, not from parsing the trial name -- the
# naming convention is not a contract.
# ------------------------------------------------------------

by_year <- same |>
  dplyr::group_by(species, scale, year) |>
  dplyr::summarise(n_folds = dplyr::n(),
                   lambda_mean = mean(lambda),
                   lambda_sd = stats::sd(lambda),
                   r_accession_mean = mean(r_accession),
                   .groups = "drop")

readr::write_csv(by_year, file.path(out_dir, "crossval_lambda_year.csv"))

cat("\n=== Lambda by year (diagnostic; the global estimate is what power uses) ===\n")
print(as.data.frame(by_year), row.names = FALSE, digits = 3)

spread <- by_year |>
  dplyr::filter(scale == "raw") |>
  dplyr::group_by(species) |>
  dplyr::summarise(gap = abs(diff(lambda_mean)),
                   pooled_sd = mean(lambda_sd), .groups = "drop")
cat("\n")
for (i in seq_len(nrow(spread))) {
  r <- spread[i, ]
  cat("  ", r$species, ": the two years differ by ", round(r$gap, 2),
      " against a within-year SD of ", round(r$pooled_sd, 2),
      if (is.finite(r$gap) && is.finite(r$pooled_sd) && r$gap < r$pooled_sd)
        " -- not a year effect worth reading"
      else " -- worth looking at",
      "\n", sep = "")
}

# ------------------------------------------------------------
# Lambda on the two scales, side by side
# ------------------------------------------------------------

scales_tbl <- summary_tbl |>
  dplyr::filter(set == "informative folds only") |>
  dplyr::select(species, scale, n_folds, lambda_mean, lambda_gm2_mean,
                lambda_sd, lambda_se_of_mean, interaction_frac,
                interaction_frac_corrected, interaction_floored,
                mean_lambda_se2, se_inflation_mean)

readr::write_csv(scales_tbl, file.path(out_dir, "crossval_lambda_scales.csv"))

cat("\n=== Does the across-fold spread survive standardising? ===\n")
cat("    interaction_frac is sd(lambda)/|mean(lambda)|, so it is scale-free\n")
cat("    and the two rows per species are directly comparable.\n\n")
print(as.data.frame(scales_tbl), row.names = FALSE, digits = 3)

cat("\n")
for (sp in unique(scales_tbl$species)) {
  a <- dplyr::filter(scales_tbl, species == sp, scale == "raw")
  b <- dplyr::filter(scales_tbl, species == sp, scale == "z")
  if (nrow(a) && nrow(b)) {
    cat("  ", sp, ": interaction_frac ", round(a$interaction_frac, 2),
        " raw vs ", round(b$interaction_frac, 2), " standardised",
        if (b$interaction_frac < 0.8 * a$interaction_frac)
          " -- a good part of the raw spread is scale heterogeneity, not interaction"
        else " -- standardising does not explain it, so the interaction is real",
        "\n", sep = "")
    cat("      after removing fold-level estimation noise: ",
        round(a$interaction_frac_corrected, 2), " (raw)",
        if (isTRUE(a$interaction_floored))
          "  [floored at zero: not distinguishable from no interaction]" else "",
        "\n", sep = "")
  }
}

cat("\n=== What this means for the validation trial ===\n")

inputs <- validation_inputs(min_partners = min_partners)
n_loc  <- validation_setting("n_locations")
p_loc  <- validation_setting("plots_per_location")
n_pool <- validation_setting("n_per_pool")

measured <- purrr::map(unique(summary_tbl$species), \(sp) {
  # The RAW scale is what the power calculation works in: the trial is sized in
  # g/m2. The standardised scale is reported alongside for the interaction
  # comparison, not substituted here.
  s <- dplyr::filter(summary_tbl, species == sp, scale == "raw",
                     set == "informative folds only")
  inp <- inputs[[sp]]
  bp  <- build_pools(inp, n_pool[[sp]], pr_quantile = pr_quantile)

  # The mean of a handful of fold slopes is itself uncertain; carry that
  # through rather than treating lambda as known.
  se_lambda <- s$lambda_sd / sqrt(s$n_folds)
  lo <- max(0, s$lambda_mean - 1.96 * se_lambda)
  hi <- s$lambda_mean + 1.96 * se_lambda

  # Named arguments throughout: contrast_power() takes only `dAs` positionally.
  # The per-pool within variance goes in as a length-2 vector so the Welch df
  # reflects any difference between the two pools.
  pw <- function(lam, int) contrast_power(
    dAs = bp$dAs,
    sigma2_within = c(bp$sigma2_within_plus, bp$sigma2_within_minus),
    n = n_pool[[sp]], sigma2_e = inp$sigma2_e, P = p_loc * n_loc,
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
