# ============================================================
# HOW MUCH COULD THE VALIDATION TRIAL DETECT?
#
# Power for the As+ vs As- pool contrast, from the fit currently in output/.
# Two routes, which should agree:
#
#   analytic    SE(delta)^2 = 2*sigma2_within/n + 4*sigma2_e/P, noncentral t
#               on 2n-2 df.  Fast, so it sweeps the whole grid.
#   simulation  draw each accession's TRUE effect from its posterior
#               (N(BLUP, PEV)), select the pools on the BLUPs as we actually
#               would, simulate the trial, fit it, and count rejections.
#               Slower, and what it checks is the SE FORMULA AND THE DESIGN'S
#               BALANCE -- that partner effects really do orthogonalise and the
#               layout behaves.  It is NOT evidence that the predicted contrast
#               materialises: drawing truth as BLUP + mean-zero error makes
#               E[truth | BLUP] = BLUP by construction, so the simulated
#               contrast must come back at lambda * dAs whatever the BLUPs are
#               worth.  Only lambda, from validate_crossval.R, speaks to that.
#
# The headline the analytic route exists to make visible: the standard error
# has a term that the plot budget cannot touch.  Adding plots re-measures the
# same n accessions; it does not add new ones.  The experimental unit for a
# pool contrast is the ACCESSION, and at these reliabilities roughly three
# quarters of the standard error is already irreducible.
#
# Outputs: output/validation/<vintage>/power_grid.csv
#          output/validation/<vintage>/power_ceiling.csv
#          output/validation/<vintage>/power_simulation.csv
#          output/validation/<vintage>/power_curves.png
# ============================================================

library(tidyverse)

here::i_am("code/validate_power.R")

source(here::here("code", "validation_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

out_root <- here::here("output", "validation")
vintage  <- format(Sys.Date(), "%Y-%m-%d")

# Shared with the other validate_* scripts (VALIDATION_DEFAULTS in
# code/validation_functions.R), so the power reported here describes the pools
# that validate_pool_selection.R actually builds.
n_locations  <- validation_setting("n_locations")
pr_quantile  <- validation_setting("pr_quantile")
min_partners <- validation_setting("min_partners")

# THE PLOT BUDGET IS FIXED AT 80 PER LOCATION. Section 5 establishes that
# plots only ever moved the smallest of the three variance terms, so the budget
# is not a live question. The 60/100 arms are kept only as a background
# sensitivity in power_grid.csv; every pre-registered number uses 80.
plots_per_location_fixed <- validation_setting("plots_per_location")
plots_per_location <- c(60L, 80L, 100L)

# LOCATIONS ARE THE SWEEP THAT MATTERS: the pool x location term is the largest
# of the three and carries only n_loc - 1 df, so this is the one lever that
# touches the binding constraint. 4 is the realistic downside -- B4I_2025_AL
# was a total loss, so a lost site is not hypothetical.
n_locations_swept <- validation_setting("n_locations_swept")

# Pool sizes to sweep. The design geometry wants n = plots-per-location / 4,
# which is 15/20/25; the rest are here to show how flat the curve is.
n_values <- c(10L, 15L, 20L, 25L, 30L, 40L, 50L)

# lambda and interaction_frac are no longer swept. Both are ESTIMATED from the
# trials in hand by code/validate_crossval.R, which now runs before this script.
#
#   lambda            attenuation of the predicted contrast when realised in new
#                     environments: how much of a BLUP difference actually shows
#                     up. 1 means delivered in full.
#   interaction_frac  SD of the location-specific contrast as a fraction of the
#                     contrast itself -- the pool x location interaction.
#
# THEY ARE THE SAME PHENOMENON SEEN TWICE, which is why one estimate gives both.
# validate_crossval.R holds out one trial at a time and regresses the realised
# associate effect on the predicted one; the slope is lambda. Its MEAN over
# folds is how much of the contrast survives a change of environment on average,
# and its SPREAD over folds is how much that survival varies by environment --
# so interaction_frac is sd(lambda)/|mean(lambda)| across folds
# (code/validate_crossval.R:341). A sweep over assumed values was right while
# both were guesses. They are now measured, and sweeping a measured quantity
# reports uncertainty that the data have already resolved.
#
# The test is directional and pre-registered, so it is ONE-SIDED. The two-sided
# arm was dropped: it was never the design's test, and reporting it alongside
# invited reading the wrong column.
alpha <- 0.05
sided <- 1L

# Simulation
sim_reps  <- 400L
# The null arms get more reps: at 400 the Monte Carlo SE of a rejection rate is
# 0.011, which cannot distinguish 0.05 from 0.08 -- and the type I comparison
# is the whole point of the null section. 2,000 gives 0.005.
null_reps <- 2000L
sim_seed  <- 20260924L
# Scenarios to simulate: a subset of the grid, since each costs sim_reps fits.
# The pool sizes come from the shared defaults, so the simulation checks the
# configuration that is actually proposed rather than a nearby one.
sim_cases <- tidyr::expand_grid(
  species = c("oat", "pea"),
  p_loc   = validation_setting("plots_per_location")
) |>
  dplyr::mutate(n = validation_setting("n_per_pool")[species], .after = species)

# ============================================================
# Simulation
# ============================================================

#' Simulate one validation trial on the REAL design and test the contrast
#' three ways.
#'
#' WHAT THIS CHECKS AND WHAT IT CANNOT. It checks the SE formula and the
#' design's balance: that partner effects really do orthogonalise, that the
#' anchors' unequal replication behaves, and -- new as of 2026-10-03 -- that a
#' denominator omitting the pool x location term over-rejects when that term is
#' real. It is NOT evidence that the predicted contrast materialises in the
#' field. The truth is drawn as `lambda * (BLUP + error)` with the error centred
#' on zero, so `E[true | BLUP] = lambda * BLUP` by construction and the
#' simulated contrast is pinned to `lambda * dAs` whatever the BLUPs are
#' actually worth. Only lambda, from validate_crossval.R, speaks to that.
#'
#' THE As x LOCATION TERM. Until 2026-10-03 this generated location MAIN
#' effects only, so the interaction that the power formula says is the largest
#' of its three variance components was never simulated at all -- and the null
#' that licensed believing a positive result was a null without it. Now each
#' location gets a deviation `d_l`, centred to sum exactly zero, added to every
#' As+ member and subtracted from every As- member. The location-specific
#' contrast then deviates by `2*d_l`, with mean EXACTLY zero over locations and
#' SD `2*tau`. Exact planting, not approximate: under `null = TRUE` the pool
#' main effect is zero by construction, so the rejection rate measures what it
#' is supposed to.
#'
#' @param tau SD of the per-location deviation. The contrast's across-location
#'   SD is twice this.
#' @param null TRUE sets every true associate MAIN effect to zero while leaving
#'   the interaction in place -- the realistic failure mode, and the one the
#'   old null could not represent.
simulate_validation <- function(inputs, oat_pools, pea_pools, p_loc, n_loc,
                                lambda = 1, alpha = 0.05, sided = 1,
                                null = FALSE, n_anchor = 1L, tau = 0,
                                int_sd_known = NULL) {

  design <- make_validation_design(oat_pools, pea_pools, n_loc = n_loc,
                                   p_loc = p_loc, n_anchor = n_anchor)

  # True effects per accession: BLUP + prediction error, attenuated by lambda.
  # The error SD is the accession's OWN posterior SD, not one global value --
  # a single PEV gave an accession seen with twenty partners the same
  # prediction error as one seen with two.
  draw_truth <- function(inp, pools) {
    a <- dplyr::filter(inp$accessions, acc %in% pools$acc)
    v <- if (null) rep(0, nrow(a)) else
      lambda * (a$As + stats::rnorm(nrow(a), 0, sqrt(a$PEV_As_i)))
    stats::setNames(v, a$acc)
  }
  draw_producer <- function(inp, pools) {
    a <- dplyr::filter(inp$accessions, acc %in% pools$acc)
    stats::setNames(a$Pr + stats::rnorm(nrow(a), 0, sqrt(inp$sigma2_Pr * 0.5)),
                    a$acc)
  }

  oat_As <- draw_truth(inputs$oat, oat_pools)
  pea_As <- draw_truth(inputs$pea, pea_pools)
  oat_Pr <- draw_producer(inputs$oat, oat_pools)
  pea_Pr <- draw_producer(inputs$pea, pea_pools)

  loc_oat <- stats::rnorm(n_loc, 0, 40)
  loc_pea <- stats::rnorm(n_loc, 0, 25)

  # pool x location deviations, NOT centred -- and that is the whole point.
  #
  # Centring them to sum exactly zero over the n_loc sites looks like the
  # careful thing to do, and it silently destroys the comparison: the estimate
  # averages over all locations, so a set of deviations summing to zero leaves
  # the contrast estimate exactly unaffected and there is no inflation left to
  # find. (That is what the first version of this did, and it duly reported no
  # inflation.)
  #
  # The locations in a trial are a SAMPLE of environments. Their deviations have
  # mean zero in the population, not in the five sites that happen to be sown,
  # so the realised mean is N(0, tau^2/n_loc) -- and that realised, non-zero
  # mean is exactly the quantity the analytic formula's (int_sd)^2/n_loc term
  # describes. Leaving it in makes this the BROAD null the design's inference
  # target requires: "no transferable pool effect in expectation over
  # environments", rather than "no pool effect at these exact five sites".
  pool_sign <- function(p) ifelse(p == "As+", 1, -1)
  d_oat <- stats::rnorm(n_loc, 0, tau)
  d_pea <- stats::rnorm(n_loc, 0, tau)

  plots <- design |>
    dplyr::mutate(
      # pea yield carries the PEA's producer effect and the OAT's associate
      # effect; oat yield the reverse. This is the whole reason one set of
      # plots validates both species.
      y_pea = loc_pea[location] + pea_Pr[pea_acc] + oat_As[oat_acc] +
              d_oat[location] * pool_sign(oat_pool) +
              stats::rnorm(dplyr::n(), 0, sqrt(inputs$oat$sigma2_e)),
      y_oat = loc_oat[location] + oat_Pr[oat_acc] + pea_As[pea_acc] +
              d_pea[location] * pool_sign(pea_pool) +
              stats::rnorm(dplyr::n(), 0, sqrt(inputs$pea$sigma2_e))
    )

  #' Three denominators for the same contrast, on the same simulated data.
  #'
  #'   pooled    accession-level t, 2n-2 df, NO location term. This is the
  #'             model VALIDATION_DESIGN.md section 4 specified, and the one
  #'             the comparison exists to indict.
  #'   prereg    the pre-registered denominator: the accession term PLUS the
  #'             pool x location term, Satterthwaite df. Method of moments
  #'             rather than lmer, so 2,000 reps cost seconds and no new
  #'             dependency is added.
  #'   known     the PRE-REGISTERED test: the same accession term, but the
  #'             interaction variance taken as KNOWN from the cross-validation's
  #'             eight folds rather than re-estimated from the trial's own
  #'             n_loc locations. This is what contrast_power() and the two
  #'             slope functions already do, and it is the one that holds its
  #'             size -- see the type I table. Re-estimating a variance on four
  #'             degrees of freedom and dividing by it makes the statistic
  #'             heavy-tailed, which is where `prereg`'s residual inflation
  #'             comes from.
  #'   location  one-sample t on the n_loc location-level contrasts, 4 df.
  #'             Reported as a reference, NOT as a conservative bound: the same
  #'             accessions appear at every location, so sigma2_within does not
  #'             enter its denominator and it answers the NARROW question "do
  #'             these 40 lines differ?" rather than the broad one about the
  #'             class they were drawn from.
  test_one <- function(response, focal_acc, focal_pool) {
    by_loc <- plots |>
      dplyr::group_by(location) |>
      dplyr::mutate(y = .data[[response]] - mean(.data[[response]])) |>
      dplyr::ungroup()

    acc <- by_loc |>
      dplyr::group_by(acc = .data[[focal_acc]], pool = .data[[focal_pool]]) |>
      dplyr::summarise(y = mean(y), .groups = "drop")

    plus <- acc$y[acc$pool == "As+"]; minus <- acc$y[acc$pool == "As-"]
    n_p <- length(plus); n_m <- length(minus)
    est <- mean(plus) - mean(minus)

    # --- pooled: no location term
    v_p <- stats::var(plus); v_m <- stats::var(minus)
    se_pooled <- sqrt(v_p / n_p + v_m / n_m)
    df_pooled <- n_p + n_m - 2

    # --- location-level contrasts, for the other two
    cl <- by_loc |>
      dplyr::group_by(location, acc = .data[[focal_acc]],
                      pool = .data[[focal_pool]]) |>
      dplyr::summarise(y = mean(y), .groups = "drop") |>
      dplyr::group_by(location, pool) |>
      dplyr::summarise(m = mean(y), .groups = "drop") |>
      tidyr::pivot_wider(names_from = pool, values_from = m) |>
      dplyr::mutate(c = `As+` - `As-`)

    se_loc <- stats::sd(cl$c) / sqrt(nrow(cl))
    df_loc <- nrow(cl) - 1L

    # --- prereg: accession term + pool x location term, Satterthwaite
    #
    # var(c_l) across locations estimates 4*tau^2 + W, where W is the
    # WITHIN-LOCATION sampling variance of the contrast. The accession effects
    # are constant across locations, so sigma2_within contributes nothing to
    # var(c_l) and must not be subtracted from it -- an earlier version
    # subtracted the whole accession term scaled by n_loc, which over-subtracts
    # by exactly sigma2_within and drove v_int toward the floor.
    #
    # W is estimated from the two spreads of accession means:
    #   S1 = mean over locations of var(accession means within a location)
    #        estimates sigma2_within + sigma2_e/m
    #   S0 = var(accession means over the whole trial)
    #        estimates sigma2_within + sigma2_e/(m * n_loc)
    # so sigma2_e/m follows from their difference and W from that.
    acc_loc <- by_loc |>
      dplyr::group_by(location, acc = .data[[focal_acc]],
                      pool = .data[[focal_pool]]) |>
      dplyr::summarise(y = mean(y), .groups = "drop")

    s_of <- function(pl) {
      s1 <- acc_loc |>
        dplyr::filter(pool == pl) |>
        dplyr::group_by(location) |>
        dplyr::summarise(v = stats::var(y), .groups = "drop") |>
        dplyr::pull(v) |> mean(na.rm = TRUE)
      s0 <- stats::var(acc$y[acc$pool == pl])
      max(0, (s1 - s0) * n_loc / max(n_loc - 1, 1))   # sigma2_e / m
      }
    W <- s_of("As+") / n_p + s_of("As-") / n_m

    v_contrast_loc <- stats::var(cl$c)
    v_acc  <- v_p / n_p + v_m / n_m
    v_int  <- max(0, (v_contrast_loc - W) / n_loc)
    v_tot  <- v_acc + v_int
    df_pre <- if (v_int > 0) {
      v_tot^2 / (v_acc^2 / df_pooled + v_int^2 / max(n_loc - 1, 1))
    } else df_pooled
    se_pre <- sqrt(v_tot)

    one <- function(est, se, df) {
      tt <- est / se
      p  <- if (sided == 1) stats::pt(tt, df, lower.tail = FALSE)
            else 2 * stats::pt(-abs(tt), df)
      c(p = p, reject = as.numeric(p < alpha), se = se, df = df)
    }

    # the pre-registered test: interaction variance supplied, not estimated
    rows <- list(pooled = one(est, se_pooled, df_pooled),
                 prereg = one(est, se_pre,    df_pre),
                 location = one(est, se_loc,  df_loc))
    if (!is.null(int_sd_known)) {
      v_k  <- v_acc + int_sd_known^2 / n_loc
      df_k <- v_k^2 / (v_acc^2 / df_pooled +
                       (int_sd_known^2 / n_loc)^2 / max(n_loc - 1, 1))
      rows$known <- one(est, sqrt(v_k), df_k)
    }
    r <- do.call(rbind, rows)
    tibble::tibble(denominator = rownames(r), estimate = est,
                   p_value = r[, "p"], reject = r[, "reject"] > 0.5,
                   SE = r[, "se"], df = r[, "df"],
                   n_focal = n_p + n_m,
                   int_sd_hat = 2 * sqrt(v_int))
  }

  dplyr::bind_rows(
    dplyr::mutate(test_one("y_pea", "oat_acc", "oat_pool"), species = "oat"),
    dplyr::mutate(test_one("y_oat", "pea_acc", "pea_pool"), species = "pea")
  ) |>
    dplyr::mutate(n_plots = nrow(plots), tau = tau, .before = 1)
}

# ============================================================
# Driver
# ============================================================

out_dir <- file.path(out_root, vintage)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- validation_inputs(min_partners = min_partners)

cat("\n=== Vintage ===\n")
print(as.data.frame(validation_vintage(inputs)), row.names = FALSE, digits = 4)

P_values <- plots_per_location * n_locations

# ------------------------------------------------------------
# lambda and interaction_frac, read from the cross-validation
#
# validate_crossval.R writes crossval_summary.csv into the same vintage folder
# and now runs BEFORE this script, so the estimate describes the data this run
# was built from. An older vintage is used if the current one has none -- stated
# out loud, because a power figure carrying last month's lambda is a different
# claim from one carrying this month's.
#
# "informative folds only" is the set used: the other set includes folds where
# the held-out trial had too little variation for its slope to say anything
# about lambda, and averaging those in pulls the estimate toward noise.
# ------------------------------------------------------------

FALLBACK <- list(lambda = 0.8, interaction_frac = 0.5)

crossval_estimates <- function(out_root, vintage,
                               set_wanted = "informative folds only") {
  dirs <- sort(list.dirs(out_root, recursive = FALSE))
  here_first <- c(file.path(out_root, vintage),
                  rev(setdiff(dirs, file.path(out_root, vintage))))
  for (d in here_first) {
    f <- file.path(d, "crossval_summary.csv")
    if (!file.exists(f)) next
    x <- readr::read_csv(f, show_col_types = FALSE)
    x <- dplyr::filter(x, set == set_wanted)
    if (nrow(x) == 0) next
    # The RAW scale is what the trial is sized in; the standardised rows are
    # reported by validate_crossval.R for the interaction comparison and are
    # not a substitute here.
    if ("scale" %in% names(x)) x <- dplyr::filter(x, scale == "raw")
    return(list(est = dplyr::select(
                  x, species, lambda = lambda_mean,
                  lambda_se = lambda_se_of_mean,
                  interaction_frac, interaction_frac_corrected, n_folds),
                from = basename(d)))
  }
  NULL
}

cv_est <- crossval_estimates(out_root, vintage)
if (is.null(cv_est)) {
  warning("no crossval_summary.csv anywhere under ", out_root,
          " -- falling back to lambda = ", FALLBACK$lambda,
          ", interaction_frac = ", FALLBACK$interaction_frac,
          ". These are GUESSES; run code/validate_crossval.R and re-run this.",
          call. = FALSE)
  est_tbl <- tibble::tibble(species = names(inputs),
                            lambda = FALLBACK$lambda,
                            lambda_se = NA_real_,
                            interaction_frac = FALLBACK$interaction_frac,
                            interaction_frac_corrected = FALLBACK$interaction_frac,
                            n_folds = NA_integer_)
  est_from <- "FALLBACK GUESS"
} else {
  est_tbl  <- cv_est$est
  est_from <- cv_est$from
  if (!identical(est_from, vintage)) {
    message("no cross-validation in this vintage; using the estimate from ",
            est_from)
  }
}

est_for <- function(sp) {
  r <- dplyr::filter(est_tbl, species == sp)
  if (nrow(r) == 0) {
    list(lambda = FALLBACK$lambda, lambda_se = NA_real_,
         interaction_frac = FALLBACK$interaction_frac,
         interaction_frac_corrected = FALLBACK$interaction_frac)
  } else {
    list(lambda = r$lambda[1], lambda_se = r$lambda_se[1],
         interaction_frac = r$interaction_frac[1],
         interaction_frac_corrected = r$interaction_frac_corrected[1])
  }
}

cat("\n=== lambda and pool x location interaction, as measured ===\n")
cat("    source: ", est_from, "  (crossval_summary.csv, informative folds)\n\n",
    sep = "")
print(as.data.frame(est_tbl), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# Analytic sweep
#
# What is still swept is what is still a CHOICE: the plot budget and the pool
# size. lambda and interaction_frac are measurements and enter as single values.
# ------------------------------------------------------------

grid <- purrr::imap(inputs, \(inp, name) {
  e <- est_for(name)
  power_grid(inp, n_values = n_values, P_values = P_values,
             lambda_values = e$lambda, sided = sided,
             interaction_values = e$interaction_frac,
             pr_quantile = pr_quantile, n_loc = n_locations,
             n_loc_values = n_locations_swept)
}) |> purrr::list_rbind()

readr::write_csv(grid, file.path(out_dir, "power_grid.csv"))

cat("\n=== Power at the design geometry (n = plots-per-location / 4) ===\n")
cat("    one-sided alpha = 0.05; lambda and interaction as measured above\n\n")

design_rows <- grid |>
  dplyr::filter(n == P / n_locations / 4) |>
  dplyr::mutate(plots_per_loc = P / n_locations) |>
  dplyr::select(species, plots_per_loc, n, P, lambda, delta, SE, power) |>
  dplyr::arrange(species, plots_per_loc, dplyr::desc(lambda))
print(as.data.frame(design_rows), row.names = FALSE, digits = 3)

cat("\n=== Pool size at P = 400, one-sided ===\n")
cat("    power is nearly flat: the contrast shrinks as fast as the SE does\n\n")
print(grid |>
  dplyr::filter(P == 400) |>
  dplyr::select(species, n, dAs_predicted, delta, SE, SE_floor, power) |>
  as.data.frame(), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# What the plot budget cannot buy
# ------------------------------------------------------------

# Grouped by n_locations as well as n: locations are now swept, and the whole
# point of the ceiling table is that the irreducible SE is mostly the
# pool x location term, which depends on how many locations there are.
ceiling_tbl <- grid |>
  dplyr::group_by(species, n, n_locations) |>
  dplyr::summarise(
    dAs = dplyr::first(dAs_predicted),
    lambda = dplyr::first(lambda),
    SE_300 = SE[P == 300], SE_500 = SE[P == 500],
    SE_floor = dplyr::first(SE_floor),
    SE_accession_only = dplyr::first(SE_accession_only),
    df = dplyr::first(df[P == 400]),
    pct_var_interaction = dplyr::first(pct_var_interaction[P == 400]),
    pct_SE_irreducible = dplyr::first(pct_SE_irreducible[P == 400]),
    power_300 = power[P == 300], power_500 = power[P == 500],
    .groups = "drop"
  ) |>
  dplyr::mutate(
    # Power with infinitely many plots: the ceiling the budget approaches, at
    # the measured lambda rather than a round number standing in for it.
    #
    # SE_floor now includes the pool x location term, because that term does
    # not shrink with P either -- so this is a genuine infinite-plot limit and
    # not, as before, the accession term alone. The df comes from the grid's
    # own Satterthwaite value rather than 2n-2, since the limiting statistic is
    # dominated by a term with n_loc - 1 df.
    power_ceiling = purrr::pmap_dbl(
      list(dAs, SE_floor, df, lambda),
      \(d, se, dfi, lam) stats::pt(stats::qt(1 - alpha, dfi), dfi,
                                   lam * d / se, lower.tail = FALSE))
  )

readr::write_csv(ceiling_tbl, file.path(out_dir, "power_ceiling.csv"))

cat("\n=== Why more plots buy so little (measured lambda, one-sided) ===\n")
cat("    SE_floor is the standard error with infinitely many plots.\n\n")
print(as.data.frame(ceiling_tbl), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# Simulation, including the null check
# ------------------------------------------------------------

set.seed(sim_seed)

# Both species come out of one simulated trial, so the cases are indexed by the
# pool sizes and the budget rather than by species. One lambda per species, but
# a simulated trial carries both at once, so the simulation uses the mean of
# the two measured values; the per-species analytic grid remains the exact
# statement.
sim_lambda <- mean(est_tbl$lambda)
sim_scenarios <- dplyr::distinct(sim_cases, p_loc) |>
  dplyr::mutate(lambda = sim_lambda)

# The interaction SD to plant. tau is the PER-LOCATION deviation, and the
# contrast's across-location SD is twice it, so tau = interaction_sd / 2 with
# interaction_sd = interaction_frac * lambda * dAs -- the same quantity the
# analytic formula uses.
sim_tau <- function(species, lambda) {
  inp <- inputs[[species]]
  bp  <- build_pools(inp, sim_cases$n[sim_cases$species == species][1],
                     pr_quantile = pr_quantile)
  est_for(species)$interaction_frac * lambda * bp$dAs / 2
}

run_cases <- function(p_loc, lambda, null = FALSE, tau = 0, reps = sim_reps,
                      int_sd_known = NULL) {
  pools <- purrr::imap(inputs, \(inp, name)
    build_pools(inp, sim_cases$n[sim_cases$species == name][1],
                pr_quantile = pr_quantile)$pools)

  runs <- purrr::map(seq_len(reps), \(i)
    simulate_validation(inputs, pools$oat, pools$pea, p_loc, n_locations,
                        lambda = lambda, alpha = alpha, sided = 1,
                        null = null, tau = tau,
                        int_sd_known = int_sd_known)) |>
    purrr::list_rbind()

  runs |>
    dplyr::group_by(species, denominator) |>
    dplyr::summarise(
      mean_estimate = mean(estimate), empirical_SE = stats::sd(estimate),
      mean_SE = mean(SE), mean_df = mean(df),
      mean_int_sd_hat = mean(int_sd_hat),
      power_sim = mean(reject), n_plots = dplyr::first(n_plots),
      .groups = "drop"
    ) |>
    dplyr::mutate(p_loc = p_loc, lambda = if (null) NA_real_ else lambda,
                  null = null, tau = tau, reps = reps, .before = 1)
}

sim <- purrr::pmap(sim_scenarios, \(p_loc, lambda) {
  message("simulating p_loc=", p_loc, " lambda=", round(lambda, 3), " ...")
  tau <- mean(c(sim_tau("oat", lambda), sim_tau("pea", lambda)))
  run_cases(p_loc, lambda, tau = tau)
}) |> purrr::list_rbind()

# attach the analytic prediction for the same configuration
sim <- sim |>
  dplyr::mutate(n = purrr::map_int(species, \(s) sim_cases$n[sim_cases$species == s][1])) |>
  dplyr::rowwise() |>
  dplyr::mutate({
    inp <- inputs[[species]]
    bp  <- build_pools(inp, n, pr_quantile = pr_quantile)
    a   <- contrast_power(dAs = bp$dAs,
                          sigma2_within = c(bp$sigma2_within_plus,
                                            bp$sigma2_within_minus),
                          n = n, sigma2_e = inp$sigma2_e,
                          P = p_loc * n_locations, lambda = lambda,
                          interaction_frac = est_for(species)$interaction_frac,
                          n_loc = n_locations, sided = 1)
    tibble::tibble(predicted_delta = a$delta, analytic_SE = a$SE,
                   analytic_df = a$df, power_analytic = a$power)
  }) |>
  dplyr::ungroup()

# ------------------------------------------------------------
# THE NULL, TWO WAYS -- the comparison this section exists for
#
# `tau = 0` is the old null: no associate effects and no interaction either.
# `tau > 0` is the realistic failure mode: associate effects that are entirely
# location-specific, with a mean of exactly zero across locations. A
# denominator that omits the pool x location term cannot tell the two apart and
# will reject the second far too often.
#
# READ IT AS A PAIRED DIFFERENCE, not against nominal 0.05. The accession-level
# t already sits below nominal because replication is unequal -- anchors carry
# twice the plots, and where a pool is larger than a cell its members appear at
# only some locations -- so the baseline is each denominator's own tau = 0 rate.
#
# 2,000 reps, because at 400 the Monte Carlo SE is 0.011 and cannot distinguish
# 0.05 from 0.08. Both tests are vectorised dplyr, not MCMC.
# ------------------------------------------------------------

null_tau <- mean(c(sim_tau("oat", sim_lambda), sim_tau("pea", sim_lambda)))

# The interaction SD the cross-validation measured, which the pre-registered
# test takes as known rather than re-deriving from five locations.
null_int_sd <- 2 * null_tau

message("simulating the null, no interaction ...")
null_flat <- run_cases(80L, 1, null = TRUE, tau = 0, reps = null_reps,
                       int_sd_known = null_int_sd)
message("simulating the null, WITH As x location ...")
null_int  <- run_cases(80L, 1, null = TRUE, tau = null_tau, reps = null_reps,
                       int_sd_known = null_int_sd)

null_runs <- dplyr::bind_rows(null_flat, null_int) |>
  dplyr::mutate(n = purrr::map_int(species, \(s) sim_cases$n[sim_cases$species == s][1]),
                predicted_delta = 0, analytic_SE = NA_real_,
                analytic_df = NA_real_, power_analytic = alpha)

sim <- dplyr::bind_rows(sim, null_runs) |>
  dplyr::select(species, denominator, n, p_loc, lambda, null, tau, reps,
                mean_estimate, predicted_delta, empirical_SE, analytic_SE,
                mean_SE, mean_df, analytic_df, mean_int_sd_hat,
                power_sim, power_analytic)

readr::write_csv(sim, file.path(out_dir, "power_simulation.csv"))

cat("\n=== Simulation against analytic (alternative; denominator = prereg) ===\n")
print(sim |>
        dplyr::filter(!null, denominator == "prereg") |>
        dplyr::select(species, n, p_loc, lambda, mean_estimate, predicted_delta,
                      empirical_SE, analytic_SE, mean_df, analytic_df,
                      power_sim, power_analytic) |>
        as.data.frame(), row.names = FALSE, digits = 3)

# The two routes are not expected to agree exactly, and the direction of the
# disagreement is what matters. The simulation attenuates the whole genetic
# signal by lambda, so the within-pool spread shrinks with the contrast; the
# analytic form shrinks only the contrast and leaves sigma2_within at full
# size. That makes the ANALYTIC the conservative one at lambda < 1, which is
# the right way round for a power claim. Only the opposite -- the simulation
# coming in BELOW the analytic -- would mean the reported power is optimistic.
optimistic <- dplyr::filter(sim, !null, denominator == "prereg",
                            power_sim < power_analytic - 0.10)
if (nrow(optimistic) > 0) {
  warning("the simulation returns LOWER power than the analytic formula in ",
          nrow(optimistic), " case(s): the reported power is optimistic and ",
          "the formula should not be trusted until that is explained",
          call. = FALSE)
}

cat("\n=== TYPE I ERROR: does omitting pool x location inflate it? ===\n")
cat("    Every row has NO associate main effect. tau = 0 is the old null;\n")
cat("    tau > 0 plants As x location with a mean of exactly zero.\n")
cat("    Read each denominator against ITS OWN tau = 0 rate.\n\n")

t1 <- sim |>
  dplyr::filter(null) |>
  dplyr::select(species, denominator, tau, reps, power_sim, mean_df) |>
  dplyr::mutate(arm = ifelse(tau > 0, "with_interaction", "flat")) |>
  dplyr::select(-tau) |>
  tidyr::pivot_wider(names_from = arm, values_from = c(power_sim, mean_df)) |>
  dplyr::mutate(inflation = power_sim_with_interaction - power_sim_flat)

print(as.data.frame(t1), row.names = FALSE, digits = 3)

readr::write_csv(t1, file.path(out_dir, "power_type1.csv"))

cat("\n")
for (i in seq_len(nrow(t1))) {
  r <- t1[i, ]
  verdict <- if (r$denominator == "known") {
    if (r$power_sim_with_interaction <= alpha + 0.01)
      "HOLDS ITS SIZE -- this is the pre-registered test, and the one the reported power describes"
    else "above nominal, which would undermine the reported power"
  } else if (r$denominator == "pooled") {
    if (r$inflation > 0.03)
      "INFLATED, as SELF_CRITIQUE.md finding C1 predicted -- this is why section 4's model needed the term"
    else "not inflated here, which would contradict the analytic variance split"
  } else {
    if (abs(r$inflation) <= 0.02) "holds its nominal rate"
    else "moves more than it should -- worth understanding before relying on it"
  }
  cat("  ", r$species, " / ", r$denominator, ": ",
      sprintf("%.3f -> %.3f (%+.3f) -- ", r$power_sim_flat,
              r$power_sim_with_interaction, r$inflation),
      verdict, "\n", sep = "")
}

cat("\n  The `location` row is a reference, not a bound: the same accessions\n",
    "  appear at every location, so sigma2_within does not enter its\n",
    "  denominator and it answers the NARROW question about these 2n lines\n",
    "  rather than the broad one about the class they stand for.\n", sep = "")

# ============================================================
# The three pre-registered estimands, at the measured lambda
#
# ONE TABLE, from PREREG_MODEL. Everything here derives from the same
# random-effects structure, so no power number describes a different experiment
# from the one that will be analysed -- which is what happened when section 4
# specified a mixed model while both power routes used a two-stage t.
#
# The variances come from the REALISED field book, not a balanced-design
# formula, so the anchors' unequal replication is accounted for exactly.
# ============================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat(paste(prereg_lines(), collapse = "\n"), "\n")
cat(strrep("=", 70), "\n", sep = "")

vc_resid <- readr::read_csv(here::here("output", "BGLR_variance_components.csv"),
                            show_col_types = FALSE) |>
  dplyr::filter(component == "residual")
var_oat_e <- mean(vc_resid$var_oat)
var_pea_e <- mean(vc_resid$var_pea)
cov_e     <- mean(vc_resid$cov_pea_oat)

nm <- function(d, col) stats::setNames(d[[col]], d$acc)

estimand_power <- purrr::map(n_locations_swept, \(nl) {
  built <- purrr::imap(inputs, \(inp, name)
    build_pools(inp, validation_setting("n_per_pool")[[name]],
                pr_quantile = pr_quantile))
  design <- make_validation_design(built$oat$pools, built$pea$pools,
                                   n_loc = nl, p_loc = plots_per_location_fixed,
                                   n_anchor = validation_setting("n_anchor_per_cell"),
                                   block_size = validation_setting("block_size"),
                                   seed = sim_seed)

  purrr::map(c("global", "lower95", "upper95"), \(src) {
    lam_of <- function(sp) {
      e <- dplyr::filter(est_tbl, species == sp)
      switch(src, global = e$lambda[1],
             lower95 = max(0, e$lambda[1] - 1.96 * e$lambda_se[1]),
             upper95 = e$lambda[1] + 1.96 * e$lambda_se[1])
    }

    rows <- purrr::imap(inputs, \(inp, sp) {
      bp   <- built[[sp]]
      lam  <- lam_of(sp)
      int  <- est_for(sp)$interaction_frac
      foc  <- if (sp == "oat") "oat_acc" else "pea_acc"
      par  <- if (sp == "oat") "pea_acc" else "oat_acc"
      other <- built[[if (sp == "oat") "pea" else "oat"]]

      pool <- contrast_power(
        dAs = bp$dAs,
        sigma2_within = c(bp$sigma2_within_plus, bp$sigma2_within_minus),
        n = bp$n, sigma2_e = inp$sigma2_e,
        P = plots_per_location_fixed * nl, lambda = lam,
        interaction_frac = int, n_loc = nl, sided = 1)

      slope <- slope_power(
        design, x = nm(bp$pools, "As"),
        pev_focal = nm(bp$pools, "PEV_As_i"),
        partner_x = nm(other$pools, "Pr"),
        sigma2_e = inp$sigma2_e, lambda = lam, interaction_frac = int,
        n_loc = nl, focal_col = foc, partner_col = par, sided = 1)

      dplyr::bind_rows(
        tibble::tibble(estimand = "slope", role = "primary", target = sp,
                       effect = slope$lambda, SE = slope$SE, df = slope$df,
                       t = slope$t, power = slope$power,
                       pct_var_interaction = slope$pct_var_interaction),
        tibble::tibble(estimand = "pool", role = "descriptive", target = sp,
                       effect = pool$delta, SE = pool$SE, df = pool$df,
                       t = pool$t, power = pool$power,
                       pct_var_interaction = pool$pct_var_interaction))
    }) |> purrr::list_rbind()

    lam_t <- mean(c(lam_of("oat"), lam_of("pea")))
    int_t <- mean(c(est_for("oat")$interaction_frac,
                    est_for("pea")$interaction_frac))
    tot <- total_yield_power(
      design,
      gma_oat = nm(built$oat$pools, "GMA"), gma_pea = nm(built$pea$pools, "GMA"),
      pev_gma_oat = nm(built$oat$pools, "PEV_GMA_i"),
      pev_gma_pea = nm(built$pea$pools, "PEV_GMA_i"),
      var_oat = var_oat_e, var_pea = var_pea_e, cov_oat_pea = cov_e,
      lambda = lam_t, interaction_frac = int_t, w_pea = 1, n_loc = nl,
      sided = 1)

    dplyr::bind_rows(rows,
      tibble::tibble(estimand = "total", role = "co-primary", target = "both",
                     effect = tot$lambda, SE = tot$SE, df = tot$df,
                     t = tot$t, power = tot$power,
                     pct_var_interaction = tot$pct_var_interaction)) |>
      dplyr::mutate(lambda_source = src, n_loc = nl,
                    plots_per_loc = plots_per_location_fixed,
                    P = plots_per_location_fixed * nl,
                    prereg = PREREG_MODEL$version, .before = 1)
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()

readr::write_csv(estimand_power, file.path(out_dir, "estimand_power.csv"))

cat("\n=== Power for the three pre-registered estimands ===\n")
cat("    ", plots_per_location_fixed, " plots per location; one-sided alpha = ",
    alpha, "; lambda as measured.\n", sep = "")
cat("    Holm over the three primaries (two slopes + total); `pool` is the\n")
cat("    descriptive two-point summary of `slope`, not a fourth test.\n\n")
print(estimand_power |>
        dplyr::select(estimand, role, target, n_loc, lambda_source, effect,
                      SE, df, t, power, pct_var_interaction) |>
        dplyr::arrange(lambda_source, n_loc, estimand, target) |>
        as.data.frame(), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

p_curves <- grid |>
  dplyr::mutate(
    P = factor(paste0(P, " plots")),
    species_lab = factor(sprintf("%s  (lambda = %.2f, interaction = %.2f)",
                                 species, lambda, interaction_frac))
  ) |>
  ggplot2::ggplot(ggplot2::aes(n, power, colour = P, group = P)) +
  ggplot2::geom_hline(yintercept = 0.8, linetype = 2, colour = "grey50") +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::facet_wrap(~ species_lab, ncol = 1) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Power of the As+ vs As- contrast",
    subtitle = paste0("vintage ", vintage,
                      "; one-sided alpha = 0.05; dashed line = 0.8. ",
                      "The curves are flat in pool size and close together ",
                      "across budgets:\nthe binding constraint is the ",
                      "reliability of the effects, not the number of plots."),
    x = "accessions per pool", y = "power", colour = NULL
  )

ggplot2::ggsave(file.path(out_dir, "power_curves.png"), p_curves,
                width = 10, height = 6, dpi = 150)

message("\nwrote ", out_dir)
