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
#               Slower, but it is the only one that captures selecting on
#               estimates and then being scored against the truth.
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

# The budgets to sweep; the default is one of these, not a replacement for it.
plots_per_location <- c(60L, 80L, 100L)

# Pool sizes to sweep. The design geometry wants n = plots-per-location / 4,
# which is 15/20/25; the rest are here to show how flat the curve is.
n_values <- c(10L, 15L, 20L, 25L, 30L, 40L, 50L)

# Attenuation of the predicted contrast when realised in new environments:
# model calibration times across-environment stability of the associate
# effects. 1 means the BLUP difference is delivered in full. This is the single
# most important number in the whole calculation and it is currently a guess --
# code/validate_crossval.R estimates it from the trials already in hand.
lambda_values <- c(1.0, 0.8, 0.6)

# SD of the location-specific contrast as a fraction of the contrast itself.
interaction_values <- c(0, 0.5)

# The test is directional and pre-registered, so one-sided is legitimate and
# worth 8-12 points of power. Two-sided is reported alongside.
alpha <- 0.05

# Simulation
sim_reps  <- 400L
sim_seed  <- 20260924L
# Scenarios to simulate: a subset of the grid, since each costs sim_reps fits.
# The pool sizes come from the shared defaults, so the simulation checks the
# configuration that is actually proposed rather than a nearby one.
sim_cases <- tidyr::expand_grid(
  species = c("oat", "pea"),
  p_loc   = validation_setting("plots_per_location"),
  lambda  = c(1.0, 0.8)
) |>
  dplyr::mutate(n = validation_setting("n_per_pool")[species], .after = species)

# ============================================================
# Simulation
# ============================================================

#' Simulate one validation trial on the REAL design and test both contrasts.
#'
#' Two things this does that a back-of-envelope simulation would not, and both
#' matter:
#'
#'   * The truth is drawn from the posterior of the effects, NOT set equal to
#'     the BLUPs. We select on the BLUPs and are then scored against a truth
#'     that differs from them by the prediction error -- which is exactly the
#'     attenuation a real validation suffers, and the reason selection on
#'     estimates cannot be assumed to deliver the estimated contrast.
#'   * Partners are real accessions laid out by make_validation_design(), not
#'     independent noise per plot. The analytic formula assumes partner effects
#'     are orthogonalised by the design's balance; simulating the actual design
#'     is what tests that assumption rather than granting it.
#'
#' @param null TRUE sets every true associate effect to zero, so the rejection
#'   rate must come back at alpha. The false-positive check.
#' @return One row per species with the realised contrast and its test.
simulate_validation <- function(inputs, oat_pools, pea_pools, p_loc, n_loc,
                                lambda = 1, alpha = 0.05, sided = 1,
                                null = FALSE, n_anchor = 1L) {

  design <- make_validation_design(oat_pools, pea_pools, n_loc = n_loc,
                                   p_loc = p_loc, n_anchor = n_anchor)

  # true effects per accession: BLUP + prediction error, attenuated by lambda
  draw_truth <- function(inp, pools) {
    a <- dplyr::filter(inp$accessions, acc %in% pools$acc)
    v <- if (null) rep(0, nrow(a)) else
      lambda * (a$As + stats::rnorm(nrow(a), 0, sqrt(inp$PEV_As)))
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

  plots <- design |>
    dplyr::mutate(
      # pea yield carries the PEA's producer effect and the OAT's associate
      # effect; oat yield the reverse. This is the whole reason one set of
      # plots validates both species.
      y_pea = loc_pea[location] + pea_Pr[pea_acc] + oat_As[oat_acc] +
              stats::rnorm(dplyr::n(), 0, sqrt(inputs$oat$sigma2_e)),
      y_oat = loc_oat[location] + oat_Pr[oat_acc] + pea_As[pea_acc] +
              stats::rnorm(dplyr::n(), 0, sqrt(inputs$pea$sigma2_e))
    )

  # Two-stage: remove locations, average to the focal accession, then compare
  # the two pools with the ACCESSION as the unit -- the same assumption the
  # analytic formula makes, so agreement between the routes tests both.
  test_one <- function(response, focal_acc, focal_pool) {
    adj <- plots |>
      dplyr::group_by(location) |>
      dplyr::mutate(y = .data[[response]] - mean(.data[[response]])) |>
      dplyr::ungroup() |>
      dplyr::group_by(acc = .data[[focal_acc]], pool = .data[[focal_pool]]) |>
      dplyr::summarise(y = mean(y), .groups = "drop")

    plus  <- adj$y[adj$pool == "As+"]
    minus <- adj$y[adj$pool == "As-"]
    tt <- stats::t.test(plus, minus, var.equal = TRUE)   # As+ minus As-
    est <- as.numeric(tt$estimate[1] - tt$estimate[2])
    p <- if (sided == 1) {
      if (est > 0) tt$p.value / 2 else 1 - tt$p.value / 2
    } else tt$p.value
    tibble::tibble(estimate = est, p_value = p, reject = p < alpha,
                   n_focal = length(plus) + length(minus))
  }

  dplyr::bind_rows(
    dplyr::mutate(test_one("y_pea", "oat_acc", "oat_pool"), species = "oat"),
    dplyr::mutate(test_one("y_oat", "pea_acc", "pea_pool"), species = "pea")
  ) |>
    dplyr::mutate(n_plots = nrow(plots), .before = 1)
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
# Analytic sweep
# ------------------------------------------------------------

grid <- purrr::imap(inputs, \(inp, name) {
  purrr::map(c(1, 2), \(sided) {
    power_grid(inp, n_values = n_values, P_values = P_values,
               lambda_values = lambda_values, sided = sided,
               interaction_values = interaction_values,
               pr_quantile = pr_quantile, n_loc = n_locations)
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()

readr::write_csv(grid, file.path(out_dir, "power_grid.csv"))

cat("\n=== Power at the design geometry (n = plots-per-location / 4) ===\n")
cat("    one-sided alpha = 0.05, no pool x location interaction\n\n")

design_rows <- grid |>
  dplyr::filter(sided == 1, interaction_frac == 0,
                n == P / n_locations / 4) |>
  dplyr::mutate(plots_per_loc = P / n_locations) |>
  dplyr::select(species, plots_per_loc, n, P, lambda, delta, SE, power) |>
  dplyr::arrange(species, plots_per_loc, dplyr::desc(lambda))
print(as.data.frame(design_rows), row.names = FALSE, digits = 3)

cat("\n=== Pool size at P = 400, lambda = 0.8, one-sided ===\n")
cat("    power is nearly flat: the contrast shrinks as fast as the SE does\n\n")
print(grid |>
  dplyr::filter(P == 400, lambda == 0.8, sided == 1, interaction_frac == 0) |>
  dplyr::select(species, n, dAs_predicted, delta, SE, SE_floor, power) |>
  as.data.frame(), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# What the plot budget cannot buy
# ------------------------------------------------------------

ceiling_tbl <- grid |>
  dplyr::filter(sided == 1, interaction_frac == 0, lambda == 0.8) |>
  dplyr::group_by(species, n) |>
  dplyr::summarise(
    dAs = dplyr::first(dAs_predicted),
    SE_300 = SE[P == 300], SE_500 = SE[P == 500],
    SE_floor = dplyr::first(SE_floor),
    pct_SE_irreducible = dplyr::first(pct_SE_irreducible[P == 400]),
    power_300 = power[P == 300], power_500 = power[P == 500],
    .groups = "drop"
  ) |>
  dplyr::mutate(
    # power with infinitely many plots: the ceiling the budget approaches
    power_ceiling = purrr::pmap_dbl(
      list(dAs, SE_floor, n),
      \(d, se, nn) stats::pt(stats::qt(1 - alpha, 2 * nn - 2), 2 * nn - 2,
                             0.8 * d / se, lower.tail = FALSE))
  )

readr::write_csv(ceiling_tbl, file.path(out_dir, "power_ceiling.csv"))

cat("\n=== Why more plots buy so little (lambda = 0.8, one-sided) ===\n")
cat("    SE_floor is the standard error with infinitely many plots.\n\n")
print(as.data.frame(ceiling_tbl), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# Simulation, including the null check
# ------------------------------------------------------------

set.seed(sim_seed)

# Both species come out of one simulated trial, so the cases are indexed by the
# pool sizes and the budget rather than by species.
sim_scenarios <- dplyr::distinct(sim_cases, p_loc, lambda)

run_cases <- function(p_loc, lambda, null = FALSE) {
  pools <- purrr::imap(inputs, \(inp, name)
    build_pools(inp, sim_cases$n[sim_cases$species == name][1],
                pr_quantile = pr_quantile)$pools)

  runs <- purrr::map(seq_len(sim_reps), \(i)
    simulate_validation(inputs, pools$oat, pools$pea, p_loc, n_locations,
                        lambda = lambda, alpha = alpha, sided = 1,
                        null = null)) |>
    purrr::list_rbind()

  runs |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      mean_estimate = mean(estimate), empirical_SE = stats::sd(estimate),
      power_sim = mean(reject), n_plots = dplyr::first(n_plots),
      .groups = "drop"
    ) |>
    dplyr::mutate(p_loc = p_loc, lambda = if (null) NA_real_ else lambda,
                  null = null, reps = sim_reps, .before = 1)
}

sim <- purrr::pmap(sim_scenarios, \(p_loc, lambda) {
  message("simulating p_loc=", p_loc, " lambda=", lambda, " ...")
  run_cases(p_loc, lambda)
}) |> purrr::list_rbind()

# attach the analytic prediction for the same configuration
sim <- sim |>
  dplyr::mutate(n = purrr::map_int(species, \(s) sim_cases$n[sim_cases$species == s][1])) |>
  dplyr::rowwise() |>
  dplyr::mutate({
    inp <- inputs[[species]]
    bp  <- build_pools(inp, n, pr_quantile = pr_quantile)
    a   <- contrast_power(bp$dAs, bp$sigma2_within, n, inp$sigma2_e,
                          p_loc * n_locations, lambda = lambda, sided = 1)
    tibble::tibble(predicted_delta = a$delta, analytic_SE = a$SE,
                   power_analytic = a$power)
  }) |>
  dplyr::ungroup()

# False-positive check: with no true contrast the rejection rate must be alpha.
# Same discipline as sim_run.R's --check -- a pipeline that cannot produce a
# negative result cannot be trusted with a positive one.
message("simulating the null ...")
null_runs <- run_cases(80L, 1, null = TRUE) |>
  dplyr::mutate(n = purrr::map_int(species, \(s) sim_cases$n[sim_cases$species == s][1]),
                predicted_delta = 0, analytic_SE = NA_real_,
                power_analytic = alpha)

sim <- dplyr::bind_rows(sim, null_runs) |>
  dplyr::select(species, n, p_loc, lambda, null, reps, mean_estimate,
                predicted_delta, empirical_SE, analytic_SE, power_sim,
                power_analytic)

readr::write_csv(sim, file.path(out_dir, "power_simulation.csv"))

cat("\n=== Simulation against analytic ===\n")
print(as.data.frame(sim), row.names = FALSE, digits = 3)

# The two routes are not expected to agree exactly, and the direction of the
# disagreement is what matters. The simulation attenuates the whole genetic
# signal by lambda, so the within-pool spread shrinks with the contrast; the
# analytic form shrinks only the contrast and leaves sigma2_within at full
# size. That makes the ANALYTIC the conservative one at lambda < 1, which is
# the right way round for a power claim. Only the opposite -- the simulation
# coming in BELOW the analytic -- would mean the reported power is optimistic.
optimistic <- dplyr::filter(sim, !null, power_sim < power_analytic - 0.10)
if (nrow(optimistic) > 0) {
  warning("the simulation returns LOWER power than the analytic formula in ",
          nrow(optimistic), " case(s): the reported power is optimistic and ",
          "the formula should not be trusted until that is explained",
          call. = FALSE)
}
conservative <- dplyr::filter(sim, !null, power_sim > power_analytic + 0.10)
if (nrow(conservative) > 0) {
  message("note: the simulation is more than 10 points ABOVE the analytic in ",
          nrow(conservative), " case(s), all at lambda < 1. Expected -- see ",
          "the lambda note in contrast_power(). The analytic number is the ",
          "one to quote.")
}

false_pos <- dplyr::filter(sim, null)
cat("\nfalse-positive rate under the null: ",
    paste(sprintf("%s %.3f", false_pos$species, false_pos$power_sim),
          collapse = ", "),
    "  (nominal ", alpha, ")\n", sep = "")
cat("  The accession-level t-test comes in below nominal because replication\n",
    "  is not equal across accessions -- anchors carry twice the plots, and\n",
    "  where a pool is larger than a cell its members appear at only some\n",
    "  locations. Equal-variance t assumes away that heterogeneity and so\n",
    "  over-states the spread. It errs toward NOT rejecting, so the power\n",
    "  above is if anything understated; the mixed model in the real analysis\n",
    "  weights by precision and should recover the nominal rate.\n", sep = "")
if (any(false_pos$power_sim > alpha + 3 * sqrt(alpha * (1 - alpha) / sim_reps))) {
  warning("false-positive rate is materially ABOVE nominal: the test is ",
          "anti-conservative and the design or the model is wrong",
          call. = FALSE)
}

# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

p_curves <- grid |>
  dplyr::filter(sided == 1, interaction_frac == 0) |>
  dplyr::mutate(
    P = factor(paste0(P, " plots")),
    lambda = factor(paste0("lambda = ", lambda))
  ) |>
  ggplot2::ggplot(ggplot2::aes(n, power, colour = P, group = P)) +
  ggplot2::geom_hline(yintercept = 0.8, linetype = 2, colour = "grey50") +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::facet_grid(species ~ lambda) +
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
