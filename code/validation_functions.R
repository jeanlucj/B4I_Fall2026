# ============================================================
# SHARED MACHINERY FOR THE VALIDATION TRIAL
#
# Selecting the As+ / As- pools, and working out what a validation trial of a
# given size could detect.  The design reasoning is in VALIDATION_DESIGN.md.
#
# Nothing here is a constant.  Reliability, PEV, the residual variances and the
# partner counts are all derived at run time from whatever fit is currently in
# output/, because three more trials are arriving and every one of those
# quantities moves when they do.  The numbers quoted in VALIDATION_DESIGN.md
# are a snapshot of one vintage, not parameters.
#
# Two ideas do the work:
#
#   * An associate effect is an effect on the PARTNER's yield, so the oat
#     contrast is read on pea yield and the pea contrast on oat yield.  Every
#     variance below is therefore paired with the partner's residual.
#
#   * The pools are compared as SAMPLES OF ACCESSIONS, not as treatments.  The
#     experimental unit for the contrast is the accession, which is why the
#     standard error has a term that no amount of plot replication touches:
#
#       SE(delta)^2 = 2 * sigma2_within / n  +  4 * sigma2_e / P
#                     \___ accessions ___/     \___ plots ___/
#
#     with sigma2_within = PEV + within-pool variance of the BLUPs.  At the
#     reliabilities this data supports, the first term is 75-80% of the total.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

# ------------------------------------------------------------
# Shared settings
#
# These live here rather than in each driver's settings block, which is the
# convention elsewhere in code/, because the three validate_* scripts must
# agree about them or the vintage is not internally consistent: power computed
# under one partner filter would describe different pools than the ones
# actually selected, and a field book built for a different pool size would
# not match either.  A driver that means to depart from a default overrides it
# explicitly and says so.
# ------------------------------------------------------------

VALIDATION_DEFAULTS <- list(
  # accessions per pool, per species; see VALIDATION_DESIGN.md section 5 for
  # why this is set by the design geometry rather than by the power curve
  n_per_pool = c(oat = 20L, pea = 30L),

  # candidates must clear this quantile of Pr among eligible accessions
  pr_quantile = 0.5,

  # largest acceptable difference in mean Pr between the two pools, g/m2
  pr_tolerance = 1.0,

  # below this many distinct partners an accession's effect is almost pure
  # shrinkage, so there is nothing in it to validate
  min_partners = 3L,

  n_locations = 5L,
  plots_per_location = 80L,

  # combinations per cell fixed across locations AND duplicated within each,
  # which is what separates plot error from combination x location
  n_anchor_per_cell = 1L,

  block_size = 10L
)

#' Read a shared default, allowing a local override.
validation_setting <- function(name, override = NULL) {
  if (!is.null(override)) return(override)
  if (!name %in% names(VALIDATION_DEFAULTS)) {
    stop("no validation default called '", name, "'", call. = FALSE)
  }
  VALIDATION_DEFAULTS[[name]]
}

# ------------------------------------------------------------
# Inputs, re-derived from the current fit
# ------------------------------------------------------------

#' Which variance component is which.
#'
#' In BGLR_multi_trait_model.R the roles are assigned per genetic term: for the
#' oat kernel the effect on oat yield is the producer effect and the effect on
#' pea yield the associate effect; for pea it is reversed.  So an accession's
#' associate variance is always paired with the residual variance of the OTHER
#' species' yield, and that pairing is what the power calculation needs.
VALIDATION_SPECIES <- list(
  oat = list(term = "G_oat", effects_file = "BGLR_oat_effects_all_seeds.csv",
             id_col = "oatAcc", pheno_col = "germplasmName",
             partner_pheno_col = "intercropGermplasmName",
             response = "pea yield", residual_col = "var_pea"),
  pea = list(term = "G_pea", effects_file = "BGLR_pea_effects_all_seeds.csv",
             id_col = "peaAcc", pheno_col = "intercropGermplasmName",
             partner_pheno_col = "germplasmName",
             response = "oat yield", residual_col = "var_oat")
)

#' Everything the pool selection and the power calculation read.
#'
#' @param out_dir Where the fit's outputs live.
#' @param pheno_file Plot table, for the partner counts.
#' @param min_partners Accessions seen with fewer distinct partners than this
#'   are dropped from the candidate set: their effects are almost pure
#'   shrinkage and there is nothing to validate.
#' @return A list with one entry per species, each carrying the accession table
#'   (BLUPs, partner counts, plots), the variance components, the derived
#'   reliability and PEV, and the partner-yield residual variance.
validation_inputs <- function(out_dir = here::here("output"),
                              pheno_file = here::here("output", "B4I_intercrop_pheno.rds"),
                              min_partners = 3L) {

  vc <- readr::read_csv(file.path(out_dir, "BGLR_variance_components.csv"),
                        show_col_types = FALSE)
  genetic  <- dplyr::filter(vc, component == "genetic")
  residual <- dplyr::filter(vc, component == "residual")

  pheno <- readRDS(pheno_file) |>
    dplyr::filter(!is.na(oat_yield), !is.na(pea_yield))

  purrr::imap(VALIDATION_SPECIES, \(sp, name) {
    gg <- dplyr::filter(genetic, term == sp$term)
    sigma2_Pr <- mean(gg$var_Pr)
    sigma2_As <- mean(gg$var_As)
    sigma2_e  <- mean(residual[[sp$residual_col]], na.rm = TRUE)

    # distinct partners and plots per accession, from the plot table
    counts <- pheno |>
      dplyr::distinct(.data[[sp$pheno_col]], .data[[sp$partner_pheno_col]]) |>
      dplyr::count(.data[[sp$pheno_col]], name = "n_partners") |>
      dplyr::rename(acc = 1) |>
      dplyr::left_join(
        dplyr::count(pheno, acc = .data[[sp$pheno_col]], name = "n_plots"),
        by = "acc"
      )

    # the BLUP is the mean over MCMC chains; averaging removes chain noise
    acc <- readr::read_csv(file.path(out_dir, sp$effects_file),
                           show_col_types = FALSE) |>
      dplyr::group_by(acc = .data[[sp$id_col]]) |>
      dplyr::summarise(Pr = mean(PrEff), As = mean(AsEff), .groups = "drop") |>
      dplyr::left_join(counts, by = "acc") |>
      dplyr::mutate(
        n_partners = tidyr::replace_na(n_partners, 0L),
        n_plots    = tidyr::replace_na(n_plots, 0L),
        eligible   = n_partners >= min_partners
      )

    # Var(true) = Var(BLUP) + E[PEV] for a conditionally unbiased predictor,
    # so the reliability of the BLUPs and the PEV both follow from the fit.
    # Clamped: a var(BLUP) above the component means the chains have not
    # settled, and a negative PEV would silently flatter every power estimate.
    rel <- function(v_blup, sigma2) min(max(v_blup / sigma2, 1e-6), 0.999)

    list(
      species      = name,
      accessions   = acc,
      n_eligible   = sum(acc$eligible),
      min_partners = min_partners,
      sigma2_Pr    = sigma2_Pr,
      sigma2_As    = sigma2_As,
      sigma2_e     = sigma2_e,
      response     = sp$response,
      rel_Pr       = rel(stats::var(acc$Pr), sigma2_Pr),
      rel_As       = rel(stats::var(acc$As), sigma2_As),
      PEV_As       = sigma2_As * (1 - rel(stats::var(acc$As), sigma2_As)),
      n_trials     = dplyr::n_distinct(pheno$studyName),
      n_plots      = nrow(pheno)
    )
  })
}

# ------------------------------------------------------------
# Pool construction
#
# Maximise mean(As+) - mean(As-) subject to: the pools disjoint, both pool
# means of Pr above a threshold, and the two pool means of Pr equal to within
# a tolerance.
#
# Selecting the top n on (Pr + As) and the top n on (Pr - As), as first
# proposed, does NOT give disjoint pools: a high-Pr accession makes both lists
# whatever its As.  In this data that was 13 of 30 shared at n = 30.
#
# The constraint is handled with a Lagrange multiplier `theta` on Pr, pushed in
# opposite directions for the two pools and bisected until the Pr difference
# crosses zero.  Disjointness is enforced by building the As+ pool first and
# excluding its members from the As- candidate set.
# ------------------------------------------------------------

#' One (theta) evaluation: the two pools and the differences they imply.
.pools_at <- function(cand, n, theta) {
  plus <- cand |>
    dplyr::mutate(score = As + theta * Pr) |>
    dplyr::slice_max(score, n = n, with_ties = FALSE)

  minus <- cand |>
    dplyr::filter(!acc %in% plus$acc) |>
    dplyr::mutate(score = As - theta * Pr) |>
    dplyr::slice_min(score, n = n, with_ties = FALSE)

  list(plus = plus, minus = minus,
       dAs = mean(plus$As) - mean(minus$As),
       dPr = mean(plus$Pr) - mean(minus$Pr))
}

#' Build the As+ / As- pools for one species.
#'
#' @param inp One element of validation_inputs().
#' @param n Accessions PER POOL (so 2n of that species enter the trial).
#' @param pr_quantile Candidates must be above this quantile of Pr among the
#'   eligible accessions.  0.5 keeps the better half: "all of these are good
#'   producers" is the claim the design has to support.
#' @param tol Largest acceptable |mean Pr difference| between pools, g/m2.
#' @param theta_max Upper end of the bisection on the Pr multiplier.
build_pools <- function(inp, n, pr_quantile = 0.5, tol = 1.0,
                        theta_max = 10, iterations = 60L) {
  cand <- dplyr::filter(inp$accessions, eligible)
  pr_min <- stats::quantile(cand$Pr, pr_quantile, names = FALSE)
  cand <- dplyr::filter(cand, Pr >= pr_min)

  if (nrow(cand) < 2 * n) {
    stop(inp$species, ": ", nrow(cand), " candidates above the Pr threshold ",
         "but ", 2 * n, " needed. Lower pr_quantile or min_partners, or use a ",
         "smaller pool.", call. = FALSE)
  }

  # theta = 0 selects on As alone; because Pr and As are negatively correlated
  # in the BLUPs that leaves the As+ pool short on Pr, so dPr starts negative
  # and rises with theta. Bisect on that.
  lo <- 0; hi <- theta_max; theta <- 0
  for (i in seq_len(iterations)) {
    theta <- (lo + hi) / 2
    if (.pools_at(cand, n, theta)$dPr < 0) lo <- theta else hi <- theta
  }
  res <- .pools_at(cand, n, theta)

  if (abs(res$dPr) > tol) {
    warning(inp$species, ": producer means differ by ", round(res$dPr, 2),
            " g/m2, outside the tolerance of ", tol,
            ". The Pr constraint could not be met at n = ", n, ".",
            call. = FALSE)
  }

  pools <- dplyr::bind_rows(
    dplyr::mutate(res$plus,  pool = "As+"),
    dplyr::mutate(res$minus, pool = "As-")
  ) |>
    dplyr::transmute(species = inp$species, pool, acc, Pr, As,
                     n_partners, n_plots) |>
    dplyr::arrange(pool, dplyr::desc(As))

  # within-pool variance of the TRUE effects: what the BLUPs still spread by
  # after selection, plus what we do not know about each one
  var_blup_within <- mean(c(stats::var(res$plus$As), stats::var(res$minus$As)))

  list(
    pools = pools, n = n, theta = theta, pr_min = pr_min,
    n_candidates = nrow(cand),
    dAs = res$dAs, dPr = res$dPr,
    mean_Pr_plus = mean(res$plus$Pr), mean_Pr_minus = mean(res$minus$Pr),
    mean_As_plus = mean(res$plus$As), mean_As_minus = mean(res$minus$As),
    var_blup_within = var_blup_within,
    sigma2_within = inp$PEV_As + var_blup_within
  )
}

#' What changed between two vintages of the pools.
#'
#' Run at every refresh.  Churn is itself a result: if a few new trials
#' reshuffle most of a pool, the effects are not stable enough to be worth
#' validating, and that is better known before seed is ordered.
pool_diff <- function(old_pools, new_pools) {
  o <- dplyr::select(old_pools, species, acc, pool_old = pool)
  n <- dplyr::select(new_pools, species, acc, pool_new = pool)

  full <- dplyr::full_join(o, n, by = c("species", "acc")) |>
    dplyr::mutate(status = dplyr::case_when(
      is.na(pool_old)          ~ "entered",
      is.na(pool_new)          ~ "left",
      pool_old != pool_new     ~ "switched pool",
      TRUE                     ~ "retained"
    ))

  summary <- full |>
    dplyr::count(species, status, name = "n") |>
    tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0L)

  list(detail = dplyr::arrange(full, species, status, acc), summary = summary)
}

# ------------------------------------------------------------
# The field design
#
# A 2x2 factorial of oat pool x pea pool, equally represented at every
# location.  Both species are validated from the same plots, because an
# associate effect is read on the PARTNER's yield: the oat contrast comes off
# pea yield and the pea contrast off oat yield, so no plot is spent on only
# one of them.
#
# Within a cell at a location, the oats and peas are paired by a permutation,
# and the permutation is redrawn at each location so an accession meets a
# different partner every time.  Averaging an accession over many partners is
# what makes its pool mean precise; it is worth more than repeating
# combinations.
#
# Balance is the load-bearing part.  Because each focal accession meets both
# partner pools equally often, partner effects are orthogonal to the focal
# contrast and drop out of it.  Without that, partner variance inflates the
# contrast's standard error -- and the analytic power formula would be wrong.
# ------------------------------------------------------------

#' Lay out the validation trial.
#'
#' @param oat_pools,pea_pools Pool tables from build_pools()$pools.
#' @param n_loc Locations.
#' @param p_loc Plots per location. Divided equally among the four cells.
#' @param n_anchor Combinations per cell that are fixed across locations AND
#'   duplicated within a location. Duplication within a location is what
#'   separates plot error from combination x location; repeating a combination
#'   only across locations confounds the two.
#' @param block_size Target plots per incomplete block within a location.
#' @return A tibble, one row per plot.
make_validation_design <- function(oat_pools, pea_pools, n_loc = 5L,
                                   p_loc = 80L, n_anchor = 1L,
                                   block_size = 10L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # The anchors' second copies come OUT of the budget, not on top of it, so
  # p_loc is what actually goes in the ground.
  cell_size <- (p_loc - 4L * n_anchor) %/% 4L
  if (cell_size < 4L) {
    stop("p_loc = ", p_loc, " with ", n_anchor, " anchor(s) per cell leaves ",
         "only ", cell_size, " plots per cell; too few to pair meaningfully",
         call. = FALSE)
  }
  if (n_anchor > cell_size) {
    stop("n_anchor = ", n_anchor, " exceeds the cell size of ", cell_size,
         call. = FALSE)
  }

  split_pool <- function(d) split(d$acc, d$pool)
  oats <- split_pool(oat_pools)
  peas <- split_pool(pea_pools)

  # Take `m` members of a pool starting at an offset, wrapping around. When the
  # pool is bigger than the cell, the offset rotates membership across
  # locations so every accession appears about equally often overall.
  take <- function(x, m, offset) x[((offset + seq_len(m) - 1L) %% length(x)) + 1L]

  cells <- tidyr::expand_grid(oat_pool = c("As+", "As-"),
                              pea_pool = c("As+", "As-"))

  design <- purrr::map(seq_len(n_loc), \(loc) {
    purrr::pmap(cells, \(oat_pool, pea_pool) {
      o_all <- oats[[oat_pool]]; p_all <- peas[[pea_pool]]
      o <- take(o_all, cell_size, (loc - 1L) * cell_size)
      p <- take(p_all, cell_size, (loc - 1L) * cell_size)

      # anchors first, identical at every location; the rest re-permuted
      anchor_idx <- seq_len(min(n_anchor, cell_size))
      rest <- setdiff(seq_len(cell_size), anchor_idx)

      pairing <- integer(cell_size)
      pairing[anchor_idx] <- anchor_idx                 # fixed pairing
      pairing[rest] <- if (length(rest) > 1) sample(rest) else rest

      tibble::tibble(
        location = loc, oat_pool = oat_pool, pea_pool = pea_pool,
        oat_acc = o, pea_acc = p[pairing],
        is_anchor = seq_len(cell_size) %in% anchor_idx
      )
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind()

  # Anchors are sown twice at each location -- the within-location duplication
  # that a clean plot-error estimate needs.
  #
  # The anchor combination is deliberately the SAME at every location, which is
  # what makes combination x location estimable as well. The price is that the
  # anchor accessions appear in twice as many plots as everyone else. That is
  # benign for the contrast: the analysis averages to one value per accession
  # before comparing pools, so an anchor contributes a single, slightly more
  # precise mean rather than extra weight.
  design <- dplyr::bind_rows(design, dplyr::filter(design, is_anchor)) |>
    dplyr::mutate(combination = paste(oat_acc, pea_acc, sep = "::"))

  # Randomise to plots, in resolvable incomplete blocks that each carry the
  # four cells in equal proportion.
  design |>
    dplyr::group_by(location) |>
    dplyr::mutate(.r = sample(dplyr::n())) |>
    dplyr::arrange(location, .r, .by_group = TRUE) |>
    dplyr::mutate(
      plot  = dplyr::row_number(),
      block = ((plot - 1L) %/% block_size) + 1L
    ) |>
    dplyr::ungroup() |>
    dplyr::select(location, block, plot, oat_pool, pea_pool,
                  oat_acc, pea_acc, combination, is_anchor)
}

#' Does the design have the balance the analysis assumes?
check_validation_design <- function(design, oat_pools, pea_pools) {
  bal <- function(d, focal, partner) {
    d |>
      dplyr::count(.data[[focal]], .data[[partner]], name = "n") |>
      tidyr::pivot_wider(names_from = dplyr::all_of(partner), values_from = n,
                         values_fill = 0L) |>
      dplyr::rename(acc = 1) |>
      dplyr::mutate(imbalance = abs(`As+` - `As-`))
  }

  oat_bal <- bal(design, "oat_acc", "pea_pool")
  pea_bal <- bal(design, "pea_acc", "oat_pool")

  list(
    plots_per_location = dplyr::count(design, location, name = "plots"),
    cells_per_location = dplyr::count(design, location, oat_pool, pea_pool,
                                      name = "plots"),
    oat_reps = dplyr::count(design, oat_acc, name = "plots"),
    pea_reps = dplyr::count(design, pea_acc, name = "plots"),
    oat_partner_balance = oat_bal,
    pea_partner_balance = pea_bal,
    max_oat_imbalance = max(oat_bal$imbalance),
    max_pea_imbalance = max(pea_bal$imbalance),
    partners_per_oat = design |> dplyr::distinct(oat_acc, pea_acc) |>
      dplyr::count(oat_acc, name = "partners"),
    replicated_combinations = design |> dplyr::count(combination, name = "n") |>
      dplyr::filter(n > 1),
    anchor_plots = sum(design$is_anchor)
  )
}

# ------------------------------------------------------------
# Power
#
# The pools are fixed sets but the accessions within them are a random sample
# of "accessions predicted to be in this class", so the accession is the
# experimental unit and enters the denominator.  df = 2n - 2.
# ------------------------------------------------------------

#' Standard error of the pool contrast.
#'
#' @param sigma2_within PEV + within-pool BLUP variance (true-effect spread
#'   among the accessions in a pool).
#' @param n Accessions per pool.
#' @param sigma2_e Plot residual variance of the PARTNER's yield.
#' @param P Total plots across all locations.
#' @param interaction_frac SD of the location-specific contrast, as a fraction
#'   of the contrast itself.  0 assumes the As effects behave the same
#'   everywhere.
#' @param delta The contrast, needed only to scale interaction_frac.
#' @param n_loc Locations.
contrast_se <- function(sigma2_within, n, sigma2_e, P,
                        interaction_frac = 0, delta = 0, n_loc = 5) {
  sqrt(2 * sigma2_within / n +
       4 * sigma2_e / P +
       (interaction_frac * delta)^2 / n_loc)
}

#' Power of the pool contrast.
#'
#' @param dAs The predicted contrast, from build_pools().
#' @param lambda Attenuation of that contrast when realised in new
#'   environments: model calibration times across-environment stability.  1 is
#'   the optimistic case in which the BLUP difference is delivered in full;
#'   estimate it with validate_crossval.R rather than assuming it.
#'
#'   Note what lambda is doing here and what it is not.  It scales the
#'   CONTRAST but leaves sigma2_within at full size, i.e. it assumes the
#'   signal shrinks while the spread among accessions within a pool does not.
#'   The alternative reading -- that the whole genetic signal attenuates, so
#'   the within-pool spread shrinks by lambda too -- gives a larger t and more
#'   power (the simulation in validate_power.R does it that way, which is why
#'   it returns 2-8 points more).  Attenuation caused by genotype x
#'   environment interaction adds variance rather than removing it, so the
#'   conservative reading used here is the more defensible one for a power
#'   claim.
#' @param sided 1 for the pre-registered directional test, 2 otherwise.
contrast_power <- function(dAs, sigma2_within, n, sigma2_e, P,
                           lambda = 1, alpha = 0.05, sided = 1,
                           interaction_frac = 0, n_loc = 5) {
  delta <- lambda * dAs
  se    <- contrast_se(sigma2_within, n, sigma2_e, P, interaction_frac,
                       delta, n_loc)
  df    <- 2 * n - 2
  ncp   <- delta / se

  power <- if (sided == 1) {
    stats::pt(stats::qt(1 - alpha, df), df, ncp, lower.tail = FALSE)
  } else {
    stats::pt(-stats::qt(1 - alpha / 2, df), df, ncp) +
      stats::pt(stats::qt(1 - alpha / 2, df), df, ncp, lower.tail = FALSE)
  }

  tibble::tibble(
    n = n, P = P, lambda = lambda, sided = sided,
    interaction_frac = interaction_frac,
    delta = delta, SE = se, df = df, t = ncp, power = power,
    # What the plot budget cannot buy: the SE that would remain with infinitely
    # many plots.  Reported two ways because they answer different questions --
    # the variance share says how the budget splits, the SE ratio says how much
    # of the standard error you are stuck with.
    SE_floor = sqrt(2 * sigma2_within / n),
    pct_var_from_accessions = 100 * (2 * sigma2_within / n) / se^2,
    pct_SE_irreducible = 100 * sqrt(2 * sigma2_within / n) / se
  )
}

#' Power over a grid of pool sizes, plot budgets and assumptions.
#'
#' Pools are rebuilt at every n, because the contrast shrinks as the pool grows
#' -- that trade-off is the whole question and must not be held fixed.
power_grid <- function(inp, n_values, P_values, lambda_values = c(1, 0.8, 0.6),
                       sided = 1, interaction_values = 0, pr_quantile = 0.5,
                       n_loc = 5) {
  purrr::map(n_values, \(n) {
    bp <- build_pools(inp, n, pr_quantile = pr_quantile)
    tidyr::expand_grid(P = P_values, lambda = lambda_values,
                       interaction_frac = interaction_values) |>
      purrr::pmap(\(P, lambda, interaction_frac) contrast_power(
        dAs = bp$dAs, sigma2_within = bp$sigma2_within, n = n,
        sigma2_e = inp$sigma2_e, P = P, lambda = lambda, sided = sided,
        interaction_frac = interaction_frac, n_loc = n_loc)) |>
      purrr::list_rbind() |>
      dplyr::mutate(species = inp$species, dAs_predicted = bp$dAs,
                    dPr = bp$dPr, sigma2_within = bp$sigma2_within,
                    .before = 1)
  }) |>
    purrr::list_rbind()
}

# ------------------------------------------------------------
# Reporting
# ------------------------------------------------------------

#' One-line summary of a vintage, for the top of every report.
validation_vintage <- function(inputs) {
  purrr::map(inputs, \(i) tibble::tibble(
    species = i$species, n_trials = i$n_trials, n_plots = i$n_plots,
    n_accessions = nrow(i$accessions), n_eligible = i$n_eligible,
    sigma2_As = i$sigma2_As, reliability_As = i$rel_As, PEV_As = i$PEV_As,
    response = i$response, sigma2_e = i$sigma2_e
  )) |>
    purrr::list_rbind()
}
