# test_validation.R
#
# The validation-trial machinery in code/validation_functions.R: pool
# construction, the field design's balance, and the power formula.
#
# These run against a synthetic accession table, not against output/, so they
# work in a fresh clone with no fitted model present. The design's balance is the
# property the whole contrast rests on -- if a focal accession does not meet both
# partner pools equally often, partner effects leak into its contrast -- so it is
# asserted exactly rather than approximately.
#
# Run: Rscript tests/test_validation.R

library(tidyverse)
here::i_am("tests/test_validation.R")
source(here::here("tests", "helper.R"))
load_code("validation_functions.R")

set.seed(5)

# A synthetic species: producer and associate effects with a known negative
# correlation, and enough partners that everyone is eligible.
make_inputs <- function(n = 120, cor_pr_as = -0.4, pev = 50, sigma2_e = 500) {
  z1 <- rnorm(n); z2 <- rnorm(n)
  z2 <- cor_pr_as * z1 + sqrt(1 - cor_pr_as^2) * residuals(lm(z2 ~ z1))
  np <- sample(3:10, n, replace = TRUE)

  # A PER-ACCESSION PEV that falls with partner count, which is the real
  # pattern: an accession grown with two partners is worse known than one grown
  # with ten. Scaled so its mean is the global `pev`, so a test that compares
  # the per-accession and global routes compares like with like.
  raw   <- 10 / np
  pev_i <- pev * raw / mean(raw)
  # Var(true) = Var(BLUP) + E[PEV], so the component variance follows from the
  # two rather than being chosen independently of them.
  sigma2_As <- stats::var(as.vector(scale(z2)) * 8) + pev

  list(
    species = "test",
    accessions = tibble::tibble(
      acc = sprintf("acc%03d", seq_len(n)),
      Pr = as.vector(scale(z1)) * 10,
      As = as.vector(scale(z2)) * 8,
      GMA = Pr + As,
      n_partners = np,
      n_plots = sample(4:12, n, replace = TRUE),
      PEV_Pr_i  = pev_i,
      PEV_As_i  = pev_i,
      PEV_GMA_i = 2 * pev_i,
      rel_As_i  = pmin(pmax(1 - pev_i / sigma2_As, 1e-6), 0.999),
      rel_GMA_i = pmin(pmax(1 - 2 * pev_i / (2 * sigma2_As), 1e-6), 0.999),
      eligible = TRUE),
    PEV_As = pev, sigma2_e = sigma2_e, sigma2_Pr = 100,
    sigma2_As = sigma2_As, sigma2_GMA = 2 * sigma2_As,
    rel_As = 0.3, response = "partner yield", n_trials = 8L, n_plots = 3181L,
    # The eligibility record validation_inputs() carries. rel_min is set below
    # every fixture accession's reliability, so `eligible = TRUE` above is
    # consistent with it; a test that needs an ineligible accession plants one
    # on a copy.
    eligibility = list(rule = "reliability", rel_min = 0.2,
                       rel_min_source = "fixture", min_partners = 3L,
                       partner_floor = 2L, n_eligible = n, n_total = n,
                       n_ref = n)
  )
}
inp <- make_inputs()

# The contract build_pools() depends on, asserted in ONE place. When
# validation_inputs() grows a column that build_pools() requires, this list and
# make_inputs() are the two things to extend -- without it, adding a required
# column breaks ~50 assertions at once with no clue which contract moved.
REQUIRED_ACC_COLS <- c("acc", "Pr", "As", "GMA", "n_partners", "n_plots",
                       "PEV_As_i", "rel_As_i", "PEV_GMA_i", "eligible")
check(all(REQUIRED_ACC_COLS %in% names(inp$accessions)),
      "the fixture supplies every column build_pools() requires")

# The candidate set build_pools() selects from, derived from its own returned
# threshold. Three places used to re-implement this filter by hand; a change to
# the filter then had to be mirrored in all three or the oracles below would
# silently compare against the wrong candidate set.
.cand_of <- function(inp, bp) {
  dplyr::filter(inp$accessions, eligible, Pr >= bp$pr_min)
}

# ------------------------------------------------------------
# 1. build_pools: disjoint, high on Pr, matched on Pr
#
# The naive rule -- top n on (Pr + As) against top n on (Pr - As) -- does NOT
# give disjoint pools, because a high-Pr accession makes both lists. That is the
# whole reason for the constrained version, so both halves are asserted.
# ------------------------------------------------------------

for (n in c(10, 20, 30)) {
  bp <- suppressWarnings(build_pools(inp, n, pr_quantile = 0.5, tol = 1.0))
  plus  <- bp$pools$acc[bp$pools$pool == "As+"]
  minus <- bp$pools$acc[bp$pools$pool == "As-"]

  check(length(plus) == n && length(minus) == n,
        sprintf("n = %d gives n per pool", n))
  check(length(intersect(plus, minus)) == 0,
        sprintf("the pools are disjoint at n = %d", n))
  # The constraint cannot always be met exactly: dPr is a STEP function of
  # theta, so the achieved gap is whatever the nearest jump leaves. What is
  # assertable is that it never does worse than selecting on As alone, and that
  # the scan really found the best |dPr| available rather than a local one.
  cand_n <- .cand_of(inp, bp)
  check(abs(bp$dPr) <= abs(.pools_at(cand_n, n, 0)$dPr) + 1e-8,
        sprintf("constraining never widens the producer gap at n = %d", n))

  # The oracle for the scan: an independent dense sweep must not find a theta
  # with a smaller |dPr|. This is what the old bisection could fail -- it
  # stopped at a sign change, which on a non-monotone step function need not be
  # the global minimum.
  sweep_dPr <- purrr::map_dbl(seq(0, 10, length.out = 1001),
                              \(th) .pools_at(cand_n, n, th)$dPr)
  check(abs(bp$dPr) <= min(abs(sweep_dPr)) + 1e-8,
        sprintf("the scan finds the global best |dPr| at n = %d", n))
  check(bp$dAs > 0, sprintf("the associate contrast is positive at n = %d", n))
  check(all(bp$pools$Pr >= bp$pr_min),
        sprintf("every member clears the Pr threshold at n = %d", n))
  check(bp$mean_Pr_plus > mean(inp$accessions$Pr) &&
          bp$mean_Pr_minus > mean(inp$accessions$Pr),
        sprintf("both pools are above the population mean Pr at n = %d", n))
}

# the naive index really does overlap -- the negative control for the above
cand <- dplyr::filter(inp$accessions, Pr >= stats::median(Pr))
naive_plus  <- dplyr::slice_max(cand, Pr + As, n = 30)$acc
naive_minus <- dplyr::slice_max(cand, Pr - As, n = 30)$acc
check(length(intersect(naive_plus, naive_minus)) > 0,
      "the naive Pr+As / Pr-As rule DOES overlap, which is why it is not used")

# ------------------------------------------------------------
# 1b. The selection index is the Lagrangian of the Pr constraint
#
# ONE index serves both pools -- rank by As + theta*Pr, take the top n and the
# bottom n -- so theta trades As extremity for Pr balance along a single axis.
# Three properties follow, and all three are independent of any particular
# number the function returns.
# ------------------------------------------------------------

cand_m <- .cand_of(inp, suppressWarnings(build_pools(inp, 20, pr_quantile = 0.5)))
theta_seq <- seq(0, 5, by = 0.05)
sweep <- purrr::map(theta_seq, \(th) .pools_at(cand_m, 20, th))
dPr_by_theta <- purrr::map_dbl(sweep, "dPr")
dAs_by_theta <- purrr::map_dbl(sweep, "dAs")

check(dPr_by_theta[1] < 0,
      "selecting on As alone leaves the As+ pool short on Pr")
check(stats::cor(theta_seq, dPr_by_theta) > 0.8,
      "dPr rises with theta in trend, which is what makes theta the right knob")

# theta = 0 is the unconstrained optimum, so nothing can beat it on contrast
check(dAs_by_theta[1] >= max(dAs_by_theta) - 1e-8,
      "theta = 0 maximises the contrast: it is the unconstrained solution")
check(stats::cor(theta_seq, dAs_by_theta) < -0.8,
      "and buying Pr balance costs contrast, monotonically in trend")

# disjointness is STRUCTURAL under one index: the top n and bottom n of a
# single ordering cannot collide. No exclusion step is needed, and this holds
# at every theta rather than only at the chosen one.
overlaps <- purrr::map_int(sweep,
  \(r) length(intersect(r$plus$acc, r$minus$acc)))
check(all(overlaps == 0),
      "top-n and bottom-n of one ordering are disjoint at every theta")

# ------------------------------------------------------------
# 1c. NEGATIVE: the two-index rule this replaced is worse
#
# It scored the As- pool by As - theta*Pr, which REWARDS high Pr there, so
# raising theta dragged both pools up rather than equalising them. Without this
# check a revert would pass the whole suite.
# ------------------------------------------------------------

two_index <- function(cand, n, theta) {
  plus <- cand |> dplyr::mutate(sc = As + theta * Pr) |>
    dplyr::slice_max(sc, n = n, with_ties = FALSE)
  minus <- cand |> dplyr::filter(!acc %in% plus$acc) |>
    dplyr::mutate(sc = As - theta * Pr) |>
    dplyr::slice_min(sc, n = n, with_ties = FALSE)
  c(dAs = mean(plus$As) - mean(minus$As),
    dPr = mean(plus$Pr) - mean(minus$Pr))
}
old_sweep <- vapply(theta_seq, \(th) two_index(cand_m, 20, th), numeric(2))
best_old <- max(old_sweep["dAs", abs(old_sweep["dPr", ]) <= 1.0])
best_new <- max(dAs_by_theta[abs(dPr_by_theta) <= 1.0])
check(best_new > best_old,
      sprintf("one index beats two at the same tolerance (%.2f vs %.2f)",
              best_new, best_old))

# the contrast must shrink as the pools grow: that is the trade-off the power
# calculation balances, and if it stopped holding the pool-size argument changes
d10 <- suppressWarnings(build_pools(inp, 10))$dAs
d30 <- suppressWarnings(build_pools(inp, 30))$dAs
check(d10 > d30, "a larger pool gives a smaller contrast")

# and when it cannot be met, the caller must be TOLD, not left to assume the
# pools are matched
w <- tryCatch({ build_pools(inp, 10, tol = 0.01); NULL },
              warning = function(w) conditionMessage(w))
check(!is.null(w) && grepl("Pr constraint could not be met", w),
      "an unmeetable producer tolerance raises a warning that says so")
check(is.null(tryCatch({ build_pools(inp, 20, tol = 5); NULL },
                       warning = function(w) conditionMessage(w))),
      "and a tolerance it can meet raises none")

# a pool larger than the candidate set must be refused, not silently truncated
check_error(build_pools(inp, 100, pr_quantile = 0.5),
            "a pool bigger than the candidate set is refused")

# Only eligible accessions may be selected. Note what this does and does not
# pin: build_pools() honours the `eligible` COLUMN, and the rule that computes
# that column lives upstream in validation_inputs(). So this assertion survives
# a change of eligibility rule unchanged, which is the point -- it separates the
# selection machinery from the policy it is handed.
inp_partial <- inp
inp_partial$accessions$eligible <- inp_partial$accessions$n_partners >= 8
bp_e <- suppressWarnings(build_pools(inp_partial, 8))
check(all(bp_e$pools$n_partners >= 8),
      "ineligible accessions are never selected")
check(nrow(bp_e$pools) == 16,
      "and an eligibility filter that bites still fills both pools")

# ------------------------------------------------------------
# 2. pool_diff: churn between vintages
# ------------------------------------------------------------

old_p <- suppressWarnings(build_pools(inp, 20))$pools
new_p <- old_p
new_p$acc[1] <- "brand_new"                       # one entered, one left
d <- pool_diff(old_p, new_p)
check(sum(d$detail$status == "entered") == 1, "an entering accession is reported")
check(sum(d$detail$status == "left") == 1, "a leaving accession is reported")
check(sum(d$detail$status == "retained") == nrow(old_p) - 1,
      "and the rest are retained")

swapped <- old_p
swapped$pool <- ifelse(swapped$pool == "As+", "As-", "As+")
check(all(pool_diff(old_p, swapped)$detail$status == "switched pool"),
      "an accession changing pool is reported as switched, not retained")

# ------------------------------------------------------------
# 3. contrast_se / contrast_power: the shape of the formula
#
# The claim the design document rests on is that a large share of the standard
# error cannot be bought with plots. That is a property of the formula, so it is
# testable without any field data.
#
# EVERY CALL GOES THROUGH PW() AND NAMES ITS ARGUMENTS. There used to be nine
# fully positional calls here relying on the order
# (dAs, sigma2_within, n, sigma2_e, P). Inserting any argument before `P` would
# have silently rebound all nine while the suite kept reporting PASS -- the
# quietest way this file could have stopped testing what it says it tests.
# ------------------------------------------------------------

PW <- function(dAs = 15, sigma2_within = 100, n = 20, sigma2_e = 500,
               P = 400, ...) {
  contrast_power(dAs = dAs, sigma2_within = sigma2_within, n = n,
                 sigma2_e = sigma2_e, P = P, ...)
}

p1 <- PW()
check(p1$power > 0 && p1$power < 1, "power is a probability")
check(p1$SE > p1$SE_floor, "the SE exceeds its infinite-plot floor")
check(p1$pct_SE_irreducible > 0 && p1$pct_SE_irreducible <= 100,
      "the irreducible share is a percentage")

# contrast_se is exactly the formula the documentation displays. An algebraic
# identity, so the tolerance is 1e-12: this is the one assertion that would
# catch a term being dropped, rescaled or double-counted.
check_near(
  contrast_se(sigma2_within = 100, n = 20, sigma2_e = 500, P = 400,
              interaction_frac = 0.4, delta = 10, n_loc = 5),
  sqrt(2 * 100 / 20 + 4 * 500 / 400 + (0.4 * 10)^2 / 5),
  tol = 1e-12, "contrast_se is exactly its documented three-term formula")

# more plots must help, but only up to the floor
se_of <- function(P) PW(P = P)$SE
check(se_of(300) > se_of(500), "more plots reduce the SE")
check(se_of(1e9) > p1$SE_floor * 0.999 && se_of(1e9) < p1$SE_floor * 1.001,
      "and with enough plots the SE converges on the floor, not on zero")

# The same limit claim with a pool x location interaction present, stated
# WITHOUT reference to SE_floor so it is true both before and after the floor is
# redefined to include that term. The interaction does not shrink with P, so the
# infinite-plot SE must stay strictly above the accession-only term.
lim_int <- PW(P = 1e12, interaction_frac = 0.5, lambda = 1)$SE
check(lim_int > sqrt(2 * 100 / 20) * 1.001,
      "with an interaction present the infinite-plot SE exceeds the accession-only term")

# more accessions reduce the part plots cannot touch
check(PW(n = 40)$SE_floor < PW(n = 20)$SE_floor,
      "more accessions per pool lowers the floor")

# a one-sided test is more powerful than two-sided for a true positive effect
check(PW(sided = 1)$power > PW(sided = 2)$power,
      "one-sided beats two-sided on a positive contrast")

# attenuation and interaction both cost power
check(PW(lambda = 0.6)$power < PW(lambda = 1.0)$power,
      "attenuation costs power")
check(PW(interaction_frac = 0.5)$power < PW(interaction_frac = 0)$power,
      "a pool x location interaction costs power")

# a zero contrast must give power at the nominal level, not zero and not one
check_near(PW(dAs = 0, sided = 1, alpha = 0.05)$power, 0.05, tol = 1e-6,
           "a zero contrast gives power equal to alpha")

# ------------------------------------------------------------
# 3b. The degrees of freedom are pinned
#
# Nothing used to assert this, so a change from 2n-2 to anything else would have
# passed the whole suite in silence. The experimental unit for a pool contrast
# is the ACCESSION, so df = 2n-2 is a claim about the design and not an
# implementation detail; it is the boundary case any later Satterthwaite
# correction must still reproduce exactly when there is no interaction term.
# ------------------------------------------------------------

check(PW(n = 20)$df == 38, "df is 2n-2 with no pool x location term")
check(PW(n = 30)$df == 58, "and it tracks n")
check(PW(n = 20, interaction_frac = 0)$df == 38,
      "an interaction of exactly zero is the same case")

# ------------------------------------------------------------
# 3c. power_grid: it rebuilds the pools at every n
#
# The trade-off between contrast and precision is the whole pool-size question,
# so holding the pools fixed across n would make the curve meaningless. Nothing
# tested this function at all.
# ------------------------------------------------------------

g <- suppressWarnings(power_grid(inp, n_values = c(10, 20),
                                 P_values = c(300, 400),
                                 lambda_values = 0.7,
                                 interaction_values = 0.5))
check(nrow(g) == 4, "power_grid returns one row per (n, P, lambda, interaction)")
check(dplyr::n_distinct(g$dAs_predicted) == 2,
      "the pools are rebuilt at each n -- the trade-off is not held fixed")
check(g$dAs_predicted[g$n == 10][1] > g$dAs_predicted[g$n == 20][1],
      "and the rebuilt contrast shrinks with n, as build_pools gives it")
# the grid must carry the sigma2_within build_pools computed, not a recycled or
# averaged stand-in: this is what catches a length-2 value silently recycling
check_near(g$sigma2_within[g$n == 20][1],
           suppressWarnings(build_pools(inp, 20))$sigma2_within,
           tol = 1e-12,
           "the grid's sigma2_within is the one build_pools computed")

# ------------------------------------------------------------
# 3c-bis. Plots per location is the coherent budget axis
#
# Crossing a TOTAL plot budget with a location count produces cells that are
# not the same design: 400 plots is 100 per site at four locations and 80 at
# five. A curve drawn against total P therefore joins points from two different
# designs, which is what made power_curves.png saw-tooth. Sweeping plots per
# location and deriving P keeps the location arms comparable.
# ------------------------------------------------------------

gp <- suppressWarnings(power_grid(inp, n_values = 20,
                                  plots_per_loc_values = c(60, 80),
                                  lambda_values = 0.8,
                                  interaction_values = 0.5,
                                  n_loc_values = c(4, 5)))
check(nrow(gp) == 4, "one row per (plots/loc, locations)")
check_near(sort(gp$P), sort(c(60*4, 60*5, 80*4, 80*5)), tol = 1e-12,
           "P is derived as plots_per_loc * n_locations")
# the two arms share a per-site effort, which the total-P parameterisation
# could not deliver
check(setequal(gp$plots_per_loc[gp$n_locations == 4],
               gp$plots_per_loc[gp$n_locations == 5]),
      "both location arms are evaluated at the same plots per location")
# and at a FIXED per-site effort, more locations must help: the interaction
# term is the one divided by n_loc, and it is the largest of the three
check(gp$power[gp$plots_per_loc == 80 & gp$n_locations == 5] >
        gp$power[gp$plots_per_loc == 80 & gp$n_locations == 4],
      "more locations raise power at the same plots per location")
check(gp$df[gp$plots_per_loc == 80 & gp$n_locations == 5] >
        gp$df[gp$plots_per_loc == 80 & gp$n_locations == 4],
      "and buy degrees of freedom for the stratum that has fewest")
# NEGATIVE: the same TOTAL budget at different location counts is NOT the same
# design, which is the whole reason for the change
g400 <- suppressWarnings(power_grid(inp, n_values = 20, P_values = 400,
                                    lambda_values = 0.8,
                                    interaction_values = 0.5,
                                    n_loc_values = c(4, 5)))
check(abs(diff(g400$power)) > 0.01,
      "400 plots at 4 locations and at 5 is NOT one design -- they differ in power")
check_error(power_grid(inp, n_values = 20, lambda_values = 0.8),
            "a grid with neither budget axis is refused")

# ------------------------------------------------------------
# 3d. The named-only barrier, and that it is still there
#
# `...` sits ahead of every argument but the first, so R cannot match them
# positionally or partially. These two assertions are what stop a future
# refactor quietly removing that barrier: without them the discipline rots and
# the nine-positional-call hazard comes straight back.
# ------------------------------------------------------------

check_error(contrast_power(15, 100, 20, 500, 400),
            "positional arguments past `dAs` are refused")
check_error(contrast_se(100, 20, 500, 400),
            "contrast_se refuses positional arguments past `sigma2_within`")
# A misspelled name would otherwise vanish into `...` and be silently ignored,
# which is the quietest failure of all.
check_error(contrast_power(dAs = 15, sigma2_withn = 100, n = 20,
                           sigma2_e = 500, P = 400),
            "a misspelled argument name is refused, not silently ignored")

# ------------------------------------------------------------
# 3e. SE_floor is the P -> infinity limit, interaction included
#
# It used to be sqrt(2*sigma2_within/n), which was the limit only because the
# interaction term defaulted to zero. Two of the three terms are free of P, so
# the floor must contain both -- otherwise power_ceiling.csv and the design
# document's "with infinitely many plots" claim mean different things.
# ------------------------------------------------------------

seI <- PW(P = 1e12, interaction_frac = 0.5, lambda = 1)
check_near(seI$SE, seI$SE_floor, tol = 1e-6 * seI$SE,
           "SE_floor is the P -> infinity limit with the interaction term too")
check(seI$SE_accession_only < seI$SE_floor - 1e-8,
      "and the accession-only term is strictly below that limit")
# more accessions lower the ACCESSION term; they do nothing to the interaction
check(PW(n = 40)$SE_accession_only < PW(n = 20)$SE_accession_only,
      "more accessions per pool lowers the accession-only floor")
check_near(PW(n = 40, interaction_frac = 0.5)$v_int,
           PW(n = 20, interaction_frac = 0.5)$v_int, tol = 1e-12,
           "but the pool x location term does not depend on n at all")
# the three shares are a decomposition, so they sum to 100
pv <- PW(interaction_frac = 0.5)
check_near(pv$pct_var_accessions + pv$pct_var_interaction + pv$pct_var_plots,
           100, tol = 1e-8, "the three variance shares decompose the total")

# ------------------------------------------------------------
# 3f. Pool-specific sigma2_within: identical SE, different df
#
# (A + B)/n is identically 2*mean/n, so splitting the within-pool variance by
# pool cannot change the standard error. It DOES change the degrees of freedom,
# through the Welch term. Both halves of that are asserted, because getting the
# first wrong would be a silent bias and getting the second wrong would be a
# silent anticonservatism.
# ------------------------------------------------------------

check_near(contrast_se(sigma2_within = c(100, 100), n = 20,
                       sigma2_e = 500, P = 400),
           contrast_se(sigma2_within = 100, n = 20, sigma2_e = 500, P = 400),
           tol = 1e-12,
           "pool-specific sigma2_within reduces to the scalar form when equal")
check_near(contrast_se(sigma2_within = c(80, 120), n = 20,
                       sigma2_e = 500, P = 400),
           contrast_se(sigma2_within = 100, n = 20, sigma2_e = 500, P = 400),
           tol = 1e-12,
           "and only the MEAN of the two enters the SE: (A+B)/n = 2*mean/n")
check_near(PW(sigma2_within = c(100, 100))$df_acc, 38, tol = 1e-10,
           "the Welch accession df is exactly 2n-2 at equal pool variances")
check(PW(sigma2_within = c(40, 160))$df_acc < 38,
      "and unequal pool variances cost accession df")

# ------------------------------------------------------------
# 3g. The Satterthwaite df, bounded by the strata it combines
#
# The pool x location term carries n_loc - 1 = 4 df and is the largest of the
# three at the measured interaction, so the effective df must sit between 4 and
# the naive 2n - 2. Nothing used to assert any of this.
# ------------------------------------------------------------

pS <- PW(interaction_frac = 0.5, lambda = 1, n_loc = 5)
check(pS$df > 4 && pS$df < 38,
      "Satterthwaite df lies between the interaction df and the naive df")
check(PW(interaction_frac = 1.0, lambda = 1, n_loc = 5)$df < pS$df,
      "a larger interaction share lowers the effective df")
check(pS$power < PW(interaction_frac = 0.5, lambda = 1,
                    df_mode = "naive")$power,
      "and the correction costs power, which is why it was worth making")
check(PW(interaction_frac = 0.5, df_mode = "naive")$df == 38,
      "df_mode = 'naive' recovers 2n-2, so the correction is detectably on")
check(PW(interaction_frac = 0.5, n_loc = 10)$df >
        PW(interaction_frac = 0.5, n_loc = 5)$df,
      "more locations buy df for the term that has the fewest")

# ------------------------------------------------------------
# 3h. The two interaction parameterisations
#
# interaction_frac is sd(lambda)/|mean(lambda)|, so interaction_frac * delta is
# the interaction SD only when `lambda` is the lambda the spread was measured
# around. The two forms must agree exactly there and diverge when they are not,
# and supplying both must be refused rather than double-counted.
# ------------------------------------------------------------

lam <- 0.83
check_near(PW(lambda = lam, interaction_frac = 0.73)$SE,
           PW(lambda = lam, interaction_sd = 0.73 * lam * 15)$SE,
           tol = 1e-10,
           "the two interaction forms agree exactly at the measured lambda")
check(PW(lambda = 0.7, interaction_sd = 0.73 * lam * 15)$power <
        PW(lambda = 0.7, interaction_frac = 0.73)$power,
      "and the absolute form is the conservative one at a different lambda")
check_error(PW(interaction_frac = 0.5, interaction_sd = 5),
            "supplying both interaction parameterisations is refused")

# ------------------------------------------------------------
# 3i. build_pools returns the candidate set it selected from
#
# The oracle this replaces was a hand-written copy of build_pools()'s own
# filter, maintained in three places in this file.
# ------------------------------------------------------------

bp_c <- suppressWarnings(build_pools(inp, 20, pr_quantile = 0.5))
check(setequal(bp_c$candidates$acc, .cand_of(inp, bp_c)$acc),
      "the returned candidate set is exactly the filter's result")
check(all(bp_c$pools$acc %in% bp_c$candidates$acc),
      "and every selected accession came from it")
check_near(bp_c$sigma2_within,
           mean(c(bp_c$sigma2_within_plus, bp_c$sigma2_within_minus)),
           tol = 1e-12,
           "the scalar sigma2_within is the mean of the two pool values")
check_near(bp_c$var_As_within_pools,
           mean(c(stats::var(bp_c$pools$As[bp_c$pools$pool == "As+"]),
                  stats::var(bp_c$pools$As[bp_c$pools$pool == "As-"]))),
           tol = 1e-12,
           "var_As_within_pools is the within-pool spread of the predictions")

# ------------------------------------------------------------
# 4. make_validation_design: the balance the contrast depends on
# ------------------------------------------------------------

oat_pools <- dplyr::mutate(suppressWarnings(build_pools(inp, 20))$pools,
                           species = "oat")
pea_pools <- dplyr::mutate(suppressWarnings(build_pools(make_inputs(130), 20))$pools,
                           species = "pea")

des <- make_validation_design(oat_pools, pea_pools, n_loc = 5L, p_loc = 80L,
                              n_anchor = 1L, block_size = 10L, seed = 1L)
chk <- check_validation_design(des, oat_pools, pea_pools)

check(all(chk$plots_per_location$plots == 80),
      "every location gets exactly the requested number of plots")
check(dplyr::n_distinct(chk$cells_per_location$plots) == 1,
      "all four cells are equally represented at every location")
check(chk$max_oat_imbalance == 0,
      "every oat meets both pea pools equally often -- the load-bearing property")
check(chk$max_pea_imbalance == 0,
      "every pea meets both oat pools equally often")
check(chk$anchor_plots == 2 * 4 * 1 * 5,
      "the anchors are duplicated within every location")
check(nrow(chk$replicated_combinations) > 0,
      "some combinations are replicated, so plot error is estimable")
check(all(sort(unique(des$oat_pool)) == c("As-", "As+")),
      "both oat pools appear")
check(all(sort(unique(des$pea_pool)) == c("As-", "As+")),
      "both pea pools appear")

# the anchor combinations must be the SAME at every location, or combination x
# location is not estimable
anchors_by_loc <- des |>
  dplyr::filter(is_anchor) |>
  dplyr::distinct(location, combination) |>
  dplyr::count(combination)
check(all(anchors_by_loc$n == 5),
      "each anchor combination appears at all five locations")

# the price of a fixed anchor: those accessions appear in more plots than the
# rest. Benign -- the analysis averages to one value per accession first -- but
# it should be visible rather than surprising.
anchor_accs <- unique(dplyr::filter(des, is_anchor)$oat_acc)
reps <- chk$oat_reps
check(all(reps$plots[reps$oat_acc %in% anchor_accs] >
            max(reps$plots[!reps$oat_acc %in% anchor_accs])),
      "anchor accessions appear in more plots than the rest, as documented")

# the anchor plots must come OUT of the budget, not be added on top
check(sum(chk$plots_per_location$plots) == 5 * 80,
      "anchors are inside the plot budget")

# a budget too small for four cells plus anchors must be refused
check_error(make_validation_design(oat_pools, pea_pools, n_loc = 5L, p_loc = 8L),
            "too small a budget is refused")

# ------------------------------------------------------------
# 5. The Pr-adjusted contrast
#
# A pool pair with a larger raw dAs can be WORSE than a smaller one when it
# also diverges in Pr in the direction the negative As-Pr correlation already
# predicts: part of its As gap is then explained by Pr rather than being
# evidence of selectable variation at fixed Pr. These pin the three changes
# the adjustment makes -- the contrast, the variance, the df -- and, first of
# all, that supplying nothing changes nothing.
# ------------------------------------------------------------

PA_ARGS <- list(sigma2_within = 100, n = 20, sigma2_e = 500, P = 400,
                interaction_frac = 0)

# --- the backwards-compatibility contract ---
#
# Not a value freeze: the point is that the unadjusted path is untouched, so
# the ~50 assertions above that describe it keep describing it. pr_adjust is
# applied only inside `if (!is.null(pr_adjust))`, so this must be EXACT, not
# merely close.
pa_off <- do.call(contrast_power, c(list(dAs = 20), PA_ARGS))
pa_nul <- do.call(contrast_power, c(list(dAs = 20), PA_ARGS,
                                    list(pr_adjust = NULL)))
num_cols <- names(pa_off)[vapply(pa_off, is.numeric, logical(1))]
check(identical(as.list(pa_off[num_cols]), as.list(pa_nul[num_cols])),
      "pr_adjust = NULL is bit-identical to omitting it")
check(pa_off$df_acc == 2 * 20 - 2 && pa_off$df == 2 * 20 - 2,
      "unadjusted df_acc is still exactly 2n - 2")
check(is.na(pa_off$b_Pr) && pa_off$var_inflation == 1 &&
        pa_off$df_cost_covariate == 0L,
      "the unadjusted row reports no inflation and no covariate df")

# --- dPr = 0: the covariate costs a df even when it corrects nothing ---
#
# The sharp case. With the pools matched on Pr there is nothing to adjust, so
# the contrast and the SE must be untouched -- but the slope was still
# estimated, so one df is still spent.
pa_zero <- do.call(contrast_power,
                   c(list(dAs = 20), PA_ARGS,
                     list(pr_adjust = list(dPr = 0, b = -0.2, Sxx = 1000))))
check_near(pa_zero$SE, pa_off$SE, 1e-12,
           "dPr = 0 leaves the SE exactly unchanged")
check_near(pa_zero$dAs_adj, 20, 1e-12,
           "dPr = 0 leaves the contrast unchanged")
check_near(pa_zero$var_inflation, 1, 1e-12, "dPr = 0 gives f = 1")
check_near(pa_zero$df_acc, 2 * 20 - 3, 1e-12,
           "dPr = 0 still costs one df: df_acc is exactly 2n - 3")

# --- the algebra, against an oracle computed by hand in the test ---
pa_in  <- list(dPr = 2.5, b = -0.25, Sxx = 1200)
pa_adj <- do.call(contrast_power, c(list(dAs = 20), PA_ARGS,
                                    list(pr_adjust = pa_in)))
f_hand <- 1 + (20 / 2) * pa_in$dPr^2 / pa_in$Sxx
se_hand <- sqrt(f_hand * (2 * 100 / 20) + f_hand * (4 * 500 / 400))
check_near(pa_adj$var_inflation, f_hand, 1e-12, "f matches 1 + (n/2) dPr^2 / Sxx")
check_near(pa_adj$SE, se_hand, 1e-12,
           "the adjusted SE is the hand-computed closed form")
check_near(pa_adj$dAs_adj, 20 - pa_in$b * pa_in$dPr, 1e-12,
           "dAs_adj = dAs - b * dPr")
# The whole reason the inflation is a FACTOR on both accession-level strata
# rather than an additive fourth term.
check_near(pa_adj$pct_var_accessions + pa_adj$pct_var_interaction +
             pa_adj$pct_var_plots, 100, 1e-10,
           "the adjusted variance still splits into three shares summing to 100")

# --- the mechanism, in the direction that matters ---
#
# b < 0 and dPr < 0 is the "wrong direction": the As+ pool is the WEAKER
# producer, so the negative As-Pr correlation already predicts part of its As
# advantage, and the adjusted contrast must shrink below the raw one.
pa_wrong <- do.call(contrast_power,
                    c(list(dAs = 20), PA_ARGS,
                      list(pr_adjust = list(dPr = -3, b = -0.25, Sxx = 1200))))
check(pa_wrong$dAs_adj < 20,
      "dPr in the direction the As-Pr correlation predicts shrinks the contrast")

# monotone in both channels, at a FIXED raw dAs
pw_at <- function(dPr, Sxx = 1200) {
  do.call(contrast_power,
          c(list(dAs = 20), PA_ARGS,
            list(pr_adjust = list(dPr = dPr, b = 0, Sxx = Sxx))))$power
}
check(pw_at(0) > pw_at(2) && pw_at(2) > pw_at(4),
      "at fixed dAs, a wider producer gap lowers power")
check(pw_at(3, Sxx = 3000) > pw_at(3, Sxx = 800),
      "more within-pool Pr spread to adjust along raises power")

# --- the negative control: the dAs-maximal pools score WORSE ---
#
# This is the test that proves the index does the job it was added for. On the
# fixture, as on the real data, theta = 0 maximises raw dAs but leaves the
# pools badly mismatched on Pr, and the adjusted score must prefer the
# Pr-matched solution despite its smaller dAs.
bp20  <- build_pools(inp, 20)
at0   <- .pools_at(bp20$candidates, 20, 0)
b_fix <- pr_adjust_slope(bp20$candidates)
score_of <- function(dAs, dPr, s2w, Sxx) {
  contrast_power(dAs = dAs, sigma2_within = s2w, n = 20,
                 sigma2_e = inp$sigma2_e, P = 400, interaction_frac = 0,
                 pr_adjust = list(dPr = dPr, b = b_fix, Sxx = Sxx))$dAs_adj
}
sxx_of <- function(p, m) 19 * (var(p$Pr) + var(m$Pr))
check(abs(at0$dPr) > abs(bp20$dPr),
      "theta = 0 really is the Pr-mismatched solution on the fixture")
check(at0$dAs > bp20$dAs, "and it really does have the larger raw dAs")
check(score_of(at0$dAs, at0$dPr,
               c(mean(at0$plus$PEV_As_i) + var(at0$plus$As),
                 mean(at0$minus$PEV_As_i) + var(at0$minus$As)),
               sxx_of(at0$plus, at0$minus)) <
        score_of(bp20$dAs, bp20$dPr,
                 c(bp20$sigma2_within_plus, bp20$sigma2_within_minus),
                 bp20$Sxx_w),
      "yet the Pr-adjusted contrast prefers the Pr-matched pools: the index re-ranks")

# --- the slope is the candidate-set one, and it is defined once ---
check_near(b_fix, unname(coef(lm(As ~ Pr, data = bp20$candidates))[2]), 1e-12,
           "pr_adjust_slope() is the As-on-Pr slope among the candidates")
check(abs(b_fix - unname(coef(lm(As ~ Pr, data = inp$accessions))[2])) > 1e-9,
      "and it differs from the eligible-set slope, so which set is used matters")
check_error(pr_adjust_slope(bp20$candidates[1:2, ]),
            "too few candidates to regress is refused")

# --- the pr_adjust list is all-or-nothing ---
check_error(contrast_se(100, n = 20, sigma2_e = 500, P = 400,
                        pr_adjust = list(dPr = 1, b = -0.2)),
            "a half-supplied pr_adjust is refused")
check_error(contrast_se(100, n = 20, sigma2_e = 500, P = 400,
                        pr_adjust = list(dPr = 1, b = -0.2, Sxx = 0)),
            "a non-positive Sxx is refused")
# A misspelling must not vanish into `...` -- the quietest failure of all.
check_error(contrast_se(100, n = 20, sigma2_e = 500, P = 400,
                        pr_adjst = list(dPr = 1, b = -0.2, Sxx = 10)),
            "a misspelled pr_adjust is refused rather than ignored")

# power_grid stays unadjusted unless asked
pg_off <- power_grid(inp, n_values = 20, plots_per_loc_values = 80,
                     lambda_values = 1, interaction_values = 0)
pg_on  <- power_grid(inp, n_values = 20, plots_per_loc_values = 80,
                     lambda_values = 1, interaction_values = 0,
                     pr_adjust = TRUE)
check(all(is.na(pg_off$b_Pr)) && !any(is.na(pg_on$b_Pr)),
      "power_grid() leaves the adjustment off by default and on when asked")

# ------------------------------------------------------------
# 6. The pool search
#
# The search exists for the ANCHORED problem, not to improve the theta rule --
# which on the real data it reproduces to four significant figures. So the
# property that matters most is that it can never do WORSE, and that is a
# consequence of seeding it with build_pools()' own solution rather than a
# hope about hill-climbing.
# ------------------------------------------------------------

sctx <- pool_search_context(bp20$candidates, n = 20, sigma2_e = inp$sigma2_e,
                            P = 400)
check(identical(sctx$acc, sort(bp20$candidates$acc)),
      "the search context is sorted by accession, so it is order-independent")
check_near(mean(sctx$As), 0, 1e-12, "and As is centred, killing the cancellation")

# --- the O(1) sufficient statistics, against a from-scratch oracle ---
#
# The one place a silent wrong answer is easy: an incremental update that
# drifts would quietly optimise the wrong objective.
fam <- .theta_family(sctx)
check(length(fam) > 1L && length(fam) < 2001L,
      "the theta scan collapses to far fewer distinct pool pairs than thetas")
set.seed(11)
incr_ok <- TRUE
for (r in seq_len(50)) {
  s <- fam[[sample.int(length(fam), 1L)]]
  from_scratch <- pool_score(s$plus, s$minus, sctx)
  shuffled <- pool_score(sample(s$plus), sample(s$minus), sctx)
  if (!isTRUE(all.equal(from_scratch$score, shuffled$score, tolerance = 1e-13))) {
    incr_ok <- FALSE
  }
}
check(incr_ok, "pool_score() is invariant to the order members are listed in")

# --- non-regression: the theorem the seeding buys ---
bp_score <- pool_score(match(bp20$pools$acc[bp20$pools$pool == "As+"], sctx$acc),
                       match(bp20$pools$acc[bp20$pools$pool == "As-"], sctx$acc),
                       sctx)
po20 <- pools_optimal(inp, 20, P = 400)
check(po20$score >= bp_score$score - 1e-12,
      "pools_optimal() is never worse than build_pools(): the non-regression theorem")
check(length(po20$plus) == 20 && length(po20$minus) == 20 &&
        length(intersect(po20$plus, po20$minus)) == 0L,
      "the search returns n per pool, disjoint")
check(all(c(po20$acc_plus, po20$acc_minus) %in% bp20$candidates$acc),
      "every selected accession is a candidate: eligible and above the Pr threshold")
check(po20$dAs > 0, "As+ really is the higher-As pool")

# --- selection is invariant to lambda and interaction_frac ---
#
# This is what makes it safe for pool selection to run without this vintage's
# cross-validation: the score is computed at lambda = 1, interaction = 0, and
# the membership it picks does not depend on either. If a future vintage makes
# the landscape sharp enough for them to matter, this fails rather than the
# pools moving in silence.
sigs <- vapply(list(c(0.5, 0), c(0.76, 0.5), c(1.0, 0.74)), function(cfg) {
  cc <- pool_search_context(bp20$candidates, n = 20, sigma2_e = inp$sigma2_e,
                            P = 400)
  ss <- pools_search(cc, starts = list(list(plus = po20$plus,
                                            minus = po20$minus)))
  paste(sort(c(cc$acc[ss$plus], cc$acc[ss$minus])), collapse = "|")
}, character(1))
check(length(unique(sigs)) == 1L,
      "the selected membership is invariant to lambda and interaction_frac")

# --- determinism, including under a reordered input table ---
check(identical(sort(po20$acc_plus), sort(pools_optimal(inp, 20, P = 400)$acc_plus)),
      "two runs of the search agree")
inp_shuf <- inp
set.seed(99)
inp_shuf$accessions <- inp$accessions[sample.int(nrow(inp$accessions)), ]
po_shuf <- pools_optimal(inp_shuf, 20, P = 400)
check(setequal(c(po_shuf$acc_plus, po_shuf$acc_minus),
               c(po20$acc_plus, po20$acc_minus)),
      "shuffling the accession table does not change the pools chosen")

# --- RNG hygiene ---
#
# validate_design.R's design_seed and validate_power.R's sim_seed must keep
# seeing exactly the stream they see today, so a search must not consume or
# reset randomness the caller can observe.
set.seed(1); want <- runif(1)
set.seed(1); invisible(pools_optimal(inp, 20, P = 400)); got <- runif(1)
check(identical(want, got), "the search leaves the caller's RNG stream untouched")
set.seed(1); invisible(pools_search(sctx, n_restarts = 3L)); got2 <- runif(1)
check(identical(want, got2), "and that holds when restarts do consume RNG")
if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
invisible(.with_seed(7, stats::runif(1)))
check(!exists(".Random.seed", envir = globalenv()),
      ".with_seed() leaves no .Random.seed behind when there was none")

# --- the context refuses an impossible problem ---
check_error(pool_search_context(bp20$candidates[1:10, ], n = 20,
                                sigma2_e = inp$sigma2_e, P = 400),
            "fewer candidates than 2n is refused")

# ------------------------------------------------------------
# 7. The anchored search
#
# Round 2 onward: hold the accessions actually found in seed storage and
# re-optimise around them, so each round asks for a handful of new names
# rather than a fresh pool to go hunting for.
#
# The feasibility arithmetic is the part most worth pinning. A guard whose
# ADVICE is wrong is worse than no guard, so the suggested allowance is parsed
# back out of the error message and checked, not merely triggered.
# ------------------------------------------------------------

cand_acc <- sort(bp20$candidates$acc)
err_of <- function(expr) {
  tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))
}

# --- the worked case, verbatim from the requirement ---
#
# 32 usable anchors, n = 20: a budget of 3 needs 2 * 17 = 34, so it is
# infeasible, and the smallest feasible budget is 20 - 32 %/% 2 = 4.
anch32 <- cand_acc[seq_len(32)]
msg32 <- err_of(pools_anchored(inp, 20, anch32, max_new_per_pool = 3,
                               reference = FALSE, P = 400))
check(!is.na(msg32), "32 anchors with a budget of 3 is refused at n = 20")
check(grepl("max_new_per_pool = 4", msg32, fixed = TRUE),
      "and the message names 4 as the smallest feasible allowance")
check(grepl("find 2 more", msg32, fixed = TRUE),
      "and offers the alternative of finding 2 more of the existing proposals")
pa32 <- pools_anchored(inp, 20, anch32, max_new_per_pool = 4,
                       reference = FALSE, P = 400)
check(pa32$n_anchors_usable == 32 && pa32$m_min == 16,
      "at the feasible budget it runs, holding 16 anchors per pool")

# --- the feasibility arithmetic, swept against a hand-written oracle ---
#
# Swept on pool_anchor_feasibility() rather than through pools_anchored(),
# which keeps a 300-cell sweep instantaneous instead of running 300 searches.
# The wiring from the arithmetic to the error is checked end to end separately,
# just above and just below.
feas_ok <- TRUE; advice_ok <- TRUE; floor_ok <- TRUE
for (nn in c(20L, 30L)) {
  for (A in seq.int(2L * nn - 14L, 2L * nn)) {
    for (mn in 0:10) {
      fs <- pool_anchor_feasibility(nn, A, mn)
      if (fs$feasible != (2L * (nn - mn) <= A)) feas_ok <- FALSE
      if (fs$max_new_min != max(0L, nn - A %/% 2L)) advice_ok <- FALSE
      # The advised allowance must actually BE feasible, and one less must not
      # be -- the property that makes it the *minimum*.
      if (!pool_anchor_feasibility(nn, A, fs$max_new_min)$feasible) floor_ok <- FALSE
      if (fs$max_new_min > 0L &&
          pool_anchor_feasibility(nn, A, fs$max_new_min - 1L)$feasible) {
        floor_ok <- FALSE
      }
    }
  }
}
check(feas_ok, "the anchor constraint is feasible exactly when 2(n - max_new) <= A")
check(advice_ok,
      "and the allowance it advises is max(0, n - A %/% 2)")
check(floor_ok,
      "the advised allowance is feasible and one smaller is not: it is the minimum")
# The integer floor is the whole point: 33 anchors support only 16 per pool.
check(pool_anchor_feasibility(20L, 33L, 4L)$feasible &&
        !pool_anchor_feasibility(20L, 33L, 3L)$feasible,
      "an odd anchor count floors rather than rounds: 33 anchors need max_new 4, not 3")

# --- HARD ORACLE: exactly one feasible membership ---
#
# With 2n anchors and no substitutions allowed, every anchor must be used and
# the membership is forced. Only the split between the pools is free.
anch_exact <- c(utils::head(bp20$candidates$acc[order(bp20$candidates$As,
                                                      decreasing = TRUE)], 20),
                utils::tail(bp20$candidates$acc[order(bp20$candidates$As,
                                                      decreasing = TRUE)], 20))
pa_exact <- pools_anchored(inp, 20, anch_exact, max_new_per_pool = 0,
                           reference = FALSE, P = 400)
check(setequal(c(pa_exact$acc_plus, pa_exact$acc_minus), anch_exact),
      "2n anchors with no substitutions returns exactly those accessions")
check(pa_exact$n_new == 0L,
      "and asks for no new accessions at all")
# The anchors were built as the extreme top and bottom on As, so the split is
# forced too: any other assignment would collapse the contrast.
check(setequal(pa_exact$acc_plus, utils::head(anch_exact, 20)),
      "and puts the high-As anchors in As+, so the split is forced as well")

# --- anchors are honoured, on the RETURNED pools ---
#
# Asserted on the output rather than trusted from the filter, the same way
# section 1 asserts build_pools()' own guarantees.
set.seed(21)
anch_ok <- TRUE; count_ok <- TRUE
for (r in seq_len(6)) {
  A <- sample(34:44, 1)
  an <- sample(cand_acc, A)
  mn <- max(0L, 20L - A %/% 2L) + sample(0:2, 1)
  if (mn > 20L) next
  pa <- pools_anchored(inp, 20, an, max_new_per_pool = mn,
                       reference = FALSE, P = 400)
  if (sum(pa$acc_plus %in% an) < 20L - mn) anch_ok <- FALSE
  if (sum(pa$acc_minus %in% an) < 20L - mn) anch_ok <- FALSE
  if (pa$n_new != sum(!c(pa$acc_plus, pa$acc_minus) %in% an)) count_ok <- FALSE
  if (length(intersect(pa$acc_plus, pa$acc_minus)) > 0L) anch_ok <- FALSE
}
check(anch_ok,
      "every anchored solution keeps at least n - max_new anchors in each pool")
check(count_ok, "and n_new counts exactly the non-anchor members")

# --- exclude is absolute, and applied before everything else ---
#
# An accession the unconstrained optimum wants must disappear when excluded,
# and the score must not rise as a result.
drop1 <- po20$acc_plus[1]
pa_ex <- pools_anchored(inp, 20, cand_acc[seq_len(40)], exclude = drop1,
                        max_new_per_pool = 4, reference = FALSE, P = 400)
check(!drop1 %in% c(pa_ex$acc_plus, pa_ex$acc_minus),
      "an excluded accession never appears in the result")
check(!drop1 %in% pa_ex$ctx$acc,
      "and it is removed from the candidate set, not merely from the pools")
check(pa_ex$score <= po20$score + 1e-12,
      "excluding a wanted accession cannot raise the score")

# --- the cost of the constraint is reported ---
pa_ref <- pools_anchored(inp, 20, cand_acc[seq_len(40)],
                         max_new_per_pool = 4, P = 400)
check(is.finite(pa_ref$score_unconstrained) &&
        pa_ref$score <= pa_ref$score_unconstrained + 1e-12,
      "the anchored score is reported against the unconstrained optimum")
check(pa_ref$score_frac_of_unconstrained <= 1 + 1e-12,
      "and as a fraction of it, so the cost of anchoring is readable")

# --- a budget wider than the pool is refused ---
check_error(pools_anchored(inp, 20, cand_acc[seq_len(40)],
                           max_new_per_pool = 21, P = 400),
            "a substitution budget wider than the pool is refused")

# ------------------------------------------------------------
# 8. The parsimony stage
#
# Maximising the score alone is the WRONG specification: the greedy optimum
# spends its whole substitution allowance whenever a further swap gains
# anything at all, which is the inverse of the protocol's purpose when every
# new accession costs a trip to seed storage. Stage B hands slots back while
# the score holds up, and both answers are reported so the trade is visible.
# ------------------------------------------------------------

anch_sl <- cand_acc[seq_len(40)]
pa_greedy <- pools_anchored(inp, 20, anch_sl, max_new_per_pool = 8,
                            parsimony = FALSE, reference = FALSE, P = 400)
pa_parsi  <- pools_anchored(inp, 20, anch_sl, max_new_per_pool = 8,
                            parsimony = TRUE, score_slack_frac = 0.005,
                            reference = FALSE, P = 400)

check(pa_parsi$n_new <= pa_parsi$n_new_greedy,
      "the parsimony stage never asks for MORE new accessions than the greedy optimum")
check(pa_parsi$n_new_greedy == pa_greedy$n_new,
      "and it reports the greedy count it started from")
check(pa_parsi$score >= pa_parsi$score_greedy * (1 - 0.005) - 1e-12,
      "the pruned solution stays within the slack of the stage-A optimum")
check(pa_parsi$score <= pa_parsi$score_greedy + 1e-12,
      "and never exceeds it: pruning trades score for fewer bags to find")
check(pa_parsi$n_pruned == pa_parsi$n_new_greedy - pa_parsi$n_new,
      "n_pruned accounts exactly for the slots handed back")

# zero slack must hand nothing back: the stage is then a no-op by construction
pa_zero_slack <- pools_anchored(inp, 20, anch_sl, max_new_per_pool = 8,
                                score_slack_frac = 0, reference = FALSE,
                                P = 400)
check(pa_zero_slack$n_new == pa_zero_slack$n_new_greedy &&
        pa_zero_slack$n_pruned == 0L,
      "with zero slack the parsimony stage is a no-op")

# a wider slack can only help, never hurt, the number of accessions to find
pa_wide <- pools_anchored(inp, 20, anch_sl, max_new_per_pool = 8,
                          score_slack_frac = 0.05, reference = FALSE, P = 400)
check(pa_wide$n_new <= pa_parsi$n_new,
      "a wider slack asks for no more new accessions than a narrow one")

# the pruned members really are anchors, and the pools stay legal
check(sum(c(pa_parsi$acc_plus, pa_parsi$acc_minus) %in% anch_sl) >=
        sum(c(pa_greedy$acc_plus, pa_greedy$acc_minus) %in% anch_sl),
      "pruning raises the number of anchors used")
check(sum(pa_parsi$acc_plus %in% anch_sl) >= 12L &&
        sum(pa_parsi$acc_minus %in% anch_sl) >= 12L &&
        length(intersect(pa_parsi$acc_plus, pa_parsi$acc_minus)) == 0L,
      "and cannot break the anchor constraint or disjointness")

# ------------------------------------------------------------
# 8b. The marginal-value curve
#
# The output that turns "is another search round worth it?" into a column.
# ------------------------------------------------------------

curve <- pool_search_curve(inp, 20, anch_sl, lambda = 0.76,
                           interaction_frac = 0.5, P = 400,
                           max_new_values = 0:8)
check(nrow(curve) == 9L, "the curve has one row per allowance asked for")
check(all(curve$max_new_min == pool_anchor_feasibility(20L, 40L, 0L)$max_new_min),
      "it records where the feasible range starts")
check(all(curve$feasible == (0:8 >= curve$max_new_min[1])),
      "and flags the infeasible allowances rather than dropping them")
feas <- dplyr::filter(curve, feasible)
check(all(is.na(dplyr::filter(curve, !feasible)$score)),
      "infeasible rows carry no score")
check(all(feas$n_new <= 2 * feas$max_new_per_pool),
      "no row asks for more new accessions than its allowance permits")
check(all(feas$score <= feas$score_unconstrained + 1e-12),
      "no anchored solution beats the unconstrained optimum")
check(all(feas$delta_power_vs_unconstrained <= 1e-12),
      "so the power delta against the ceiling is never positive")
check(all(diff(cummax(feas$score)) >= -1e-12),
      "the best score attainable is non-decreasing in the allowance")

# ------------------------------------------------------------
# 9. What the seed store says, and what the fit can use
#
# A round's input is typed by a person reading bag labels, and two things
# stand between that and the search: the bag carries an ORIGINAL germplasm
# name while the fit knows a synthetic analysis name, and a found accession
# may not be eligible at all.
#
# The accounting must be exhaustive. An accession the user believes they
# reported as found, which quietly never reached the search, is the worst
# failure this step can have -- so every input name gets exactly one row.
# ------------------------------------------------------------

# Injectable fixtures, so none of this needs output/.
map_fix <- tibble::tibble(
  germplasmName = c("entryA", "entryB", "entryC", "acc001"),
  analysis_name = c("acc005", "acc005", "acc006", "acc001"))
collapsed_fix <- tibble::tibble(
  species = "test", analysis_name = c("acc005", "acc006"),
  n_entries = c(2L, 1L), entries = c("entryA; entryB", "entryC"))

# --- the three routes ---
rsn <- resolve_seed_names(c("acc010", "entryA", "entryC", "nope"),
                          known = inp$accessions$acc, map = map_fix)
check(nrow(rsn) == 4L, "one row per input name")
check(rsn$route[rsn$input_name == "acc010"] == "identity" &&
        rsn$resolved_acc[rsn$input_name == "acc010"] == "acc010",
      "a name the fit knows is taken as is")
check(rsn$route[rsn$input_name == "entryA"] == "member_entry" &&
        rsn$resolved_acc[rsn$input_name == "entryA"] == "acc005",
      "a bag labelled with a pooled entry resolves to its analysis name")
check(rsn$route[rsn$input_name == "nope"] == "unresolved" &&
        is.na(rsn$resolved_acc[rsn$input_name == "nope"]),
      "an unknown name resolves to nothing rather than being dropped")

# A name that is BOTH an analysis name and a germplasm name must be taken as
# is, not remapped -- the `<seed_parent>` curation case, where the analysis
# name is a real line name.
rsn_both <- resolve_seed_names("acc001", known = inp$accessions$acc,
                               map = map_fix)
check(rsn_both$route == "identity" && rsn_both$resolved_acc == "acc001",
      "a name that is both an analysis name and an entry is taken as is")

# duplicates collapse rather than multiplying rows
check(nrow(resolve_seed_names(c("acc010", "acc010", "entryA"),
                              known = inp$accessions$acc, map = map_fix)) == 2L,
      "repeated input names are reported once")

check_error(resolve_seed_names("acc010", known = inp$accessions$acc,
                               map = tibble::tibble(x = 1)),
            "a map without the required columns is refused")

# --- one planted case per bucket, and the accounting is exhaustive ---
# One fixture used for BOTH the classification and the search below, so the
# two agree about which accessions exist and which are eligible. Classifying
# against one population and searching another is a mistake this test made
# once already.
inp_cls <- inp
inp_cls$accessions$eligible[2] <- FALSE   # ineligible on reliability
inp_cls$accessions$rel_As_i[2] <- 0.1
inelig  <- inp_cls$accessions$acc[2]

pc_cls   <- .pool_candidates(inp_cls, 0.5)
cand_cls <- sort(pc_cls$cand$acc)
usable   <- cand_cls[1]
below_pr <- dplyr::filter(inp_cls$accessions, eligible, Pr < pc_cls$pr_min)$acc[1]

probe <- c(usable, below_pr, inelig, "entryA", "typo999")
cls <- classify_search_input(inp_cls, probe, map = map_fix)

check(nrow(cls) == length(unique(probe)) && !anyDuplicated(cls$input_name),
      "every input name is accounted for exactly once")
check(sum(table(cls$bucket)) == length(unique(probe)),
      "the bucket counts sum to the number of inputs")
check(cls$bucket[cls$input_name == usable] == "in_candidates",
      "a usable accession lands in in_candidates")
check(cls$bucket[cls$input_name == below_pr] ==
        "eligible_but_below_Pr_threshold",
      "an eligible accession below the Pr threshold is bucketed, not admitted")
check(cls$bucket[cls$input_name == inelig] == "ineligible",
      "an ineligible accession is bucketed")
check(cls$bucket[cls$input_name == "typo999"] == "not_in_accession_table",
      "an unknown name is bucketed as a name problem")
check(grepl("below the threshold", cls$reason[cls$input_name == below_pr]) &&
        cls$pr_shortfall[cls$input_name == below_pr] > 0,
      "and the Pr shortfall is reported in g/m2, so relaxing stays deliberate")
check(grepl("reliability", cls$reason[cls$input_name == inelig]),
      "the ineligibility reason names which floor was missed")

# --- THE CONSEQUENCE: non-candidates cannot act as anchors ---
#
# The classifier's verdict has to bind. A found-but-unusable accession that
# still counted toward the anchor total would make an infeasible problem look
# feasible, which is the failure the bucket accounting exists to prevent.
unusable <- cls$input_name[cls$bucket != "in_candidates"]
pa_mixed <- pools_anchored(inp_cls, 20, c(cand_cls[seq_len(36)], unusable),
                           max_new_per_pool = 2, reference = FALSE, P = 400)
check(pa_mixed$n_anchors_usable == 36L,
      "only in_candidates accessions count toward the anchor total")
check(!any(unusable %in% c(pa_mixed$acc_plus, pa_mixed$acc_minus)),
      "and none of them can reach the pools, by name or by resolution")
check(!any(cls$resolved_acc[cls$bucket != "in_candidates"] %in%
             c(pa_mixed$acc_plus, pa_mixed$acc_minus), na.rm = TRUE),
      "nor can the accessions they resolve to")

# --- the forward direction: which entries satisfy a proposal ---
ce <- collapsed_entries(c("acc005", "acc010"), "test", collapsed = collapsed_fix)
check(ce$is_collapsed[ce$acc == "acc005"] &&
        ce$n_entries[ce$acc == "acc005"] == 2L &&
        ce$entries[ce$acc == "acc005"] == "entryA; entryB",
      "a collapsed proposal says which member entries would satisfy it")
check(!ce$is_collapsed[ce$acc == "acc010"] &&
        ce$n_entries[ce$acc == "acc010"] == 1L &&
        ce$entries[ce$acc == "acc010"] == "acc010",
      "a plain accession reports itself, so a caller can join unconditionally")

# ------------------------------------------------------------
# 1f. THE ROUND PROTOCOL
#
# The ledger is hand-typed primary data, so the properties worth asserting are
# the ones that protect what a person wrote: a typo must not be read as "not
# found", a re-run must not advance the round or rewrite an observation, and a
# status must not disagree with the exclusions derived from it.
# ------------------------------------------------------------

led_file <- tempfile(fileext = ".csv")
on.exit(unlink(led_file), add = TRUE)

check(identical(read_seed_ledger(tempfile(fileext = ".csv")),
                .empty_seed_ledger()),
      "a missing ledger reads as an empty one, so round 1 needs no setup")

# --- a typo in `status` must stop, not be read as "not found" ---
readr::write_csv(tibble::tibble(species = "test", acc = "acc001",
                                round_proposed = 1L, pool_proposed = "As+",
                                status = "fuond"), led_file)
check_error(read_seed_ledger(led_file),
            "an unrecognised status is refused rather than treated as missing")

# --- a duplicated accession would double-count one bag ---
readr::write_csv(tibble::tibble(species = "test", acc = c("acc001", "acc001"),
                                round_proposed = 1L, pool_proposed = "As+",
                                status = c("found", "missing")), led_file)
check_error(read_seed_ledger(led_file),
            "one accession listed twice is refused")

# --- blank status means unsearched, and case does not matter ---
readr::write_csv(tibble::tibble(species = "test",
                                acc = c("acc001", "acc002", "acc003"),
                                round_proposed = 1L, pool_proposed = "As+",
                                status = c("Found", NA, "MISSING_RETRY")),
                 led_file)
led <- read_seed_ledger(led_file)
check(identical(led$status, c("found", "unsearched", "missing_retry")),
      "status is case-insensitive and a blank reads as unsearched")

# --- the three statuses mean three different things to the search ---
st <- ledger_state(
  tibble::tibble(species = "test", acc = paste0("acc", 1:4),
                 round_proposed = 1L, pool_proposed = "As+",
                 status = c("found", "missing", "missing_retry",
                            "unsearched")),
  "test")
check(identical(st$anchors, "acc1"), "only `found` anchors the search")
check(identical(st$exclude, "acc2"), "only `missing` is excluded")
# missing_retry is the one that is easy to get wrong in either direction: it
# must not anchor (we do not have it) and must not exclude (we may yet).
check(!("acc3" %in% c(st$anchors, st$exclude)),
      "missing_retry neither anchors nor excludes")

# --- the round is open while anything it asked for is unsearched ---
check(seed_ledger_round(.empty_seed_ledger()) == 1L,
      "an empty ledger is round 1")
l1 <- tibble::tibble(species = "test", acc = paste0("acc", 1:3),
                     round_proposed = 1L, pool_proposed = "As+",
                     status = "unsearched")
check(seed_ledger_round(l1) == 1L,
      "re-running before anything is searched stays in the same round")
check(seed_ledger_round(dplyr::mutate(l1, status = "found")) == 2L,
      "the round advances once the request is answered")
l2 <- dplyr::bind_rows(dplyr::mutate(l1, status = "found"),
                       tibble::tibble(species = "test", acc = "acc9",
                                      round_proposed = 2L,
                                      pool_proposed = "As-",
                                      status = "unsearched"))
check(seed_ledger_round(l2) == 2L,
      "an outstanding round-2 request keeps the count at 2 across re-runs")

# --- append-only: what a person typed comes back untouched ---
typed <- tibble::tibble(
  species = "test", acc = c("acc001", "acc002"), round_proposed = 1L,
  pool_proposed = c("As+", "As-"), status = c("found", "missing"),
  searched_on = "2026-11-01", sowable_entry = c("entryA", NA),
  note = c("bottom shelf", NA))
new_pools <- tibble::tibble(species = "test",
                            acc = c("acc001", "acc003"),
                            pool = c("As-", "As+"))
upd <- update_seed_ledger(typed, new_pools, round = 2L)
back <- dplyr::filter(upd$ledger, acc %in% typed$acc) |>
  dplyr::arrange(acc)
check(identical(as.data.frame(back), as.data.frame(dplyr::arrange(typed, acc))),
      "a searched row is returned byte for byte, pool_proposed included")
check(nrow(upd$added) == 1L && upd$added$acc == "acc003" &&
        upd$added$round_proposed == 2L &&
        upd$added$status == "unsearched",
      "only the newly proposed accession is appended")

# --- an unsearched row follows the optimiser between pools ---
unsearched <- tibble::tibble(
  species = "test", acc = "acc001", round_proposed = 1L,
  pool_proposed = "As+", status = "unsearched", searched_on = NA_character_,
  sowable_entry = NA_character_, note = NA_character_)
u2 <- update_seed_ledger(unsearched,
                         tibble::tibble(species = "test", acc = "acc001",
                                        pool = "As-"), round = 2L)
check(u2$ledger$pool_proposed == "As-" && u2$n_repooled == 1L,
      "an unsearched row is re-pooled, so the request stays current")
check(u2$ledger$round_proposed == 1L,
      "round_proposed stays at the round it was FIRST asked for")

# ------------------------------------------------------------
# The search request
# ------------------------------------------------------------

res_u <- pools_optimal(inp, 10, pr_quantile = 0.5, tol = Inf)

# --- every open slot is listed, and only open slots ---
alt <- pool_slot_alternates(res_u, n_alternates = 2L)
check(dplyr::n_distinct(alt$slot) == 20L,
      "with no anchors every one of the 2n slots needs a bag")
check(sum(alt$is_first_choice) == 20L &&
        all(alt$acc[alt$is_first_choice] %in%
              c(res_u$acc_plus, res_u$acc_minus)),
      "each slot's first choice is the accession the search actually chose")

# --- THE PROPERTY THE GREEDY ASSIGNMENT EXISTS FOR ---
# Scored slot by slot the best replacement barely depends on which member it
# replaces, so the naive version offers one accession to many slots and a
# searcher who finds it can still only use it once.
check(!anyDuplicated(alt$acc),
      "no accession is offered to two slots, so the lists cannot collide")
check(all(alt$score_cost[!alt$is_first_choice] >= 0) &&
        all(alt$score_cost[alt$is_first_choice] == 0),
      "score_cost is the non-negative drop in t from taking an alternate")

# --- an alternate really does score what it claims ---
one <- dplyr::filter(alt, !is_first_choice) |> dplyr::slice(1)
sl  <- dplyr::filter(alt, slot == one$slot, is_first_choice)
side <- if (sl$pool == "As+") "plus" else "minus"
memb <- res_u[[side]]
memb[match(sl$acc, res_u$ctx$acc[memb])] <- match(one$acc, res_u$ctx$acc)
oracle <- if (side == "plus") pool_score(memb, res_u$minus, res_u$ctx) else
  pool_score(res_u$plus, memb, res_u$ctx)
check_near(one$score, oracle$score, tol = 1e-10,
           "an alternate's score equals the from-scratch score of that swap")

# --- ranks are ordered by score within a slot ---
check(all(alt |> dplyr::group_by(slot) |>
            dplyr::summarise(ok = !is.unsorted(rev(score)),
                             .groups = "drop") |> dplyr::pull(ok)),
      "ranks within a slot run from best to worst")

# --- a slot held by an anchor needs no search ---
anchors_all <- c(res_u$acc_plus, res_u$acc_minus)
res_a <- pools_anchored(inp, 10, anchors_all, max_new_per_pool = 3L,
                        tol = Inf, reference = FALSE)
check(nrow(pool_slot_alternates(res_a)) == 0L,
      "when every member is in hand there is nothing to go and find")

# THE END-TO-END ORACLE: every proposed accession found must reproduce the
# unconstrained proposal exactly and ask for nothing.
check(setequal(c(res_a$acc_plus, res_a$acc_minus), anchors_all) &&
        res_a$n_new == 0L,
      "a round in which everything was found returns the same pools, unchanged")

req_empty <- search_request(res_a, inp, round = 2L)
check(nrow(req_empty) == 0L, "and an empty search request")
md <- search_request_md(req_empty, vintage = "2026-11-01",
                        ledger_path = "x.csv", round = 2L)
check(any(grepl("round 2", md)) && any(grepl("already in hand", md)),
      "the empty request still records which round it closed")

# --- the request carries what a person at a freezer needs ---
req <- search_request(res_u, inp, round = 1L, n_alternates = 2L,
                      collapsed = collapsed_fix)
check(nrow(req) == nrow(alt) && all(req$species == "test") &&
        all(req$round == 1L),
      "one request row per offered accession")
check(!any(is.na(req$Pr)) && !any(is.na(req$As)),
      "every offered accession is described by the fit's own numbers")
check(all(req$n_entries[req$acc %in% collapsed_fix$analysis_name] > 1L) ||
        !any(req$acc %in% collapsed_fix$analysis_name),
      "a collapsed name in the request carries its member entries")
check(all(req$sowable_entries[!req$acc %in% collapsed_fix$analysis_name] ==
            req$acc[!req$acc %in% collapsed_fix$analysis_name]),
      "a plain accession is sowable as itself")

# --- the anchored request asks for exactly the new accessions ---
keep <- utils::head(anchors_all, 14L)
res_p <- pools_anchored(inp, 10, keep, max_new_per_pool = 3L, tol = Inf,
                        reference = FALSE)
req_p <- search_request(res_p, inp, round = 2L)
firsts <- dplyr::filter(req_p, is_first_choice)
check(nrow(firsts) == res_p$n_new &&
        setequal(firsts$acc, res_p$acc_new),
      "the request's first choices are exactly the accessions reported as new")
check(all(!firsts$acc %in% keep),
      "nothing already in hand is sent back to the freezer")

# --- a confirmed-missing accession is never offered, not even as an alternate ---
gone <- res_u$acc_plus[1]
res_x <- pools_anchored(inp, 10, setdiff(anchors_all, gone),
                        exclude = gone, max_new_per_pool = 3L, tol = Inf,
                        reference = FALSE)
check(!gone %in% search_request(res_x, inp, round = 2L)$acc,
      "an accession recorded as missing cannot come back as an alternate")

# ------------------------------------------------------------
# Is the store itself biased in As or Pr?
#
# Reported every round, never acted on: if storage is short of (say) old
# low-As lines, anchoring entrenches that, which is a caveat on the As- pool's
# interpretation rather than a bug in the search.
# ------------------------------------------------------------

led_bal <- tibble::tibble(
  species = "test", acc = inp$accessions$acc[1:20], round_proposed = 1L,
  pool_proposed = "As+",
  status = rep(c("found", "missing", "missing_retry", "unsearched"), each = 5))
cls_bal <- classify_search_input(inp, led_bal$acc, map = map_fix)
bal <- found_missing_balance(led_bal, cls_bal)
check(nrow(bal) == 2L && setequal(bal$variable, c("As", "Pr")),
      "the balance is reported for both effects")
check(all(bal$n_found == 5L) && all(bal$n_missing == 5L),
      "missing_retry is left out: an unfinished search is not evidence of absence")
check_near(bal$mean_found[bal$variable == "As"],
           mean(inp$accessions$As[1:5]), tol = 1e-10,
           "the means come from the fit, not from the ledger")

# ------------------------------------------------------------
# 1g. READING THE ADOPTED POOLS BACK
#
# validate_power.R's estimand section reports the pools that will actually be
# sown, while power_grid() re-derives membership because it sweeps n. The risk
# this guards is the two describing different accessions, which is exactly
# what happened before search-v3 made the adopted membership differ from
# build_pools()'s.
# ------------------------------------------------------------

pools_file <- tempfile(fileext = ".csv")
on.exit(unlink(pools_file), add = TRUE)

bp10 <- build_pools(inp, 10, pr_quantile = 0.5)
readr::write_csv(bp10$pools, pools_file)
pfc <- pools_from_csv(inp, pools_file, n = 10)

# The oracle: reading back the pools build_pools() just wrote must reproduce
# every number it reported. That is what makes .pool_stats() one definition
# rather than two that agree today.
check(setequal(pfc$pools$acc, bp10$pools$acc) &&
        identical(pfc$pools$pool, bp10$pools$pool),
      "the adopted membership is read back unchanged")
for (f in c("dAs", "dPr", "sigma2_within", "sigma2_within_plus",
            "sigma2_within_minus", "PEV_plus", "PEV_minus", "Sxx_w",
            "b_As_on_Pr", "var_As_within_pools", "mean_As_plus",
            "mean_Pr_minus", "pr_min")) {
  check_near(pfc[[f]], bp10[[f]], tol = 1e-12,
             paste0("pools_from_csv() reproduces ", f))
}
check(pfc$n == 10L && pfc$n_candidates == bp10$n_candidates,
      "and the pool size and candidate count")

# --- the numbers come from `inp`, not from the file ---
# A pools.csv is a membership list. If the fit moves, the power must move with
# it, so an edited As column in the file must change nothing.
readr::write_csv(dplyr::mutate(bp10$pools, As = As + 1000), pools_file)
check_near(pools_from_csv(inp, pools_file, n = 10)$dAs, bp10$dAs, tol = 1e-12,
           "effect sizes are recomputed from the fit, not read from the file")

# --- a pool member the fit no longer carries is a human decision ---
readr::write_csv(dplyr::mutate(bp10$pools,
                               acc = dplyr::if_else(dplyr::row_number() == 1L,
                                                    "no-such-accession", acc)),
                 pools_file)
check_error(pools_from_csv(inp, pools_file, n = 10),
            "an adopted accession missing from the fit stops the run")

# --- the wrong pool size is refused rather than silently reported ---
readr::write_csv(bp10$pools, pools_file)
check_error(pools_from_csv(inp, pools_file, n = 12),
            "a pools.csv at a different n than asked for is refused")

# --- the fresh-clone fallback ---
fb <- pools_from_csv(inp, tempfile(fileext = ".csv"), n = 10)
check(setequal(fb$pools$acc, bp10$pools$acc),
      "with no pools.csv it falls back to re-deriving them")

finish("validation tests")
