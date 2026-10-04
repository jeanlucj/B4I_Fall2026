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
    rel_As = 0.3, response = "partner yield", n_trials = 8L, n_plots = 3181L
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

finish("validation tests")
