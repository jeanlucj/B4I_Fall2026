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
  list(
    species = "test",
    accessions = tibble::tibble(
      acc = sprintf("acc%03d", seq_len(n)),
      Pr = as.vector(scale(z1)) * 10,
      As = as.vector(scale(z2)) * 8,
      n_partners = sample(3:10, n, replace = TRUE),
      n_plots = sample(4:12, n, replace = TRUE),
      eligible = TRUE),
    PEV_As = pev, sigma2_e = sigma2_e, sigma2_Pr = 100,
    rel_As = 0.3, response = "partner yield", n_trials = 6L, n_plots = 2000L
  )
}
inp <- make_inputs()

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
  # theta, so the bisection lands on a jump and the achieved gap can be
  # whatever the jump leaves. What is assertable is that it never does worse
  # than selecting on As alone, and that it sits at the zero crossing.
  cand_n <- dplyr::filter(inp$accessions, eligible, Pr >= bp$pr_min)
  check(abs(bp$dPr) <= abs(.pools_at(cand_n, n, 0)$dPr) + 1e-8,
        sprintf("constraining never widens the producer gap at n = %d", n))
  check(.pools_at(cand_n, n, bp$theta * 0.5)$dPr <= 0 &&
          .pools_at(cand_n, n, bp$theta * 2)$dPr >= 0,
        sprintf("the solution brackets the zero crossing at n = %d", n))
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

# What the bisection assumes about theta. dPr is a step function and is NOT
# locally monotone -- it reverses between adjacent steps -- so the bisection
# finds *a* zero crossing rather than the global best |dPr|. The trend it relies
# on is what is assertable, and pinning it here is also a record that the
# stronger claim is false.
cand_m <- dplyr::filter(inp$accessions, eligible,
                        Pr >= stats::quantile(Pr, 0.5, names = FALSE))
theta_seq <- seq(0, 5, by = 0.25)
dPr_by_theta <- purrr::map_dbl(theta_seq, \(th) .pools_at(cand_m, 20, th)$dPr)
check(dPr_by_theta[1] < 0,
      "selecting on As alone leaves the As+ pool short on Pr")
check(dPr_by_theta[length(dPr_by_theta)] > 0,
      "and a large theta overshoots, so a crossing exists to bisect on")
check(stats::cor(theta_seq, dPr_by_theta) > 0.8,
      "dPr rises with theta in trend, which is what the bisection needs")
check(any(diff(dPr_by_theta) < 0),
      "but NOT monotonically -- so the crossing found is local, not optimal")

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

# only eligible accessions may be selected
inp_partial <- inp
inp_partial$accessions$eligible <- inp_partial$accessions$n_partners >= 8
bp_e <- suppressWarnings(build_pools(inp_partial, 8))
check(all(bp_e$pools$n_partners >= 8),
      "ineligible accessions are never selected")

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
# 3. contrast_power: the shape of the formula
#
# The claim the design document rests on is that a large share of the standard
# error cannot be bought with plots. That is a property of the formula, so it is
# testable without any field data.
# ------------------------------------------------------------

pw <- function(...) contrast_power(dAs = 15, sigma2_within = 100, n = 20,
                                  sigma2_e = 500, P = 400, ...)
p1 <- pw()
check(p1$power > 0 && p1$power < 1, "power is a probability")
check(p1$SE > p1$SE_floor, "the SE exceeds its infinite-plot floor")
check(p1$pct_SE_irreducible > 0 && p1$pct_SE_irreducible <= 100,
      "the irreducible share is a percentage")

# more plots must help, but only up to the floor
se_of <- function(P) contrast_power(15, 100, 20, 500, P)$SE
check(se_of(300) > se_of(500), "more plots reduce the SE")
check(se_of(1e9) > p1$SE_floor * 0.999 && se_of(1e9) < p1$SE_floor * 1.001,
      "and with enough plots the SE converges on the floor, not on zero")

# more accessions reduce the part plots cannot touch
check(contrast_power(15, 100, 40, 500, 400)$SE_floor <
        contrast_power(15, 100, 20, 500, 400)$SE_floor,
      "more accessions per pool lowers the floor")

# a one-sided test is more powerful than two-sided for a true positive effect
check(contrast_power(15, 100, 20, 500, 400, sided = 1)$power >
        contrast_power(15, 100, 20, 500, 400, sided = 2)$power,
      "one-sided beats two-sided on a positive contrast")

# attenuation and interaction both cost power
check(contrast_power(15, 100, 20, 500, 400, lambda = 0.6)$power <
        contrast_power(15, 100, 20, 500, 400, lambda = 1.0)$power,
      "attenuation costs power")
check(contrast_power(15, 100, 20, 500, 400, interaction_frac = 0.5)$power <
        contrast_power(15, 100, 20, 500, 400, interaction_frac = 0)$power,
      "a pool x location interaction costs power")

# a zero contrast must give power at the nominal level, not zero and not one
p0 <- contrast_power(0, 100, 20, 500, 400, sided = 1, alpha = 0.05)
check_near(p0$power, 0.05, tol = 1e-6,
           "a zero contrast gives power equal to alpha")

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
