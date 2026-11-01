# ============================================================
# BUILD THE As+ / As- POOLS FOR THE VALIDATION TRIAL
#
# Two pools per species, both drawn from high-producer accessions and
# contrasted only on their associate effects.  Finding accessions with a LOW
# producer effect is not interesting -- poor lines are easy to find -- so the
# question the validation asks is whether, among good producers, we can tell
# apart the ones that help their partner from the ones that hurt it.
#
# The rule is a constrained optimisation: maximise the Pr-adjusted contrast
# mean(As+) - mean(As-) - b * dPr subject to the pools being disjoint, both
# pool means of Pr sitting above a threshold, and the two pool means of Pr
# being equal to within a tolerance.  See build_pools() and THE POOL SEARCH in
# code/validation_functions.R for why the obvious rule -- top n on (Pr + As)
# against top n on (Pr - As) -- does not work.
#
# TWO MODES, BECAUSE SEED HAS TO BE FOUND BEFORE IT CAN BE SOWN.
#
#   unconstrained  the optimum over all candidates.  This is round 1: the
#                  proposal you take to the seed store the first time.
#   anchored       hold the accessions already found in storage, allow at most
#                  --max-new substitutions per pool, and ask for a short list
#                  of replacements for the ones that were missing.
#
# The protocol between them is data/seed_search_ledger.csv, which this script
# appends proposals to and a person fills the `status` column of.  --search
# auto picks anchored as soon as that ledger has a `found` row, and says which
# it chose and why: a refresh must never switch mode silently.
#
# Usage:
#   Rscript code/validate_pool_selection.R                     # auto
#   Rscript code/validate_pool_selection.R --search unconstrained
#   Rscript code/validate_pool_selection.R --max-new 5 --round 2
#   Rscript code/validate_pool_selection.R --no-curve          # skip the sweep
#
# Pool membership is provisional until seed must be ordered.  Re-run this after
# every new trial is loaded; the diff against the previous vintage is written
# out, and churn in it is a result in its own right.
#
# Outputs, all under output/validation/<vintage>/:
#   pools.csv                 both species, both pools -- the adopted design
#   pool_summary.csv          what each pool achieved, and at what constraint
#   search_request.csv/.md    one row per bag to find, with ranked alternates
#   search_curve.csv          what each extra substitution would buy
#   search_input_buckets.csv  every name the ledger supplied, accounted for
#   found_missing_balance.csv is the store biased in As or Pr?
#   pool_diff.csv             vs the previous vintage
#   pools_Pr_vs_As.png
# and data/seed_search_ledger.csv is appended to in place.
# ============================================================

library(tidyverse)

here::i_am("code/validate_pool_selection.R")

source(here::here("code", "validation_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

out_root <- here::here("output", "validation")

# Vintage folder. Dated, so a refresh never overwrites the pools a decision
# was taken on.
vintage <- format(Sys.Date(), "%Y-%m-%d")

# The settings that the three validate_* scripts must agree about live in
# VALIDATION_DEFAULTS, in code/validation_functions.R: power computed under a
# different partner filter would describe different pools than the ones
# selected here. Override a default only deliberately.
n_per_pool   <- validation_setting("n_per_pool")
pr_quantile  <- validation_setting("pr_quantile")
pr_tolerance <- validation_setting("pr_tolerance")
min_partners <- validation_setting("min_partners")
n_locations  <- validation_setting("n_locations")
plots_per_location <- validation_setting("plots_per_location")
score_slack_frac   <- validation_setting("score_slack_frac")

# The number of plots the reported power is computed at. ONE value for both
# species, because a plot is an oat x pea combination and the same 400 plots
# serve both -- which is why the drivers pass P explicitly rather than taking
# the search functions' 4 * n * n_loc convenience default (that default is the
# oat geometry and overstates pea's plot count by half).
P_total <- plots_per_location * n_locations

search_mode_req <- arg_value("--search", "auto")
if (!search_mode_req %in% c("auto", "unconstrained", "anchored")) {
  stop("--search must be auto, unconstrained or anchored; got ",
       search_mode_req, call. = FALSE)
}
ledger_path <- arg_value("--ledger", here::here("data", "seed_search_ledger.csv"))
max_new_per_pool <- as.integer(arg_value("--max-new",
                                         validation_setting("max_new_per_pool")))
want_curve <- !("--no-curve" %in% args)

# ============================================================
# Driver
# ============================================================

out_dir <- file.path(out_root, vintage)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# min_partners is now the REFERENCE count from which the reliability threshold
# is derived (the first quartile of rel_As among accessions sitting exactly at
# it), not the threshold itself. partner_floor is the hard identifiability
# guard that still applies.
inputs <- validation_inputs(min_partners = min_partners)

cat("\n=== Vintage ===\n")
vint <- validation_vintage(inputs)
print(as.data.frame(vint), row.names = FALSE, digits = 4)

cat("\nReliability is what limits this experiment. PEV enters the contrast's\n",
    "standard error divided by the pool size, and nothing about the field\n",
    "design can reduce it -- only better estimates can.\n", sep = "")

# ------------------------------------------------------------
# lambda, for the power REPORTED beside the proposal
#
# Not for the selection: the search maximises a lambda-free score and the
# chosen membership is invariant to lambda over the whole plausible range, so
# a stale lambda misstates the printed power and provably not the accessions.
# Said out loud below, because "the power moved" and "the design moved" are
# very different pieces of news.
# ------------------------------------------------------------

cv_est <- crossval_estimates(out_root, vintage)
if (is.null(cv_est)) {
  warning("no crossval_summary.csv anywhere under ", out_root,
          " -- reporting power at lambda = ", FALLBACK$lambda,
          ", interaction_frac = ", FALLBACK$interaction_frac,
          ". These are GUESSES; the proposed accessions do not depend on ",
          "them, but the power printed beside them does.", call. = FALSE)
  est_tbl <- tibble::tibble(species = names(inputs),
                            lambda = FALLBACK$lambda, lambda_se = NA_real_,
                            interaction_frac = FALLBACK$interaction_frac)
  lambda_from <- "FALLBACK GUESS"
} else {
  est_tbl     <- cv_est$est
  lambda_from <- cv_est$from
}
est_for <- function(sp) {
  r <- dplyr::filter(est_tbl, species == sp)
  if (nrow(r) == 0) {
    list(lambda = FALLBACK$lambda, lambda_se = NA_real_,
         interaction_frac = FALLBACK$interaction_frac)
  } else {
    list(lambda = r$lambda[1], lambda_se = r$lambda_se[1],
         interaction_frac = r$interaction_frac[1])
  }
}

cat("\n=== lambda, for the reported power only ===\n")
cat("    source: ", lambda_from,
    if (identical(lambda_from, vintage)) " (this vintage)"
    else "  <- NOT this vintage",
    "\n    lambda scales the contrast, so it moves every power below. It does\n",
    "    not affect WHICH accessions are selected.\n", sep = "")
print(as.data.frame(est_tbl), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# What the seed store has told us so far
# ------------------------------------------------------------

ledger <- read_seed_ledger(ledger_path)
round_no <- as.integer(arg_value("--round", seed_ledger_round(ledger)))

cat("\n=== Seed-search ledger ===\n")
cat(basename(ledger_path), ": ", nrow(ledger), " row(s), ",
    sum(ledger$status != "unsearched"), " searched. Round ", round_no, ".\n",
    sep = "")
if (nrow(ledger) > 0) {
  print(as.data.frame(dplyr::count(ledger, species, status)), row.names = FALSE)
}

# Every name the ledger supplied, accounted for in exactly one bucket. The
# failure this guards against is an accession reported as found that quietly
# never reached the search.
buckets <- purrr::imap(inputs, \(inp, name) {
  nm <- dplyr::filter(ledger, species == name)$acc
  if (length(nm) == 0L) return(NULL)
  classify_search_input(inp, nm, pr_quantile = pr_quantile)
}) |> purrr::compact()

if (length(buckets) > 0) {
  buckets_tbl <- purrr::list_rbind(buckets)
  readr::write_csv(buckets_tbl, file.path(out_dir, "search_input_buckets.csv"))
  cat("\nnames the ledger supplied, by bucket:\n")
  print(as.data.frame(dplyr::count(buckets_tbl, species, bucket)),
        row.names = FALSE)
  unusable <- dplyr::filter(buckets_tbl, bucket != "in_candidates")
  if (nrow(unusable) > 0) {
    cat("\nNOT usable as anchors -- reported, never quietly admitted.",
        "\nThe Pr threshold is part of the hypothesis (the pools are",
        "\ncontrasted among GOOD producers), so lowering pr_quantile to let",
        "\none of these in has to be a deliberate, recorded act:\n")
    print(as.data.frame(dplyr::select(unusable, species, input_name,
                                      resolved_acc, bucket, reason)),
          row.names = FALSE)
  }
} else {
  buckets_tbl <- NULL
}

# Mode, decided per species and announced with its reason.
usable_anchors <- purrr::imap(inputs, \(inp, name) {
  st <- ledger_state(ledger, name)
  b  <- buckets[[name]]
  anch <- if (is.null(b)) character(0) else {
    dplyr::filter(b, input_name %in% st$anchors,
                  bucket == "in_candidates")$resolved_acc
  }
  # `exclude` applies to every resolvable missing name, candidate or not, and
  # is cumulative over all rounds: an accession confirmed absent is never
  # proposed again.
  excl <- if (is.null(b)) character(0) else {
    stats::na.omit(dplyr::filter(b, input_name %in% st$exclude)$resolved_acc)
  }
  list(anchors = unique(anch), exclude = unique(as.character(excl)),
       state = st)
})

modes <- purrr::imap_chr(inputs, \(inp, name) {
  n_anch <- length(usable_anchors[[name]]$anchors)
  if (search_mode_req != "auto") return(search_mode_req)
  if (n_anch > 0L) "anchored" else "unconstrained"
})

cat("\n=== Search mode ===\n")
for (name in names(inputs)) {
  n_anch <- length(usable_anchors[[name]]$anchors)
  cat(name, ": ", modes[[name]], "  (--search ", search_mode_req, "; ",
      n_anch, " usable anchor(s), ",
      length(usable_anchors[[name]]$exclude), " excluded)\n", sep = "")
}
if (any(modes == "anchored" &
        purrr::map_int(usable_anchors, \(u) length(u$anchors)) == 0L)) {
  stop("--search anchored was asked for, but a species has no usable anchor. ",
       "Nothing to anchor on: run --search unconstrained, record the results ",
       "in ", basename(ledger_path), ", and come back.", call. = FALSE)
}

# ------------------------------------------------------------
# What each extra substitution would buy
#
# Written BEFORE the pools are built, because the curve is computable even
# when the chosen allowance is infeasible -- and that is exactly the run in
# which you need it.
# ------------------------------------------------------------

curve <- NULL
if (want_curve && any(modes == "anchored")) {
  curve <- purrr::imap(inputs, \(inp, name) {
    if (modes[[name]] != "anchored") return(NULL)
    u <- usable_anchors[[name]]
    message("--- ", name, ": pricing the substitution allowance ---")
    pool_search_curve(inp, n_per_pool[[name]], u$anchors,
                      exclude = u$exclude, pr_quantile = pr_quantile,
                      tol = pr_tolerance, lambda = est_for(name)$lambda,
                      interaction_frac = est_for(name)$interaction_frac,
                      P = P_total, n_loc = n_locations,
                      score_slack_frac = score_slack_frac)
  }) |> purrr::compact() |> purrr::list_rbind()

  readr::write_csv(curve, file.path(out_dir, "search_curve.csv"))
  cat("\n=== What each extra substitution buys ===\n")
  cat("THE SMALLEST FEASIBLE ALLOWANCE IS USUALLY THE WORST PLACE TO SIT: at\n",
      "the minimum every usable anchor must be used, so the search cannot\n",
      "drop a badly-estimated one. One or two extra units normally recover it\n",
      "in full, for the same number of bags to find.\n", sep = "")
  print(as.data.frame(dplyr::select(curve, species, max_new_per_pool, feasible,
                                    n_new, n_new_greedy, power,
                                    power_unconstrained, dPr, at_min_budget)),
        row.names = FALSE, digits = 3)
}

# ------------------------------------------------------------
# Is the allowance even feasible?
#
# Checked for EVERY species before any is built, and by arithmetic rather
# than by letting the first search throw.  Otherwise a run that fails on oat
# never computes pea, and you learn the two answers one re-run at a time --
# with each answer being "go and find more seed", the slowest possible way to
# discover that you have to do it twice.
#
# Never auto-relaxed: widening the allowance spends somebody's afternoon, so
# it is reported and stopped on, exactly as agreed.
# ------------------------------------------------------------

infeasible <- purrr::imap(inputs, \(inp, name) {
  if (modes[[name]] != "anchored") return(NULL)
  n  <- n_per_pool[[name]]
  # A found accession that is not a candidate cannot be used, and
  # usable_anchors already restricts to bucket == "in_candidates"; `exclude`
  # cannot remove one, since an accession has a single status.
  fs <- pool_anchor_feasibility(n, length(usable_anchors[[name]]$anchors),
                                max_new_per_pool)
  if (fs$feasible) return(NULL)
  tibble::tibble(species = name, n_per_pool = n, n_usable_anchors = fs$n_usable,
                 max_new_per_pool = max_new_per_pool, anchors_required = fs$required,
                 shortfall = fs$shortfall, min_feasible_max_new = fs$max_new_min,
                 new_accessions_at_min = 2L * fs$max_new_min)
}) |> purrr::compact() |> purrr::list_rbind()

if (nrow(infeasible) > 0) {
  cat("\n=== INFEASIBLE at --max-new ", max_new_per_pool, " ===\n", sep = "")
  print(as.data.frame(infeasible), row.names = FALSE)
  cat("\nEach pool must hold at least n - max_new_per_pool anchors, and the two\n",
      "pools are disjoint, so the allowance needs 2 * (n - max_new) usable\n",
      "anchors. Two ways out, both yours to choose:\n",
      "  - find `shortfall` more of the accessions already proposed, or\n",
      "  - re-run with --max-new `min_feasible_max_new`, which asks for\n",
      "    `new_accessions_at_min` new bags instead.\n",
      "search_curve.csv prices the second option; note that the minimum\n",
      "feasible allowance is usually the worst place to sit, and one unit\n",
      "above it normally costs nothing extra to find.\n", sep = "")
  stop("anchor constraint infeasible for: ",
       paste(infeasible$species, collapse = ", "),
       " -- see the table above.", call. = FALSE)
}

# ------------------------------------------------------------
# Build
# ------------------------------------------------------------

built <- purrr::imap(inputs, \(inp, name) {
  n <- n_per_pool[[name]]
  u <- usable_anchors[[name]]
  message("\n--- ", name, ": ", n, " accessions per pool (", modes[[name]],
          ") ---")

  res <- if (modes[[name]] == "unconstrained") {
    pools_optimal(inp, n, pr_quantile = pr_quantile, tol = pr_tolerance,
                  P = P_total, n_loc = n_locations)
  } else {
    # Never auto-relaxed. An infeasible allowance is reported with the
    # smallest one that would work, and the run stops -- quietly widening the
    # allowance would spend somebody's afternoon without being asked.
    pools_anchored(inp, n, u$anchors, exclude = u$exclude,
                   max_new_per_pool = max_new_per_pool,
                   pr_quantile = pr_quantile, tol = pr_tolerance,
                   P = P_total, n_loc = n_locations,
                   score_slack_frac = score_slack_frac)
  }
  bp <- res$build_pools

  e <- inp$eligibility
  cat("\n=== ", toupper(name), " ===\n", sep = "")
  # The rule is printed from inp$eligibility rather than described here, so the
  # console, the figure subtitle and the generated document cannot drift apart.
  cat("eligible accessions: ", e$n_eligible, " of ", e$n_total, "  [rule: ",
      if (e$rule == "reliability")
        sprintf("rel_As >= %.3f (%s, from %d accessions at exactly %d partners) AND >= %d partners",
                e$rel_min, e$rel_min_source, e$n_ref, e$min_partners, e$partner_floor)
      else sprintf(">= %d partners", e$min_partners),
      "]\n", sep = "")
  cat("candidates above the Pr threshold (", round(bp$pr_min, 2), "): ",
      bp$n_candidates,
      if (length(u$exclude) > 0)
        paste0(" (", length(u$exclude), " excluded as confirmed missing)")
      else "", "\n", sep = "")
  cat("associate contrast  : ", round(res$dAs, 2), " g/m2 on ", inp$response,
      "\n", sep = "")
  cat("Pr-adjusted contrast: ", round(res$dAs_adj, 2),
      " g/m2  (b = ", round(res$ctx$b, 4), " As per unit Pr)\n", sep = "")
  cat("producer difference : ", round(res$dPr, 3), " g/m2 (tolerance ",
      pr_tolerance, ")", if (isTRUE(res$dpr_tol_relaxed))
        "  <- TOLERANCE RELAXED, see dpr_tol_relaxed" else "", "\n", sep = "")
  # The SELECTION score, not the trial's t: it is the noncentrality at
  # lambda = 1 with no pool x location interaction, which is what makes
  # membership independent of the cross-validation. The t the trial will
  # actually see is in pool_summary.csv beside power_pr_adj.
  cat("selection score     : ", round(res$score, 4),
      " (t at lambda = 1, no interaction)\n", sep = "")
  cat("  reached from theta = ", signif(res$theta, 4), " (", res$seed_kind,
      ") in ", res$n_moves, " move(s)",
      if (res$n_moves == 0L)
        " -- the theta rule was already optimal, which is the usual outcome"
      else "", "\n", sep = "")
  cat("mean Pr, As+ / As-  : ",
      round(mean(res$pools$Pr[res$pools$pool == "As+"]), 2), " / ",
      round(mean(res$pools$Pr[res$pools$pool == "As-"]), 2),
      "   (population mean ", round(mean(inp$accessions$Pr), 2), ")\n", sep = "")
  cat("within-pool SD of true effects: ",
      round(sqrt(mean(c(res$sigma2_within_plus, res$sigma2_within_minus))), 2),
      "\n", sep = "")
  # The PEV above is MEASURED from the posterior draws, per accession. The old
  # route backed it out of Var(true) = Var(BLUP) + E[PEV], an identity this fit
  # does not satisfy, and was 2.3x larger -- which inflated this SD and
  # understated power. Both are reported in vintage.csv as PEV_As (backed out)
  # and PEV_As_measured, with their ratio.
  cat("  (backed-out global PEV would have been ", round(inp$PEV_As, 1),
      ", a factor of ", round(inp$PEV_ratio_backed_out_to_measured, 2),
      " larger -- see vintage.csv)\n", sep = "")

  if (modes[[name]] == "anchored") {
    cat("anchors: ", res$n_anchors_used, " of ", res$n_anchors_usable,
        " usable used; ", res$n_new, " new accession(s) to find",
        if (res$n_new_greedy > res$n_new)
          paste0(" (the unpruned optimum asked for ", res$n_new_greedy,
                 " to gain ", round(res$score_greedy - res$score, 4), " in t)")
        else "", "\n", sep = "")
    cat("allowance: ", res$max_new_per_pool, " per pool (minimum feasible ",
        res$max_new_min, ")",
        if (isTRUE(res$at_min_budget))
          "  <- AT THE MINIMUM: every usable anchor must be kept, which costs power"
        else "", "\n", sep = "")
    cat("t against the unconstrained optimum: ", round(res$score, 4), " vs ",
        round(res$score_unconstrained, 4), " = ",
        round(100 * res$score_frac_of_unconstrained, 2), "%  (regret ",
        round(res$score_unconstrained - res$score, 4), " in t)\n", sep = "")
  }

  stopifnot(
    "pools are not disjoint" =
      length(intersect(res$pools$acc[res$pools$pool == "As+"],
                       res$pools$acc[res$pools$pool == "As-"])) == 0,
    # Both halves of the eligibility rule, asserted on the SELECTED accessions
    # rather than trusted from the candidate filter. The partner floor is the
    # one that matters: an accession grown with a single partner has its
    # producer and associate effects perfectly aliased, and the pools exist to
    # separate exactly those two.
    "a pool member falls below the reliability threshold" =
      (inp$eligibility$rule != "reliability") ||
        all(res$pools$rel_As_i >= inp$eligibility$rel_min),
    "a pool member falls below the partner floor" =
      all(res$pools$n_partners >= inp$eligibility$partner_floor),
    "a pool member has producer and associate effects aliased" =
      all(res$pools$n_partners >= 2L),
    # An accession confirmed absent from storage must not reappear.
    "a pool member was recorded as missing from seed storage" =
      !any(res$pools$acc %in% u$exclude)
  )
  res
})

pools <- purrr::map(built, "pools") |> purrr::list_rbind()

# ------------------------------------------------------------
# Power, and what stopping now costs
#
# THREE NUMBERS AGAINST PATH DEPENDENCE, printed every round including when
# they are zero, so a non-zero value is noticed rather than discovered:
#
#   regret                    the unconstrained optimum RECOMPUTED from this
#                             fit, minus what the anchors allow.  What the
#                             searching done so far has cost.
#   score_if_remaining_found  the same anchored search with the not-yet-
#                             searched accessions assumed found.  What
#                             stopping now costs -- the number that decides
#                             whether another round is worth it.
#   search_curve.csv          what one more substitution would buy.
# ------------------------------------------------------------

power_of <- function(res, inp, name) {
  e <- est_for(name)
  contrast_power(
    dAs = res$dAs,
    sigma2_within = c(res$sigma2_within_plus, res$sigma2_within_minus),
    n = res$ctx$n, sigma2_e = inp$sigma2_e, P = P_total,
    lambda = e$lambda, interaction_frac = e$interaction_frac,
    n_loc = n_locations, sided = 1,
    pr_adjust = list(dPr = res$dPr, b = res$ctx$b, Sxx = res$Sxx))
}

stop_cost <- purrr::imap(inputs, \(inp, name) {
  u <- usable_anchors[[name]]
  if (modes[[name]] != "anchored" || length(u$state$unsearched) == 0L) {
    return(tibble::tibble(species = name,
                          n_unsearched = length(u$state$unsearched),
                          score_if_remaining_found = built[[name]]$score))
  }
  b <- buckets[[name]]
  more <- dplyr::filter(b, input_name %in% u$state$unsearched,
                        bucket == "in_candidates")$resolved_acc
  hyp <- try(pools_anchored(inp, n_per_pool[[name]],
                            unique(c(u$anchors, more)), exclude = u$exclude,
                            max_new_per_pool = max_new_per_pool,
                            pr_quantile = pr_quantile, tol = pr_tolerance,
                            P = P_total, n_loc = n_locations,
                            score_slack_frac = score_slack_frac,
                            reference = FALSE), silent = TRUE)
  tibble::tibble(species = name, n_unsearched = length(u$state$unsearched),
                 score_if_remaining_found =
                   if (inherits(hyp, "try-error")) NA_real_ else hyp$score)
}) |> purrr::list_rbind()

pool_summary <- purrr::imap(built, \(res, name) {
  bp  <- res$build_pools
  pw  <- power_of(res, inputs[[name]], name)
  e   <- est_for(name)
  pl  <- res$pools$pool == "As+"
  anch_switch <- {
    prev <- dplyr::filter(ledger, species == name)
    j <- dplyr::inner_join(dplyr::select(res$pools, acc, pool),
                           dplyr::select(prev, acc, pool_proposed),
                           by = "acc")
    sum(!is.na(j$pool_proposed) & j$pool_proposed != j$pool)
  }
  ref_power <- if (modes[[name]] == "anchored" &&
                   !is.na(res$score_unconstrained)) {
    r <- pools_optimal(inputs[[name]], res$ctx$n, pr_quantile = pr_quantile,
                       tol = pr_tolerance, P = P_total, n_loc = n_locations)
    power_of(r, inputs[[name]], name)$power
  } else pw$power

  tibble::tibble(
    species = name, n_per_pool = res$ctx$n, n_candidates = bp$n_candidates,
    pr_threshold = bp$pr_min, theta = res$theta,
    pool_index = POOL_INDEX_VERSION,
    search_mode = res$search_mode, round = round_no,
    n_moves = res$n_moves, seed_kind = res$seed_kind,
    dAs = res$dAs, dPr = res$dPr, dAs_adj = res$dAs_adj,
    b_As_on_Pr = res$ctx$b,
    mean_As_plus = mean(res$pools$As[pl]),
    mean_As_minus = mean(res$pools$As[!pl]),
    mean_Pr_plus = mean(res$pools$Pr[pl]),
    mean_Pr_minus = mean(res$pools$Pr[!pl]),
    sigma2_within = mean(c(res$sigma2_within_plus, res$sigma2_within_minus)),
    sigma2_within_plus = res$sigma2_within_plus,
    sigma2_within_minus = res$sigma2_within_minus,
    PEV_plus = mean(res$pools$PEV_As_i[pl]),
    PEV_minus = mean(res$pools$PEV_As_i[!pl]),
    var_As_within_pools = mean(c(stats::var(res$pools$As[pl]),
                                 stats::var(res$pools$As[!pl]))),
    PEV = inputs[[name]]$PEV_As,
    PEV_measured = inputs[[name]]$PEV_As_measured,
    reliability_As = inputs[[name]]$rel_As,
    eligibility_rule = inputs[[name]]$eligibility$rule,
    rel_min = inputs[[name]]$eligibility$rel_min,
    # The constraint and what it cost.
    n_anchors_available = length(usable_anchors[[name]]$anchors),
    n_anchors_usable = res$n_anchors_usable,
    n_anchors_used = res$n_anchors_used %||% 0L,
    n_excluded = length(usable_anchors[[name]]$exclude),
    max_new_per_pool = res$max_new_per_pool %||% NA_integer_,
    max_new_min = res$max_new_min %||% NA_integer_,
    at_min_budget = res$at_min_budget %||% NA,
    n_new = res$n_new, n_new_greedy = res$n_new_greedy %||% res$n_new,
    n_pruned = res$n_pruned %||% 0L,
    anchor_pool_switches = anch_switch,
    dpr_tol_relaxed = res$dpr_tol_relaxed %||% FALSE,
    score = res$score,
    score_unconstrained = res$score_unconstrained %||% res$score,
    pct_of_unconstrained = 100 * (res$score_frac_of_unconstrained %||% 1),
    regret = (res$score_unconstrained %||% res$score) - res$score,
    # Power is REPORTED, not selected on: see the lambda block above.
    power_pr_adj = pw$power, power_unconstrained = ref_power,
    # Two t's, deliberately both here: `score` is the lambda-free selection
    # objective and `t_reported` is what the trial will see at the measured
    # lambda and interaction. They differ by about a factor of two, and
    # mistaking one for the other would misstate the experiment.
    t_reported = pw$t, SE = pw$SE, df = pw$df,
    lambda = e$lambda, lambda_se = e$lambda_se,
    interaction_frac = e$interaction_frac,
    lambda_source_vintage = lambda_from
  )
}) |>
  purrr::list_rbind() |>
  dplyr::left_join(stop_cost, by = "species")

readr::write_csv(pools, file.path(out_dir, "pools.csv"))
readr::write_csv(pool_summary, file.path(out_dir, "pool_summary.csv"))
readr::write_csv(vint, file.path(out_dir, "vintage.csv"))

cat("\n=== Pool summary ===\n")
print(as.data.frame(dplyr::select(pool_summary, species, search_mode, round,
                                  theta, n_moves, dAs, dPr, dAs_adj, score,
                                  power_pr_adj, lambda, lambda_se)),
      row.names = FALSE, digits = 3)
cat("\nThe power above carries lambda's own standard error: +/- 1.96 se on\n",
    "lambda is roughly +/- 15 power points, so differences of a few tenths\n",
    "of a point between proposals are below the noise floor. The search is a\n",
    "screen against bad choices, not a fine-grained preference order.\n",
    sep = "")

cat("\n=== The constraint, and what stopping now costs ===\n")
print(as.data.frame(dplyr::select(pool_summary, species, n_anchors_usable,
                                  n_anchors_used, max_new_per_pool, n_new,
                                  n_new_greedy, n_pruned, score,
                                  score_unconstrained, regret,
                                  score_if_remaining_found, n_unsearched,
                                  anchor_pool_switches, dpr_tol_relaxed)),
      row.names = FALSE, digits = 4)
if (any(pool_summary$anchor_pool_switches > 0)) {
  cat("\nAnchors are ONE global list the search may reassign between pools, so\n",
      "an accession found for As+ can end up in As-. That is deliberate -- a\n",
      "found bag is found, and the As+/As- label is a property of the design,\n",
      "not of the seed lot -- but see previous_pool in search_request.csv.\n",
      sep = "")
}

# PREREG TRIPWIRE. Pr-adjusted power as a SELECTION index needs no model
# version bump, and an unconstrained proposal is Pr-matched by construction.
# An anchored one need not be: when |dPr| exceeds the tolerance the `pool`
# estimand should condition on Pr, which is a change to the pre-registered
# analysis and a decision for a person.
if (any(abs(pool_summary$dPr) > pr_tolerance)) {
  warning("|dPr| exceeds the tolerance of ", pr_tolerance,
          " for: ", paste(pool_summary$species[abs(pool_summary$dPr) >
                                               pr_tolerance], collapse = ", "),
          ". The pools are no longer Pr-matched by construction, so adopting ",
          "this proposal means the `pool` estimand must condition on Pr -- a ",
          "PREREG_MODEL version bump, not a code change. See prereg_lines().",
          call. = FALSE)
}

# ------------------------------------------------------------
# Is the seed store itself biased?
# ------------------------------------------------------------

if (!is.null(buckets_tbl) && any(ledger$status %in% c("found", "missing"))) {
  balance <- found_missing_balance(ledger, buckets_tbl)
  readr::write_csv(balance, file.path(out_dir, "found_missing_balance.csv"))
  cat("\n=== Found vs missing, in the effects themselves ===\n")
  cat("If storage is systematically short of (say) old low-As lines, then\n",
      "anchoring entrenches that, and the As- pool ends up describing the\n",
      "store's history rather than the trait. A caveat on the interpretation\n",
      "if it is real, not a bug in the search -- so it is reported, never\n",
      "acted on automatically.\n", sep = "")
  print(as.data.frame(balance), row.names = FALSE, digits = 3)
}

# ------------------------------------------------------------
# The deliverable: what to go and find
# ------------------------------------------------------------

request <- purrr::imap(built, \(res, name) {
  search_request(res, inputs[[name]], round = round_no, ledger = ledger,
                 n_alternates = 3L, tol = pr_tolerance)
}) |> purrr::list_rbind()

readr::write_csv(request, file.path(out_dir, "search_request.csv"))
writeLines(search_request_md(request, vintage = vintage, round = round_no,
                             ledger_path = basename(ledger_path)),
           file.path(out_dir, "search_request.md"))

firsts <- dplyr::filter(request, is_first_choice)
cat("\n=== Search request ===\n")
if (nrow(firsts) == 0) {
  cat("Nothing to find: every pool member is already in hand.\n")
} else {
  cat(nrow(firsts), " slot(s) to fill, each with up to 3 ranked alternates.\n",
      "`score_cost` is the drop in t from taking an alternate instead -- a\n",
      "few thousandths is a genuine tie at this reliability.\n", sep = "")
  print(as.data.frame(dplyr::select(firsts, species, pool, slot, acc, Pr, As,
                                    n_entries, previous_status)),
        row.names = FALSE, digits = 3)
}

# ------------------------------------------------------------
# Record the proposal in the ledger
#
# Append-only: rows a person has already filled in come back untouched, and
# only pool_proposed on a still-unsearched row is refreshed.
# ------------------------------------------------------------

upd <- update_seed_ledger(ledger, pools, round = round_no)
write_seed_ledger(upd$ledger, ledger_path)
cat("\n", basename(ledger_path), ": ", nrow(upd$added), " row(s) added, ",
    upd$n_repooled, " re-pooled, ", nrow(upd$ledger), " total.",
    if (nrow(upd$added) > 0)
      "\nFill in `status` (found / missing / missing_retry) and `searched_on`, then re-run."
    else "", "\n", sep = "")

# ------------------------------------------------------------
# What changed since last time
#
# The pools are only worth validating if they are stable. If a couple of new
# trials reshuffle most of a pool, that says the effects are not settled, and
# it is much better to know that before seed is ordered than after.
# ------------------------------------------------------------

previous <- list.dirs(out_root, recursive = FALSE) |>
  setdiff(out_dir) |>
  (\(d) d[file.exists(file.path(d, "pools.csv"))])() |>
  sort()

if (length(previous) > 0) {
  prev_dir <- utils::tail(previous, 1)
  prev <- readr::read_csv(file.path(prev_dir, "pools.csv"), show_col_types = FALSE)
  dif <- pool_diff(prev, pools)

  readr::write_csv(dif$detail, file.path(out_dir, "pool_diff.csv"))

  cat("\n=== Change since ", basename(prev_dir), " ===\n", sep = "")
  print(as.data.frame(dif$summary), row.names = FALSE)

  churn <- dif$detail |>
    dplyr::filter(status != "retained") |>
    dplyr::count(species, name = "n_changed")
  retained <- dif$detail |>
    dplyr::filter(status == "retained") |>
    dplyr::count(species, name = "n_retained")

  dplyr::full_join(retained, churn, by = "species") |>
    dplyr::mutate(dplyr::across(c(n_retained, n_changed), \(x) tidyr::replace_na(x, 0L)),
                  pct_retained = round(100 * n_retained / (n_retained + n_changed))) |>
    as.data.frame() |>
    print(row.names = FALSE)

  # WHAT THIS IS AND IS NOT. Churn is reported because it is cheap and because
  # somebody ordering seed will want to know what moved. It is NOT a criterion:
  # the analysis method is still being settled, so churn between vintages
  # mostly measures changes to the method rather than instability in the
  # effects, and nothing in this chain gates on it. The earlier 70%-retention
  # guidance has been removed for that reason.
  #
  # The pool_index stamp says whether the two vintages are even comparable.
  prev_sum <- file.path(prev_dir, "pool_summary.csv")
  prev_idx <- if (file.exists(prev_sum)) {
    ps <- readr::read_csv(prev_sum, show_col_types = FALSE)
    if ("pool_index" %in% names(ps)) unique(ps$pool_index)[1] else NA_character_
  } else NA_character_

  cat("\nselection rule: ", POOL_INDEX_VERSION, " now, ",
      if (is.na(prev_idx)) "unstamped" else prev_idx, " then",
      if (!identical(prev_idx, POOL_INDEX_VERSION))
        " -- so this diff reflects the rule change as well as the new data."
      else " -- same rule, so this diff is the new data alone.",
      "\n", sep = "")
} else {
  cat("\nNo previous vintage to compare with; this is the baseline.\n")
}

# ------------------------------------------------------------
# Figure: where the pools sit in the Pr x As plane
# ------------------------------------------------------------

plot_data <- purrr::imap(inputs, \(inp, name) {
  inp$accessions |>
    dplyr::mutate(species = name) |>
    dplyr::left_join(dplyr::select(pools, species, acc, pool),
                     by = c("species", "acc")) |>
    dplyr::mutate(pool = tidyr::replace_na(pool, "not selected"))
}) |> purrr::list_rbind()

p <- plot_data |>
  dplyr::mutate(pool = factor(pool, levels = c("As+", "As-", "not selected"))) |>
  ggplot2::ggplot(ggplot2::aes(Pr, As, colour = pool, alpha = pool, size = pool)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  ggplot2::geom_point() +
  ggplot2::scale_colour_manual(values = c("As+" = "#1b7837", "As-" = "#762a83",
                                          "not selected" = "grey60")) +
  ggplot2::scale_alpha_manual(values = c("As+" = 1, "As-" = 1, "not selected" = 0.35),
                              guide = "none") +
  ggplot2::scale_size_manual(values = c("As+" = 2.2, "As-" = 2.2, "not selected" = 1.2),
                             guide = "none") +
  ggplot2::facet_wrap(~ species, scales = "free") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Validation pools: high producers, contrasted for associate effect",
    # Stated from inp$eligibility, not hard-coded: the caption is one of the
    # four places that used to describe the eligibility rule in its own words.
    subtitle = paste0(
      "vintage ", vintage, "; ", paste(unique(pool_summary$search_mode),
                                       collapse = "/"),
      " search; candidates restricted to ",
      if (inputs[[1]]$eligibility$rule == "reliability")
        sprintf("associate reliability >= %.3f and >= %d partners",
                inputs[[1]]$eligibility$rel_min,
                inputs[[1]]$eligibility$partner_floor)
      else sprintf(">= %d distinct partners", inputs[[1]]$eligibility$min_partners),
      ", and Pr above the ", round(100 * pr_quantile), "th percentile"),
    x = "producer effect (own yield, g/m2)",
    y = "associate effect (partner yield, g/m2)", colour = NULL
  )

ggplot2::ggsave(file.path(out_dir, "pools_Pr_vs_As.png"), p,
                width = 10, height = 5, dpi = 150)

message("\nwrote ", out_dir)
