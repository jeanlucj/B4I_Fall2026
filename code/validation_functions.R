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
#     standard error has terms that no amount of plot replication touches:
#
#       SE(delta)^2 = 2*sigma2_within/n + (int_sd)^2/n_loc + 4*sigma2_e/P
#                     \__ accessions __/   \_ As x loc _/    \_ plots _/
#
#     with sigma2_within = PEV + within-pool variance of the BLUPs.  TWO of the
#     three are free of P, so roughly 90% of the standard error cannot be
#     bought with plots -- and at the measured interaction most of that is the
#     MIDDLE term, which only more locations reduce.  It also carries just
#     n_loc - 1 degrees of freedom, so it drives the effective df as well as
#     the variance; see contrast_se().
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

# b4i_fit_frame() and apply_trial_qc() live in dge_ige_functions.R, and
# validation_inputs() must use exactly the filter chain the fit used -- that
# shared chain is the fix for SELF_CRITIQUE.md finding A. Sourced here rather
# than left to each driver, so no driver can forget it and fall back to reading
# the raw phenotype file. Guarded, because several drivers source both.
if (!exists("b4i_fit_frame", mode = "function")) {
  source(here::here("code", "dge_ige_functions.R"))
}

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

  # How many members of each pool an ANCHORED search may draw from outside the
  # accessions already found in seed storage -- i.e. how many new bags you are
  # willing to go and look for, per pool, per round.  3 is the user's choice
  # and is tight: it needs 2 * (n - 3) anchors, which is 85% of oat's 40 found
  # but 90% of pea's 60, so pea will often need more.  pool_search_curve()
  # prices each extra unit.
  max_new_per_pool = 3L,

  # How much of the optimal score the parsimony stage may give back to ask for
  # fewer new accessions.  0.5% of t is about 0.3 power points, far below the
  # +/- 15 that lambda's own standard error carries -- see THE PARSIMONY STAGE.
  score_slack_frac = 0.005,

  n_locations = 5L,
  # Locations are swept, and this is the sweep that matters: the pool x location
  # term is the largest of the three variance components and carries only
  # n_loc - 1 degrees of freedom, so it is the only lever that touches the
  # binding constraint. 4 is the realistic downside -- B4I_2025_AL was a total
  # loss, so a lost site is not hypothetical.
  n_locations_swept = c(4L, 5L),

  # Plots per location. Fixed, not swept: section 5 of VALIDATION_DESIGN.md
  # establishes that plots only ever moved the smallest of the three terms, so
  # the budget is no longer a live question.
  plots_per_location = 80L,

  # combinations per cell fixed across locations AND duplicated within each,
  # which is what separates plot error from combination x location
  n_anchor_per_cell = 1L,

  block_size = 10L
)

# ------------------------------------------------------------
# THE PRE-REGISTERED ANALYSIS
#
# One random-effects structure, three instantiations. Everything that computes
# power derives from this object, so a change to the analysis cannot leave a
# power claim describing a different experiment -- which is what happened when
# section 4 specified a mixed model while both power routes used a two-stage
# t-test, and neither carried the pool x location term that the power formula
# says dominates.
#
# Held as DATA rather than prose so the generated block in VALIDATION_DESIGN.md
# and the field book both print the same thing.
# ------------------------------------------------------------

PREREG_MODEL <- list(
  version = "prereg-v1",
  dated   = "2026-10-03",

  base = c("y ~ location + (1|location:block) + <FIXED> + <AsxE>",
           "      + (1|acc_focal) + (1|acc_partner) + (1|combination) + e"),

  estimands = list(
    slope = list(
      role     = "PRIMARY, one per species",
      response = "the PARTNER's yield",
      fixed    = "x_focal + x_partner   (predicted associate effects, centred)",
      AsxE     = "(0 + x_focal | location) + (0 + x_partner | location)",
      test     = "slope on x_focal > 0, one-sided alpha = 0.05",
      why      = paste(
        "The slope IS lambda measured in new data, on the same scale as the",
        "cross-validation's interaction_frac -- so the trial measures the",
        "quantity the whole design is conditioned on, and the random slope by",
        "location measures the interaction the power calculation assumes.",
        "It also does not depend on where the pool boundary fell.")),

    total = list(
      role     = "CO-PRIMARY",
      response = "oat_yield + w_pea * pea_yield   (w_pea = 1, the physical total)",
      fixed    = "x_total = GMA_oat + w_pea * GMA_pea, centred",
      AsxE     = "(0 + x_total | location)",
      test     = "slope on x_total > 0, one-sided alpha = 0.05",
      why      = paste(
        "Total productivity is the breeding objective, and it costs nothing:",
        "both yields are already recorded on every plot. GMA = Pr + As is",
        "exactly the per-accession contribution to a plot total, so no new",
        "quantity is estimated. Note it brings a different RESPONSE, not a",
        "different predictor -- because the pools are matched on Pr, dGMA is",
        "almost exactly dAs.")),

    pool = list(
      role     = "DESCRIPTIVE, not a separate test",
      response = "the PARTNER's yield",
      fixed    = "pool_focal + pool_partner + pool_focal:pool_partner",
      AsxE     = "(1|location:pool_focal) + (1|location:pool_partner)",
      test     = "the pool_focal contrast, reported with its interval",
      why      = paste(
        "The pre-registered two-point summary of `slope`, on the same data.",
        "Reported because it is what a breeder reads, and counted as a",
        "separate test would double-penalise: the two statistics are nearly",
        "equivalent (t 3.54 against 3.64 in a worked case), because selection",
        "has already removed the within-pool spread in the predictor that a",
        "regression would otherwise exploit."))
  ),

  multiplicity = paste(
    "HOLM over the THREE primaries: the two per-species slopes and the",
    "total-yield slope. Directions are pre-registered, so every test is",
    "one-sided at alpha = 0.05. The pool contrasts are descriptive summaries",
    "of the first two and are reported unadjusted; counting them as separate",
    "tests would penalise reporting two views of the same regression.",
    "This settles open decision 4."),

  secondary = c(
    "pool x pool interaction",
    "specific-combination variance from the anchors -- estimable almost only from the four anchor combinations, so expect it near the boundary",
    "the realised associate effect regressed on the focal accession's phenology BLUE from PRIOR trials (not from this trial: in-trial phenology is a MEDIATOR of the associate effect, so adjusting for it would bias the slope toward zero)",
    "an economically weighted total, w_pea != 1",
    "accession-level correlation between predicted and realised, for comparison with the cross-validation folds"),

  not_done = c(
    "a third pool selected on a phenotypic proxy -- considered and declined; it would cost about 25% more plots",
    "monoculture checks, so no land-equivalent ratio",
    "extra plots replicating specific combinations beyond the four anchors -- to be argued at the field-design stage via n_anchor_per_cell")
)

#' The pre-registered analysis, as printable lines.
prereg_lines <- function(m = PREREG_MODEL) {
  out <- c(sprintf("Analysis %s, fixed %s.", m$version, m$dated), "",
           "One random-effects structure:", paste0("    ", m$base), "")
  for (nm in names(m$estimands)) {
    e <- m$estimands[[nm]]
    out <- c(out, sprintf("%s  [%s]", nm, e$role),
             sprintf("    response : %s", e$response),
             sprintf("    <FIXED>  : %s", e$fixed),
             sprintf("    <AsxE>   : %s", e$AsxE),
             sprintf("    test     : %s", e$test), "")
  }
  c(out, "Multiplicity:", paste0("    ", m$multiplicity))
}

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
             pev_file = "BGLR_oat_pev.csv",
             id_col = "oatAcc", pheno_col = "germplasmName",
             partner_pheno_col = "intercropGermplasmName",
             response = "pea yield", residual_col = "var_pea"),
  pea = list(term = "G_pea", effects_file = "BGLR_pea_effects_all_seeds.csv",
             pev_file = "BGLR_pea_pev.csv",
             id_col = "peaAcc", pheno_col = "intercropGermplasmName",
             partner_pheno_col = "germplasmName",
             response = "oat yield", residual_col = "var_oat")
)

#' Everything the pool selection and the power calculation read.
#'
#' THE PLOT TABLE MUST BE THE ONE THE FIT SAW. This function used to read the
#' raw phenotype file with no trial QC and report `n_trials` from it, while the
#' fit applied QC and a hard-coded whitelist. The result was a vintage that
#' advertised nine trials and 3,567 plots behind BLUPs that saw five and 1,985,
#' and partner counts drawn from four trials the model never used -- so the
#' eligibility filter admitted exactly the accessions it existed to exclude.
#' See SELF_CRITIQUE.md finding A. Now both go through `b4i_fit_frame()`, and
#' the fit writes out what it saw so the two can be checked rather than assumed
#' equal.
#'
#' @param out_dir Where the fit's outputs live.
#' @param pheno_file Plot table. Filtered through `b4i_fit_frame()`, so trial
#'   QC and the monoculture drop are applied exactly as the fit applied them.
#' @param eligibility `"reliability"` (default) keeps accessions whose own
#'   associate effect is estimated well enough to be worth validating;
#'   `"partners"` is the older count rule, kept so an earlier vintage can be
#'   reproduced.
#' @param rel_min Reliability threshold. `NULL` derives it from the data as the
#'   first quartile of `rel_As_i` among accessions with exactly
#'   `min_partners` partners -- a deliberately permissive bar, set by the
#'   worst-estimated accessions the old count rule admitted.
#' @param min_partners Reference partner count for deriving `rel_min`, and the
#'   threshold itself under `eligibility = "partners"`.
#' @param partner_floor Hard connectivity floor, applied under BOTH rules.
#'   Reliability measures posterior precision, not identifiability: an
#'   accession grown with a single partner has its producer and associate
#'   effects perfectly aliased and can still score a respectable reliability by
#'   borrowing from well-genotyped relatives. The pools exist to separate
#'   exactly those two effects, so such an accession must not enter them
#'   whatever its reliability.
#' @return A list with one entry per species: the accession table (BLUPs,
#'   per-accession PEV and reliability, partner counts, plots), the variance
#'   components, the global reliability and PEV, the partner-yield residual
#'   variance, and the fit's provenance.
validation_inputs <- function(out_dir = here::here("output"),
                              pheno_file = here::here("output", "B4I_intercrop_pheno.rds"),
                              eligibility = c("reliability", "partners"),
                              rel_min = NULL,
                              min_partners = 3L,
                              partner_floor = 2L) {

  eligibility <- match.arg(eligibility)

  vc <- readr::read_csv(file.path(out_dir, "BGLR_variance_components.csv"),
                        show_col_types = FALSE)
  genetic  <- dplyr::filter(vc, component == "genetic")
  residual <- dplyr::filter(vc, component == "residual")

  # The same filter chain the fit ran, so partner counts describe the plots the
  # model actually used.
  pheno <- b4i_fit_frame(pheno_file = pheno_file,
                         qc_file = file.path(out_dir, "trial_qc.csv"),
                         quiet = TRUE)

  # WHAT THE FIT SAW, read from the fit rather than inferred. A missing file is
  # a hard stop, deliberately unlike apply_trial_qc()'s lenient missing-file
  # policy: falling back to the phenotype table's trial count is precisely the
  # bug this guard exists to prevent.
  fit_trials_file <- file.path(out_dir, "BGLR_fit_trials.csv")
  if (!file.exists(fit_trials_file)) {
    stop("no record of which trials the fit used: ", fit_trials_file,
         " is missing.\n  Re-run code/BGLR_multi_trait_model.R. Falling back ",
         "to the phenotype table's trial count is what produced ",
         "SELF_CRITIQUE.md finding A, so this is refused rather than guessed.",
         call. = FALSE)
  }
  fit_trials   <- readr::read_csv(fit_trials_file, show_col_types = FALSE)
  pheno_trials <- sort(unique(as.character(pheno$studyName)))

  if (!setequal(fit_trials$studyName, pheno_trials)) {
    stop("the fit used a different trial set from the QC-filtered plot table.",
         "\n  in the fit only:   ",
         paste(setdiff(fit_trials$studyName, pheno_trials), collapse = ", "),
         "\n  in the table only: ",
         paste(setdiff(pheno_trials, fit_trials$studyName), collapse = ", "),
         "\n  This is SELF_CRITIQUE.md finding A. Re-run ",
         "code/BGLR_multi_trait_model.R so the two agree.", call. = FALSE)
  }
  if (sum(fit_trials$n_plots) != nrow(pheno)) {
    stop("the fit saw ", sum(fit_trials$n_plots), " plots but the same filter ",
         "chain yields ", nrow(pheno),
         ". b4i_fit_frame() is supposed to be the one chain both use.",
         call. = FALSE)
  }

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
      dplyr::mutate(GMA = Pr + As) |>
      dplyr::left_join(counts, by = "acc") |>
      dplyr::mutate(
        n_partners = tidyr::replace_na(n_partners, 0L),
        n_plots    = tidyr::replace_na(n_plots, 0L)
      )

    # PER-ACCESSION prediction error variance, from the streamed MCMC draws
    # (code/BGLR_multi_trait_model.R). Averaged over chains, matching how Pr
    # and As are averaged above. Without this every accession gets the same
    # PEV, so one grown with twenty partners and one grown with a single
    # partner are treated as equally well known.
    pev_file <- file.path(out_dir, sp$pev_file)
    if (!file.exists(pev_file)) {
      stop("no per-accession PEV at ", pev_file,
           ".\n  Re-run code/BGLR_multi_trait_model.R, which streams the ",
           "coefficient draws and writes it.", call. = FALSE)
    }
    pev <- readr::read_csv(pev_file, show_col_types = FALSE) |>
      dplyr::group_by(acc) |>
      dplyr::summarise(PEV_Pr_i  = mean(PEV_Pr),
                       PEV_As_i  = mean(PEV_As),
                       PEV_GMA_i = mean(PEV_GMA), .groups = "drop")

    missing_pev <- setdiff(acc$acc, pev$acc)
    if (length(missing_pev) > 0) {
      stop(length(missing_pev), " accession(s) have an effect but no PEV, e.g. ",
           paste(utils::head(missing_pev, 5), collapse = ", "),
           ".\n  The effects and the draws come from the same fit, so this ",
           "means one of the two files is stale.", call. = FALSE)
    }

    sigma2_GMA <- sigma2_Pr + sigma2_As + 2 * mean(gg$cov_PrAs)

    # Var(true) = Var(BLUP) + E[PEV] for a conditionally unbiased predictor,
    # so the reliability of the BLUPs and the PEV both follow from the fit.
    # Clamped: a var(BLUP) above the component means the chains have not
    # settled, and a negative PEV would silently flatter every power estimate.
    rel <- function(v_blup, sigma2) min(max(v_blup / sigma2, 1e-6), 0.999)
    # The per-accession version is clamped ELEMENTWISE for the same reason: a
    # short chain can put one accession's PEV above the component variance,
    # which would otherwise give it a negative reliability.
    rel_i <- function(pev_i, sigma2) pmin(pmax(1 - pev_i / sigma2, 1e-6), 0.999)

    acc <- acc |>
      dplyr::left_join(pev, by = "acc") |>
      dplyr::mutate(rel_As_i  = rel_i(PEV_As_i,  sigma2_As),
                    rel_GMA_i = rel_i(PEV_GMA_i, sigma2_GMA))

    # THE ELIGIBILITY THRESHOLD, derived rather than chosen. The old rule kept
    # anything with at least `min_partners` partners, so the first quartile of
    # reliability among accessions sitting exactly at that count is the bar the
    # old rule was already willing to accept -- a permissive threshold, now
    # stated as a reliability instead of a proxy for one.
    ref <- dplyr::filter(acc, n_partners == min_partners)
    if (nrow(ref) < 4L) {
      ref <- dplyr::filter(acc, n_partners <= min_partners, n_partners > 0L)
      message("  ", name, ": only ", sum(acc$n_partners == min_partners),
              " accession(s) have exactly ", min_partners,
              " partners; deriving rel_min from the ", nrow(ref),
              " with at most that many")
    }
    rel_min_used <- if (!is.null(rel_min)) rel_min else
      stats::quantile(ref$rel_As_i, 0.25, names = FALSE)

    acc <- dplyr::mutate(acc, eligible = switch(
      eligibility,
      reliability = rel_As_i >= rel_min_used & n_partners >= partner_floor,
      partners    = n_partners >= min_partners
    ))

    list(
      species      = name,
      accessions   = acc,
      n_eligible   = sum(acc$eligible),
      # How eligibility was decided, carried so the console output, the figure
      # subtitle and the generated document all state the same rule instead of
      # each hard-coding its own description of it.
      eligibility  = list(
        rule          = eligibility,
        rel_min       = rel_min_used,
        rel_min_source = if (is.null(rel_min)) "derived" else "supplied",
        min_partners  = min_partners,
        partner_floor = partner_floor,
        n_eligible    = sum(acc$eligible),
        n_total       = nrow(acc),
        n_ref         = nrow(ref)
      ),
      # retained under its old name: validate_pool_selection.R asserts on it
      min_partners = min_partners,
      sigma2_Pr    = sigma2_Pr,
      sigma2_As    = sigma2_As,
      sigma2_GMA   = sigma2_GMA,
      sigma2_e     = sigma2_e,
      response     = sp$response,
      # GLOBAL reliability and PEV. Kept, with their meanings unchanged, so
      # every existing consumer keeps working; the per-accession versions are
      # additions in `accessions`, suffixed `_i`.
      rel_Pr       = rel(stats::var(acc$Pr), sigma2_Pr),
      rel_As       = rel(stats::var(acc$As), sigma2_As),
      rel_GMA      = rel(stats::var(acc$GMA), sigma2_GMA),
      # THE OLD PEV, BACKED OUT OF AN IDENTITY THAT DOES NOT HOLD IN THIS FIT.
      # It assumes Var(true) = Var(BLUP) + E[PEV] and solves for the second
      # term, giving sigma2_As - var(BLUP). Measured against the posterior
      # itself that is badly wrong here: on the 8-trial fit the identity is out
      # by a factor of 0.6 (oat: var(BLUP) 46.9 + E[PEV] 49.0 = 95.9 against a
      # component of 160.9, with mean(diag(G)) = 1.00, so the GRM scale is not
      # the explanation). Heavy shrinkage plus an inverse-Wishart prior on the
      # component will do that: the component absorbs variance the shrunken
      # BLUPs never express.
      #
      # So this number is 114 where the posterior says 49 -- 2.3x too large --
      # and it inflated sigma2_within and UNDERSTATED power. It is retained
      # only so the two routes can be compared in the vintage; nothing should
      # compute with it. build_pools() uses the measured per-accession PEV.
      PEV_As       = sigma2_As * (1 - rel(stats::var(acc$As), sigma2_As)),
      PEV_As_measured = mean(acc$PEV_As_i),
      # Ratio of the two. Far from 1 means the fitted component and the
      # posterior disagree about how much is known, which is worth seeing on
      # every vintage rather than discovering once.
      PEV_ratio_backed_out_to_measured =
        (sigma2_As * (1 - rel(stats::var(acc$As), sigma2_As))) /
          mean(acc$PEV_As_i),
      # Provenance of the fit these numbers describe, read from the fit itself.
      fit_trials   = sort(fit_trials$studyName),
      n_trials     = nrow(fit_trials),
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
# The constraint is handled with a Lagrange multiplier `theta` on Pr, and the
# key point is that ONE index serves both pools: rank the candidates by
# As + theta * Pr and take the top n and the bottom n.  Both selections then
# move along a single axis, so theta buys Pr balance efficiently, and
# disjointness is structural -- top n and bottom n of one ordering cannot
# collide once there are at least 2n candidates, which build_pools() already
# requires.
#
# WHY NOT A SEPARATE INDEX PER POOL, which is what this did until 2026-10-03.
# It scored the As+ pool by As + theta * Pr and the As- pool by As - theta * Pr.
# The second term REWARDS high Pr in the As- pool, so raising theta dragged BOTH
# pools toward high Pr instead of equalising them.  Balance was reached only
# indirectly and expensively -- oat needed theta = 0.61, spending a great deal of
# As extremity to get there.  The single index reaches the same balance at
# theta = 0.19.  Measured on the 2026-10-03 data at P = 400, both targeting
# dPr = 0:
#
#     oat n=20   dAs 15.24 -> 19.59   dPr -0.50 -> +0.23   power 0.691 -> 0.779
#     pea n=30   dAs 18.01 -> 25.07   dPr -0.66 -> -0.13   power 0.656 -> 0.839
#
# More power than any affordable increase in plots or pool size would buy, and
# the producer balance gets TIGHTER at the same time.  See docs/B4I_followups.md
# item 13.
#
# THE SEARCH IS A SCAN, NOT A BISECTION, for a related reason.  dPr(theta) is a
# step function and is not monotone -- a dense scan finds three sign changes for
# oat -- so bisecting for a sign change can converge on the wrong crossing.  A
# scan costs nothing (5,000 theta in 0.13 s through .pool_scan) and finds the
# global minimum of |dPr| rather than a local crossing.
# ------------------------------------------------------------

# Bump this whenever the selection rule changes, so a diff against an older
# vintage can be recognised as incomparable. This is PROVENANCE, not a
# criterion: pool churn is not a go/no-go for this project while the analysis
# method is still being settled, because churn then measures the method rather
# than the effects.
#   single-v1  2026-10-03  one index for both pools, dense scan on theta
#   single-v2  2026-10-03  eligibility is a per-accession reliability threshold
#                          plus a hard partner floor, replacing the partner
#                          count; and PEV is measured from the posterior draws
#                          rather than backed out of a variance identity
#   search-v3  2026-11-01  the membership is the optimum of the Pr-ADJUSTED
#                          contrast, dAs - b * dPr, found by steepest ascent
#                          seeded from the theta family, rather than the best
#                          point of that family alone; and the search can be
#                          ANCHORED to the accessions already found in seed
#                          storage.  On the 2026-10-31 vintage the
#                          unconstrained answer moves 3 of 40 oat members and
#                          2 of 60 pea, gaining 0.0015 in t -- so most of a
#                          diff across this bump is still data, not rule.
#
# search_mode is a SEPARATE column in pool_summary.csv rather than part of
# this stamp, because an anchored and an unconstrained run under the same rule
# ARE comparable: the stamp answers "were these built the same way", and the
# constraint's cost is exactly what search_curve.csv measures.
POOL_INDEX_VERSION <- "search-v3"

#' Scan many theta cheaply: the two summaries only, no tibbles.
#'
#' Base R on plain vectors.  The dplyr path in .pools_at() costs 1.7 ms a call,
#' which a few thousand theta would make felt inside validate_crossval.R's fold
#' loop; this is ~57x faster, so the dense scan is free.
#'
#' @return A two-column matrix, dAs and dPr, one row per theta.
.pool_scan <- function(As, Pr, n, thetas) {
  t(vapply(thetas, function(th) {
    o <- order(As + th * Pr, decreasing = TRUE)
    p <- o[seq_len(n)]
    m <- o[seq.int(length(o) - n + 1L, length(o))]
    c(dAs = mean(As[p]) - mean(As[m]), dPr = mean(Pr[p]) - mean(Pr[m]))
  }, numeric(2)))
}

#' One (theta) evaluation: the two pools and the differences they imply.
.pools_at <- function(cand, n, theta) {
  o     <- order(cand$As + theta * cand$Pr, decreasing = TRUE)
  plus  <- cand[utils::head(o, n), , drop = FALSE]
  minus <- cand[utils::tail(o, n), , drop = FALSE]

  list(plus = plus, minus = minus,
       dAs = mean(plus$As) - mean(minus$As),
       dPr = mean(plus$Pr) - mean(minus$Pr))
}

#' The candidate set a selection is made from: eligible, then above the Pr
#' threshold.
#'
#' Factored out of build_pools() so the searches and the search-input
#' classifier cannot drift from the filter the selection actually used -- the
#' kind of duplication that makes an oracle silently compare against the wrong
#' population.
.pool_candidates <- function(inp, pr_quantile = 0.5) {
  elig   <- dplyr::filter(inp$accessions, eligible)
  pr_min <- stats::quantile(elig$Pr, pr_quantile, names = FALSE)
  list(cand = dplyr::filter(elig, Pr >= pr_min),
       eligible = elig, pr_min = pr_min)
}

#' The As-on-Pr slope the Pr adjustment uses.
#'
#' THE ONE DEFINITION of `b`, so the document and the optimiser cannot drift.
#'
#' It takes the POST-THRESHOLD candidate set -- eligible AND above the Pr
#' quantile, i.e. `build_pools()$candidates` -- not the eligible set and not
#' the within-pool rows. The three differ materially: on the 2026-10-31
#' vintage the candidate slope is -0.171 (oat) / -0.443 (pea) while the
#' eligible-set slope is -0.272 / -0.525, and the within-pool slope is
#' selection-biased by the `As + theta * Pr >= cut` truncation. The candidate
#' set is the right one because it is the population the selection actually
#' chooses from.
pr_adjust_slope <- function(cand) {
  if (!all(c("As", "Pr") %in% names(cand))) {
    stop("pr_adjust_slope() needs columns As and Pr", call. = FALSE)
  }
  if (nrow(cand) < 3L) {
    stop("pr_adjust_slope() needs at least 3 candidates to regress; got ",
         nrow(cand), call. = FALSE)
  }
  unname(stats::coef(stats::lm(As ~ Pr, data = cand))[2])
}

#' Build the As+ / As- pools for one species.
#'
#' @param inp One element of validation_inputs().
#' @param n Accessions PER POOL (so 2n of that species enter the trial).
#' @param pr_quantile Candidates must be above this quantile of Pr among the
#'   eligible accessions.  0.5 keeps the better half: "all of these are good
#'   producers" is the claim the design has to support.
#' @param tol Largest acceptable |mean Pr difference| between pools, g/m2.
#' @param theta_max Upper end of the scan over the Pr multiplier.
#' @param n_grid Points in the coarse scan.  A refinement pass around the winner
#'   follows, so this sets where the search looks rather than how precisely.
build_pools <- function(inp, n, pr_quantile = 0.5, tol = 1.0,
                        theta_max = 10, n_grid = 2001L) {
  pc     <- .pool_candidates(inp, pr_quantile)
  cand   <- pc$cand
  pr_min <- pc$pr_min

  if (nrow(cand) < 2 * n) {
    stop(inp$species, ": ", nrow(cand), " candidates above the Pr threshold ",
         "but ", 2 * n, " needed. Lower pr_quantile or min_partners, or use a ",
         "smaller pool.", call. = FALSE)
  }

  # theta = 0 selects on As alone; because Pr and As are negatively correlated
  # in the BLUPs that leaves the As+ pool short on Pr, so dPr starts negative
  # and rises with theta in TREND.
  #
  # IN TREND, not monotonically: dPr is a step function and adjacent steps can
  # reverse, so a bisection for a sign change can converge on a local crossing
  # rather than the best one.  Scan instead, and take the global minimum of
  # |dPr|.  n_sign_changes is reported so the non-monotonicity is a number
  # rather than a comment.
  #
  # TIES ARE BROKEN BY LARGER dAs.  A step function makes exact ties in |dPr|
  # common -- a whole interval of theta gives the identical pair of pools -- and
  # without the tie-break the pick among them is whichever the grid happened to
  # land on first.  Among theta within `eps` of the best |dPr|, take the one
  # with the largest contrast.
  pick <- function(thetas) {
    sc  <- .pool_scan(cand$As, cand$Pr, n, thetas)
    eps <- max(1e-9, 1e-6 * max(abs(sc[, "dPr"])))
    ok  <- which(abs(sc[, "dPr"]) <= min(abs(sc[, "dPr"])) + eps)
    thetas[ok[which.max(sc[ok, "dAs"])]]
  }

  coarse <- seq(0, theta_max, length.out = n_grid)
  theta  <- pick(coarse)

  # Refine around the winner so the answer does not depend on grid resolution.
  step  <- theta_max / (n_grid - 1)
  theta <- pick(seq(max(0, theta - step), min(theta_max, theta + step),
                    length.out = 201L))

  res <- .pools_at(cand, n, theta)

  dPr_grid      <- .pool_scan(cand$As, cand$Pr, n, coarse)[, "dPr"]
  n_sign_changes <- sum(diff(sign(dPr_grid)) != 0)

  if (abs(res$dPr) > tol) {
    warning(inp$species, ": producer means differ by ", round(res$dPr, 2),
            " g/m2, outside the tolerance of ", tol,
            ". The Pr constraint could not be met at n = ", n, ".",
            call. = FALSE)
  }

  pools <- .pools_from_acc(inp, res$plus$acc, res$minus$acc)

  # Within-pool variance of the TRUE effects: what the BLUPs still spread by
  # after selection, plus what we do not know about each one. Computed PER POOL
  # and then averaged, because the two pools need not be equally well estimated
  # -- the As+ extreme and the As- extreme are different parts of the
  # distribution, and with a per-accession PEV they can differ materially.
  #
  # The scalar `sigma2_within` is retained as their mean. It is what
  # contrast_se() has always taken and what pool_summary.csv and power_grid.csv
  # carry as a single column, so adding the per-pool values rather than
  # replacing the scalar keeps every existing consumer working. The SE is
  # unaffected by the split anyway -- (A + B)/n is identically 2*mean/n -- but
  # the DEGREES OF FREEDOM are not, via the Welch term in contrast_power().
  c(.pool_stats(pools, n = n, cand = cand, pr_min = pr_min),
    list(theta = theta, n_sign_changes = n_sign_changes))
}

#' Everything a pair of pools implies, from their membership alone.
#'
#' ONE DEFINITION of the pool summaries, shared by build_pools() and by
#' pools_from_csv().  The reason it is factored out rather than recomputed is
#' that validate_power.R's estimand section reports the ADOPTED pools while
#' power_grid() re-derives them: a second copy of this arithmetic is how
#' estimand_power.csv and pool_summary.csv would come to describe the same
#' accessions with different numbers.
#'
#' @param pools The .pools_from_acc() table -- species, pool, acc and the
#'   per-accession columns.  Order within a pool does not matter: everything
#'   here is a mean or a variance.
.pool_stats <- function(pools, ..., n, cand, pr_min) {
  if (...length() > 0L) {
    stop("'.pool_stats()' takes only `pools` positionally; name the rest",
         call. = FALSE)
  }
  pl <- pools$pool == "As+"
  stopifnot(sum(pl) == n, sum(!pl) == n)

  # Within-pool variance of the TRUE effects: what the BLUPs still spread by
  # after selection, plus what we do not know about each one. Computed PER POOL
  # and then averaged, because the two pools need not be equally well estimated
  # -- the As+ extreme and the As- extreme are different parts of the
  # distribution, and with a per-accession PEV they can differ materially.
  #
  # The scalar `sigma2_within` is retained as their mean. It is what
  # contrast_se() has always taken and what pool_summary.csv and power_grid.csv
  # carry as a single column, so adding the per-pool values rather than
  # replacing the scalar keeps every existing consumer working. The SE is
  # unaffected by the split anyway -- (A + B)/n is identically 2*mean/n -- but
  # the DEGREES OF FREEDOM are not, via the Welch term in contrast_power().
  pev_plus  <- mean(pools$PEV_As_i[pl])
  pev_minus <- mean(pools$PEV_As_i[!pl])
  var_blup_plus  <- stats::var(pools$As[pl])
  var_blup_minus <- stats::var(pools$As[!pl])
  var_blup_within <- mean(c(var_blup_plus, var_blup_minus))
  sigma2_within_plus  <- pev_plus  + var_blup_plus
  sigma2_within_minus <- pev_minus + var_blup_minus

  # What the Pr adjustment needs. Sxx_w is the within-pool sum of squares of Pr
  # pooled over both pools -- the denominator of the ANCOVA inflation.
  var_Pr_plus  <- stats::var(pools$Pr[pl])
  var_Pr_minus <- stats::var(pools$Pr[!pl])

  list(
    pools = pools, n = n, pr_min = pr_min,
    # Stamped into pool_summary.csv so a later run can tell whether the previous
    # vintage's pools are comparable with its own. Churn between vintages built
    # by different rules measures the RULE rather than anything about the data,
    # which is the whole reason the stamp exists.
    pool_index = POOL_INDEX_VERSION,
    n_candidates = nrow(cand),
    dAs = mean(pools$As[pl]) - mean(pools$As[!pl]),
    dPr = mean(pools$Pr[pl]) - mean(pools$Pr[!pl]),
    mean_Pr_plus = mean(pools$Pr[pl]), mean_Pr_minus = mean(pools$Pr[!pl]),
    mean_As_plus = mean(pools$As[pl]), mean_As_minus = mean(pools$As[!pl]),
    var_blup_within = var_blup_within,
    var_blup_plus = var_blup_plus, var_blup_minus = var_blup_minus,
    PEV_plus = pev_plus, PEV_minus = pev_minus,
    sigma2_within_plus = sigma2_within_plus,
    sigma2_within_minus = sigma2_within_minus,
    sigma2_within = mean(c(sigma2_within_plus, sigma2_within_minus)),
    # Spread of the PREDICTED effect among the 2n selected accessions, which is
    # what the continuous-slope estimand regresses on: Var(x) over both pools is
    # (dAs/2)^2 + var_As_within_pools.
    var_As_within_pools = var_blup_within,
    var_Pr_plus = var_Pr_plus, var_Pr_minus = var_Pr_minus,
    Sxx_w = (n - 1) * (var_Pr_plus + var_Pr_minus),
    b_As_on_Pr = pr_adjust_slope(cand),
    # The candidate set this selection was made from. Returned so a caller --
    # or a test -- can assert against it instead of re-deriving the filter,
    # which three places in tests/test_validation.R used to do by hand.
    candidates = cand
  )
}

#' The ADOPTED pools, read back from pools.csv.
#'
#' WHY THIS EXISTS.  code/validate_power.R has two jobs that pull in opposite
#' directions.  power_grid() SWEEPS n, and the adopted pools exist at one n
#' only, so it must keep re-deriving membership with build_pools() -- that
#' table answers "what would a design of size n achieve".  The estimand
#' section describes the trial that will actually be sown, so it must read the
#' pools that were actually adopted, which since search-v3 are the search's
#' and not build_pools'.  Without this the field book and estimand_power.csv
#' described different accessions.
#'
#' Everything numeric is RECOMPUTED from `inp`, so a stale pools.csv cannot
#' smuggle last vintage's effect sizes into this vintage's power.  Only the
#' membership is taken from the file.
#'
#' Falls back to build_pools() with a message when the file is absent, which
#' is the fresh-clone and test case.
#'
#' @param theta Carried through from pool_summary.csv when available; it is
#'   descriptive only.
pools_from_csv <- function(inp, pools_file, ..., n = NULL, pr_quantile = 0.5) {
  if (...length() > 0L) {
    stop("pools_from_csv() takes inp and pools_file positionally; name the ",
         "rest", call. = FALSE)
  }
  if (!file.exists(pools_file)) {
    message("no pools.csv at ", pools_file,
            " -- re-deriving the pools with build_pools()")
    return(build_pools(inp, n, pr_quantile = pr_quantile))
  }
  tbl <- readr::read_csv(pools_file, show_col_types = FALSE) |>
    dplyr::filter(species == inp$species)
  if (nrow(tbl) == 0L) {
    stop(pools_file, " has no rows for species ", inp$species, call. = FALSE)
  }

  # A pools.csv naming an accession the fit no longer carries means the FIT
  # moved under the adopted pools. That is a human decision -- re-run the pool
  # selection, or deliberately keep the old design -- so stop pointedly rather
  # than silently dropping a pool member and reporting power for 19 of 20.
  gone <- setdiff(tbl$acc, inp$accessions$acc)
  if (length(gone) > 0L) {
    stop(pools_file, ": ", length(gone), " adopted ", inp$species,
         " accession(s) are not in this fit -- ",
         paste(utils::head(gone, 5), collapse = ", "),
         if (length(gone) > 5L) ", ..." else "",
         ".
  The fit has moved under the adopted pools. Re-run ",
         "code/validate_pool_selection.R to adopt a new design, or point ",
         "this at the vintage the pools belong to.", call. = FALSE)
  }

  acc_plus  <- tbl$acc[tbl$pool == "As+"]
  acc_minus <- tbl$acc[tbl$pool == "As-"]
  n_found   <- length(acc_plus)
  if (!is.null(n) && n_found != n) {
    stop(pools_file, ": ", inp$species, " has ", n_found,
         " accessions per pool but n = ", n, " was asked for.", call. = FALSE)
  }
  pc <- .pool_candidates(inp, pr_quantile)
  c(.pool_stats(.pools_from_acc(inp, acc_plus, acc_minus), n = n_found,
                cand = pc$cand, pr_min = pc$pr_min),
    list(theta = NA_real_, n_sign_changes = NA_integer_, from_file = pools_file))
}

#' What changed between two vintages of the pools.
#'
# ------------------------------------------------------------
# THE POOL SEARCH
#
# build_pools() finds the pools by a dense scan over one Lagrange multiplier
# (theta) on a single index, As + theta * Pr.  That is a one-parameter family,
# and on the 2026-10-31 vintage it already sits at the optimum of the
# Pr-adjusted score to four significant figures: a steepest-ascent climb from
# it moves 3 of 40 oat members and gains 0.0015 in t.  So this search is NOT
# here to improve the unconstrained proposal, and must not be sold as doing
# so.  It is here because the ANCHORED problem -- "keep at least this many of
# the accessions we actually found in seed storage" -- has no one-parameter
# family to scan, and because an explicit optimum is what the constraint's
# cost is measured against.
#
# THE OBJECTIVE IS THE NONCENTRALITY t, NOT POWER.  Power saturates -- pea is
# already 0.89 and would press against 1 at a larger pool -- and a flat
# objective gives a hill-climb nothing to climb.  t is strictly monotone in
# power at fixed df, so maximising it maximises power while staying
# well-conditioned.
#
# IT IS COMPUTED AT lambda = 1 AND interaction_frac = 0, which makes selection
# a pure function of the fit and independent of validate_crossval.R.  That is
# not an approximation adopted for convenience: the selected membership is
# invariant to both across the whole plausible range (identical accessions for
# lambda in 0.5..1.0 crossed with interaction_frac in 0..0.74), because the
# landscape is flat.  A test pins the invariance, so a future vintage that
# breaks it fails loudly instead of moving the pools in silence.  The power
# REPORTED beside a proposal is still the Pr-adjusted power at the measured
# lambda, via contrast_power().  Selection and reporting are deliberately
# different quantities.
# ------------------------------------------------------------

#' Everything the search needs about one species, in the shape it needs it.
#'
#' Plain numeric vectors rather than a tibble: the climb evaluates thousands
#' of candidate moves, and the dplyr path costs ~1.7 ms a call (the same
#' reason .pool_scan() works on vectors).
#'
#' TWO DELIBERATE CHOICES HERE.
#'
#' Rows are sorted by `acc`, so the search is order-independent BY
#' CONSTRUCTION rather than in expectation: "the lowest candidate index" is
#' then exactly "the lexicographically first accession", which makes the
#' tie-break deterministic however `inp$accessions` happened to be ordered.
#'
#' As and Pr are CENTRED on the candidate set.  Every quantity the score uses
#' is a difference of means or a variance, so centring changes nothing, but it
#' removes the cancellation in the computational variance formula
#' (Q - S^2/m)/(m-1) that the O(1) updates rely on.
pool_search_context <- function(cand, ..., n, sigma2_e, P, b = NULL,
                                anchors = character(0),
                                exclude = character(0)) {
  if (...length() > 0L) {
    stop("pool_search_context() takes only `cand` positionally; name every ",
         "other argument", call. = FALSE)
  }
  if (length(exclude) > 0L) cand <- dplyr::filter(cand, !acc %in% exclude)
  cand <- dplyr::arrange(cand, acc)
  if (nrow(cand) < 2L * n) {
    stop(nrow(cand), " candidates but ", 2L * n, " needed for pools of ", n,
         call. = FALSE)
  }
  list(acc = cand$acc,
       As  = cand$As - mean(cand$As),
       Pr  = cand$Pr - mean(cand$Pr),
       PEV = cand$PEV_As_i,
       is_anchor = cand$acc %in% anchors,
       n = n, sigma2_e = sigma2_e, P = P,
       b = b %||% pr_adjust_slope(cand),
       N = nrow(cand), candidates = cand)
}

#' The six running sums one pool's score depends on.
#'
#' A swap updates each with one subtraction and one addition, which is what
#' makes a full sweep of the move set affordable.
.pool_sums <- function(ctx, idx) {
  c(m    = length(idx),
    S_As = sum(ctx$As[idx]),  Q_As = sum(ctx$As[idx]^2),
    S_Pr = sum(ctx$Pr[idx]),  Q_Pr = sum(ctx$Pr[idx]^2),
    S_PEV = sum(ctx$PEV[idx]),
    a    = sum(ctx$is_anchor[idx]))
}

#' One pool's sums after losing `out` and gaining `ins`.
#'
#' Either argument may be a vector, so one call produces a whole family of
#' candidate moves. The O(1) update -- one subtraction and one addition per
#' sum -- is what makes a full sweep of the move set affordable, and having
#' ONE definition of it is what stops the climb and the prune from drifting.
.sums_bump <- function(s, out, ins, ctx) {
  list(m     = s[["m"]],
       S_As  = s[["S_As"]]  - ctx$As[out]    + ctx$As[ins],
       Q_As  = s[["Q_As"]]  - ctx$As[out]^2  + ctx$As[ins]^2,
       S_Pr  = s[["S_Pr"]]  - ctx$Pr[out]    + ctx$Pr[ins],
       Q_Pr  = s[["Q_Pr"]]  - ctx$Pr[out]^2  + ctx$Pr[ins]^2,
       S_PEV = s[["S_PEV"]] - ctx$PEV[out]   + ctx$PEV[ins],
       a     = s[["a"]] - ctx$is_anchor[out] + ctx$is_anchor[ins])
}

#' Score from the two pools' sufficient statistics.
#'
#' Every argument may be a vector of the same length, so one call scores a
#' whole family of candidate moves. This is the only place the objective is
#' written down.
.score_from_sums <- function(sp, sm, ctx) {
  n <- ctx$n
  mean_As_p <- sp[["S_As"]] / sp[["m"]]
  var_As_p  <- (sp[["Q_As"]] - sp[["S_As"]]^2 / sp[["m"]]) / (sp[["m"]] - 1)
  mean_Pr_p <- sp[["S_Pr"]] / sp[["m"]]
  var_Pr_p  <- (sp[["Q_Pr"]] - sp[["S_Pr"]]^2 / sp[["m"]]) / (sp[["m"]] - 1)
  mean_As_m <- sm[["S_As"]] / sm[["m"]]
  var_As_m  <- (sm[["Q_As"]] - sm[["S_As"]]^2 / sm[["m"]]) / (sm[["m"]] - 1)
  mean_Pr_m <- sm[["S_Pr"]] / sm[["m"]]
  var_Pr_m  <- (sm[["Q_Pr"]] - sm[["S_Pr"]]^2 / sm[["m"]]) / (sm[["m"]] - 1)

  s2w_p <- sp[["S_PEV"]] / sp[["m"]] + var_As_p
  s2w_m <- sm[["S_PEV"]] / sm[["m"]] + var_As_m

  dAs <- mean_As_p - mean_As_m
  dPr <- mean_Pr_p - mean_Pr_m
  Sxx <- (n - 1) * (var_Pr_p + var_Pr_m)
  f   <- 1 + (n / 2) * dPr^2 / Sxx
  dAs_adj <- dAs - ctx$b * dPr
  v   <- f * (s2w_p + s2w_m) / n + f * 4 * ctx$sigma2_e / ctx$P

  list(score = dAs_adj / sqrt(v), se = sqrt(v),
       dAs = dAs, dPr = dPr, dAs_adj = dAs_adj, f = f, Sxx = Sxx,
       sigma2_within_plus = s2w_p, sigma2_within_minus = s2w_m,
       a_plus = sp[["a"]], a_minus = sm[["a"]])
}

#' Score one membership, from scratch.
#'
#' The public entry point, and the oracle the O(1) updates are tested against.
pool_score <- function(plus, minus, ctx) {
  .score_from_sums(.pool_sums(ctx, plus), .pool_sums(ctx, minus), ctx)
}

#' Run `expr` under `seed` without disturbing the caller's RNG stream.
#'
#' The climb consumes no randomness, but restarts do, and validate_design.R's
#' design_seed and validate_power.R's sim_seed must keep seeing exactly the
#' stream they see today.
.with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv())
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv())) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

#' The distinct pool pairs the theta scan can reach.
#'
#' dPr(theta) is a step function, so a few thousand theta collapse to far
#' fewer memberships -- 71 distinct pairs for oat, 83 for pea. Those are the
#' search's seeds: already near-optimal, already spread out, and free, since
#' .pool_scan() covers the grid in about 0.13 s.
.theta_family <- function(ctx, theta_max = 10, n_grid = 2001L) {
  n <- ctx$n; N <- ctx$N
  out <- list(); seen <- character(0)
  for (th in seq(0, theta_max, length.out = n_grid)) {
    o     <- order(ctx$As + th * ctx$Pr, decreasing = TRUE)
    plus  <- o[seq_len(n)]
    minus <- o[seq.int(N - n + 1L, N)]
    key   <- paste(c(sort(plus), -1L, sort(minus)), collapse = ",")
    if (key %in% seen) next
    seen <- c(seen, key)
    out[[length(out) + 1L]] <- list(plus = plus, minus = minus, theta = th)
  }
  out
}

#' Steepest-ascent exchange from one starting membership.
#'
#' STEEPEST ASCENT, not first improvement: the sweep is cheap, and taking the
#' best move in the whole neighbourhood makes the result independent of the
#' order candidates are enumerated in, which first-improvement is not. No RNG
#' is consumed.
#'
#' The move set is three kinds, and the third is NOT reachable by composing
#' the first two -- both legs would have to pass an accession through the
#' reservoir, which it never enters:
#'
#'   1. a member of As+ for an unused candidate
#'   2. a member of As- for an unused candidate
#'   3. a member of As+ for a member of As-
#'
#' A move is accepted only if it improves the score by more than
#' `min_gain_frac` of it, so the climb cannot chase floating-point noise and
#' is guaranteed to terminate.
#'
#' @param m_min Minimum members of each pool that must come from
#'   `ctx$is_anchor`. 0 leaves the search unconstrained.
#' @param tol Largest acceptable |dPr|. Moves that breach it are rejected, so
#'   feasibility is preserved by construction rather than repaired afterwards.
.pools_climb <- function(plus, minus, ctx, ..., m_min = 0L, tol = Inf,
                         min_gain_frac = 1e-4, max_moves = 1000L) {
  if (...length() > 0L) {
    stop("'.pools_climb()' takes plus, minus and ctx positionally; name the ",
         "rest", call. = FALSE)
  }
  n <- ctx$n
  anch <- ctx$is_anchor
  sp <- .pool_sums(ctx, plus); sm <- .pool_sums(ctx, minus)
  cur <- .score_from_sums(sp, sm, ctx)
  n_moves <- 0L

  bump <- function(s, out, ins) .sums_bump(s, out, ins, ctx)

  repeat {
    best <- list(gain = min_gain_frac * abs(cur$score), apply = NULL)
    pool_of <- list(plus = plus, minus = minus)
    reservoir <- setdiff(seq_len(ctx$N), c(plus, minus))

    # Moves 1 and 2, vectorised over the reservoir.
    if (length(reservoir) > 0L) {
      for (side in c("plus", "minus")) {
        members <- pool_of[[side]]
        for (j in seq_len(n)) {
          out <- members[j]
          if (side == "plus") {
            cand_p <- bump(sp, out, reservoir)
            cand_m <- as.list(sm)
          } else {
            cand_p <- as.list(sp)
            cand_m <- bump(sm, out, reservoir)
          }
          sc <- .score_from_sums(cand_p, cand_m, ctx)
          ok <- cand_p$a >= m_min & cand_m$a >= m_min & abs(sc$dPr) <= tol
          if (!any(ok)) next
          gain <- ifelse(ok, sc$score - cur$score, -Inf)
          k <- which.max(gain)
          if (gain[k] > best$gain) {
            best <- list(gain = gain[k],
                         apply = list(side = side, j = j, inn = reservoir[k],
                                      cross = FALSE))
          }
        }
      }
    }

    # Move 3, vectorised over the As- members.
    for (j in seq_len(n)) {
      out_p <- plus[j]
      cand_p <- bump(sp, out_p, minus)
      cand_m <- bump(sm, minus, out_p)
      sc <- .score_from_sums(cand_p, cand_m, ctx)
      ok <- cand_p$a >= m_min & cand_m$a >= m_min & abs(sc$dPr) <= tol
      if (!any(ok)) next
      gain <- ifelse(ok, sc$score - cur$score, -Inf)
      k <- which.max(gain)
      if (gain[k] > best$gain) {
        best <- list(gain = gain[k],
                     apply = list(side = "plus", j = j, inn = minus[k],
                                  cross = TRUE, k = k))
      }
    }

    if (is.null(best$apply)) break
    mv <- best$apply
    if (isTRUE(mv$cross)) {
      out_p <- plus[mv$j]
      plus[mv$j]  <- mv$inn
      minus[mv$k] <- out_p
    } else if (mv$side == "plus") {
      plus[mv$j] <- mv$inn
    } else {
      minus[mv$j] <- mv$inn
    }
    sp <- .pool_sums(ctx, plus); sm <- .pool_sums(ctx, minus)
    cur <- .score_from_sums(sp, sm, ctx)
    n_moves <- n_moves + 1L
    if (n_moves >= max_moves) {
      warning("pool climb hit max_moves = ", max_moves,
              "; returning the best membership found so far", call. = FALSE)
      break
    }
  }

  # Feasible by construction, but an indexing bug would be silent, so assert.
  stopifnot(length(plus) == n, length(minus) == n,
            length(intersect(plus, minus)) == 0L,
            sum(anch[plus]) >= m_min, sum(anch[minus]) >= m_min)

  c(list(plus = plus, minus = minus, n_moves = n_moves), cur)
}

#' Search for the best pools, optionally anchored.
#'
#' Seeded from build_pools()'s own solution and from the theta family, which
#' is what makes `pools_optimal()$score >= pool_score(build_pools())` a
#' property of the construction rather than a hope.
#'
#' THE WINNING SEED IS REPORTED, as `theta` and `seed_kind`.  Every start is
#' reached from some point of the one-parameter theta family -- the
#' unconstrained family for an unanchored search, `.anchored_start()`'s for an
#' anchored one -- so `theta` remains defined in both modes and keeps its old
#' meaning in pool_summary.csv: the Pr multiplier whose ranking the adopted
#' membership came from.  Read it with `n_moves`: `n_moves == 0` says the theta
#' rule was already optimal and the climb changed nothing.
#'
#' @param starts Memberships to climb from in addition to the theta family.
#'   Each may carry `theta` and `kind`, which are passed through to the result.
pools_search <- function(ctx, ..., m_min = 0L, tol = Inf, n_seeds = 5L,
                         n_restarts = 0L, seed = 20261101L,
                         starts = list(), min_gain_frac = 1e-4) {
  if (...length() > 0L) {
    stop("pools_search() takes only `ctx` positionally; name the rest",
         call. = FALSE)
  }
  fam <- .theta_family(ctx)
  scored <- purrr::map_dbl(fam, \(s) {
    sc <- pool_score(s$plus, s$minus, ctx)
    if (sc$a_plus >= m_min && sc$a_minus >= m_min && abs(sc$dPr) <= tol) {
      sc$score
    } else -Inf
  })

  # theta = 0 always, then the best of the rest: the family is already
  # near-optimal, so these beat random starts outright.
  ord <- order(scored, decreasing = TRUE)
  pick <- unique(c(1L, utils::head(ord[is.finite(scored[ord])], n_seeds)))
  seeds <- c(purrr::map(starts, \(s) {
               s$kind <- s$kind %||% "supplied"
               s$theta <- s$theta %||% NA_real_
               s
             }),
             purrr::map(fam[pick], \(s) c(s[c("plus", "minus", "theta")],
                                          list(kind = "theta_family"))))
  seeds <- purrr::keep(seeds, \(s) {
    sc <- pool_score(s$plus, s$minus, ctx)
    sc$a_plus >= m_min && sc$a_minus >= m_min && abs(sc$dPr) <= tol
  })
  if (length(seeds) == 0L) return(NULL)

  if (n_restarts > 0L) {
    seeds <- c(seeds, .with_seed(seed, purrr::map(seq_len(n_restarts), \(i) {
      s <- seeds[[1L]]
      s$kind <- "restart"; s$theta <- NA_real_
      k <- sample.int(min(6L, ctx$n), 1L) + 1L
      res <- setdiff(seq_len(ctx$N), c(s$plus, s$minus))
      if (length(res) < 2L * k) return(s)
      s$plus[sample.int(ctx$n, k)]  <- sample(res, k)
      res <- setdiff(res, s$plus)
      s$minus[sample.int(ctx$n, k)] <- sample(res, k)
      s
    })))
  }

  res <- purrr::map(seeds, \(s) {
    r <- .pools_climb(s$plus, s$minus, ctx, m_min = m_min, tol = tol,
                      min_gain_frac = min_gain_frac)
    if (is.null(r)) return(NULL)
    c(r, list(theta = s$theta, seed_kind = s$kind))
  })
  res <- purrr::keep(res, \(r) !is.null(r))
  res[[which.max(purrr::map_dbl(res, \(r) r$score))]]
}

#' The pools.csv schema, from two lists of accession names.
#'
#' ONE DEFINITION of the adopted-pools table, used by build_pools() and by
#' both searches, because code/validate_design.R and code/validation_report.R
#' read pools.csv by column name and a second hand-written copy of this
#' transmute is exactly how they would come to disagree.
.pools_from_acc <- function(inp, acc_plus, acc_minus) {
  acc_tbl <- inp$accessions
  stopifnot(!anyDuplicated(acc_tbl$acc),
            length(intersect(acc_plus, acc_minus)) == 0L,
            all(c(acc_plus, acc_minus) %in% acc_tbl$acc))
  dplyr::bind_rows(
    dplyr::mutate(acc_tbl[match(acc_plus,  acc_tbl$acc), ], pool = "As+"),
    dplyr::mutate(acc_tbl[match(acc_minus, acc_tbl$acc), ], pool = "As-")
  ) |>
    dplyr::transmute(species = inp$species, pool, acc, Pr, As, GMA,
                     n_partners, n_plots, PEV_As_i, rel_As_i, PEV_GMA_i) |>
    dplyr::arrange(pool, dplyr::desc(As))
}

#' The unconstrained proposal: round 1.
#'
#' Reproduces build_pools() on the current vintage to four significant
#' figures (see the section note above). It exists to produce the round-1
#' membership, to supply the reference optimum the anchored search's cost is
#' measured against, and as a regression guard on the theta rule.
#'
#' @param P Total plots in the trial. The default 4 * n * n_loc is a
#'   CONVENIENCE for tests and happens to equal the real design for oat
#'   (4 * 20 * 5 = 400 = 80 plots x 5 locations); it overstates pea's count by
#'   half, because a plot is an oat x pea combination and the same 400 plots
#'   serve both species. Every driver passes P explicitly for that reason, and
#'   P does affect the ranking -- it is the plot-error term in the score's
#'   denominator.
pools_optimal <- function(inp, n, ..., pr_quantile = 0.5, tol = Inf,
                          b = NULL, n_seeds = 5L, n_restarts = 0L,
                          seed = 20261101L, P = NULL, n_loc = 5L) {
  if (...length() > 0L) {
    stop("pools_optimal() takes inp and n positionally; name the rest",
         call. = FALSE)
  }
  bp  <- build_pools(inp, n, pr_quantile = pr_quantile)
  ctx <- pool_search_context(bp$candidates, n = n, sigma2_e = inp$sigma2_e,
                             P = P %||% (4L * n * n_loc), b = b)
  seed_start <- list(list(
    plus  = match(bp$pools$acc[bp$pools$pool == "As+"], ctx$acc),
    minus = match(bp$pools$acc[bp$pools$pool == "As-"], ctx$acc),
    theta = bp$theta, kind = "build_pools"))
  out <- pools_search(ctx, tol = tol, n_seeds = n_seeds,
                      n_restarts = n_restarts, seed = seed,
                      starts = seed_start)
  c(out, list(ctx = ctx, build_pools = bp,
              acc_plus = ctx$acc[out$plus], acc_minus = ctx$acc[out$minus],
              pools = .pools_from_acc(inp, ctx$acc[out$plus],
                                      ctx$acc[out$minus]),
              search_mode = "unconstrained",
              n_anchors_usable = 0L, n_new = 2L * n))
}

#' Can the anchor constraint be satisfied, and if not, what would it take?
#'
#' THE ARITHMETIC. The two pools are disjoint, so their anchor counts are
#' distinct and sum to at most the number of usable anchors `A`. Each pool
#' needs `m = n - max_new_per_pool` of them, so the constraint is
#' `2 * m <= A`.
#'
#' The floor in `max_new_min` is INTEGER, not `n - A/2`, because a per-pool
#' count is a whole number: with A = 33 and n = 20 each pool still needs 17,
#' and 2 * 17 = 34 > 33.
#'
#' Separated from pools_anchored() so the arithmetic has one definition and
#' can be checked without running a search.
#'
#' @param n_usable_anchors Anchors that are actually candidates -- a found
#'   accession that is ineligible or below the Pr threshold cannot be used, so
#'   it does not count here.
pool_anchor_feasibility <- function(n, n_usable_anchors, max_new_per_pool) {
  n <- as.integer(n)
  A <- as.integer(n_usable_anchors)
  m <- n - as.integer(max_new_per_pool)
  list(feasible    = 2L * m <= A,
       m_min       = m,
       required    = 2L * m,
       n_usable    = A,
       max_new_min = max(0L, n - A %/% 2L),
       shortfall   = max(0L, 2L * m - A))
}

#' A feasible anchored starting membership, built from one theta.
#'
#' The theta index generalises to the anchored problem: ONE ordering both
#' splits the anchors between the pools and fills the remaining slots from the
#' reservoir, so the whole start is built by a single sort. Scanning theta then
#' gives a family of feasible starts to climb from, exactly as .pool_scan()
#' does for the unconstrained case.
#'
#' This is not a formality. A naive start -- rank the anchors on As alone and
#' fill the rest by As -- produces an INFEASIBLE membership whenever the
#' anchors are concentrated at one end, because the two halves can collide.
#' Returns NULL when this theta cannot produce a feasible start.
#'
#' @param m Minimum anchors required in each pool (n - max_new_per_pool).
.anchored_start <- function(ctx, theta, m) {
  n   <- ctx$n
  idx <- ctx$As + theta * ctx$Pr
  anch <- which(ctx$is_anchor)
  if (length(anch) < 2L * m) return(NULL)

  if (m > 0L) {
    ao    <- anch[order(idx[anch], decreasing = TRUE)]
    plus  <- ao[seq_len(m)]
    minus <- ao[seq.int(length(ao) - m + 1L, length(ao))]
  } else {
    plus <- integer(0); minus <- integer(0)
  }

  k <- n - m
  if (k > 0L) {
    rest <- setdiff(seq_len(ctx$N), c(plus, minus))
    if (length(rest) < 2L * k) return(NULL)
    ro    <- rest[order(idx[rest], decreasing = TRUE)]
    plus  <- c(plus,  ro[seq_len(k)])
    minus <- c(minus, ro[seq.int(length(ro) - k + 1L, length(ro))])
  }
  if (length(intersect(plus, minus)) > 0L) return(NULL)
  list(plus = plus, minus = minus, theta = theta)
}

# ------------------------------------------------------------
# THE PARSIMONY STAGE
#
# WHY THIS EXISTS. Maximising the score alone is the WRONG SPECIFICATION for
# this problem, and measurably so. The greedy anchored optimum will spend its
# entire substitution allowance whenever a further swap gains anything at all:
# prototyping found it asking for 11 new accessions instead of 8 to gain 0.03
# power points. Finding an accession in seed storage is costly manual labour,
# so that trade is the exact inverse of what the protocol is for.
#
# So the search runs in two stages. Stage A climbs to the optimum. Stage B
# then hands slots BACK -- replacing non-anchor members with anchors already
# in hand -- for as long as the score stays within `score_slack_frac` of the
# optimum. Both answers are reported, so the choice is visible rather than
# made silently: "find 4 and gain nothing, or find 11 and gain 0.03 points".
#
# WHY THE DEFAULT SLACK IS 0.005. It is 0.5% of t, which is roughly 0.3 power
# points. On the 2026-10-31 vintage lambda is 0.759 +/- 0.198 for oat (26%
# relative) and 0.784 +/- 0.140 for pea, so the contrast itself is uncertain
# to about +/- 15 power points; associate reliability is 0.29 and 0.36. A
# difference of 0.3 power points is therefore far below the sampling error of
# the quantity being predicted, and is not worth an afternoon in a cold room.
# Raising the slack trades more of a difference you cannot measure for fewer
# bags to find.
#
# Stage B's moves are always anchor-feasible by construction -- they only ever
# RAISE a pool's anchor count -- so the constraint cannot be broken here.
# ------------------------------------------------------------

#' Hand substitution slots back while the score holds up.
#'
#' @param score_star The stage-A optimum, which the slack is measured against.
#' @return The pruned membership, or the input unchanged if nothing qualified.
.pool_prune <- function(plus, minus, ctx, ..., score_star,
                        score_slack_frac = 0.005, tol = Inf) {
  if (...length() > 0L) {
    stop("'.pool_prune()' takes plus, minus and ctx positionally; name the ",
         "rest", call. = FALSE)
  }
  floor_score <- score_star * (1 - score_slack_frac)
  n_steps <- 0L

  repeat {
    sp <- .pool_sums(ctx, plus); sm <- .pool_sums(ctx, minus)
    spare <- setdiff(which(ctx$is_anchor), c(plus, minus))
    if (length(spare) == 0L) break

    best <- NULL
    for (side in c("plus", "minus")) {
      members <- if (side == "plus") plus else minus
      for (j in which(!ctx$is_anchor[members])) {
        out <- members[j]
        if (side == "plus") {
          cand_p <- .sums_bump(sp, out, spare, ctx); cand_m <- as.list(sm)
        } else {
          cand_p <- as.list(sp); cand_m <- .sums_bump(sm, out, spare, ctx)
        }
        sc <- .score_from_sums(cand_p, cand_m, ctx)
        ok <- sc$score >= floor_score & abs(sc$dPr) <= tol
        if (!any(ok)) next
        # Smallest loss = highest resulting score among those that clear the
        # floor, so the cheapest slot is always handed back first.
        k <- which.max(ifelse(ok, sc$score, -Inf))
        if (is.null(best) || sc$score[k] > best$score) {
          best <- list(score = sc$score[k], side = side, j = j,
                       inn = spare[k])
        }
      }
    }
    if (is.null(best)) break
    if (best$side == "plus") plus[best$j] <- best$inn
    else                     minus[best$j] <- best$inn
    n_steps <- n_steps + 1L
  }

  stopifnot(length(intersect(plus, minus)) == 0L)
  list(plus = plus, minus = minus, n_pruned = n_steps)
}

#' The anchored proposal: round 2 onward.
#'
#' Holds the accessions already found in seed storage and re-optimises around
#' them, so each round asks for only a handful of new names rather than a new
#' pool to go hunting for.
#'
#' ANCHORS ARE ONE GLOBAL LIST, which the optimiser may reassign between the
#' pools, rather than per-pool quotas. The cost being minimised is
#' pool-agnostic -- a found bag is found -- and the As+/As- label is a derived
#' property of the design rather than a property of the seed lot, so forbidding
#' reassignment would prune the feasible set for no scientific gain. It also
#' collapses feasibility to one inequality instead of two coupled ones. The
#' price is that an anchor can change pool between rounds, so the caller is
#' told which ones did.
#'
#' POOL LABELS ARE ASSIGNED BY OUTCOME, never inherited from the previous
#' round: As+ is whichever pool ends up with the higher mean As.
#'
#' @param anchors Accessions confirmed present in seed storage.
#' @param exclude Accessions confirmed ABSENT. Removed from the candidate set
#'   before anything else -- before the feasibility count and before the
#'   As-on-Pr slope -- so the slope describes the set actually selected from.
#' @param max_new_per_pool How many members of each pool may come from outside
#'   `anchors`, i.e. how many new accessions you are willing to go and find.
pools_anchored <- function(inp, n, anchors, ..., exclude = character(0),
                           max_new_per_pool = 3L, pr_quantile = 0.5,
                           tol = Inf, n_seeds = 5L, n_restarts = 0L,
                           seed = 20261101L, P = NULL, n_loc = 5L, b = NULL,
                           parsimony = TRUE, score_slack_frac = 0.005,
                           reference = TRUE) {
  if (...length() > 0L) {
    stop("pools_anchored() takes inp, n and anchors positionally; name the ",
         "rest", call. = FALSE)
  }
  max_new_per_pool <- as.integer(max_new_per_pool)
  if (max_new_per_pool < 0L || max_new_per_pool > n) {
    stop("max_new_per_pool must be between 0 and n = ", n, call. = FALSE)
  }
  anchors <- unique(anchors)

  bp  <- build_pools(inp, n, pr_quantile = pr_quantile)
  P   <- P %||% (4L * n * n_loc)
  ctx <- pool_search_context(bp$candidates, n = n, sigma2_e = inp$sigma2_e,
                             P = P, b = b, anchors = anchors,
                             exclude = exclude)

  A   <- sum(ctx$is_anchor)
  fs  <- pool_anchor_feasibility(n, A, max_new_per_pool)
  m   <- fs$m_min
  max_new_min <- fs$max_new_min
  if (!fs$feasible) {
    stop(inp$species, ": ", A, " of ", length(anchors), " anchor accession(s) ",
         "are usable candidates, but max_new_per_pool = ", max_new_per_pool,
         " requires 2 * (", n, " - ", max_new_per_pool, ") = ", fs$required,
         " of them.\n  The smallest feasible allowance is max_new_per_pool = ",
         max_new_min, ", which asks you to find ", 2L * max_new_min,
         " new accession(s).\n  Alternatively, find ", fs$shortfall,
         " more of the accessions already proposed.", call. = FALSE)
  }

  starts <- purrr::compact(purrr::map(seq(0, 10, length.out = 201L),
                                      \(th) .anchored_start(ctx, th, m)))
  keys <- purrr::map_chr(starts, \(s)
    paste(c(sort(s$plus), -1L, sort(s$minus)), collapse = ","))
  starts <- starts[!duplicated(keys)]

  out <- pools_search(ctx, m_min = m, tol = tol, n_seeds = n_seeds,
                      n_restarts = n_restarts, seed = seed, starts = starts)

  # A SECOND FEASIBILITY CONSTRAINT, which anchor feasibility does not imply:
  # the anchors can make |dPr| <= tol unreachable. Refusing would leave no
  # proposal at all, so proceed and say so -- loudly, and as a returned FLAG
  # rather than only a warning, because a warning in a long console log is
  # exactly what gets missed.
  relaxed <- FALSE
  if (is.null(out) && is.finite(tol)) {
    warning(inp$species, ": no anchored membership reaches |dPr| <= ", tol,
            ". Proceeding with the Pr tolerance relaxed; check dPr on the ",
            "result before adopting it.", call. = FALSE)
    relaxed <- TRUE
    out <- pools_search(ctx, m_min = m, tol = Inf, n_seeds = n_seeds,
                        n_restarts = n_restarts, seed = seed, starts = starts)
  }
  if (is.null(out)) {
    stop(inp$species, ": the anchor constraint is satisfiable by the counts (",
         A, " usable anchors for 2 * ", m, " required) but no feasible ",
         "starting membership could be built. This is a bug in ",
         ".anchored_start(), not a property of the data.", call. = FALSE)
  }

  # STAGE A's answer, kept so the trade is reportable rather than silent.
  greedy <- out
  n_new_greedy <- sum(!ctx$is_anchor[out$plus]) +
                  sum(!ctx$is_anchor[out$minus])

  # STAGE B: hand substitution slots back while the score holds up.
  n_pruned <- 0L
  if (isTRUE(parsimony)) {
    pr <- .pool_prune(out$plus, out$minus, ctx, score_star = out$score,
                      score_slack_frac = score_slack_frac,
                      tol = if (relaxed) Inf else tol)
    n_pruned <- pr$n_pruned
    if (n_pruned > 0L) {
      out <- c(list(plus = pr$plus, minus = pr$minus, n_moves = out$n_moves,
                    theta = out$theta, seed_kind = out$seed_kind),
               .score_from_sums(.pool_sums(ctx, pr$plus),
                                .pool_sums(ctx, pr$minus), ctx))
    }
  }

  # Labels by outcome, so a round never inherits the previous round's pools.
  if (out$dAs < 0) {
    sw <- out$plus; out$plus <- out$minus; out$minus <- sw
    out <- c(list(plus = out$plus, minus = out$minus, n_moves = out$n_moves,
                  theta = out$theta, seed_kind = out$seed_kind),
             .score_from_sums(.pool_sums(ctx, out$plus),
                              .pool_sums(ctx, out$minus), ctx))
  }
  stopifnot(out$dAs > 0)

  ref <- if (isTRUE(reference)) {
    pools_optimal(inp, n, pr_quantile = pr_quantile, tol = tol, b = b,
                  n_seeds = n_seeds, seed = seed, P = P, n_loc = n_loc)
  } else NULL

  new_plus  <- !ctx$is_anchor[out$plus]
  new_minus <- !ctx$is_anchor[out$minus]

  c(out,
    list(ctx = ctx, build_pools = bp,
         acc_plus = ctx$acc[out$plus], acc_minus = ctx$acc[out$minus],
         pools = .pools_from_acc(inp, ctx$acc[out$plus], ctx$acc[out$minus]),
         search_mode = "anchored",
         max_new_per_pool = max_new_per_pool, m_min = m,
         max_new_min = max_new_min,
         # OPERATING AT THE MINIMUM FEASIBLE ALLOWANCE IS THE EXPENSIVE
         # CHOICE, and it is easy to land on by accident since the minimum is
         # what the infeasibility error advises. At the minimum every usable
         # anchor must be used, so the search cannot drop a poorly-estimated
         # one: measured cost 1.3 (oat) to 1.9 (pea) power points, recovered
         # in full by one or two extra units. The caller is told so it can say
         # so.
         at_min_budget = identical(max_new_per_pool, max_new_min),
         n_anchors_supplied = length(anchors), n_anchors_usable = A,
         n_anchors_used = sum(!new_plus) + sum(!new_minus),
         n_new = sum(new_plus) + sum(new_minus),
         n_new_plus = sum(new_plus), n_new_minus = sum(new_minus),
         acc_new = c(ctx$acc[out$plus][new_plus], ctx$acc[out$minus][new_minus]),
         # Both stages reported: "find n_new and score, or find n_new_greedy
         # and score_greedy". The whole point of the parsimony stage is that
         # this choice is visible.
         n_new_greedy = n_new_greedy, score_greedy = greedy$score,
         n_pruned = n_pruned, score_slack_frac = score_slack_frac,
         dpr_tol_relaxed = relaxed,
         score_unconstrained = if (is.null(ref)) NA_real_ else ref$score,
         score_frac_of_unconstrained =
           if (is.null(ref)) NA_real_ else out$score / ref$score))
}

# ------------------------------------------------------------
# WHAT THE SEED STORE SAYS, AND WHAT THE FIT CAN USE
#
# A round's input is a list of accessions found or missing, typed by a person
# reading bag labels. Two things stand between that list and the search.
#
# FIRST, NAMES. Curation pools entries it finds not to segregate under a
# single synthetic ANALYSIS NAME -- IL18-735_IL17-10306_no_cross stands for 13
# original entries (VALIDATION_DESIGN.md section 7c). The seed store does not
# hold a bag labelled with the synthetic name; it holds bags labelled with the
# ORIGINAL germplasm names. So the operationally important direction is the
# reverse of the one section 7c first suggests: a name read off a bag has to
# be resolved FORWARD to the analysis name the fit knows. Resolving it does
# not change the design -- it only recognises that seed which is present was
# reported as absent. Conversely, a proposal naming a collapsed accession
# should tell whoever is searching which member entries would satisfy it, so
# the request reads "find any one of these 13".
#
# SECOND, ELIGIBILITY. A found accession is not automatically usable: it may
# be below the reliability floor, below the partner floor, or below the Pr
# threshold. Those are reported, never quietly admitted -- see
# classify_search_input().
# ------------------------------------------------------------

#' Resolve names read off seed bags to the accession names the fit knows.
#'
#' @param names Names as typed, e.g. from a seed-store inventory.
#' @param known The accession names the fit carries (`inp$accessions$acc`).
#' @param map Two columns, `germplasmName` and `analysis_name`. Injectable so
#'   this is testable without `output/`; read from
#'   `output/<species>_analysis_names.csv` when NULL.
#' @return One row per input name, with the route taken. `identity` means the
#'   name was already an accession the fit knows; `member_entry` means it was
#'   an original entry pooled under an analysis name.
resolve_seed_names <- function(names, ..., known, map = NULL,
                               species = NULL) {
  if (...length() > 0L) {
    stop("resolve_seed_names() takes only `names` positionally; name the rest",
         call. = FALSE)
  }
  if (is.null(map)) {
    if (is.null(species)) {
      stop("resolve_seed_names() needs `map` or `species`", call. = FALSE)
    }
    f <- here::here("output", paste0(species, "_analysis_names.csv"))
    if (!file.exists(f)) {
      stop("no name map at ", f, ". Run code/curate_", species,
           "_accessions.R, or pass `map` directly.", call. = FALSE)
    }
    map <- readr::read_csv(f, show_col_types = FALSE)
  }
  if (!all(c("germplasmName", "analysis_name") %in% names(map))) {
    stop("`map` needs columns germplasmName and analysis_name", call. = FALSE)
  }
  lookup <- stats::setNames(map$analysis_name, map$germplasmName)

  tibble::tibble(input_name = unique(as.character(names))) |>
    dplyr::mutate(
      mapped = unname(lookup[input_name]),
      # Order matters: a name the fit already knows is taken AS IS, because an
      # analysis name can coincide with a germplasm name (the `<seed_parent>`
      # case) and remapping it would be wrong.
      route = dplyr::case_when(
        input_name %in% known                      ~ "identity",
        !is.na(mapped) & mapped %in% known         ~ "member_entry",
        TRUE                                       ~ "unresolved"),
      resolved_acc = dplyr::case_when(
        route == "identity"     ~ input_name,
        route == "member_entry" ~ mapped,
        TRUE                    ~ NA_character_)) |>
    dplyr::select(input_name, resolved_acc, route)
}

#' Account for every name a round's search results supplied.
#'
#' Exactly one row per unique input name, so nothing can be silently dropped.
#' That is the point: an accession the user believes they reported as found,
#' which quietly never reached the search, is the worst failure this step can
#' have.
#'
#' THE Pr THRESHOLD IS NOT RELAXED to admit a found accession that sits below
#' it. The threshold is part of the hypothesis -- the pools are contrasted
#' "among good producers" -- so admitting a sub-threshold line would change
#' what the trial tests. It is reported with the shortfall in g/m2 instead, so
#' lowering `pr_quantile` stays a deliberate, recorded act.
classify_search_input <- function(inp, names, ..., pr_quantile = 0.5,
                                  map = NULL) {
  if (...length() > 0L) {
    stop("classify_search_input() takes inp and names positionally; name the ",
         "rest", call. = FALSE)
  }
  pc   <- .pool_candidates(inp, pr_quantile)
  acc  <- inp$accessions
  el   <- inp$eligibility

  res <- resolve_seed_names(names, known = acc$acc, map = map,
                            species = inp$species)

  res |>
    dplyr::left_join(
      dplyr::select(acc, resolved_acc = acc, Pr, As, rel_As_i, n_partners,
                    eligible),
      by = "resolved_acc") |>
    dplyr::mutate(
      below_rel     = !is.na(rel_As_i) & rel_As_i < el$rel_min,
      below_partner = !is.na(n_partners) & n_partners < el$partner_floor,
      pr_shortfall  = dplyr::if_else(!is.na(Pr) & Pr < pc$pr_min,
                                     pc$pr_min - Pr, NA_real_),
      bucket = dplyr::case_when(
        route == "unresolved"                 ~ "not_in_accession_table",
        resolved_acc %in% pc$cand$acc         ~ "in_candidates",
        !eligible                             ~ "ineligible",
        TRUE                                  ~ "eligible_but_below_Pr_threshold"),
      reason = dplyr::case_when(
        bucket == "not_in_accession_table" ~
          "no such accession in the fit, and not an entry pooled under one: check the name",
        bucket == "in_candidates" ~ NA_character_,
        bucket == "eligible_but_below_Pr_threshold" ~
          sprintf("Pr %.2f is %.2f g/m2 below the threshold of %.2f",
                  Pr, pr_shortfall, pc$pr_min),
        below_rel & below_partner ~
          sprintf("associate reliability %.3f < %.3f and %d partner(s) < %d",
                  rel_As_i, el$rel_min, n_partners, el$partner_floor),
        below_rel ~ sprintf("associate reliability %.3f < %.3f",
                            rel_As_i, el$rel_min),
        below_partner ~ sprintf("%d partner(s) < %d",
                                n_partners, el$partner_floor),
        TRUE ~ "ineligible"),
      species = inp$species, .before = 1) |>
    dplyr::select(species, input_name, resolved_acc, route, bucket, reason,
                  Pr, As, rel_As_i, n_partners, pr_shortfall)
}

#' Which member entries would satisfy a proposed accession.
#'
#' The forward direction of the name problem: a collapsed analysis name has no
#' bag of its own, so a search request naming one must say which original
#' entries count. Accessions that are not collapsed come back with
#' `n_entries = 1` and themselves as the entry, so a caller can join on this
#' unconditionally.
#'
#' @param collapsed Injectable `output/summary_collapsed_names.csv`, so this is
#'   testable without `output/`.
collapsed_entries <- function(accs, species, collapsed = NULL) {
  if (is.null(collapsed)) {
    f <- here::here("output", "summary_collapsed_names.csv")
    collapsed <- if (file.exists(f)) {
      readr::read_csv(f, show_col_types = FALSE)
    } else {
      tibble::tibble(species = character(), analysis_name = character(),
                     n_entries = integer(), entries = character())
    }
  }
  tibble::tibble(acc = unique(as.character(accs))) |>
    dplyr::left_join(
      collapsed |>
        dplyr::filter(species == !!species) |>
        dplyr::select(acc = analysis_name, n_entries, entries),
      by = "acc") |>
    dplyr::mutate(
      is_collapsed = !is.na(n_entries),
      n_entries = dplyr::coalesce(as.integer(n_entries), 1L),
      entries = dplyr::coalesce(entries, acc))
}

#' What each extra substitution buys: the marginal-value curve.
#'
#' THE DECISION THIS ANSWERS is "is another search round worth it?", which
#' otherwise has to be probed for by hand. One anchored run per allowance,
#' from the smallest feasible one up, each about a tenth of a second, gives
#' the whole trade-off as one monotone column.
#'
#' Worth knowing before reading it: the SMALLEST FEASIBLE allowance is
#' usually the worst place to operate. At the minimum every usable anchor must
#' be used, so the search has no freedom to drop a poorly-estimated one, and
#' on the 2026-10-31 vintage that costs 1.3 (oat) to 1.9 (pea) power points
#' which one or two extra units of allowance recover in full.
#'
#' @param lambda,interaction_frac Used only for the REPORTED power column, not
#'   for selection -- the chosen membership is invariant to both.
#' @return One row per allowance. `power` is the Pr-adjusted power at the
#'   supplied lambda; `power_unconstrained` is the same for the unanchored
#'   optimum, which is the ceiling the constraint is measured against.
pool_search_curve <- function(inp, n, anchors, ..., exclude = character(0),
                              lambda = 1, interaction_frac = 0,
                              max_new_values = NULL, pr_quantile = 0.5,
                              tol = Inf, P = NULL, n_loc = 5L, b = NULL,
                              parsimony = TRUE, score_slack_frac = 0.005,
                              seed = 20261101L) {
  if (...length() > 0L) {
    stop("pool_search_curve() takes inp, n and anchors positionally; name ",
         "the rest", call. = FALSE)
  }
  P   <- P %||% (4L * n * n_loc)
  ref <- pools_optimal(inp, n, pr_quantile = pr_quantile, tol = tol, b = b,
                       seed = seed, P = P, n_loc = n_loc)
  ctx0 <- pool_search_context(ref$build_pools$candidates, n = n,
                              sigma2_e = inp$sigma2_e, P = P, b = b,
                              anchors = anchors, exclude = exclude)
  A   <- sum(ctx0$is_anchor)
  fs0 <- pool_anchor_feasibility(n, A, 0L)

  pw <- function(s, ctx) {
    contrast_power(dAs = s$dAs,
                   sigma2_within = c(s$sigma2_within_plus,
                                     s$sigma2_within_minus),
                   n = n, sigma2_e = inp$sigma2_e, P = P, lambda = lambda,
                   interaction_frac = interaction_frac, n_loc = n_loc,
                   pr_adjust = list(dPr = s$dPr, b = ctx$b, Sxx = s$Sxx))$power
  }
  power_ref <- pw(ref, ref$ctx)

  # DEFAULT RANGE. The curve is flat almost immediately: on the 2026-10-31
  # vintage, oat reaches the unconstrained ceiling one unit above the minimum
  # feasible allowance and never asks for more new accessions however wide the
  # allowance gets, and pea is at the ceiling from the minimum. Sweeping to
  # n costs ~30 s (oat) to ~70 s (pea) to show a flat line, so the default
  # stops five units past the minimum; pass max_new_values to go further.
  vals <- max_new_values %||%
    seq.int(fs0$max_new_min, min(n, fs0$max_new_min + 5L))
  purrr::map(vals, \(mn) {
    fs <- pool_anchor_feasibility(n, A, mn)
    if (!fs$feasible) {
      return(tibble::tibble(
        species = inp$species, max_new_per_pool = as.integer(mn),
        feasible = FALSE, n_new = NA_integer_, n_new_greedy = NA_integer_,
        score = NA_real_, power = NA_real_, dAs = NA_real_, dPr = NA_real_,
        dAs_adj = NA_real_, at_min_budget = NA, dpr_tol_relaxed = NA))
    }
    r <- pools_anchored(inp, n, anchors, exclude = exclude,
                        max_new_per_pool = mn, pr_quantile = pr_quantile,
                        tol = tol, seed = seed, P = P, n_loc = n_loc, b = b,
                        parsimony = parsimony,
                        score_slack_frac = score_slack_frac,
                        reference = FALSE)
    tibble::tibble(
      species = inp$species, max_new_per_pool = as.integer(mn),
      feasible = TRUE, n_new = r$n_new, n_new_greedy = r$n_new_greedy,
      score = r$score, power = pw(r, r$ctx),
      dAs = r$dAs, dPr = r$dPr, dAs_adj = r$dAs_adj,
      at_min_budget = r$at_min_budget,
      dpr_tol_relaxed = r$dpr_tol_relaxed)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      n_anchors_usable = A, n_anchors_supplied = length(unique(anchors)),
      max_new_min = fs0$max_new_min,
      score_unconstrained = ref$score, power_unconstrained = power_ref,
      delta_power_vs_unconstrained = power - power_ref,
      score_frac_of_unconstrained = score / ref$score,
      .after = feasible)
}

# ------------------------------------------------------------
# lambda AND THE POOL x LOCATION INTERACTION, AS MEASURED
#
# code/validate_crossval.R writes crossval_summary.csv into a vintage folder,
# and every power figure in this project is conditional on the lambda it
# reports.  An older vintage is used when the current one has none -- stated
# out loud by the caller, because a power figure carrying last month's lambda
# is a different claim from one carrying this month's.
#
# "informative folds only" is the set used: the other set includes folds where
# the held-out trial had too little variation for its slope to say anything
# about lambda, and averaging those in pulls the estimate toward noise.
#
# WHAT LAMBDA DOES AND DOES NOT AFFECT.  It scales the contrast, so it moves
# every reported power.  It does NOT affect which accessions are selected:
# the pool search maximises a lambda-free score, and the chosen membership is
# invariant to lambda across the whole plausible range (see THE POOL SEARCH).
# So a stale lambda misstates the number printed beside a proposal and
# provably not the proposal itself.
# ------------------------------------------------------------

FALLBACK <- list(lambda = 0.8, interaction_frac = 0.5)

#' lambda and interaction_frac from the most recent cross-validation.
#'
#' @return `est` (one row per species) and `from` (the vintage it came from),
#'   or NULL when no vintage anywhere has a crossval_summary.csv.
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

# ------------------------------------------------------------
# THE ROUND PROTOCOL
#
# Round 1 proposes the unconstrained optimum.  Somebody then walks to the seed
# store with search_request.csv, and records what they found.  Round 2 holds
# the found accessions and re-optimises around them, asking for a short list
# of replacements for the ones that were missing.
#
# THE LEDGER IS PRIMARY DATA AND LIVES IN data/, NOT IN A VINTAGE FOLDER.
# Everything under output/validation/<vintage>/ is derived and regenerated on
# demand; the ledger is hand-typed and cannot be, so it follows
# data/trial_qc_manual.csv -- a small CSV under version control that a person
# edits and a script only ever reads and appends to.
#
# Three statuses, each with an operational meaning:
#
#   found          -> an ANCHOR.  The search may keep it.
#   missing        -> EXCLUDED, cumulatively across rounds, so an accession
#                     confirmed absent is never proposed again.
#   missing_retry  -> neither.  "I did not look hard enough" or "it may be in
#                     the other freezer": not usable as an anchor this round,
#                     but still a candidate the search may propose again.
#
# `exclude` is DERIVED from `status`, never typed, so the two cannot disagree.
# ------------------------------------------------------------

SEED_LEDGER_STATUSES <- c("unsearched", "found", "missing", "missing_retry")

#' An empty ledger, which is also the schema.
#'
#' `species`, `acc`, `round_proposed` and `pool_proposed` are written by
#' code/validate_pool_selection.R; `status`, `searched_on`, `sowable_entry` and
#' `note` are written by whoever did the searching.
.empty_seed_ledger <- function() {
  tibble::tibble(species = character(), acc = character(),
                 round_proposed = integer(), pool_proposed = character(),
                 status = character(), searched_on = character(),
                 sowable_entry = character(), note = character())
}

#' Read the seed-search ledger, or return an empty one.
#'
#' Every column is read as text and then coerced, because readr's type
#' guessing on a mostly-empty hand-edited column is a source of surprises --
#' an all-blank `searched_on` guesses logical, and then one typed date fails.
#'
#' Two things are checked rather than assumed, because both would corrupt the
#' search silently: an unrecognised `status` (a typo would be read as
#' "not found", quietly dropping an anchor) and a duplicated species x acc
#' (which would double-count one bag).
read_seed_ledger <- function(path = NULL) {
  path <- path %||% here::here("data", "seed_search_ledger.csv")
  if (!file.exists(path)) return(.empty_seed_ledger())

  x <- readr::read_csv(path, col_types = readr::cols(.default = "c"),
                       progress = FALSE)
  need <- c("species", "acc", "status")
  miss <- setdiff(need, names(x))
  if (length(miss) > 0L) {
    stop(path, " is missing the column(s) ", paste(miss, collapse = ", "),
         ". The ledger needs at least species, acc and status.", call. = FALSE)
  }
  for (col in setdiff(names(.empty_seed_ledger()), names(x))) {
    x[[col]] <- NA_character_
  }
  x <- x |>
    dplyr::mutate(
      species = trimws(species), acc = trimws(acc),
      status = tolower(trimws(status)),
      status = dplyr::if_else(is.na(status) | status == "",
                              "unsearched", status),
      round_proposed = suppressWarnings(as.integer(round_proposed))) |>
    dplyr::select(dplyr::all_of(names(.empty_seed_ledger())),
                  dplyr::everything())

  bad <- setdiff(unique(x$status), SEED_LEDGER_STATUSES)
  if (length(bad) > 0L) {
    stop(path, ": unrecognised status ", paste0('"', bad, '"', collapse = ", "),
         ". Allowed: ", paste(SEED_LEDGER_STATUSES, collapse = ", "),
         ".", call. = FALSE)
  }
  dup <- x |> dplyr::count(species, acc) |> dplyr::filter(n > 1L)
  if (nrow(dup) > 0L) {
    stop(path, ": ", nrow(dup), " accession(s) appear more than once -- ",
         paste(paste0(dup$species, "/", dup$acc), collapse = ", "),
         ". One row per accession; edit the existing row rather than adding ",
         "another.", call. = FALSE)
  }
  x
}

#' What the ledger says about one species, in the terms the search uses.
#'
#' @return `anchors` (found), `exclude` (missing, cumulative over all rounds),
#'   `retry` and `unsearched` -- neither anchored nor excluded -- and the
#'   round number the ledger has reached.
ledger_state <- function(ledger, species) {
  sp <- dplyr::filter(ledger, species == !!species)
  by <- \(st) sp$acc[sp$status %in% st]
  list(anchors = by("found"), exclude = by("missing"),
       retry = by("missing_retry"), unsearched = by("unsearched"),
       n_rows = nrow(sp),
       n_searched = sum(sp$status != "unsearched"),
       last_round = if (nrow(sp) == 0L) 0L else
         max(sp$round_proposed, 0L, na.rm = TRUE))
}

#' The round this proposal belongs to.
#'
#' A ROUND IS OPEN WHILE ANY ACCESSION IT ASKED FOR IS STILL UNSEARCHED, and a
#' re-run inside an open round is the same round, not the next one.  Counting
#' instead on "has anything ever been searched" advances on every re-run once
#' the first round is answered, so the ledger ends up claiming searches that
#' never happened.
seed_ledger_round <- function(ledger) {
  if (nrow(ledger) == 0L) return(1L)
  last <- max(ledger$round_proposed, 0L, na.rm = TRUE)
  open <- any(ledger$round_proposed == last & ledger$status == "unsearched",
              na.rm = TRUE)
  if (open) max(last, 1L) else last + 1L
}

#' Add this round's proposal to the ledger, without touching what was typed.
#'
#' APPEND-ONLY BY CONSTRUCTION.  A row whose `status` is not "unsearched"
#' carries a person's observation and is returned byte for byte; the only
#' in-place change is `pool_proposed` on a still-unsearched row, so the
#' request stays current if the optimiser moved the accession to the other
#' pool.  `round_proposed` is left at the round the accession was FIRST asked
#' for, which is what makes "this has been outstanding since round 1" legible.
#'
#' @return `ledger` (merged), `added` (the new rows) and `n_repooled`.
update_seed_ledger <- function(ledger, pools, ..., round) {
  if (...length() > 0L) {
    stop("update_seed_ledger() takes ledger and pools positionally; name the ",
         "rest", call. = FALSE)
  }
  cols <- names(.empty_seed_ledger())
  prop <- pools |>
    dplyr::transmute(species, acc, pool_now = pool,
                     round_proposed = as.integer(round))

  kept <- ledger |>
    dplyr::left_join(dplyr::select(prop, species, acc, pool_now),
                     by = c("species", "acc")) |>
    dplyr::mutate(
      repooled = status == "unsearched" & !is.na(pool_now) &
                 !identical(pool_proposed, pool_now) &
                 (is.na(pool_proposed) | pool_proposed != pool_now),
      pool_proposed = dplyr::if_else(status == "unsearched" & !is.na(pool_now),
                                     pool_now, pool_proposed))
  n_repooled <- sum(kept$repooled, na.rm = TRUE)
  kept <- dplyr::select(kept, -pool_now, -repooled)

  added <- prop |>
    dplyr::anti_join(ledger, by = c("species", "acc")) |>
    dplyr::transmute(species, acc, round_proposed,
                     pool_proposed = pool_now, status = "unsearched",
                     searched_on = NA_character_,
                     sowable_entry = NA_character_, note = NA_character_)

  out <- dplyr::bind_rows(kept, added) |>
    dplyr::arrange(species, round_proposed, acc)

  # The append-only promise, asserted rather than trusted: every row a person
  # has touched comes back unchanged.
  typed <- dplyr::filter(ledger, status != "unsearched")
  if (nrow(typed) > 0L) {
    back <- dplyr::semi_join(out, typed, by = c("species", "acc")) |>
      dplyr::arrange(species, acc)
    stopifnot(identical(
      as.data.frame(back[, cols]),
      as.data.frame(dplyr::arrange(typed, species, acc)[, cols])))
  }
  list(ledger = out, added = added, n_repooled = n_repooled)
}

write_seed_ledger <- function(ledger, path = NULL) {
  path <- path %||% here::here("data", "seed_search_ledger.csv")
  readr::write_csv(ledger, path, na = "")
  invisible(path)
}

#' Is what the store happens to hold biased in As or Pr?
#'
#' THE RISK THIS MEASURES.  Anchoring entrenches whatever the seed store has.
#' If storage is systematically short of, say, old low-As lines, then every
#' round after the first inherits that, and the As- pool ends up describing
#' the store's history rather than the trait.  It is a caveat on the
#' interpretation, not a bug in the search, so it is REPORTED every round and
#' never acted on automatically.
#'
#' `missing_retry` is deliberately left out: it records a search that was not
#' finished, which is evidence about the searcher rather than about the store.
#'
#' @param cls A classify_search_input() table whose `input_name` matches the
#'   ledger's `acc`; As and Pr come from the fit, not from the ledger.
found_missing_balance <- function(ledger, cls) {
  d <- ledger |>
    dplyr::filter(status %in% c("found", "missing")) |>
    dplyr::inner_join(dplyr::select(cls, species, acc = input_name, Pr, As),
                      by = c("species", "acc")) |>
    dplyr::filter(!is.na(Pr), !is.na(As))

  tidyr::expand_grid(species = sort(unique(d$species)),
                     variable = c("As", "Pr")) |>
    purrr::pmap(\(species, variable) {
      sp <- dplyr::filter(d, species == !!species)
      f  <- sp[[variable]][sp$status == "found"]
      m  <- sp[[variable]][sp$status == "missing"]
      tt <- if (length(f) >= 2L && length(m) >= 2L) {
        stats::t.test(f, m)
      } else NULL
      tibble::tibble(
        species = species, variable = variable,
        n_found = length(f), n_missing = length(m),
        mean_found = if (length(f)) mean(f) else NA_real_,
        mean_missing = if (length(m)) mean(m) else NA_real_,
        diff = if (length(f) && length(m)) mean(f) - mean(m) else NA_real_,
        t = if (is.null(tt)) NA_real_ else unname(tt$statistic),
        df = if (is.null(tt)) NA_real_ else unname(tt$parameter),
        p_value = if (is.null(tt)) NA_real_ else tt$p.value)
    }) |>
    purrr::list_rbind()
}

#' Near-equivalent replacements for each slot that still needs a bag.
#'
#' WHY ALTERNATES AT ALL.  The cost this protocol minimises is TRIPS to the
#' seed store, and a slot whose single named accession turns out to be missing
#' costs a whole extra round.  Two or three ranked choices per slot means one
#' trip usually resolves the slot, and the ranking says what the substitution
#' costs: `score_cost` is the DROP in the noncentrality t, reported positive,
#' which on this vintage runs a few hundredths out of ~5 -- within or near the
#' parsimony slack, and far below the error on lambda.
#'
#' THE ALTERNATES ARE DISJOINT ACROSS SLOTS, and that is not a detail.  Scored
#' slot by slot they are not: the best replacement for a pool member barely
#' depends on WHICH member it replaces, so the naive version hands the same
#' three names to every open slot -- 17 slots each, on the first live run --
#' and a searcher who finds one of them can still only use it once.  So the
#' lists are assigned greedily, best (slot, candidate) pair first, with each
#' candidate offered to exactly one slot.  Every slot's list is then genuinely
#' its own, and a searcher can work down the sheet without collisions.
#'
#' A consequence worth stating: when the reservoir cannot supply
#' `n_alternates` for every slot -- round 1 asks for 2n slots at once, and pea
#' has 134 reservoir candidates for 60 slots -- the later slots get shorter
#' lists rather than borrowed ones.
#'
#' Only slots held by a NON-ANCHOR are listed: a slot already filled from
#' stores needs no search.  Each evaluation is the same O(1) sums update the
#' climb uses, so the whole table costs one sweep of the move set.
#'
#' A reservoir accession that is itself an anchor can appear as an alternate,
#' flagged `in_hand`.  That is not a contradiction with the parsimony stage:
#' the prune only takes such a swap when it costs less than the slack, and
#' here it is offered with its price attached.
pool_slot_alternates <- function(res, ..., n_alternates = 3L, tol = Inf,
                                 m_min = NULL) {
  if (...length() > 0L) {
    stop("pool_slot_alternates() takes only `res` positionally; name the rest",
         call. = FALSE)
  }
  ctx   <- res$ctx
  m_min <- m_min %||% res$m_min %||% 0L
  sp    <- .pool_sums(ctx, res$plus)
  sm    <- .pool_sums(ctx, res$minus)
  cur   <- .score_from_sums(sp, sm, ctx)
  reservoir <- setdiff(seq_len(ctx$N), c(res$plus, res$minus))

  empty <- tibble::tibble(pool = character(), slot = integer(),
                          slot_acc = character(), rank_in_slot = integer(),
                          acc = character(), is_first_choice = logical(),
                          in_hand = logical(), score = numeric(),
                          score_cost = numeric())

  # The open slots, in a fixed order: As+ then As-, by position.
  open <- list()
  for (side in c("plus", "minus")) {
    members <- res[[side]]
    for (j in seq_along(members)) {
      if (ctx$is_anchor[members[j]]) next
      open[[length(open) + 1L]] <- list(
        side = side, j = j, idx = members[j],
        pool = if (side == "plus") "As+" else "As-")
    }
  }
  if (length(open) == 0L) return(empty)

  firsts <- purrr::imap(open, \(o, i) tibble::tibble(
    pool = o$pool, slot = i, slot_acc = ctx$acc[o$idx], rank_in_slot = 1L,
    acc = ctx$acc[o$idx], is_first_choice = TRUE, in_hand = FALSE,
    score = cur$score, score_cost = 0)) |> purrr::list_rbind()

  if (length(reservoir) == 0L || n_alternates <= 0L) return(firsts)

  # One row per open slot, one column per reservoir candidate. -Inf marks a
  # swap that would break the anchor count or the Pr tolerance.
  S <- vapply(open, \(o) {
    if (o$side == "plus") {
      cp <- .sums_bump(sp, o$idx, reservoir, ctx); cm <- as.list(sm)
    } else {
      cp <- as.list(sp); cm <- .sums_bump(sm, o$idx, reservoir, ctx)
    }
    sc <- .score_from_sums(cp, cm, ctx)
    ok <- cp$a >= m_min & cm$a >= m_min & abs(sc$dPr) <= tol
    ifelse(ok, sc$score, -Inf)
  }, numeric(length(reservoir))) |>
    matrix(nrow = length(reservoir), ncol = length(open))

  # Greedy assignment, best pair first. Deterministic: ctx is sorted by acc,
  # so order()'s tie-break is "the lexicographically first accession".
  room  <- rep(as.integer(n_alternates), length(open))
  taken <- logical(length(reservoir))
  picks <- vector("list", sum(room))
  np <- 0L
  for (e in order(S, decreasing = TRUE)) {
    if (!is.finite(S[[e]])) break
    if (all(room == 0L)) break
    k <- ((e - 1L) %% nrow(S)) + 1L
    i <- ((e - 1L) %/% nrow(S)) + 1L
    if (room[i] == 0L || taken[k]) next
    room[i]  <- room[i] - 1L
    taken[k] <- TRUE
    np <- np + 1L
    picks[[np]] <- c(slot = i, cand = k, score = S[[e]])
  }
  if (np == 0L) return(firsts)

  alts <- purrr::list_rbind(purrr::map(picks[seq_len(np)], \(x)
    tibble::tibble(slot = as.integer(x[["slot"]]),
                   cand = as.integer(x[["cand"]]),
                   score = x[["score"]]))) |>
    dplyr::arrange(slot, dplyr::desc(score)) |>
    dplyr::group_by(slot) |>
    dplyr::mutate(rank_in_slot = dplyr::row_number() + 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      pool = purrr::map_chr(slot, \(i) open[[i]]$pool),
      slot_acc = purrr::map_chr(slot, \(i) ctx$acc[open[[i]]$idx]),
      acc = ctx$acc[reservoir[cand]],
      is_first_choice = FALSE,
      in_hand = ctx$is_anchor[reservoir[cand]],
      score_cost = cur$score - score) |>
    dplyr::select(dplyr::all_of(names(empty)))

  dplyr::bind_rows(firsts, alts) |>
    dplyr::arrange(slot, rank_in_slot)
}

#' The deliverable: one row per bag to go and find.
#'
#' Built from the adopted proposal, so it cannot name an accession the design
#' does not use.  Carries what a person standing at a freezer needs: the name,
#' which pool it is for, the alternates that would do instead, and -- for a
#' collapsed analysis name, which has no bag of its own -- the original
#' entries any one of which satisfies it.
#'
#' @param ledger Used only for `previous_pool` and `previous_status`, so the
#'   request can say "you were asked for this in round 1 and marked it
#'   missing_retry".
search_request <- function(res, inp, ..., round, ledger = NULL,
                           collapsed = NULL, n_alternates = 3L, tol = Inf) {
  if (...length() > 0L) {
    stop("search_request() takes res and inp positionally; name the rest",
         call. = FALSE)
  }
  # Captured before the pipeline: `round` becomes a COLUMN below, and inside
  # dplyr::mutate() the column would win over the argument.
  rnd <- as.integer(round)
  alt <- pool_slot_alternates(res, n_alternates = n_alternates, tol = tol)
  if (nrow(alt) == 0L) {
    return(tibble::tibble(species = character(), round = integer(),
                          pool = character(), slot = integer(),
                          slot_acc = character(), rank_in_slot = integer(),
                          acc = character(), is_first_choice = logical(),
                          in_hand = logical(), Pr = numeric(), As = numeric(),
                          rel_As_i = numeric(), n_partners = integer(),
                          score = numeric(), score_cost = numeric(),
                          n_entries = integer(), sowable_entries = character(),
                          previous_pool = character(),
                          previous_status = character(),
                          n_slots_offering = integer(), note = character()))
  }

  ent <- collapsed_entries(unique(alt$acc), inp$species,
                           collapsed = collapsed) |>
    dplyr::select(acc, n_entries, sowable_entries = entries, is_collapsed)
  prev <- (ledger %||% .empty_seed_ledger()) |>
    dplyr::filter(species == inp$species) |>
    dplyr::select(acc, previous_pool = pool_proposed,
                  previous_status = status, previous_round = round_proposed)

  alt |>
    dplyr::left_join(
      dplyr::select(inp$accessions, acc, Pr, As, rel_As_i, n_partners),
      by = "acc") |>
    dplyr::left_join(ent, by = "acc") |>
    dplyr::left_join(prev, by = "acc") |>
    dplyr::mutate(species = inp$species, round = rnd, .before = 1) |>
    dplyr::group_by(acc) |>
    dplyr::mutate(n_slots_offering = dplyr::n_distinct(slot)) |>
    dplyr::ungroup() |>
    dplyr::mutate(note = purrr::pmap_chr(
      list(in_hand, is_collapsed, n_entries, previous_status, previous_round),
      \(ih, ic, ne, pstat, pround) {
        bits <- c(
          if (isTRUE(ih)) "already found in storage; no search needed",
          if (isTRUE(ic)) paste0("collapsed name: any one of its ", ne,
                                 " entries will do"),
          if (!is.na(pstat) && pstat == "missing_retry")
            paste0("marked missing_retry in round ", pround),
          # Only worth saying once it is OLD news: on the first proposal every
          # row is outstanding since round 1, which is noise.
          if (!is.na(pstat) && pstat == "unsearched" && !is.na(pround) &&
              pround < rnd)
            paste0("outstanding since round ", pround))
        if (length(bits) == 0L) NA_character_ else
          paste(bits, collapse = "; ")
      })) |>
    (\(x) { stopifnot("an accession is offered to more than one slot" =
                       all(x$n_slots_offering == 1L)); x })() |>
    dplyr::arrange(pool, slot, rank_in_slot) |>
    dplyr::select(species, round, pool, slot, slot_acc, rank_in_slot, acc,
                  is_first_choice, in_hand, Pr, As, rel_As_i, n_partners,
                  score, score_cost, n_entries, sowable_entries,
                  previous_pool, previous_status, n_slots_offering, note)
}

#' The same request as a checklist to carry to the seed store.
#'
#' The CSV is the record; this is the working copy, which is why it is a list
#' of tick-boxes grouped by slot rather than a table.  Dropping it would cost
#' nothing analytically -- said here so that if this file ever needs trimming,
#' this is the first thing to go.
search_request_md <- function(req, ..., vintage, ledger_path, round = NULL) {
  if (...length() > 0L) {
    stop("search_request_md() takes only `req` positionally; name the rest",
         call. = FALSE)
  }
  firsts <- dplyr::filter(req, is_first_choice)
  per_sp <- firsts |> dplyr::count(species, name = "n")
  # `round` is passed rather than read off the table because an EMPTY request
  # -- the good outcome, everything already in hand -- has no rows to read it
  # from, and "round -" in the heading of the one document somebody files is
  # a poor way to record a success.
  rnd <- round %||% (if (nrow(req) == 0L) NA_integer_ else req$round[1])
  if (nrow(req) == 0L) {
    return(c(paste0("# Seed search request -- round ", rnd), "",
             paste0("Vintage ", vintage,
                    ". Nothing to find: every pool member is already in ",
                    "hand."), "",
             paste0("The pools in pools.csv are sowable as they stand.")))
  }
  head <- c(
    paste0("# Seed search request -- round ", rnd),
    "",
    paste0("Vintage ", vintage, ". ", nrow(firsts), " slot(s) to fill",
           if (nrow(per_sp) > 0L)
             paste0(" (", paste(per_sp$species, per_sp$n, collapse = ", "),
                    ")")
           else "", "."),
    "",
    paste0("Find ONE accession per slot: the first choice, or the first ",
           "alternate you can lay hands on."),
    paste0("`cost` is the drop in the contrast's noncentrality t if you take ",
           "that alternate instead of the first choice -- a few thousandths ",
           "is far below the uncertainty in lambda, so a near tie really is ",
           "a tie."),
    "",
    paste0("Record what happened in `", ledger_path,
           "`: set `status` to found, missing or missing_retry, fill in ",
           "`searched_on`, and for a collapsed name put the entry you ",
           "actually have in `sowable_entry`."),
    "")

  body <- req |>
    dplyr::group_split(species, pool) |>
    purrr::map(\(g) {
      slots <- purrr::map(dplyr::group_split(g, slot), \(sl) c(
        paste0("**Slot ", sl$slot[1], "**",
               if (!is.na(sl$note[1])) paste0("  _", sl$note[1], "_") else ""),
        purrr::pmap_chr(
          list(sl$acc, sl$Pr, sl$As, sl$rank_in_slot, sl$score_cost,
               sl$n_entries, sl$sowable_entries, sl$in_hand),
          \(acc, Pr, As, rk, cost, ne, ents, ih) paste0(
            "- [ ] ", if (rk == 1L) "**" else "", acc,
            if (rk == 1L) "**" else "",
            sprintf("  (Pr %.1f, As %.1f", Pr, As),
            if (rk == 1L) ", first choice" else sprintf(", cost %.4f", cost),
            if (isTRUE(ih)) ", IN HAND" else "", ")",
            if (!is.na(ne) && ne > 1L)
              paste0("\n      any one of: ", ents) else "")),
        ""))
      c(paste0("## ", g$species[1], " -- ", g$pool[1]), "",
        unlist(slots, use.names = FALSE))
    }) |>
    unlist(use.names = FALSE)
  c(head, body)
}

#' Run at every refresh, and REPORTED rather than judged: whoever orders seed
#' needs to know what moved. Churn is deliberately not a criterion -- see
#' VALIDATION_DESIGN.md section 7. The year-3 trial is a validation trial
#' regardless, so a reshuffled pool is not a reason to postpone it, and while
#' the method is still being settled churn between vintages largely measures
#' changes to the METHOD (compare the pool_index stamp before reading a diff
#' as instability).
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
  take <- function(x, m, offset) {
    if (m <= 0L) return(x[0])
    x[((offset + seq_len(m) - 1L) %% length(x)) + 1L]
  }

  # The anchors must sit OUTSIDE that rotation. If they are taken from the
  # rotated window, the first member of the window differs at every location
  # whenever the cell is smaller than the pool -- which is the normal case,
  # because the anchors' second copies come out of the budget -- and then the
  # anchor combination is not repeated across locations at all, so combination
  # x location is not estimable. So: the first `anchor_n` members of each pool
  # are fixed, and only the remaining slots rotate.
  anchor_n <- min(as.integer(n_anchor), cell_size)
  pick <- function(x, m, offset) {
    if (length(x) <= anchor_n) {
      stop("a pool of ", length(x), " cannot supply ", anchor_n,
           " anchor(s) and still rotate; use a larger pool or fewer anchors",
           call. = FALSE)
    }
    c(x[seq_len(anchor_n)], take(x[-seq_len(anchor_n)], m - anchor_n, offset))
  }

  cells <- tidyr::expand_grid(oat_pool = c("As+", "As-"),
                              pea_pool = c("As+", "As-"))

  design <- purrr::map(seq_len(n_loc), \(loc) {
    purrr::pmap(cells, \(oat_pool, pea_pool) {
      o_all <- oats[[oat_pool]]; p_all <- peas[[pea_pool]]
      rot <- (loc - 1L) * (cell_size - anchor_n)
      o <- pick(o_all, cell_size, rot)
      p <- pick(p_all, cell_size, rot)

      # anchors first, identical at every location; the rest re-permuted
      anchor_idx <- seq_len(anchor_n)
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
# The pools are fixed sets but the accessions within them are a random sample of
# "accessions predicted to be in this class", so the accession is the
# experimental unit and enters the denominator.
#
# THREE STRATA, THREE DEGREES OF FREEDOM. The contrast's variance has an
# accession term (2n - 2 df, Welch-adjusted when the two pools differ), a
# pool x location term (n_loc - 1 df) and a plot term (effectively P df). The
# pool x location term is the LARGEST of the three at the measured interaction,
# and it has four degrees of freedom. Using 2n - 2 for the whole statistic, as
# this did until 2026-10-03, overstates power by 1-5 points -- more when the
# interaction's share is larger. The effective df is Satterthwaite's.
#
# ARGUMENTS AFTER THE FIRST ARE NAMED-ONLY, enforced by placing `...` ahead of
# them: R will not match an argument after `...` positionally or partially.
# This is deliberate. tests/test_validation.R had nine fully positional calls
# relying on the order (dAs, sigma2_within, n, sigma2_e, P), so inserting any
# argument before `P` would have silently rebound all nine while the suite kept
# reporting PASS. The `...length()` guard also catches a MISSPELLED name, which
# would otherwise vanish into `...` and be ignored -- the quietest failure of
# all.
# ------------------------------------------------------------

# ------------------------------------------------------------
# THE Pr-ADJUSTED CONTRAST
#
# The hypothesis the trial powers is "for accessions of equal Pr, there is
# variation in As on which we can select".  Strictly that is the hypothesis
# that the genetic correlation between As and Pr differs from -1; in practice
# it is the hypothesis that the pools diverge in As with Pr carried as a
# COVARIATE.  So a pool pair with a larger raw dAs can be WORSE than a smaller
# one if it also diverges in Pr in the direction the negative As-Pr
# correlation already predicts: part of its As gap is then explained by Pr
# rather than being evidence of selectable variation at fixed Pr.
#
# Accession-level ANCOVA on the 2n selected accessions gives three changes:
#
#   dAs_adj = dAs - b * dPr                      the adjusted contrast
#   f       = 1 + (n/2) * dPr^2 / Sxx            variance inflation, applied
#                                                  to v_acc AND v_plt
#   df_acc  = df_acc * (2n - 3) / (2n - 2)       one df to the slope
#
# with Sxx = (n-1) * (var_Pr_plus + var_Pr_minus), the within-pool sum of
# squares of Pr pooled over both pools.
#
# `b` IS THE CANDIDATE-SET SLOPE, not the within-pool one.  Membership is
# decided by `As + theta * Pr >= cut`, which truncates the joint distribution
# and induces a spurious negative As-Pr relation INSIDE each pool; regressing
# there would adjust using an artefact of the selection.  pr_adjust_slope()
# is the single definition, and it takes the post-threshold candidate set.
#
# WHAT IS DELIBERATELY NOT DONE.  Conditioning on Pr also reduces the
# accession-level residual to sigma2_within - b^2 * var_Pr_within.  That
# reduction is NOT applied.  It would raise power, and the same conservative
# choice is already made for `lambda` (see its note under contrast_power():
# lambda scales the contrast but leaves sigma2_within at full size).  A power
# claim should not be the optimistic reading of either.
# ------------------------------------------------------------

#' Validate a `pr_adjust` list and return it with the inflation factor.
#'
#' All three elements or none: a half-supplied adjustment is a silent wrong
#' answer, so it is one validation rather than three independent defaults.
.check_pr_adjust <- function(pr_adjust, n) {
  if (!is.list(pr_adjust)) {
    stop("pr_adjust must be a list of dPr, b and Sxx, or NULL", call. = FALSE)
  }
  need <- c("dPr", "b", "Sxx")
  miss <- setdiff(need, names(pr_adjust))
  if (length(miss) > 0L) {
    stop("pr_adjust is missing: ", paste(miss, collapse = ", "),
         ". All of dPr, b and Sxx are required, or pass NULL for the ",
         "unadjusted contrast.", call. = FALSE)
  }
  for (k in need) {
    v <- pr_adjust[[k]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop("pr_adjust$", k, " must be a single finite number", call. = FALSE)
    }
  }
  if (pr_adjust$Sxx <= 0) {
    stop("pr_adjust$Sxx must be positive: it is the within-pool sum of ",
         "squares of Pr, and a non-positive value means the pools carry no ",
         "Pr variation to adjust along.", call. = FALSE)
  }
  pr_adjust$f <- 1 + (n / 2) * pr_adjust$dPr^2 / pr_adjust$Sxx
  pr_adjust
}

#' Standard error of the pool contrast.
#'
#' @param sigma2_within PEV + within-pool BLUP variance (true-effect spread
#'   among the accessions in a pool). Length 1 for both pools, or length 2 as
#'   `c(plus, minus)`. The SE is identical either way -- `(A + B)/n` is
#'   `2*mean/n` -- but the Welch df in `contrast_power()` is not.
#' @param n Accessions per pool.
#' @param sigma2_e Plot residual variance of the PARTNER's yield.
#' @param P Total plots across all locations.
#' @param interaction_frac SD of the location-specific contrast, as a FRACTION
#'   of the contrast itself. This is the form the cross-validation measures:
#'   `sd(lambda)/|mean(lambda)|` across folds.
#' @param interaction_sd SD of the location-specific contrast in ABSOLUTE
#'   units, i.e. g/m2. Supply this instead of `interaction_frac` when `delta`
#'   is NOT the contrast the spread was measured around -- the across-location
#'   spread is a property of the environments, so it must not be rescaled by a
#'   lambda chosen afterwards. Exactly one of the two may be given.
#' @param delta The contrast, needed only to scale `interaction_frac`.
#' @param n_loc Locations.
#' @param df_e_lost Model df taken out of the plot stratum. Power is
#'   insensitive to it -- the plot term's contribution to the Satterthwaite
#'   denominator is under 1% -- so it is an argument with a default rather than
#'   a question worth time.
#' @param pr_adjust NULL (the default) for the unadjusted contrast, or a list
#'   of `dPr`, `b` and `Sxx` for the Pr-ADJUSTED one -- the pool contrast with
#'   producer effect carried as a covariate. See the algebra note above
#'   `.check_pr_adjust()`. NULL leaves every expression below bit-identical to
#'   the unadjusted form, which is the backwards-compatibility contract.
#' @param .components TRUE returns the variance split and the df alongside the
#'   SE, so `contrast_power()` does not recompute them. Duplicating the split
#'   in two functions is how an SE and its df drift apart.
contrast_se <- function(sigma2_within, ..., n, sigma2_e, P,
                        interaction_frac = NULL, interaction_sd = NULL,
                        delta = 0, n_loc = 5L, df_e_lost = 4L,
                        pr_adjust = NULL,
                        .components = FALSE) {
  if (...length() > 0L) {
    stop("contrast_se() takes only `sigma2_within` positionally; ",
         ...length(), " extra positional argument(s) given. Name every ",
         "argument: contrast_se(sigma2_within = , n = , sigma2_e = , P = , ...)",
         call. = FALSE)
  }
  if (!is.null(interaction_frac) && !is.null(interaction_sd)) {
    stop("supply interaction_frac OR interaction_sd, not both: they are two ",
         "parameterisations of the same term and would be double-counted",
         call. = FALSE)
  }
  if (!length(sigma2_within) %in% c(1L, 2L)) {
    stop("sigma2_within must be length 1 (both pools) or 2 (plus, minus)",
         call. = FALSE)
  }

  int_sd <- if (!is.null(interaction_sd)) interaction_sd
            else if (!is.null(interaction_frac)) interaction_frac * delta
            else 0

  # The ANCOVA inflation, as a FACTOR on both accession-level strata rather
  # than an additive term on v_acc alone. Both are of the form
  # (2/n) * (a per-accession variance) -- note
  # v_plt = 4*sigma2_e/P = (2/n) * (2n*sigma2_e/P) -- so the plot noise
  # propagates through the estimated slope exactly as the accession noise
  # does. Writing it as a factor also keeps the variance split at THREE
  # components, so pct_var_accessions + pct_var_interaction + pct_var_plots
  # still sums to 100; an additive fourth term would break that.
  f <- if (is.null(pr_adjust)) 1 else .check_pr_adjust(pr_adjust, n)$f

  v_acc <- sum(sigma2_within) / n * (if (length(sigma2_within) == 1L) 2 else 1)
  v_int <- int_sd^2 / n_loc
  v_plt <- 4 * sigma2_e / P
  # Applied only when adjusting, so the unadjusted path performs no extra
  # floating-point operation at all and stays bit-identical.
  if (!is.null(pr_adjust)) {
    v_acc <- f * v_acc
    v_plt <- f * v_plt
  }
  v     <- v_acc + v_int + v_plt

  if (!.components) return(sqrt(v))

  # Welch for the accession stratum: with unequal pool variances the two halves
  # of the contrast carry different weight, and this collapses to exactly
  # 2n - 2 when they are equal.
  a <- (if (length(sigma2_within) == 1L) sigma2_within else sigma2_within[1]) / n
  b <- (if (length(sigma2_within) == 1L) sigma2_within else sigma2_within[2]) / n
  df_acc <- (a + b)^2 / (a^2 / (n - 1) + b^2 / (n - 1))

  # Estimating the Pr slope costs one df. Applied as a MULTIPLIER on the Welch
  # df_acc rather than by substituting 2n - 3, so unequal pool variances still
  # propagate; at equal variances it returns exactly 2n - 3.
  if (!is.null(pr_adjust)) df_acc <- df_acc * (2 * n - 3) / (2 * n - 2)

  df <- if (v_int > 0) {
    v^2 / (v_acc^2 / df_acc + v_int^2 / max(n_loc - 1, 1) +
           v_plt^2 / max(P - df_e_lost, 1))
  } else {
    # With no interaction term the statistic is the ordinary two-sample
    # contrast and its df is 2n - 2 EXACTLY, which is the boundary case
    # tests/test_validation.R pins. Routing it through Satterthwaite would
    # return 2n - 2 only up to floating point, so it is returned directly.
    df_acc
  }

  list(se = sqrt(v), v_acc = v_acc, v_int = v_int, v_plt = v_plt,
       df_acc = df_acc, df = df, int_sd = int_sd,
       var_inflation = f,
       df_cost_covariate = if (is.null(pr_adjust)) 0L else 1L,
       interaction_mode = if (!is.null(interaction_sd)) "sd"
                          else if (!is.null(interaction_frac)) "frac" else "none")
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
#'
#'   A WARNING ABOUT COMBINING lambda WITH interaction_frac. interaction_frac
#'   is `sd(lambda)/|mean(lambda)|`, so `interaction_frac * delta` is the
#'   interaction SD only when `lambda` is the same lambda the spread was
#'   measured around. Pass `interaction_sd` instead whenever it is not.
#' @param sided 1 for the pre-registered directional test, 2 otherwise.
#' @param df_mode "satterthwaite" (default) or "naive" for the old 2n - 2,
#'   which is kept so a report can show the two side by side.
contrast_power <- function(dAs, ..., sigma2_within, n, sigma2_e, P,
                           lambda = 1, alpha = 0.05, sided = 1,
                           interaction_frac = NULL, interaction_sd = NULL,
                           n_loc = 5L, df_e_lost = 4L,
                           pr_adjust = NULL,
                           df_mode = c("satterthwaite", "naive")) {
  if (...length() > 0L) {
    stop("contrast_power() takes only `dAs` positionally; ", ...length(),
         " extra positional argument(s) given. Name every argument: ",
         "contrast_power(dAs = , sigma2_within = , n = , sigma2_e = , P = , ...)",
         call. = FALSE)
  }
  df_mode <- match.arg(df_mode)

  # The adjustment moves the CONTRAST as well as its variance, and the
  # interaction term scales off the contrast, so dAs_adj has to be formed
  # before contrast_se() is called rather than after.
  dAs_adj <- if (is.null(pr_adjust)) dAs else dAs - pr_adjust$b * pr_adjust$dPr

  delta <- lambda * dAs_adj
  cmp   <- contrast_se(sigma2_within = sigma2_within, n = n,
                       sigma2_e = sigma2_e, P = P,
                       interaction_frac = interaction_frac,
                       interaction_sd = interaction_sd,
                       delta = delta, n_loc = n_loc, df_e_lost = df_e_lost,
                       pr_adjust = pr_adjust,
                       .components = TRUE)
  se  <- cmp$se
  df  <- if (df_mode == "naive") 2 * n - 2 else cmp$df
  ncp <- delta / se

  power <- if (sided == 1) {
    stats::pt(stats::qt(1 - alpha, df), df, ncp, lower.tail = FALSE)
  } else {
    stats::pt(-stats::qt(1 - alpha / 2, df), df, ncp) +
      stats::pt(stats::qt(1 - alpha / 2, df), df, ncp, lower.tail = FALSE)
  }

  # SE_floor is the standard error with INFINITELY MANY PLOTS, which is what
  # the design document means by it and what power_ceiling.csv reports. Two of
  # the three terms are free of P, so the floor is sqrt(v_acc + v_int). It used
  # to be sqrt(v_acc) alone, which was the P -> infinity limit only because the
  # interaction term defaulted to zero. The old quantity survives as
  # SE_accession_only, because "more accessions lower the floor" is true of the
  # accession term and false of the interaction term -- and that distinction is
  # the point: at the measured interaction ~90% of the SE is irreducible by
  # plots, and most of the irreducible part is pool x location, so the one
  # remaining lever is MORE LOCATIONS.
  se_floor <- sqrt(cmp$v_acc + cmp$v_int)

  tibble::tibble(
    n = n, P = P, lambda = lambda, sided = sided,
    interaction_frac = interaction_frac %||% NA_real_,
    interaction_sd = cmp$int_sd,
    interaction_mode = cmp$interaction_mode,
    delta = delta, SE = se,
    v_acc = cmp$v_acc, v_int = cmp$v_int, v_plt = cmp$v_plt,
    pct_var_accessions  = 100 * cmp$v_acc / se^2,
    pct_var_interaction = 100 * cmp$v_int / se^2,
    pct_var_plots       = 100 * cmp$v_plt / se^2,
    df = df, df_acc = cmp$df_acc, df_naive = 2 * n - 2, df_mode = df_mode,
    t = ncp, power = power,
    # Appended, never inserted: every consumer selects by name, so growing the
    # tibble is safe while reordering it would not be.
    dAs_adj = dAs_adj,
    b_Pr   = if (is.null(pr_adjust)) NA_real_ else pr_adjust$b,
    Sxx_Pr = if (is.null(pr_adjust)) NA_real_ else pr_adjust$Sxx,
    var_inflation = cmp$var_inflation,
    df_cost_covariate = cmp$df_cost_covariate,
    SE_floor = se_floor,
    SE_accession_only = sqrt(cmp$v_acc),
    # Kept under its old name: the share of the SE the plot budget cannot buy.
    pct_var_from_accessions = 100 * cmp$v_acc / se^2,
    pct_SE_irreducible = 100 * se_floor / se
  )
}

# ------------------------------------------------------------
# The continuous estimands
#
# The pool contrast is a two-point summary of a regression. The regression
# itself is the better statement of the same evidence: its slope IS lambda
# measured in new data, on the same scale as the cross-validation's
# interaction_frac, so the trial measures the quantity the whole design is
# conditioned on. It also does not depend on where the pool boundary fell.
#
# WHAT IT DOES NOT BUY IS POWER. Pools are the extremes of the predicted
# distribution, so selection has already removed most of the within-pool spread
# in the predictor that a regression would otherwise exploit: Var(x) over the 2n
# selected accessions is (dAs/2)^2 + var_As_within_pools, and the second term is
# a few percent of the first. The slope is adopted for interpretability, and
# tests/test_validation.R pins the near-equivalence so this cannot quietly be
# oversold.
#
# THE VARIANCE IS COMPUTED FROM THE REALISED DESIGN, not from a closed form that
# assumes balance. For a slope b = Sxy/Sxx with plot-level errors and
# accession-level conditional errors shared across an accession's plots,
#
#   Var(b) = [ sum_m PEV_m * S_m^2  +  sigma2_e * Sxx ] / Sxx^2
#
# where S_m is the sum of the centred predictor over accession m's plots. That
# is an exact sandwich for the design `make_validation_design()` produced,
# including the anchors' unequal replication, which no balanced formula covers.
# ------------------------------------------------------------

#' Power for the slope of realised on predicted associate effect.
#'
#' @param design A field book from `make_validation_design()`.
#' @param x Named numeric: the predicted associate effect per FOCAL accession.
#' @param pev_focal Named numeric: per-accession PEV of the focal effect, i.e.
#'   the conditional variance of the truth given the prediction.
#' @param pev_partner,partner_x The same for the partner species, whose own
#'   producer effect is in the response. Supply both or neither.
#' @param sigma2_e Plot residual variance of the response.
#' @param lambda The true slope under the alternative. The null is zero.
#' @param interaction_frac SD of the location-specific slope as a fraction of
#'   the slope itself -- the same quantity the cross-validation measures.
#' @param focal_col,partner_col Which columns of `design` carry the two.
slope_power <- function(design, ..., x, pev_focal, sigma2_e, lambda,
                        interaction_frac = 0, n_loc = NULL,
                        focal_col = "oat_acc", partner_col = "pea_acc",
                        pev_partner = NULL, partner_x = NULL,
                        alpha = 0.05, sided = 1L) {
  if (...length() > 0L) {
    stop("slope_power() takes only `design` positionally; name the rest",
         call. = FALSE)
  }
  n_loc <- n_loc %||% dplyr::n_distinct(design$location)

  f <- as.character(design[[focal_col]])
  xi <- unname(x[f])
  if (anyNA(xi)) {
    stop("every focal accession in the design needs a predicted effect; ",
         sum(is.na(xi)), " plot(s) have none", call. = FALSE)
  }
  # The partner's own producer effect sits in the response. Including it in the
  # predictor is what the pre-registered model does (`x_F + x_Pa`), so it is
  # orthogonalised rather than left to inflate the residual.
  if (!is.null(partner_x)) {
    pj <- unname(partner_x[as.character(design[[partner_col]])])
    xi <- xi - stats::lm(xi ~ pj)$fitted.values + mean(xi)
  }

  xc  <- xi - mean(xi)
  Sxx <- sum(xc^2)

  # accession-level conditional errors, shared across an accession's plots
  S_focal <- tapply(xc, f, sum)
  v_focal <- sum(unname(pev_focal[names(S_focal)]) * S_focal^2)
  v_partner <- if (!is.null(pev_partner)) {
    g <- as.character(design[[partner_col]])
    S_p <- tapply(xc, g, sum)
    sum(unname(pev_partner[names(S_p)]) * S_p^2)
  } else 0

  v_slope_acc <- (v_focal + v_partner) / Sxx^2
  v_slope_plt <- sigma2_e / Sxx
  # The slope varies by location: that IS the As x environment interaction, and
  # on this scale its SD is interaction_frac * lambda, by the definition of
  # interaction_frac as sd(lambda)/|mean(lambda)|.
  v_slope_int <- (interaction_frac * lambda)^2 / n_loc

  v  <- v_slope_acc + v_slope_plt + v_slope_int
  se <- sqrt(v)

  n_acc <- dplyr::n_distinct(f)
  df <- if (v_slope_int > 0) {
    v^2 / (v_slope_acc^2 / max(n_acc - 2, 1) +
           v_slope_int^2 / max(n_loc - 1, 1) +
           v_slope_plt^2 / max(nrow(design) - n_acc, 1))
  } else max(n_acc - 2, 1)

  ncp   <- lambda / se
  power <- if (sided == 1) {
    stats::pt(stats::qt(1 - alpha, df), df, ncp, lower.tail = FALSE)
  } else {
    stats::pt(-stats::qt(1 - alpha / 2, df), df, ncp) +
      stats::pt(stats::qt(1 - alpha / 2, df), df, ncp, lower.tail = FALSE)
  }

  tibble::tibble(
    n_plots = nrow(design), n_accessions = n_acc, n_loc = n_loc,
    lambda = lambda, var_x = Sxx / (nrow(design) - 1), Sxx = Sxx,
    SE = se, v_acc = v_slope_acc, v_plt = v_slope_plt, v_int = v_slope_int,
    pct_var_interaction = 100 * v_slope_int / v,
    df = df, t = ncp, power = power,
    interaction_frac = interaction_frac
  )
}

#' Power for the total-yield slope: observed total against predicted total.
#'
#' The co-primary. `GMA = Pr + As` is exactly the per-accession contribution to
#' the plot total, so the predicted total for a combination is
#' `GMA_oat + w_pea * GMA_pea` and no new quantity has to be estimated.
#'
#' TWO THINGS ABOUT IT THAT ARE EASY TO MISREAD.
#'
#' Because the pools are matched on Pr, `dGMA` is almost exactly `dAs` -- so
#' this estimand brings a different RESPONSE, not a different predictor. Its
#' value is that total productivity is the breeding objective, and that it costs
#' nothing: both yields are already measured on every plot.
#'
#' And the 2x2 factorial makes the two species' contrasts ADD in the predicted
#' total, so the predictor spans roughly twice what either species alone does.
#' That extra spread is what offsets the larger plot variance of a sum.
#'
#' @param w_pea Weight on pea yield. 1 is the physical total and is the
#'   pre-registered co-primary; an economic weighting is a secondary analysis,
#'   kept separate so a physical claim is not read as an economic one.
total_yield_power <- function(design, ..., gma_oat, gma_pea,
                              pev_gma_oat, pev_gma_pea,
                              var_oat, var_pea, cov_oat_pea,
                              lambda, interaction_frac = 0, w_pea = 1,
                              n_loc = NULL, alpha = 0.05, sided = 1L) {
  if (...length() > 0L) {
    stop("total_yield_power() takes only `design` positionally; name the rest",
         call. = FALSE)
  }
  n_loc <- n_loc %||% dplyr::n_distinct(design$location)

  o <- as.character(design$oat_acc)
  p <- as.character(design$pea_acc)
  xi <- unname(gma_oat[o]) + w_pea * unname(gma_pea[p])
  if (anyNA(xi)) stop("a plot has no predicted total", call. = FALSE)

  xc  <- xi - mean(xi)
  Sxx <- sum(xc^2)

  S_o <- tapply(xc, o, sum); S_p <- tapply(xc, p, sum)
  v_acc <- (sum(unname(pev_gma_oat[names(S_o)]) * S_o^2) +
            w_pea^2 * sum(unname(pev_gma_pea[names(S_p)]) * S_p^2)) / Sxx^2

  # Residual variance of the SUM carries the within-plot covariance, which is
  # negative here: oat and pea compete, so the total is less variable than the
  # two yields separately would suggest.
  sigma2_total <- var_oat + w_pea^2 * var_pea + 2 * w_pea * cov_oat_pea
  v_plt <- sigma2_total / Sxx
  v_int <- (interaction_frac * lambda)^2 / n_loc

  v  <- v_acc + v_plt + v_int
  se <- sqrt(v)
  n_acc <- dplyr::n_distinct(o) + dplyr::n_distinct(p)
  df <- if (v_int > 0) {
    v^2 / (v_acc^2 / max(n_acc - 2, 1) + v_int^2 / max(n_loc - 1, 1) +
           v_plt^2 / max(nrow(design) - n_acc, 1))
  } else max(n_acc - 2, 1)

  ncp   <- lambda / se
  power <- if (sided == 1) {
    stats::pt(stats::qt(1 - alpha, df), df, ncp, lower.tail = FALSE)
  } else {
    stats::pt(-stats::qt(1 - alpha / 2, df), df, ncp) +
      stats::pt(stats::qt(1 - alpha / 2, df), df, ncp, lower.tail = FALSE)
  }

  tibble::tibble(
    n_plots = nrow(design), n_accessions = n_acc, n_loc = n_loc,
    w_pea = w_pea, lambda = lambda, sigma2_total = sigma2_total,
    var_x = Sxx / (nrow(design) - 1), Sxx = Sxx,
    SE = se, v_acc = v_acc, v_plt = v_plt, v_int = v_int,
    pct_var_interaction = 100 * v_int / v,
    df = df, t = ncp, power = power, interaction_frac = interaction_frac
  )
}

#' Power over a grid of pool sizes, plot budgets and locations.
#'
#' Pools are rebuilt at every n, because the contrast shrinks as the pool grows
#' -- that trade-off is the whole question and must not be held fixed.
#'
#' `n_loc_values` is swept because it is the only lever that touches the
#' pool x location term, which is the largest of the three variance components
#' at the measured interaction. The plot budget, by contrast, is settled: it
#' only ever moved the smallest term.
#'
#' The per-pool `sigma2_within` is passed through as a length-2 vector, so the
#' Welch df reflects any difference between the two pools. The scalar is still
#' reported as a column, which is what pool_summary.csv and the generated
#' tables carry.
power_grid <- function(inp, n_values, P_values = NULL,
                       lambda_values = c(1, 0.8, 0.6),
                       sided = 1, interaction_values = 0, pr_quantile = 0.5,
                       n_loc = 5L, n_loc_values = NULL,
                       plots_per_loc_values = NULL,
                       interaction_mode = c("frac", "sd"),
                       pr_adjust = FALSE,
                       df_mode = "satterthwaite") {
  interaction_mode <- match.arg(interaction_mode)
  n_loc_values <- n_loc_values %||% n_loc

  # PLOTS PER LOCATION IS THE BETTER AXIS once locations are swept too.
  # Crossing a TOTAL plot budget with a location count gives cells that are not
  # comparable: 400 plots is 100 per site at four locations and 80 at five, so
  # the two location arms never share a per-site effort and a curve drawn
  # against total P zig-zags between two different designs. Sweeping plots per
  # location and deriving P = ppl * n_loc keeps the arms comparable and asks the
  # question the design actually faces -- another site, or more plots at each.
  #
  # P_values is kept for callers that genuinely mean a fixed total.
  if (is.null(plots_per_loc_values) && is.null(P_values)) {
    stop("power_grid() needs plots_per_loc_values or P_values", call. = FALSE)
  }

  budget <- if (!is.null(plots_per_loc_values)) {
    tidyr::expand_grid(plots_per_loc = plots_per_loc_values,
                       n_locations = n_loc_values) |>
      dplyr::mutate(P = plots_per_loc * n_locations)
  } else {
    tidyr::expand_grid(P = P_values, n_locations = n_loc_values) |>
      dplyr::mutate(plots_per_loc = P / n_locations)
  }

  purrr::map(n_values, \(n) {
    bp <- build_pools(inp, n, pr_quantile = pr_quantile)
    s2w <- c(bp$sigma2_within_plus, bp$sigma2_within_minus)
    # Default FALSE: power_grid.csv is a long-standing output and every pin on
    # it describes the unadjusted contrast. The adjustment is opt-in.
    pa <- if (isTRUE(pr_adjust)) {
      list(dPr = bp$dPr, b = bp$b_As_on_Pr, Sxx = bp$Sxx_w)
    } else NULL
    tidyr::expand_grid(budget, lambda = lambda_values,
                       interaction = interaction_values) |>
      purrr::pmap(\(P, n_locations, plots_per_loc, lambda, interaction) {
        args <- list(dAs = bp$dAs, sigma2_within = s2w, n = n,
                     sigma2_e = inp$sigma2_e, P = P, lambda = lambda,
                     sided = sided, n_loc = n_locations, pr_adjust = pa,
                     df_mode = df_mode)
        args[[if (interaction_mode == "sd") "interaction_sd"
              else "interaction_frac"]] <- interaction
        do.call(contrast_power, args) |>
          dplyr::mutate(n_locations = n_locations,
                        plots_per_loc = plots_per_loc, .after = P)
      }) |>
      purrr::list_rbind() |>
      dplyr::mutate(species = inp$species, dAs_predicted = bp$dAs,
                    dPr = bp$dPr, sigma2_within = bp$sigma2_within,
                    sigma2_within_plus = bp$sigma2_within_plus,
                    sigma2_within_minus = bp$sigma2_within_minus,
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
    eligibility_rule = i$eligibility$rule, rel_min = i$eligibility$rel_min,
    sigma2_As = i$sigma2_As, reliability_As = i$rel_As, PEV_As = i$PEV_As,
    # the measured posterior variance, and how far the backed-out one is from it
    PEV_As_measured = i$PEV_As_measured,
    PEV_ratio = i$PEV_ratio_backed_out_to_measured,
    median_rel_As_i = stats::median(i$accessions$rel_As_i),
    response = i$response, sigma2_e = i$sigma2_e
  )) |>
    purrr::list_rbind()
}
