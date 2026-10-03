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
POOL_INDEX_VERSION <- "single-v2"

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

  pools <- dplyr::bind_rows(
    dplyr::mutate(res$plus,  pool = "As+"),
    dplyr::mutate(res$minus, pool = "As-")
  ) |>
    dplyr::transmute(species = inp$species, pool, acc, Pr, As, GMA,
                     n_partners, n_plots, PEV_As_i, rel_As_i, PEV_GMA_i) |>
    dplyr::arrange(pool, dplyr::desc(As))

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
  pev_plus  <- mean(res$plus$PEV_As_i)
  pev_minus <- mean(res$minus$PEV_As_i)
  var_blup_plus  <- stats::var(res$plus$As)
  var_blup_minus <- stats::var(res$minus$As)
  var_blup_within <- mean(c(var_blup_plus, var_blup_minus))
  sigma2_within_plus  <- pev_plus  + var_blup_plus
  sigma2_within_minus <- pev_minus + var_blup_minus

  list(
    pools = pools, n = n, theta = theta, pr_min = pr_min,
    # Stamped into pool_summary.csv so a later run can tell whether the previous
    # vintage's pools are comparable with its own. Churn between vintages built
    # by different rules measures the RULE, not the stability of the effect
    # estimates, and the 70% retention guidance is about the latter.
    pool_index = POOL_INDEX_VERSION,
    n_candidates = nrow(cand), n_sign_changes = n_sign_changes,
    dAs = res$dAs, dPr = res$dPr,
    mean_Pr_plus = mean(res$plus$Pr), mean_Pr_minus = mean(res$minus$Pr),
    mean_As_plus = mean(res$plus$As), mean_As_minus = mean(res$minus$As),
    var_blup_within = var_blup_within,
    var_blup_plus = var_blup_plus, var_blup_minus = var_blup_minus,
    PEV_plus = pev_plus, PEV_minus = pev_minus,
    sigma2_within_plus = sigma2_within_plus,
    sigma2_within_minus = sigma2_within_minus,
    sigma2_within = mean(c(sigma2_within_plus, sigma2_within_minus)),
    # Spread of the PREDICTED effect among the 2n selected accessions, which is
    # what the continuous-slope estimand regresses on: Var(x) over both pools is
    # (dAs/2)^2 + var_As_within_pools. Returned here because only build_pools()
    # knows which accessions were selected.
    var_As_within_pools = var_blup_within,
    # The candidate set this selection was made from. Returned so a caller --
    # or a test -- can assert against it instead of re-deriving the filter,
    # which three places in tests/test_validation.R used to do by hand.
    candidates = cand
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
#' @param .components TRUE returns the variance split and the df alongside the
#'   SE, so `contrast_power()` does not recompute them. Duplicating the split
#'   in two functions is how an SE and its df drift apart.
contrast_se <- function(sigma2_within, ..., n, sigma2_e, P,
                        interaction_frac = NULL, interaction_sd = NULL,
                        delta = 0, n_loc = 5L, df_e_lost = 4L,
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

  v_acc <- sum(sigma2_within) / n * (if (length(sigma2_within) == 1L) 2 else 1)
  v_int <- int_sd^2 / n_loc
  v_plt <- 4 * sigma2_e / P
  v     <- v_acc + v_int + v_plt

  if (!.components) return(sqrt(v))

  # Welch for the accession stratum: with unequal pool variances the two halves
  # of the contrast carry different weight, and this collapses to exactly
  # 2n - 2 when they are equal.
  a <- (if (length(sigma2_within) == 1L) sigma2_within else sigma2_within[1]) / n
  b <- (if (length(sigma2_within) == 1L) sigma2_within else sigma2_within[2]) / n
  df_acc <- (a + b)^2 / (a^2 / (n - 1) + b^2 / (n - 1))

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
                           df_mode = c("satterthwaite", "naive")) {
  if (...length() > 0L) {
    stop("contrast_power() takes only `dAs` positionally; ", ...length(),
         " extra positional argument(s) given. Name every argument: ",
         "contrast_power(dAs = , sigma2_within = , n = , sigma2_e = , P = , ...)",
         call. = FALSE)
  }
  df_mode <- match.arg(df_mode)

  delta <- lambda * dAs
  cmp   <- contrast_se(sigma2_within = sigma2_within, n = n,
                       sigma2_e = sigma2_e, P = P,
                       interaction_frac = interaction_frac,
                       interaction_sd = interaction_sd,
                       delta = delta, n_loc = n_loc, df_e_lost = df_e_lost,
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
power_grid <- function(inp, n_values, P_values, lambda_values = c(1, 0.8, 0.6),
                       sided = 1, interaction_values = 0, pr_quantile = 0.5,
                       n_loc = 5L, n_loc_values = NULL,
                       interaction_mode = c("frac", "sd"),
                       df_mode = "satterthwaite") {
  interaction_mode <- match.arg(interaction_mode)
  n_loc_values <- n_loc_values %||% n_loc

  purrr::map(n_values, \(n) {
    bp <- build_pools(inp, n, pr_quantile = pr_quantile)
    s2w <- c(bp$sigma2_within_plus, bp$sigma2_within_minus)
    tidyr::expand_grid(P = P_values, lambda = lambda_values,
                       interaction = interaction_values,
                       n_locations = n_loc_values) |>
      purrr::pmap(\(P, lambda, interaction, n_locations) {
        args <- list(dAs = bp$dAs, sigma2_within = s2w, n = n,
                     sigma2_e = inp$sigma2_e, P = P, lambda = lambda,
                     sided = sided, n_loc = n_locations, df_mode = df_mode)
        args[[if (interaction_mode == "sd") "interaction_sd"
              else "interaction_frac"]] <- interaction
        do.call(contrast_power, args) |>
          dplyr::mutate(n_locations = n_locations, .after = P)
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
