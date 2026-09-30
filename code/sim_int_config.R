# ============================================================
# THE INTERACTION-FOCUSED DESIGN
#
# A second, smaller simulation whose only job is to locate the boundary found by
# docs/BGLR_vs_MegaLMM.md: MegaLMM recovers the oat x pea interaction better than
# a Kronecker kernel does, but only above a density threshold, and only when the
# interaction is low rank. The main design could not say where either boundary
# sits, for two separate reasons.
#
# WHY A SEPARATE DESIGN RATHER THAN MORE OF THE MAIN ONE.
#
# The main grid's MegaLMM sweep is a D-OPTIMAL FRACTION. That was the right
# choice for its question -- estimate every main effect and two-way interaction
# of seven factors at a fifth of the cost -- but it means most scenarios were
# never given most settings. Any analysis that conditions on a particular
# configuration (say K = 5 with the main effect pinned) is left with 40% of the
# scenarios, and the survivors are unbalanced across exactly the axis the
# question is about. A conditional question needs a design where the condition
# is present everywhere, so this one is a FULL FACTORIAL.
#
# WHAT CHANGED, AND WHY
#
#   sparsity      1.6 / 4.8 / 9 / 16 / 48 %. The 9% level is new: the switch lies
#                 between 4.8% (MegaLMM behind by 0.24) and 16% (ahead by 0.25),
#                 and four levels cannot say where.
#   n_factors     1 / 3 / 5. Rank 3 is new, for the same reason on the other
#                 axis: at 16% observed, rank 1 wins by 0.37 and rank 5 by 0.03,
#                 so the rank boundary is somewhere between and unlocated.
#   fixed_main_effect  KEPT AND FULLY CROSSED. This is the open question rather
#                 than a setting to fix. Pinning helps GMA when sparse and costs
#                 interaction recovery everywhere, so there is no single rule
#                 that is right for both responses -- and nobody knows which way
#                 it goes at 9%. Crossing it lets the simulation answer that
#                 instead of assuming it.
#   K             fixed at 5. The main design's ANOVA found 5 better than 10 on
#                 every response, consistently though not by much.
#   eigen_variance  fixed at 0.75. Found inert: partial omega squared was
#                 NEGATIVE on all four responses, i.e. it explained less than its
#                 one degree of freedom would by chance.
#   environment   one / ten_gxe only. `ten_stable` was indistinguishable from
#                 `one` (0.735 vs 0.740 on r_addsurf_oat), so the level that costs
#                 accuracy is GxE, not the number of environments. Dropping the
#                 middle level halves the grid at no cost to the question.
#   additive      not fitted. See run_scenario_bglr(models = ).
#   interaction   no "none" level: r_int is undefined without an interaction, and
#                 this design is about the interaction.
#
# 120 data scenarios, 2 MegaLMM runs each (pinned and not), fully crossed. No
# fraction, so every combination of every factor exists.
# ============================================================

SIM_INT_LEVELS <- list(
  n_acc             = c(200L, 400L),
  sparsity          = c(0.016, 0.048, 0.09, 0.16, 0.48),
  n_factors         = c(1L, 3L, 5L),
  interaction_pct   = c(0.10, 0.20),
  n_envs            = c(1L, 10L)
)

# n_envs = 10 always carries GxE here; the stable case was the one that behaved
# like a single environment, so it is the level that adds nothing.
SIM_INT_GXE <- c(`1` = 1.0, `10` = 0.6)

SIM_INT_MEGALMM <- list(
  K                 = 5L,
  eigen_variance    = 0.75,
  fixed_main_effect = c(FALSE, TRUE)     # the open question, fully crossed
)

SIM_INT_BGLR_MODELS <- "dge_ige"

# Deliberately far above the main design's range rather than a day later than
# it. SIM_BASE_SEED is 20260921 and its seeds run to about 20262100 at ten
# replicates, so a base 9 apart would have made 591 of 600 seeds collide. The two
# simulations would then have drawn from identical RNG streams for unrelated
# scenarios -- harmless while their caches and scenario names stay separate, but
# a trap for anyone who later pooled them, and free to avoid.
SIM_INT_BASE_SEED <- 20270000L
SIM_INT_REPS      <- 3L

#' The full factorial, one row per data scenario.
#'
#' Seeds follow the main design's rule -- replicate is the SLOW index -- so
#' adding replicates later reuses what is already cached rather than renumbering
#' it. See code/sim_config.R.
sim_int_grid <- function(levels = SIM_INT_LEVELS, n_reps = SIM_INT_REPS) {
  sc <- tidyr::expand_grid(!!!levels) |>
    dplyr::mutate(
      gxe_cor = unname(SIM_INT_GXE[as.character(n_envs)]),
      scenario = sprintf("int_n%d_sp%03d_f%d_i%02d_e%02d_g%03d",
                         n_acc, round(sparsity * 1000), n_factors,
                         round(interaction_pct * 100), n_envs,
                         round(gxe_cor * 100))
    ) |>
    dplyr::arrange(n_acc, sparsity, n_factors, interaction_pct, n_envs)

  n_sc <- nrow(sc)
  tidyr::expand_grid(rep = seq_len(n_reps), sc) |>
    dplyr::mutate(seed = SIM_INT_BASE_SEED + (rep - 1L) * n_sc +
                    rep(seq_len(n_sc), times = n_reps)) |>
    dplyr::relocate(scenario, rep, seed)
}

#' Every MegaLMM run: each data scenario crossed with each setting.
sim_int_runs <- function(scenarios = sim_int_grid(),
                         mm = SIM_INT_MEGALMM) {
  tidyr::expand_grid(dplyr::select(scenarios, scenario, rep, seed),
                     tidyr::expand_grid(!!!mm)) |>
    dplyr::left_join(scenarios, by = c("scenario", "rep", "seed"))
}
