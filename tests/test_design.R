# test_design.R
#
# The simulation's design: the collapses that stop duplicate scenarios, the
# composite recoding that makes a two-way model estimable, the D-optimal
# fraction, and the cache key.
#
# The headline test here is a NEGATIVE one. A main-effects-plus-two-way model is
# singular on the FULL grid, because two axes are nested rather than crossed.
# That is easy to mistake for a consequence of fractioning, so the test asserts
# the singularity directly -- if someone "fixes" the recoding away, this fails.
#
# Run: Rscript tests/test_design.R

library(tidyverse)
here::i_am("tests/test_design.R")
source(here::here("tests", "helper.R"))
load_code("sim_config.R")

# ------------------------------------------------------------
# 1. sim_grid: the collapses, and unique names
# ------------------------------------------------------------

g <- sim_grid(SIM_LEVELS, 1)

check(nrow(g) == dplyr::n_distinct(g$scenario),
      "every scenario has a distinct name")
check(nrow(g) == dplyr::n_distinct(g$seed), "and a distinct seed")

# interaction_pct scales nothing at n_factors = 0, so those cells must appear
# once, not once per level
check(all(g$interaction_pct[g$n_factors == 0] == 0),
      "interaction_pct is pinned to 0 when there is no interaction")
check(dplyr::n_distinct(dplyr::filter(g, n_factors == 0)$scenario) ==
        nrow(dplyr::filter(g, n_factors == 0)),
      "and no duplicate no-interaction scenarios survive")

# gxe_cor has nothing to vary across with one environment
check(all(g$gxe_cor[g$n_envs == 1] == 1),
      "gxe_cor is pinned to 1 when there is one environment")
check(dplyr::n_distinct(dplyr::filter(g, n_envs == 1)$scenario) ==
        nrow(dplyr::filter(g, n_envs == 1)),
      "and no duplicate single-environment scenarios survive")

# the expected size, derived rather than hard-coded
n_int <- length(SIM_LEVELS$n_factors[SIM_LEVELS$n_factors > 0]) *
           length(SIM_LEVELS$interaction_pct) + 1L          # +1 for "none"
n_env <- 1L + (length(SIM_LEVELS$n_envs) - 1L) * length(SIM_LEVELS$gxe_cor)
check(nrow(g) == length(SIM_LEVELS$n_acc) * length(SIM_LEVELS$sparsity) *
        n_int * n_env,
      "the collapsed grid is the size the collapses imply")

# sparsity is written in tenths of a percent, so 1.6% and 2% cannot collide
check(any(stringr::str_detect(g$scenario, "sp016")),
      "the scenario name records sparsity finely enough for 1.6%")
check(!any(stringr::str_detect(g$scenario, "sp02_")),
      "and does not round 1.6% to a whole percent")

# ------------------------------------------------------------
# 2. The nesting really does make the raw coding singular
#
# This is the finding the composite recoding exists for. If it ever stops being
# true, the recoding is unnecessary and should go.
# ------------------------------------------------------------

raw <- tidyr::expand_grid(!!!SIM_LEVELS) |>
  dplyr::mutate(
    interaction_pct = dplyr::if_else(n_factors == 0, 0, interaction_pct),
    gxe_cor = dplyr::if_else(n_envs == 1, 1, gxe_cor)) |>
  dplyr::distinct() |>
  dplyr::mutate(dplyr::across(dplyr::everything(), factor))

X_raw <- stats::model.matrix(~ .^2, data = as.data.frame(raw))
check(qr(X_raw)$rank < ncol(X_raw),
      "a two-way model on the RAW coding is singular, even on the full grid")

cand <- sim_candidates()
X_comp <- stats::model.matrix(
  ~ .^2,
  data = as.data.frame(lapply(
    cand[c("n_acc", "sparsity", "interaction", "environment",
           "K", "eigen_variance", "fixed_main_effect")], factor)))
check(qr(X_comp)$rank == ncol(X_comp),
      "the composite recoding is full rank, which is why it is used")

# ------------------------------------------------------------
# 3. The composite maps are complete and consistent
#
# Every composite level must map back to something the generator can run, and
# every combination the generator would accept must have a composite level.
# ------------------------------------------------------------

check(nrow(SIM_INTERACTION_MAP) == n_int,
      "the interaction map has one level per surviving combination")
check(nrow(SIM_ENVIRONMENT_MAP) == n_env,
      "the environment map has one level per surviving combination")
check(all(SIM_INTERACTION_MAP$n_factors %in% SIM_LEVELS$n_factors),
      "every mapped n_factors is a level of the design")
check(all(SIM_ENVIRONMENT_MAP$n_envs %in% SIM_LEVELS$n_envs),
      "every mapped n_envs is a level of the design")
check(all(SIM_INTERACTION_MAP$interaction_pct[SIM_INTERACTION_MAP$n_factors == 0] == 0),
      "the no-interaction level carries zero interaction variance")
check(all(SIM_ENVIRONMENT_MAP$gxe_cor[SIM_ENVIRONMENT_MAP$n_envs == 1] == 1),
      "the single-environment level carries no GxE")

check(nrow(cand) == length(SIM_LEVELS$n_acc) * length(SIM_LEVELS$sparsity) *
        n_int * n_env * length(SIM_MEGALMM_LEVELS$K) *
        length(SIM_MEGALMM_LEVELS$eigen_variance) *
        length(SIM_MEGALMM_LEVELS$fixed_main_effect),
      "the candidate set is the full crossing of design and MegaLMM levels")

# ------------------------------------------------------------
# 4. sim_design: the fraction estimates what it claims to
# ------------------------------------------------------------

d <- sim_design(n_runs = SIM_DESIGN_RUNS)
check(nrow(d$runs) == SIM_DESIGN_RUNS, "the fraction has the requested size")
check(d$model_rank == d$model_terms,
      "and is full rank for main effects plus all two-way interactions")
check(nrow(d$scenarios) == nrow(g),
      "it still touches every data scenario, so DGE-IGE loses no coverage")

# every MegaLMM setting must appear, or its main effect is not estimable
mm_seen <- dplyr::distinct(d$runs, K, eigen_variance, fixed_main_effect)
check(nrow(mm_seen) == length(SIM_MEGALMM_LEVELS$K) *
        length(SIM_MEGALMM_LEVELS$eigen_variance) *
        length(SIM_MEGALMM_LEVELS$fixed_main_effect),
      "every MegaLMM setting appears somewhere in the fraction")

# the runs must carry what the generator needs, joined from the composites
for (col in c("n_factors", "interaction_pct", "n_envs", "gxe_cor",
              "scenario", "seed", "rep")) {
  check(col %in% names(d$runs), paste("the run list carries", col))
}
check(!anyNA(d$runs$scenario), "every run is matched to a scenario")
check(!anyNA(d$runs$seed), "every run has a seed")

# reproducible: the same seed must give the same design
d2 <- sim_design(n_runs = SIM_DESIGN_RUNS)
check(identical(d$runs$scenario, d2$runs$scenario),
      "the design is reproducible from its seed")

# a fraction too small to estimate the model must be REFUSED, with a message
# that names the number of parameters rather than AlgDesign's "columns in
# expanded X"
check_error(sim_design(n_runs = 60),
            "a fraction below the parameter count is refused")
err <- tryCatch(sim_design(n_runs = 60), error = conditionMessage)
check(grepl("parameters", err) && grepl("82", err),
      "and the message says how many parameters there are")

# --full must return the whole candidate set
full <- sim_design(full = TRUE)
check(nrow(full$runs) == nrow(cand), "full = TRUE runs every candidate")
check(full$model_rank == full$model_terms, "and is full rank")

# replicates multiply scenarios, not the run list's distinct settings
d3 <- sim_design(n_runs = SIM_DESIGN_RUNS, n_reps = 2)
check(nrow(d3$scenarios) == 2 * nrow(g), "two replicates double the scenarios")
check(dplyr::n_distinct(d3$scenarios$seed) == nrow(d3$scenarios),
      "and every replicate gets its own seed")

# ------------------------------------------------------------
# 5. --reps is ADDITIVE
#
# The load-bearing property of the seed numbering, and the one that is easy to
# lose. `rep` is the SLOW index, so a scenario's seed does not depend on how
# many replicates were requested. Get this wrong -- number the rows of
# expand_grid(scenarios, rep), which varies rep fastest -- and replicate 1 of
# --reps 5 is a different seed from replicate 1 of --reps 1. Adding replicates
# to a finished grid then silently refits all of it and leaves the old cache
# behind as orphans that still glob into the combined CSV, both calling
# themselves replicate 1 and indistinguishable once there.
#
# This is not hypothetical: it happened on the September 2026 run and put 138
# duplicated rows in the results.
# ------------------------------------------------------------

rep1_of_1 <- dplyr::arrange(sim_design(n_runs = SIM_DESIGN_RUNS, n_reps = 1)$scenarios,
                            scenario)
for (n in c(2L, 5L)) {
  d_n <- sim_design(n_runs = SIM_DESIGN_RUNS, n_reps = n)
  rep1 <- dplyr::arrange(dplyr::filter(d_n$scenarios, rep == 1), scenario)

  check(identical(rep1$scenario, rep1_of_1$scenario),
        sprintf("n_reps = %d covers the same scenarios in replicate 1", n))
  check(identical(rep1$seed, rep1_of_1$seed),
        sprintf("and gives them the SAME seeds as n_reps = 1 -- so --reps %d is additive", n))
  check(length(intersect(dplyr::filter(d_n$scenarios, rep > 1)$seed,
                         rep1_of_1$seed)) == 0,
        sprintf("while replicates 2..%d are entirely new seeds", n))
  check(dplyr::n_distinct(d_n$scenarios$seed) == nrow(d_n$scenarios),
        sprintf("every seed is distinct at n_reps = %d", n))
}

# The seed must still move when the GRID changes, which is the reason it is in
# the cache key at all: a renumbered grid must MISS its old files rather than
# reuse them under a seed it would never have assigned.
smaller <- modifyList(SIM_LEVELS, list(sparsity = SIM_LEVELS$sparsity[-1]))
d_small <- sim_design(levels = smaller, n_runs = SIM_DESIGN_RUNS)
shared <- dplyr::inner_join(
  dplyr::select(d_small$scenarios, scenario, seed),
  dplyr::select(rep1_of_1, scenario, seed_full = seed), by = "scenario")
check(nrow(shared) > 0 && any(shared$seed != shared$seed_full),
      "dropping a level renumbers seeds, so stale cache files miss rather than being reused")

finish("design tests")
