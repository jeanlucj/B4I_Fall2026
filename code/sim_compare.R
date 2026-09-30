# ============================================================
# BGLR versus MegaLMM, at MegaLMM's best settings
#
#   Rscript code/sim_compare.R
#
# Two questions, and they have different answers:
#   GMA          -- does analysing both yields jointly (BGLR Multitrait) beat
#                   two single-orientation factor models?
#   INTERACTION  -- can the factor model decompose oat x pea specific
#                   combining ability better than a Kronecker kernel can?
#
# THE COMPARISON IS PAIRED. Every model is fitted to the SAME simulated data --
# same scenario, same replicate, same seed, same train/held split -- so the unit
# of analysis is the within-scenario difference, not the group mean. Comparing
# marginal means would throw away the pairing and drown the contrast in the
# scenario-to-scenario variation that sparsity alone drives from r = 0.15 to
# r = 0.91.
#
# "MegaLMM AT ITS BEST" IS ESTIMATED, NOT SUBSETTED. The obvious approach --
# keep only the rows with the preferred settings -- fails here, because the run
# list is a D-optimal FRACTION: most scenarios were never given the preferred
# combination at all. Requiring K = 5 together with the sparsity-dependent
# pinning rule keeps 192 of 480 scenario-replicates, and the survivors are
# unbalanced across exactly the generative axis the comparison is about
# (46% of the rank-1 cells against 33% of the rank-5 ones). So the settings go
# into the model as factors and the preferred configuration is read off with
# emmeans, which uses all 480 and is balanced by construction. The subset is
# still computed, as a check that the two agree.
#
# Replicate 1 is excluded throughout: it holds two sets of cache files under
# different seeds. See code/sim_filter_results.R.
# ============================================================

library(tidyverse)

here::i_am("code/sim_compare.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
out_dir  <- here::here("output")
res_path <- arg_value("--results",
                      file.path(out_dir, "simulation_results_from_ceres_filtered.csv"))

GEN  <- c("interaction", "environment")
DAT  <- c("n_acc", "sparsity")
SET  <- c("K", "eigen_variance", "fixed_main_effect")

fz  <- function(r) atanh(pmax(pmin(r, 1 - 1e-8), -1 + 1e-8))
ifz <- tanh

# The rule under test: K = 5 always, pin the main effect only where the ANOVA
# showed it helps -- at 1.6% and 4.8% observed.
PREFER_K   <- 5
prefer_pin <- function(sparsity) sparsity <= 0.048

# Results written before the rename carry r_gma_oat / r_gma_pea for what is now
# r_addsurf_oat / r_addsurf_pea. Same quantity -- the additive part of that
# trait's yield surface -- renamed because `r_gma_oat` and `r_oat_gma` are
# different things and two words in a different order was too fine a distinction
# to hang that on. Renaming on read keeps the September Ceres results usable.
read_sim_results <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  old <- c(r_addsurf_oat = "r_gma_oat", r_addsurf_pea = "r_gma_pea")
  present <- old[old %in% names(x) & !names(old) %in% names(x)]
  if (length(present) > 0) {
    message("renaming ", paste(present, collapse = ", "),
            " from a pre-rename results file")
    x <- dplyr::rename(x, !!!present)
  }
  x
}

raw <- read_sim_results(res_path) |> dplyr::filter(rep > 1)

# NOTE ON THE INPUT. These read the September 2026 Ceres results, which used the
# pre-October scoring scheme: 20% of the observed cells held out, so the models
# were fitted at 0.8x the labelled sparsity and the per-cell metrics were scored
# on the held-out cells. Read every sparsity level as 0.8x its label. Results
# generated from now on use the never-observed cells and are fitted at the
# labelled density; they also name the fit statistic r_fit_ rather than r_obs_.
# See SIMULATION_GLOSSARY.md.

RESP <- c("r_addsurf_oat", "r_addsurf_pea", "r_int_oat", "r_int_pea")

long <- raw |>
  dplyr::select(model, scenario, rep, dplyr::all_of(c(GEN, DAT, SET)),
                dplyr::all_of(RESP)) |>
  tidyr::pivot_longer(dplyr::all_of(RESP), names_to = "response", values_to = "r") |>
  dplyr::mutate(z = fz(r))

# The two BGLR comparators, one value per scenario-replicate-response.
#
# `additive` has no interaction term, so its r_int is NA by construction rather
# than poor -- for the interaction responses the only BGLR comparator that
# exists is dge_ige. That is not a gap in the data, it is what "additive" means.
bglr <- long |>
  dplyr::filter(model %in% c("additive", "dge_ige"), !is.na(z)) |>
  dplyr::select(scenario, rep, response, comparator = model, z_bglr = z)

mm <- long |>
  dplyr::filter(model == "megalmm", !is.na(z)) |>
  dplyr::select(scenario, rep, response, dplyr::all_of(c(GEN, DAT, SET)),
                z_mm = z)

paired <- dplyr::inner_join(mm, bglr, by = c("scenario", "rep", "response"),
                            relationship = "many-to-many") |>
  dplyr::mutate(d = z_mm - z_bglr,
                unit = paste(scenario, rep, sep = "|"),
                dplyr::across(dplyr::all_of(c(GEN, DAT, SET)), factor))

message(dplyr::n_distinct(paired$unit), " scenario-replicates, ",
        nrow(paired), " paired comparisons")

# ------------------------------------------------------------
# 1. The headline: MegaLMM at the preferred settings, by subsetting
#
# Straightforward and easy to check, but on 40% of the data.
# ------------------------------------------------------------

subset_rule <- paired |>
  dplyr::filter(K == PREFER_K,
                fixed_main_effect == prefer_pin(as.numeric(as.character(sparsity)))) |>
  # where both eigen_variance levels survive, average them: the ANOVA found the
  # lever inert (negative partial omega squared), so this is noise reduction
  # rather than a choice
  dplyr::group_by(unit, response, comparator,
                  dplyr::across(dplyr::all_of(c(GEN, DAT)))) |>
  dplyr::summarise(d = mean(d), n_ev = dplyr::n(), .groups = "drop")

paired_test <- function(d) {
  if (length(d) < 3 || stats::sd(d) == 0) {
    return(tibble::tibble(n = length(d), mean_d = mean(d), lo = NA_real_,
                          hi = NA_real_, p = NA_real_, win = mean(d > 0)))
  }
  t <- stats::t.test(d)
  tibble::tibble(n = length(d), mean_d = mean(d),
                 lo = t$conf.int[1], hi = t$conf.int[2],
                 p = t$p.value, win = mean(d > 0))
}

# Differences are on the z scale; shown on the r scale as the gap at the
# comparator's own mean, which is where the comparison actually sits.
to_r <- function(mean_d, z_ref) ifz(z_ref + mean_d) - ifz(z_ref)

ref <- bglr |> dplyr::group_by(response, comparator) |>
  dplyr::summarise(z_ref = mean(z_bglr), .groups = "drop")

overall <- subset_rule |>
  dplyr::group_by(response, comparator) |>
  dplyr::reframe(paired_test(d)) |>
  dplyr::left_join(ref, by = c("response", "comparator")) |>
  dplyr::mutate(dr = to_r(mean_d, z_ref), .after = mean_d)

cat("\n", strrep("=", 78),
    "\n1. MegaLMM (K=5, pin when sparse) minus BGLR -- subset, paired\n",
    strrep("=", 78), "\n", sep = "")
print(overall |> dplyr::transmute(
  response, comparator, n,
  `mean d (z)` = round(mean_d, 3), `gap (r)` = round(dr, 3),
  `95% CI (z)` = paste0("[", round(lo, 3), ", ", round(hi, 3), "]"),
  p = format.pval(p, digits = 2, eps = 1e-16),
  `MegaLMM wins` = paste0(round(100 * win), "%")), n = 20, width = 200)

# ------------------------------------------------------------
# 2. The same contrast estimated from ALL the data
#
# Settings as factors, scenario-replicate as a random intercept because the two
# MegaLMM runs of one scenario share a BGLR value and are not independent.
# ------------------------------------------------------------

model_based <- purrr::map(RESP, function(resp) {
  purrr::map(c("additive", "dge_ige"), function(cmp) {
    d <- dplyr::filter(paired, response == resp, comparator == cmp)
    if (nrow(d) < 50) return(NULL)
    d <- droplevels(d)
    if (dplyr::n_distinct(d$K) < 2) return(NULL)
    fit <- lme4::lmer(
      stats::as.formula(paste("d ~", paste(c(SET, GEN, DAT), collapse = " + "),
                              "+ (1|unit)")),
      data = d, REML = TRUE,
      control = lme4::lmerControl(calc.derivs = FALSE))
    # the preferred configuration, averaged over eigen_variance, and averaged
    # over the generative factors at equal weight
    em <- suppressMessages(emmeans::emmeans(
      fit, specs = ~ K + fixed_main_effect + sparsity, data = d))
    s <- tibble::as_tibble(summary(em)) |>
      dplyr::mutate(dplyr::across(c(K, sparsity), ~ as.numeric(as.character(.x))),
                    fixed_main_effect = as.logical(as.character(fixed_main_effect))) |>
      dplyr::filter(K == PREFER_K, fixed_main_effect == prefer_pin(sparsity))
    tibble::tibble(response = resp, comparator = cmp,
                   mean_d = mean(s$emmean),
                   se = sqrt(sum(s$SE^2)) / nrow(s),
                   n_obs = nrow(d), n_units = dplyr::n_distinct(d$unit))
  }) |> purrr::compact() |> purrr::list_rbind()
}) |> purrr::list_rbind() |>
  dplyr::left_join(ref, by = c("response", "comparator")) |>
  dplyr::mutate(dr = to_r(mean_d, z_ref))

cat("\n", strrep("=", 78),
    "\n2. The same contrast from ALL 480 scenario-replicates (model-based)\n",
    strrep("=", 78), "\n", sep = "")
print(model_based |> dplyr::transmute(
  response, comparator, n_units, n_obs,
  `mean d (z)` = round(mean_d, 3), `gap (r)` = round(dr, 3),
  `SE` = round(se, 3)), n = 20)

cat("\nsubset vs model-based, mean d on the z scale:\n")
print(dplyr::left_join(
  dplyr::select(overall, response, comparator, subset = mean_d),
  dplyr::select(model_based, response, comparator, model_based = mean_d),
  by = c("response", "comparator")) |>
  dplyr::mutate(dplyr::across(c(subset, model_based), ~ round(.x, 3)),
                diff = round(model_based - subset, 3)), n = 20)

# ------------------------------------------------------------
# 3. Does the answer depend on the generative model?
#
# The hypothesis under test: MegaLMM should do relatively better when the
# interaction is simple (rank 1), because a low-rank surface is what a factor
# model is built to find.
# ------------------------------------------------------------

by_factor <- function(fac) {
  subset_rule |>
    dplyr::group_by(response, comparator, level = .data[[fac]]) |>
    dplyr::reframe(paired_test(d)) |>
    dplyr::left_join(ref, by = c("response", "comparator")) |>
    dplyr::mutate(factor = fac, dr = to_r(mean_d, z_ref), .before = level)
}

by_gen <- purrr::map(c(GEN, DAT), by_factor) |> purrr::list_rbind()

cat("\n", strrep("=", 78),
    "\n3. The contrast, broken down by the generative model and the data\n",
    strrep("=", 78), "\n", sep = "")
for (resp in RESP) {
  cat("\n--", resp, "--\n")
  print(by_gen |> dplyr::filter(response == resp) |>
          dplyr::transmute(comparator, factor, level, n,
                           `gap (r)` = round(dr, 3),
                           p = format.pval(p, digits = 2, eps = 1e-16),
                           wins = paste0(round(100 * win), "%")),
        n = 40, width = 200)
}

readr::write_csv(dplyr::bind_rows(
  dplyr::mutate(overall, factor = "overall", level = "all"),
  by_gen), file.path(out_dir, "simulation_compare.csv"))
readr::write_csv(model_based, file.path(out_dir, "simulation_compare_modelbased.csv"))

# ------------------------------------------------------------
# 4. How much is the "best shot" worth, and how much of it is selection?
#
# Three quantities. If (a) and (b) are close the rule is not buying much and
# the question of bias is moot; the distance from (c) is what a genuinely
# adaptive, per-cell choice would add, and is an upper bound on what any
# selection rule could achieve.
# ------------------------------------------------------------

cmp_best <- paired |>
  dplyr::group_by(unit, response, comparator) |>
  dplyr::summarise(
    rule = mean(d[K == PREFER_K &
                    fixed_main_effect ==
                      prefer_pin(as.numeric(as.character(sparsity)))]),
    average = mean(d),
    oracle = max(d),
    .groups = "drop") |>
  dplyr::group_by(response, comparator) |>
  dplyr::summarise(n_with_rule = sum(!is.nan(rule)),
                   dplyr::across(c(rule, average, oracle),
                                 ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  dplyr::left_join(ref, by = c("response", "comparator")) |>
  dplyr::mutate(dplyr::across(c(rule, average, oracle), ~ to_r(.x, z_ref)))

cat("\n", strrep("=", 78),
    "\n4. What the setting rule is worth, on the r scale\n",
    strrep("=", 78), "\n", sep = "")
print(cmp_best |> dplyr::transmute(
  response, comparator, n_with_rule,
  `rule` = round(rule, 3), `average of all settings` = round(average, 3),
  `per-cell oracle` = round(oracle, 3),
  `rule - average` = round(rule - average, 3),
  `oracle - rule` = round(oracle - rule, 3)), n = 20, width = 200)
cat("\n`oracle` takes the best setting per scenario KNOWING the answer, so it is\n",
    "not achievable -- it is the ceiling any selection rule is working towards.\n")
readr::write_csv(cmp_best, file.path(out_dir, "simulation_compare_settings.csv"))

# ------------------------------------------------------------
# 5. The question that decides whether any of this is usable
#
# MegaLMM's interaction advantage turns out to depend on TWO things at once --
# a low-rank interaction, and enough observed cells. Whether those two
# requirements are separable or must both hold is what says whether the result
# transfers to a real, sparse experiment, so it gets its own table rather than
# being left to two marginal breakdowns.
# ------------------------------------------------------------

arch_sparsity <- subset_rule |>
  dplyr::filter(stringr::str_starts(response, "r_int"), interaction != "none") |>
  dplyr::mutate(architecture = paste0(
    stringr::str_match(as.character(interaction), "^f(\\d+)")[, 2], " factor(s)")) |>
  dplyr::group_by(response, architecture, sparsity) |>
  dplyr::reframe(paired_test(d)) |>
  dplyr::left_join(dplyr::filter(ref, comparator == "dge_ige") |>
                     dplyr::select(response, z_ref), by = "response") |>
  dplyr::mutate(dr = to_r(mean_d, z_ref))

cat("\n", strrep("=", 78),
    "\n5. Interaction recovery: MegaLMM minus dge_ige, by rank and density\n",
    strrep("=", 78), "\n", sep = "")
for (resp in c("r_int_oat", "r_int_pea")) {
  cat("\n--", resp, "-- gap on the r scale (positive = MegaLMM better)\n")
  print(arch_sparsity |> dplyr::filter(response == resp) |>
          dplyr::transmute(architecture, sparsity, n,
                           `gap (r)` = round(dr, 3),
                           wins = paste0(round(100 * win), "%"),
                           p = format.pval(p, digits = 2, eps = 1e-16)),
        n = 20, width = 200)
}
readr::write_csv(arch_sparsity, file.path(out_dir, "simulation_compare_arch_sparsity.csv"))
