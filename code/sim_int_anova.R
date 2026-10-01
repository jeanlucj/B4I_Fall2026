# ============================================================
# SPLIT-PLOT ANOVA OF THE sim_int INTERACTION ACCURACIES
#
#   Rscript code/sim_int_anova.R
#   Rscript code/sim_int_anova.R --response r_int_pea
#   Rscript code/sim_int_anova.R --paired output/simulation_int_paired.csv
#
# Asks which design factors decide whether MegaLMM or DGE-IGE recovers the
# oat x pea interaction better, and -- the question that matters for reporting
# -- whether the model x sparsity and model x n_factors interactions can be read
# averaged over everything else, or whether some third factor moves them.
#
# WHY A SPLIT PLOT. The two models are not two runs. They are two analyses of
# the SAME simulated dataset, so they are paired, not crossed. Treating `model`
# as a crossed factor would double the apparent replication and shrink every
# standard error. The dataset is the whole plot and the model is the subplot, and
# the two strata are obtained exactly -- the design is a complete factorial, so
# this is orthogonal -- by analysing per-dataset means and differences:
#
#   mean_z = (z_mm + z_bglr)/2  -> BETWEEN-DATASET stratum. Main effects and
#                                 interactions of everything except `model`.
#   diff_z =  z_mm - z_bglr     -> WITHIN-DATASET stratum. Its INTERCEPT is the
#                                 `model` main effect, and a term A in this
#                                 stratum IS the `model` x A interaction.
#
# THE DECODER RING for the within-dataset stratum, which is the whole point:
#
#   intercept        =  model
#   A                =  model x A            (A a whole-plot factor)
#   A:B              =  model x A x B
#
# So the three-way interactions this script was written to find -- does n_acc,
# or interaction_pct, or n_envs, or fixed_main_effect move the model x sparsity
# gap? -- are read off the TWO-WAY rows of the subplot table. `sparsity:n_acc`
# there is model x sparsity x n_acc. No three-way terms need to be fitted.
#
# FISHER z, not r. The responses run from 0 to 0.95 and r is bounded and
# heteroscedastic. atanh(r) has variance ~1/(n-3) whatever r is. Everything is
# fitted on z and reported back on the r scale, where it means something.
#
# FULL FACTORS, not 1 df each. An earlier pass entered sparsity and n_factors as
# continuous 1-df predictors. Measured: 82% of the whole-plot residual sum of
# squares was then lack of fit, because MegaLMM's accuracy saturates -- 0.00 at
# 1.6% of combinations observed to 0.88 at 48%. That inflates the error mean
# square about fivefold, making every F conservative and every omega^2 an
# understatement. Here both are factors (4 df and 2 df), so the residual is
# dataset-to-dataset noise and curvature is accounted for.
#
# EFFECT SIZE FIRST. With 1,200 datasets almost everything clears p < 0.05, so
# the p-value is not the interesting column. Reported instead:
#   partial omega^2  variance explained, bias-corrected, and negative when a
#                    term explains less than its degrees of freedom would by
#                    chance
#   dr               the spread the term produces on the r scale, from estimated
#                    marginal means. The practical column: a term can be
#                    overwhelmingly significant and still move r by 0.01.
# ============================================================

library(tidyverse)

here::i_am("code/sim_int_anova.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

RESPONSE  <- arg_value("--response", "r_int_oat")
out_dir   <- here::here("output", "simulation_int")
paired_in <- arg_value("--paired",
                       file.path(out_dir, "simulation_int_paired.csv"))

# The two factors the conclusions are stated in terms of, and the four the
# conclusions are averaged over. The split is what the "do I need to worry"
# section below tests.
FOCAL    <- c("sparsity", "n_factors")
NUISANCE <- c("n_acc", "interaction_pct", "n_envs", "fixed_main_effect")
WHOLE    <- c(FOCAL, NUISANCE)

ifz <- tanh

# ------------------------------------------------------------
# Data: one row per dataset, carrying both strata
# ------------------------------------------------------------

d <- readr::read_csv(paired_in, show_col_types = FALSE) |>
  dplyr::filter(response == RESPONSE, !is.na(z_mm), !is.na(z_bglr)) |>
  dplyr::mutate(mean_z = (z_mm + z_bglr) / 2, diff_z = d) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(WHOLE), factor))

if (nrow(d) == 0) stop("no rows for response ", RESPONSE, " in ", paired_in)
message(nrow(d), " datasets, response ", RESPONSE)

rhs   <- paste0("(", paste(WHOLE, collapse = " + "), ")^2")
fit_w <- stats::lm(stats::as.formula(paste("mean_z ~", rhs)), data = d)
fit_s <- stats::lm(stats::as.formula(paste("diff_z ~", rhs)), data = d)

# ------------------------------------------------------------
# The term table
# ------------------------------------------------------------

# Type II sums of squares. The design is a complete factorial so Type I, II and
# III agree; Type II is used because it does not depend on the order terms were
# typed, which matters if the design is ever run as a fraction.
tab <- function(fit, stratum) {
  a    <- car::Anova(fit, type = "II")
  df_e <- fit$df.residual
  ms_e <- stats::sigma(fit)^2
  tibble::as_tibble(a, rownames = "term") |>
    dplyr::filter(term != "Residuals") |>
    dplyr::rename(SS = `Sum Sq`, df = Df, F = `F value`, p = `Pr(>F)`) |>
    dplyr::mutate(
      stratum = stratum, MS = SS / df,
      # partial omega^2: bias-corrected share of variance, honest about terms
      # that explain less than chance
      omega2_p = (SS - df * ms_e) / (SS - df * ms_e + (df_e + 1) * ms_e),
      .before = SS)
}

# Practical magnitude, on the r scale.
#
# For a MAIN effect, the spread of its marginal means is what it does. For an
# INTERACTION it is not: the cell means contain the two main effects as well,
# so a weak interaction can print a larger spread than the strong main effect
# inside it. What the interaction itself does is the cell means with their
# additive part swept out, shown as a shift around the grand mean.
dr_of <- function(fit, term, stratum, base) {
  vars <- strsplit(term, ":", fixed = TRUE)[[1]]
  if (!all(vars %in% names(d))) return(NA_real_)
  em <- try(suppressMessages(emmeans::emmeans(fit, specs = vars, data = d)),
            silent = TRUE)
  if (inherits(em, "try-error")) return(NA_real_)
  s <- summary(em)
  if (length(vars) == 1) {
    rng <- range(s$emmean)
  } else {
    M <- tapply(s$emmean, list(s[[vars[1]]], s[[vars[2]]]), identity)
    M <- M - rowMeans(M)[row(M)] - colMeans(M)[col(M)] + mean(M)
    rng <- range(M) + if (stratum == "whole") base else 0
  }
  if (stratum == "whole") {
    diff(ifz(rng))
  } else {
    # A subplot term is a spread of MegaLMM-minus-DGE-IGE differences on the z
    # scale. Shown as the width of the model gap it implies on the r scale,
    # taken around the grand mean so it is comparable with the other column.
    diff(range(ifz(base + rng / 2) - ifz(base - rng / 2)))
  }
}

base <- mean(d$mean_z)
out  <- dplyr::bind_rows(tab(fit_w, "between-dataset (whole plot)"),
                         tab(fit_s, "within-dataset (subplot)"))
out$dr <- purrr::map2_dbl(out$term, out$stratum, \(tm, st)
  if (startsWith(st, "between")) dr_of(fit_w, tm, "whole", base)
  else dr_of(fit_s, tm, "sub", base))

# The `model` main effect is the subplot INTERCEPT, which car::Anova omits.
ic  <- summary(fit_s)$coefficients["(Intercept)", ]
out <- dplyr::bind_rows(
  tibble::tibble(term = "model (megalmm vs dge_ige)",
                 stratum = "within-dataset (subplot)", df = 1,
                 MS = NA_real_, omega2_p = NA_real_, SS = NA_real_,
                 F = ic[["t value"]]^2, p = ic[["Pr(>|t|)"]],
                 dr = abs(ifz(mean(d$z_mm)) - ifz(mean(d$z_bglr)))),
  out) |>
  # What each subplot row means once decoded, so the CSV can be read alone.
  dplyr::mutate(
    reads_as = dplyr::case_when(
      stratum == "between-dataset (whole plot)" ~ term,
      term == "model (megalmm vs dge_ige)"      ~ "model",
      TRUE                                      ~ paste0("model x ", term)),
    response = RESPONSE, .before = 1) |>
  dplyr::relocate(reads_as, .after = term)

readr::write_csv(out, file.path(
  out_dir, paste0("simulation_int_anova_", RESPONSE, ".csv")))

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

fmt   <- function(x, k = 3) formatC(x, format = "f", digits = k)
stars <- function(p) dplyr::case_when(is.na(p) ~ "", p < 1e-10 ~ "****",
                                      p < 1e-4 ~ "***", p < 0.01 ~ "**",
                                      p < 0.05 ~ "*", TRUE ~ "")

for (st in c("within-dataset (subplot)", "between-dataset (whole plot)")) {
  cat("\n", strrep("=", 78), "\n", st, "\n",
      if (startsWith(st, "within"))
        "a row `A` here is the model x A interaction; `A:B` is model x A x B\n"
      else "accuracy averaged over the two models\n",
      strrep("=", 78), "\n", sep = "")
  out |>
    dplyr::filter(stratum == st) |>
    dplyr::arrange(dplyr::desc(dplyr::coalesce(omega2_p, Inf))) |>
    dplyr::transmute(reads_as, df, F = fmt(F, 1), sig = stars(p),
                     omega2_p = fmt(omega2_p), dr = fmt(dr)) |>
    print(n = 30)
}

# ------------------------------------------------------------
# Expected mean squares
#
# The split plot is what makes this worth writing down: the two strata have
# DIFFERENT error terms, so a whole-plot factor and a model x factor interaction
# are not tested against the same thing.
#
# Model for one dataset i and model j:
#
#   z_ij = mu + (whole-plot fixed effects)_i + (model effect)_j
#              + (model x whole-plot)_ij + u_i + e_ij
#
#   u_i  ~ N(0, s2_data)  dataset-to-dataset: a different simulated panel under
#                         the same settings
#   e_ij ~ N(0, s2_sub)   model-to-model within one dataset, i.e. sampler noise
#                         and whatever one fitter does that the other does not
#
# Analysing the per-dataset mean and difference gives the two strata exactly:
#
#   mean_z_i = mu + (whole-plot)_i + u_i + e_bar_i          Var = s2_data + s2_sub/2
#   diff_z_i = (model) + (model x whole)_i + (e_i1 - e_i2)  Var = 2 s2_sub
#
# so, with a whole-plot factor A of a levels and n datasets per level:
#
#   BETWEEN-DATASET STRATUM (fitted on mean_z)
#     A                 df a-1   E(MS) = s2_data + s2_sub/2 + n/(a-1) SUM alpha^2
#     whole-plot error  df N-p   E(MS) = s2_data + s2_sub/2
#
#   WITHIN-DATASET STRATUM (fitted on diff_z)
#     model             df 1     E(MS) = 2 s2_sub + N delta^2
#     model x A         df a-1   E(MS) = 2 s2_sub + 2n/(a-1) SUM (delta.alpha)^2
#     subplot error     df N-p   E(MS) = 2 s2_sub
#
# Each stratum is tested against its own residual, which is what the two lm()
# fits do, and the components come out of the two error mean squares:
#
#     s2_sub  = MS_subplot_error / 2
#     s2_data = MS_wholeplot_error - MS_subplot_error / 4
# ------------------------------------------------------------

ms_w  <- stats::sigma(fit_w)^2
ms_s  <- stats::sigma(fit_s)^2
v_sub <- ms_s / 2
v_dat <- ms_w - v_sub / 2

ems <- tibble::tibble(
  response = RESPONSE, n_datasets = nrow(d),
  df_whole = fit_w$df.residual, df_sub = fit_s$df.residual,
  ms_whole_error = round(ms_w, 5), ms_whole_expects = "s2_data + s2_sub/2",
  ms_sub_error = round(ms_s, 5),   ms_sub_expects = "2 s2_sub",
  s2_data = round(v_dat, 5), s2_sub = round(v_sub, 5),
  data_share_of_noise = paste0(round(100 * v_dat / (v_dat + v_sub)), "%"))

cat("\n", strrep("=", 78), "\nExpected mean squares and variance components\n",
    strrep("=", 78), "\n", sep = "")
print(ems, width = 200)
readr::write_csv(ems, file.path(
  out_dir, paste0("simulation_int_anova_ems_", RESPONSE, ".csv")))

# ------------------------------------------------------------
# Can model x sparsity and model x n_factors be read averaged over the rest?
#
# The question is whether any nuisance factor moves those two interactions,
# which by the decoder ring is the subplot `focal:nuisance` rows -- the
# three-way model x focal x nuisance terms. A term is only a problem if it is
# large relative to the interaction it would be disturbing, so each one is
# reported as a fraction of its parent.
# ------------------------------------------------------------

sub <- dplyr::filter(out, stratum == "within-dataset (subplot)")
parent_of <- function(term) {
  v <- strsplit(term, ":", fixed = TRUE)[[1]]
  f <- intersect(v, FOCAL)
  if (length(f) == 1) f else NA_character_
}

worry <- sub |>
  dplyr::filter(grepl(":", term)) |>
  dplyr::mutate(focal = purrr::map_chr(term, parent_of)) |>
  dplyr::filter(!is.na(focal),
                purrr::map_lgl(term, \(t)
                  any(strsplit(t, ":", fixed = TRUE)[[1]] %in% NUISANCE))) |>
  dplyr::left_join(
    sub |> dplyr::filter(term %in% FOCAL) |>
      dplyr::select(focal = term, parent_omega2 = omega2_p, parent_dr = dr),
    by = "focal") |>
  dplyr::transmute(
    disturbs = paste0("model x ", focal), by = sub("^.*:|:.*$", "", term),
    term = paste0("model x ", term), df, F = round(F, 1), p,
    omega2_p = round(omega2_p, 3), dr = round(dr, 3),
    share_of_parent_omega2 = round(omega2_p / parent_omega2, 3),
    share_of_parent_dr = round(dr / parent_dr, 3),
    verdict = dplyr::case_when(
      p > 0.05                   ~ "ignore: not detectable",
      dr < 0.05                  ~ "ignore: detectable but moves r < 0.05",
      share_of_parent_dr < 0.25  ~ "note: real but small next to its parent",
      TRUE                       ~ "MIND THIS: materially reshapes its parent")) |>
  dplyr::arrange(disturbs, dplyr::desc(dr))

cat("\n", strrep("=", 78),
    "\nDoes anything disturb model x sparsity or model x n_factors?\n",
    strrep("=", 78), "\n", sep = "")
print(worry, n = 30, width = 220)
readr::write_csv(worry, file.path(
  out_dir, paste0("simulation_int_anova_threeway_", RESPONSE, ".csv")))

# ------------------------------------------------------------
# The same question answered as a picture rather than a test: where does the
# model gap cross zero? If the crossover does not move across the levels of a
# factor, averaging over it is safe whatever the F says.
#
# Two views, because the two focal factors are the two axes the conclusions are
# stated on. Each is faceted by every OTHER design factor, including the other
# focal one, so a crossover that moves shows up as a sign change moving along a
# row. `n_plots` is on the sparsity view because sparsity and n_acc are not
# independent axes in practice -- both buy plots, and sparsity * n_acc^2 is the
# plot count the budget actually pays for.
# ------------------------------------------------------------

gap_table <- function(columns) {
  facets <- setdiff(WHOLE, columns)
  rows <- purrr::map(c("(all)", facets), function(fc) {
    g <- if (fc == "(all)") d else dplyr::group_by(d, .data[[fc]])
    g |>
      dplyr::group_by(.data[[columns]], .add = TRUE) |>
      dplyr::summarise(n = dplyr::n(),
                       r_megalmm = mean(ifz(z_mm)),
                       r_dge_ige = mean(ifz(z_bglr)),
                       .groups = "drop") |>
      dplyr::mutate(facet = fc,
                    level = if (fc == "(all)") "(all)"
                            else as.character(.data[[fc]]),
                    gap_r = round(r_megalmm - r_dge_ige, 3),
                    favours = dplyr::if_else(gap_r > 0, "megalmm", "dge_ige"))
  }) |> purrr::list_rbind()

  out <- rows |>
    dplyr::transmute(response = RESPONSE, facet, level,
                     !!columns := .data[[columns]], n_per_cell = n,
                     r_megalmm = round(r_megalmm, 3),
                     r_dge_ige = round(r_dge_ige, 3), gap_r, favours)
  if (columns == "sparsity") {
    # plots evaluated, where the facet pins n_acc; NA where it is averaged over
    out <- out |>
      dplyr::mutate(n_plots = dplyr::if_else(
        facet == "n_acc",
        round(as.numeric(as.character(sparsity)) *
                suppressWarnings(as.numeric(level))^2),
        NA_real_), .after = n_per_cell)
  }
  out
}

for (cols in c("sparsity", "n_factors")) {
  g <- gap_table(cols)
  f <- file.path(out_dir, paste0("simulation_int_gap_by_", cols, "_",
                                 RESPONSE, ".csv"))
  readr::write_csv(g, f)
  cat("\n", strrep("=", 78),
      "\nMegaLMM minus DGE-IGE on the r scale, by ", cols,
      " within each other factor",
      "\n(the crossover is where the sign flips)\n",
      strrep("=", 78), "\n", sep = "")
  g |>
    dplyr::select(facet, level, dplyr::all_of(cols), gap_r) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(cols), values_from = gap_r,
                       names_prefix = paste0(substr(cols, 1, 2), "=")) |>
    print(n = 30, width = 200)
}

message("\nwrote the term table, the EMS, the three-way summary and the two gap\n",
        "tables to output/simulation_int/")
