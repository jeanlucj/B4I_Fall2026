# ============================================================
# ANOVA OF THE MegaLMM SIMULATION RESULTS
#
#   Rscript code/sim_anova.R
#   Rscript code/sim_anova.R --results output/..._filtered.csv
#
# Asks which of the three axes -- what was simulated, how much data there was,
# and how MegaLMM was configured -- move its accuracy, and by how much.
#
# THREE THINGS THIS DOES DIFFERENTLY FROM A PLAIN aov() ON r.
#
# 1. FISHER z, not r. The responses run from about 0 to 0.95 and r is both
#    bounded and heteroscedastic -- its sampling variance collapses as |r| -> 1,
#    so a plain ANOVA gives the high cells too much weight and the residuals
#    fan. atanh(r) has variance ~1/(n-3) whatever r is. Everything is fitted on
#    z and reported back on the r scale, where it means something.
#
# 2. SPLIT PLOT, because Eta and U are not two runs. They are two readouts of
#    the SAME MegaLMM fit -- Eta_mean the predicted phenotype, U the genetic
#    value without the per-column intercept -- so they are paired, not
#    independent. Treating them as a crossed factor would double the apparent
#    replication and shrink every standard error in the table.
#
#    The run is the whole plot and the readout is the subplot, so the two
#    strata are obtained exactly, by analysing per-run means and differences:
#
#      mean_z = (z_Eta + z_U)/2   -> the whole-plot stratum. Main effects and
#                                    interactions of everything except readout.
#      diff_z =  z_Eta - z_U      -> the subplot stratum. Its INTERCEPT is the
#                                    readout main effect, and a main effect of
#                                    factor A here IS the readout x A
#                                    interaction.
#
#    Each run contributes one mean and one difference, so this is orthogonal
#    and stays exact under the unbalanced D-optimal fraction, which
#    aov(... + Error(run)) would not.
#
# 3. EFFECT SIZE FIRST. With ~800 runs almost everything clears p < 0.05, so
#    the p-value is not the interesting column. Reported instead:
#      partial omega^2  -- variance explained, bias-corrected, and negative
#                          when a term explains less than its degrees of
#                          freedom would by chance
#      dr               -- the spread the term produces on the r scale, from
#                          estimated marginal means. This is the practical
#                          column: a term can be overwhelmingly significant and
#                          still move r by 0.01.
#
# ONLY REPS 2-5 ARE USED. Replicate 1 contains two sets of cache files under
# different seeds, both calling themselves rep 1, because --reps was added to a
# finished grid; see code/sim_filter_results.R. Reps 2-5 are one coherent run.
# ============================================================

library(tidyverse)

here::i_am("code/sim_anova.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
out_dir <- here::here("output")
res_path <- arg_value("--results",
                      file.path(out_dir, "simulation_results_from_ceres_filtered.csv"))

# NOTE ON THE INPUT. These read the September 2026 Ceres results, which used the
# pre-October scoring scheme: 20% of the observed cells held out, so the models
# were fitted at 0.8x the labelled sparsity and the per-cell metrics were scored
# on the held-out cells. Read every sparsity level as 0.8x its label. Results
# generated from now on use the never-observed cells and are fitted at the
# labelled density; they also name the fit statistic r_fit_ rather than r_obs_.
# See SIMULATION_GLOSSARY.md.

RESPONSES <- c("r_addsurf_oat", "r_addsurf_pea", "r_int_oat", "r_int_pea")
GEN <- c("interaction", "environment")                    # the generative model
DAT <- c("n_acc", "sparsity")                             # how much data
ANA <- c("K", "eigen_variance", "fixed_main_effect")      # the analysis model
WHOLE <- c(GEN, DAT, ANA)

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------

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

raw <- read_sim_results(res_path) |>
  dplyr::filter(stringr::str_starts(model, "megalmm"), rep > 1)

stopifnot(all(RESPONSES %in% names(raw)))

runs <- raw |>
  dplyr::mutate(readout = dplyr::if_else(model == "megalmm", "Eta", "U")) |>
  dplyr::select(scenario, rep, dplyr::all_of(WHOLE), readout,
                dplyr::all_of(RESPONSES))

# Every run must have exactly one Eta and one U, or the pairing is not a pairing
chk <- dplyr::count(runs, scenario, rep, dplyr::across(dplyr::all_of(WHOLE)),
                    name = "n")
if (any(chk$n != 2)) {
  stop(sum(chk$n != 2), " run(s) do not have exactly one Eta and one U row. ",
       "Filter the results first (code/sim_filter_results.R).", call. = FALSE)
}
message(nrow(chk), " MegaLMM runs, ", nrow(runs), " rows, reps ",
        paste(sort(unique(runs$rep)), collapse = "/"))

# Fisher z. r is capped just inside +-1 so a degenerate cell cannot give Inf.
fz  <- function(r) atanh(pmax(pmin(r, 1 - 1e-8), -1 + 1e-8))
ifz <- tanh

paired <- runs |>
  tidyr::pivot_longer(dplyr::all_of(RESPONSES), names_to = "response",
                      values_to = "r") |>
  dplyr::mutate(z = fz(r)) |>
  dplyr::select(-r) |>
  tidyr::pivot_wider(names_from = readout, values_from = z) |>
  dplyr::filter(!is.na(Eta), !is.na(U)) |>
  dplyr::mutate(mean_z = (Eta + U) / 2, diff_z = Eta - U)

as_factors <- function(d) {
  d |> dplyr::mutate(dplyr::across(dplyr::all_of(WHOLE), factor))
}

# ------------------------------------------------------------
# The two models
#
# Whole plot: main effects of all seven factors, plus the two interaction
# families asked for -- generative x analysis and data x analysis. NOT
# generative x data, and not analysis x analysis: those are separate questions
# and spending degrees of freedom on them is not free under a fraction.
#
# Subplot: main effects only. Each one is a readout x factor interaction, and
# the intercept is the readout main effect. Going further would be asking for
# three-way interactions.
# ------------------------------------------------------------

two_way <- function(a, b) as.vector(outer(a, b, paste, sep = ":"))

# --all-two-way also fits generative x data.
#
# Leaving it out is defensible -- it is a question about the SIMULATION rather
# than about MegaLMM -- but it is not free: those terms are real and large
# (F = 16.6 on r_addsurf_oat), so omitting them pushes their sum of squares into the
# whole-plot error and inflates it by about 50%. Every test in the table is then
# conservative. Reported both ways so the cost is visible rather than implicit.
all_two_way <- "--all-two-way" %in% args
extra <- if (all_two_way) two_way(GEN, DAT) else character(0)
suffix <- if (all_two_way) "_all2w" else ""

rhs_whole <- paste(c(WHOLE, two_way(GEN, ANA), two_way(DAT, ANA), extra),
                   collapse = " + ")
rhs_sub   <- paste(WHOLE, collapse = " + ")

# ------------------------------------------------------------
# One response
# ------------------------------------------------------------

analyse_one <- function(resp) {
  d <- dplyr::filter(paired, response == resp)

  # r_int does not exist where there is no interaction to recover: the true
  # surface is exactly zero, so the correlation is undefined, not small. Drop
  # the level rather than the rows' worth of information it does not carry.
  none_rows <- dplyr::filter(d, interaction == "none")
  dropped_none <- nrow(none_rows) > 0 && all(is.na(none_rows$mean_z))
  if (dropped_none) d <- dplyr::filter(d, interaction != "none")
  d <- dplyr::filter(d, !is.na(mean_z), !is.na(diff_z))
  d <- as_factors(d) |> dplyr::mutate(interaction = droplevels(interaction))

  fit_w <- stats::lm(stats::as.formula(paste("mean_z ~", rhs_whole)), data = d)
  fit_s <- stats::lm(stats::as.formula(paste("diff_z ~", rhs_sub)), data = d)

  # Type II: the design is unbalanced (D-optimal fraction), so sequential sums
  # of squares would depend on the order terms were typed.
  tab <- function(fit, stratum) {
    a <- car::Anova(fit, type = "II")
    ss_e <- stats::sigma(fit)^2 * fit$df.residual
    df_e <- fit$df.residual
    ms_e <- ss_e / df_e
    tibble::as_tibble(a, rownames = "term") |>
      dplyr::filter(term != "Residuals") |>
      dplyr::rename(SS = `Sum Sq`, df = Df, F = `F value`, p = `Pr(>F)`) |>
      dplyr::mutate(
        stratum = stratum, MS = SS / df,
        # partial omega^2: bias-corrected share of variance, and honest about
        # terms that explain less than chance
        omega2_p = (SS - df * ms_e) / (SS - df * ms_e + (df_e + 1) * ms_e),
        .before = SS)
  }

  out <- dplyr::bind_rows(tab(fit_w, "between-run (whole plot)"),
                          tab(fit_s, "within-run (subplot)"))

  # The readout main effect is the subplot INTERCEPT, which car::Anova omits.
  ci <- stats::confint(fit_s)["(Intercept)", ]
  s_sum <- summary(fit_s)$coefficients["(Intercept)", ]
  out <- dplyr::bind_rows(
    tibble::tibble(term = "readout (Eta vs U)", stratum = "within-run (subplot)",
                   df = 1, MS = NA_real_, omega2_p = NA_real_, SS = NA_real_,
                   F = s_sum[["t value"]]^2, p = s_sum[["Pr(>|t|)"]]),
    out)

  # Practical magnitude: the spread each term produces on the r scale.
  # For a MAIN effect, the spread of its marginal means is what it does.
  #
  # For an INTERACTION, the spread of the cell means is NOT: it contains the two
  # main effects as well, which is why a weak interaction can otherwise print a
  # larger dr than the strong main effect inside it. What the interaction itself
  # does is the cell means with their additive part removed -- exactly the
  # interaction_part() the rest of the project uses -- so that is what is
  # measured, and it is expressed as a shift on the r scale around the grand
  # mean.
  dr <- function(fit, term, data, base) {
    vars <- strsplit(term, ":", fixed = TRUE)[[1]]
    if (!all(vars %in% names(data))) return(NA_real_)
    em <- try(suppressMessages(
      emmeans::emmeans(fit, specs = vars, data = data)), silent = TRUE)
    if (inherits(em, "try-error")) return(NA_real_)
    s <- summary(em)
    if (length(vars) == 1) {
      rng <- range(s$emmean)
    } else {
      M <- tapply(s$emmean, list(s[[vars[1]]], s[[vars[2]]]), identity)
      M <- M - rowMeans(M)[row(M)] - colMeans(M)[col(M)] + mean(M)
      rng <- range(M) + base          # a deviation, shown around the grand mean
    }
    diff(ifz(rng))
  }
  # whole-plot terms are on the mean-z scale, so back-transform directly
  out$dr <- purrr::map2_dbl(out$term, out$stratum, function(tm, st) {
    if (tm == "readout (Eta vs U)") {
      return(abs(ifz(mean(d$Eta)) - ifz(mean(d$U))))
    }
    if (st == "between-run (whole plot)") {
      dr(fit_w, tm, d, mean(d$mean_z))
    } else {
      # A subplot term is a spread of Eta-minus-U DIFFERENCES on the z scale.
      # Shown as the width of the Eta/U gap it implies, on the r scale, taken
      # around the grand mean so the number is comparable with the others.
      base <- mean(d$mean_z)
      g <- function(dz) ifz(base + dz / 2) - ifz(base - dz / 2)
      vars <- strsplit(tm, ":", fixed = TRUE)[[1]]
      em <- try(suppressMessages(
        emmeans::emmeans(fit_s, specs = vars, data = d)), silent = TRUE)
      if (inherits(em, "try-error")) NA_real_ else
        diff(range(g(summary(em)$emmean)))
    }
  })

  # Variance components, from the two error terms:
  #   whole-plot MS_error estimates  sigma2_run + sigma2_sub/2
  #   subplot   MS_error estimates  2 * sigma2_sub
  ms_w <- stats::sigma(fit_w)^2
  ms_s <- stats::sigma(fit_s)^2
  v_sub <- ms_s / 2
  v_run <- ms_w - v_sub / 2

  list(
    response = resp, dropped_none = dropped_none,
    n_runs = nrow(d), df_w = fit_w$df.residual, df_s = fit_s$df.residual,
    table = dplyr::mutate(out, response = resp, .before = 1),
    var_run = v_run, var_sub = v_sub, ms_whole = ms_w, ms_sub = ms_s,
    fit_w = fit_w, fit_s = fit_s, data = d
  )
}

results <- purrr::map(RESPONSES, analyse_one) |> purrr::set_names(RESPONSES)

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

anova_all <- purrr::map(results, "table") |> purrr::list_rbind()

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
stars <- function(p) dplyr::case_when(is.na(p) ~ "", p < 1e-10 ~ "****",
                                      p < 1e-4 ~ "***", p < 0.01 ~ "**",
                                      p < 0.05 ~ "*", TRUE ~ "")

for (resp in RESPONSES) {
  x <- results[[resp]]
  cat("\n", strrep("=", 78), "\n", resp, "  --  ", x$n_runs, " runs",
      if (x$dropped_none) "  (no-interaction scenarios excluded: r_int undefined there)" else "",
      "\n", strrep("=", 78), "\n", sep = "")
  tb <- x$table |>
    dplyr::mutate(sig = stars(p)) |>
    dplyr::arrange(stratum, dplyr::desc(omega2_p)) |>
    dplyr::transmute(stratum, term, df,
                     F = fmt(F, 1), p = format.pval(p, digits = 2, eps = 1e-16),
                     sig, `omega2_p` = fmt(omega2_p), `dr` = fmt(dr))
  print(tb, n = 60, width = 200)
  cat(sprintf("\nerror MS: whole plot %.5f (df %d), subplot %.5f (df %d)\n",
              x$ms_whole, x$df_w, x$ms_sub, x$df_s))
  cat(sprintf("variance components on the z scale: run-to-run %.5f, readout-to-readout %.5f\n",
              x$var_run, x$var_sub))
}

readr::write_csv(anova_all, file.path(out_dir, paste0("simulation_anova", suffix, ".csv")))
message("\nwrote ", file.path(out_dir, paste0("simulation_anova", suffix, ".csv")))

# ------------------------------------------------------------
# What the 4 df of `interaction` actually say
#
# interaction_pct and n_factors are NESTED, not crossed: there is no
# "0% with 5 factors" cell, which is why sim_config.R folds them into one
# composite factor. Once the no-interaction level is gone the remaining four
# levels ARE a 2x2, so size and architecture separate cleanly -- but only
# there. This is the decomposition, on the whole-plot stratum.
# ------------------------------------------------------------

cat("\n", strrep("=", 78), "\ninteraction: size vs architecture, within the 2x2\n",
    strrep("=", 78), "\n", sep = "")
decomp <- purrr::map(RESPONSES, function(resp) {
  # size and architecture read back out of the composite level name, since the
  # underlying columns are not carried through
  d <- results[[resp]]$data |>
    dplyr::mutate(
      int_arch = factor(stringr::str_match(as.character(interaction), "^f(\\d+)")[, 2]),
      int_size = factor(stringr::str_match(as.character(interaction), "_i(\\d+)$")[, 2]))
  if (dplyr::n_distinct(d$int_size) < 2 || dplyr::n_distinct(d$int_arch) < 2)
    return(NULL)
  f <- stats::lm(stats::as.formula(paste(
    "mean_z ~ int_size * int_arch +",
    paste(c(setdiff(WHOLE, "interaction"),
            two_way(GEN[GEN != "interaction"], ANA), two_way(DAT, ANA),
            if (all_two_way) two_way(GEN[GEN != "interaction"], DAT) else NULL),
          collapse = " + "))), data = d)
  a <- car::Anova(f, type = "II")
  ms_e <- stats::sigma(f)^2
  tibble::as_tibble(a, rownames = "term") |>
    dplyr::filter(term %in% c("int_size", "int_arch", "int_size:int_arch")) |>
    dplyr::transmute(response = resp, term, df = Df,
                     F = fmt(`F value`, 1),
                     p = format.pval(`Pr(>F)`, digits = 2, eps = 1e-16),
                     omega2_p = fmt((`Sum Sq` - Df * ms_e) /
                                      (`Sum Sq` - Df * ms_e +
                                         (f$df.residual + 1) * ms_e)))
}) |> purrr::list_rbind()
print(decomp, n = 20)
readr::write_csv(decomp, file.path(out_dir, paste0("simulation_anova_interaction", suffix, ".csv")))
message("wrote ", file.path(out_dir, paste0("simulation_anova_interaction", suffix, ".csv")))

# ------------------------------------------------------------
# Expected mean squares
#
# The split plot is what makes this worth writing down: the two strata have
# DIFFERENT error terms, so a whole-plot factor and a readout x factor
# interaction are not tested against the same thing.
#
# Model for one run i and readout j:
#
#   z_ij = mu + (whole-plot fixed effects)_i + (readout effects)_j
#              + (readout x whole-plot)_ij + u_i + e_ij
#
#   u_i  ~ N(0, s2_run)   run-to-run: one MegaLMM fit versus another under the
#                         same settings, i.e. sampler and simulation noise
#   e_ij ~ N(0, s2_sub)   readout-to-readout within a fit
#
# Analysing the per-run mean and difference gives the two strata exactly:
#
#   mean_z_i = mu + (whole-plot)_i + u_i + e_bar_i    Var = s2_run + s2_sub/2
#   diff_z_i = (readout)   + (readout x whole)_i + (e_i1 - e_i2)   Var = 2 s2_sub
#
# so, with a fixed factor A of a levels and n runs per level:
#
#   WHOLE-PLOT STRATUM (fitted on mean_z)
#     A                  df a-1        E(MS) = s2_run + s2_sub/2 + n/(a-1) * SUM alpha^2
#     whole-plot error   df N-p        E(MS) = s2_run + s2_sub/2
#
#   SUBPLOT STRATUM (fitted on diff_z)
#     readout            df 1          E(MS) = 2 s2_sub + N * delta^2
#     readout x A        df a-1        E(MS) = 2 s2_sub + 2n/(a-1) * SUM (delta.alpha)^2
#     subplot error      df N-p        E(MS) = 2 s2_sub
#
# Both strata are therefore tested against their own residual, which is what the
# two lm() fits do -- and the components come out of the two error mean squares:
#
#     s2_sub = MS_subplot_error / 2
#     s2_run = MS_wholeplot_error - MS_subplot_error / 4
#
# The coefficients above are the balanced-design ones. This design is a
# D-optimal FRACTION and so is not balanced, which is why the tests use Type II
# sums of squares rather than these coefficients; the EMS is here to show what
# is being tested against what, not to compute the F ratios.
# ------------------------------------------------------------

cat("\n", strrep("=", 78), "\nExpected mean squares, and the variance components\n",
    strrep("=", 78), "\n", sep = "")
ems <- purrr::map(RESPONSES, function(resp) {
  x <- results[[resp]]
  tibble::tibble(
    response = resp,
    `MS whole-plot error` = round(x$ms_whole, 5),
    `= s2_run + s2_sub/2` = "",
    `MS subplot error` = round(x$ms_sub, 5),
    `= 2 s2_sub` = "",
    `s2_run` = round(x$var_run, 5),
    `s2_sub` = round(x$var_sub, 5),
    `run share of noise` = paste0(round(100 * x$var_run /
                                          (x$var_run + x$var_sub)), "%"))
}) |> purrr::list_rbind()
print(ems, width = 200)
readr::write_csv(ems, file.path(out_dir, paste0("simulation_anova_ems", suffix, ".csv")))

# ------------------------------------------------------------
# Marginal means: which setting to actually use
# ------------------------------------------------------------

cat("\n", strrep("=", 78), "\nMarginal means on the r scale -- what to do\n",
    strrep("=", 78), "\n", sep = "")
marg <- purrr::map(RESPONSES, function(resp) {
  x <- results[[resp]]
  one <- function(spec) {
    em <- try(suppressMessages(emmeans::emmeans(x$fit_w, specs = spec,
                                                data = x$data)), silent = TRUE)
    if (inherits(em, "try-error")) return(NULL)
    s <- summary(em)
    tibble::tibble(response = resp, factor = paste(spec, collapse = ":"),
                   level = apply(s[spec], 1, paste, collapse = " / "),
                   r = round(tanh(s$emmean), 3))
  }
  dplyr::bind_rows(one("sparsity"), one("K"), one("eigen_variance"),
                   one("fixed_main_effect"), one("environment"),
                   one(c("sparsity", "fixed_main_effect")))
}) |> purrr::list_rbind()
print(dplyr::filter(marg, !stringr::str_detect(factor, ":")) |>
        tidyr::pivot_wider(names_from = response, values_from = r), n = 40)
cat("\nthe one interaction that matters -- sparsity x fixed_main_effect:\n")
print(dplyr::filter(marg, factor == "sparsity:fixed_main_effect") |>
        tidyr::pivot_wider(names_from = response, values_from = r), n = 20)
readr::write_csv(marg, file.path(out_dir, paste0("simulation_anova_means", suffix, ".csv")))
message("\nwrote the EMS and marginal-mean tables")

# ------------------------------------------------------------
# Two cross-tabulations worth having as numbers rather than as effect sizes
#
# partial omega^2 says how much of the explainable variation a term owns; it
# says nothing about what the cells actually look like. For the two terms whose
# mechanism matters most -- the readout x sparsity interaction and the
# architecture x size structure of the simulated interaction -- the cell means
# are the useful output.
#
# These are OBSERVED means on the r scale, not marginal means from the model:
# both layouts are near enough balanced that the two agree, and an observed mean
# needs no explaining.
# ------------------------------------------------------------

cat("\n", strrep("=", 78), "\nEta vs U across sparsity: mean r_addsurf\n",
    strrep("=", 78), "\n", sep = "")
eta_u <- raw |>
  dplyr::mutate(readout = dplyr::if_else(model == "megalmm", "Eta", "U")) |>
  dplyr::group_by(sparsity, readout) |>
  dplyr::summarise(n = dplyr::n(),
                   r_addsurf_oat = mean(r_addsurf_oat, na.rm = TRUE),
                   r_addsurf_pea = mean(r_addsurf_pea, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = readout,
                     values_from = c(n, r_addsurf_oat, r_addsurf_pea)) |>
  dplyr::mutate(gap_oat = r_addsurf_oat_Eta - r_addsurf_oat_U,
                gap_pea = r_addsurf_pea_Eta - r_addsurf_pea_U) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric) & !dplyr::starts_with("n_"),
                              ~ round(.x, 3)))
print(eta_u, width = 200)
cat("\nThe gap WIDENS with data rather than closing: the per-column intercept U\n",
    "discards is the part more data estimates best. U is not a noisy Eta, it is\n",
    "missing a component.\n")
readr::write_csv(eta_u, file.path(out_dir, paste0("simulation_anova_eta_vs_u", suffix, ".csv")))

cat("\n", strrep("=", 78),
    "\nInteraction architecture x size: mean r_int\n", strrep("=", 78), "\n", sep = "")
arch_size <- raw |>
  dplyr::filter(interaction != "none") |>
  dplyr::mutate(
    architecture = paste0(stringr::str_match(interaction, "^f(\\d+)")[, 2], " factor(s)"),
    size = paste0(as.integer(stringr::str_match(interaction, "_i(\\d+)$")[, 2]), "%")) |>
  dplyr::group_by(architecture, size) |>
  dplyr::summarise(n = dplyr::n(),
                   r_int_oat = round(mean(r_int_oat, na.rm = TRUE), 3),
                   r_int_pea = round(mean(r_int_pea, na.rm = TRUE), 3),
                   .groups = "drop")
print(arch_size)
cat("\nmargins:\n")
print(dplyr::group_by(arch_size, architecture) |>
        dplyr::summarise(dplyr::across(c(r_int_oat, r_int_pea),
                                       ~ round(mean(.x), 3)), .groups = "drop"))
print(dplyr::group_by(arch_size, size) |>
        dplyr::summarise(dplyr::across(c(r_int_oat, r_int_pea),
                                       ~ round(mean(.x), 3)), .groups = "drop"))
cat("\nRank costs roughly twice what variance share buys, and the two are\n",
    "additive -- which matters because a programme can influence how much\n",
    "specific combining ability there is far more easily than its rank.\n")
readr::write_csv(arch_size, file.path(out_dir, paste0("simulation_anova_arch_size", suffix, ".csv")))
