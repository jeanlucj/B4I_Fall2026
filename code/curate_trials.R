# ============================================================
# TRIAL-LEVEL QUALITY CONTROL: WHICH TRIALS ARE WORTH FITTING?
#
#   Rscript code/curate_trials.R                  # diagnose every trial
#   Rscript code/curate_trials.R --traits oat_yield,pea_yield,oat_biomass
#   Rscript code/curate_trials.R --min-h2 0.08 --max-cv 0.6
#
# Accession curation asks whether two entries are the same line. Nothing until
# now asked whether a TRIAL's data are worth using at all. B4I_2025_AL is the
# case in point: it returned a full set of 392 plots, so every count-based check
# passes, and it was a crop failure. Oat yield there averages 8.3 against 383 at
# IL, with a coefficient of variation of 0.99 -- the plot-to-plot spread is as
# large as the mean -- and the genotypic variance comes back at exactly zero.
# There is no signal in it to fit.
#
# WHAT IT FITS. One simple mixed model per trial per trait, with no GRM:
#
#   trait ~ (1 | germplasmName) + (1 | intercropGermplasmName) [+ (1 | block)]
#
# Kinship is deliberately absent. The question here is whether the DATA carry
# signal, and borrowing strength through a GRM would partly manufacture the
# answer -- a trial with no information of its own can still look respectable
# once relatives are allowed to speak for it. Each trial is judged alone.
#
# THE BLOCK TERM IS CONDITIONAL, and has to be. Four of the six trials loaded
# today (AL, IA, ND, NY) carry a single block, which is the problem recorded in
# docs/B4I_followups.md section 4. A `(1 | block)` on one level is not a variance
# and lme4 refuses the fit outright, so the term is included only where there are
# at least two levels. A fixed formula here would report "model failed" for four
# good trials and tell you nothing about any of them.
#
# THE DIAGNOSTICS, per trial x trait:
#   n_plots, pct_missing   did the trial return data at all
#   mean, sd, cv           cv = sd/mean. A crop failure shows up here first:
#                          the spread stops being about genotype and starts
#                          being about which plots happened to produce anything
#   rel_mean               the trial's mean over the median across trials for
#                          that trait. Absolute floors cannot work across traits
#                          with different units, so the comparison is internal
#   repeat_geno            var(genotype) / (var(genotype) + var(residual)).
#                          Plot-basis repeatability. The direct question: can
#                          this trial tell two accessions apart?
#   repeat_partner         the same for the intercrop partner
#   singular               lme4 drove a variance component to the boundary,
#                          which for the genotype term means exactly no signal
#
# THE RULE. Each trait trips a flag per threshold it violates. A trait with two
# or more flags, or with zero genotypic variance, FAILS. One flag is REVIEW.
# A trial fails if any of its traits fails. Failing trials are dropped
# downstream; REVIEW trials are kept and flagged, because a judgement that fine
# belongs to a person and not to a threshold.
#
# OVERRIDING IT. The thresholds are a starting point, not an authority. Write
# data/trial_qc_manual.csv with columns studyName, keep and (optionally) note,
# and it wins over everything computed here -- the script reports every override
# loudly rather than quietly honouring it. That file is versioned, because it
# records a human decision about the data and should outlive any single run.
# Then re-run from this step:
#
#   Rscript code/validate_refresh.R --from qc
#
# Needs lme4.
#
# Outputs: output/trial_diagnostics.csv   one row per trial x trait
#          output/trial_qc.csv            one row per trial: status, keep, reason
# ============================================================

library(tidyverse)

here::i_am("code/curate_trials.R")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

# ------------------------------------------------------------
# Configuration
#
# TRIAL_QC_TRAITS is the line to edit when the trials start carrying more than
# the two yields. Anything numeric in output/B4I_intercrop_pheno.rds can go here;
# a trait absent from a given trial is reported as such rather than failing it.
# ------------------------------------------------------------

TRIAL_QC_TRAITS <- strsplit(arg_value("--traits", "oat_yield,pea_yield"),
                            ",")[[1]]

TRIAL_QC <- list(
  min_repeat   = as.numeric(arg_value("--min-h2", "0.05")),
  max_cv       = as.numeric(arg_value("--max-cv", "0.75")),
  min_rel_mean = as.numeric(arg_value("--min-rel-mean", "0.20")),
  min_plots    = as.integer(arg_value("--min-plots", "50")),
  max_missing  = as.numeric(arg_value("--max-missing", "0.25"))
)

pheno_file  <- here::here("output", "B4I_intercrop_pheno.rds")
manual_file <- here::here("data", "trial_qc_manual.csv")
out_dir     <- here::here("output")

if (!file.exists(pheno_file)) {
  stop("no phenotype file at ", pheno_file,
       ". Run code/assemble_B4I_phenotypes.R first.", call. = FALSE)
}
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("curate_trials.R needs the lme4 package", call. = FALSE)
}

pheno <- readRDS(pheno_file)
have <- intersect(TRIAL_QC_TRAITS, names(pheno))
if (length(have) == 0) {
  stop("none of the requested traits are in the phenotype file: ",
       paste(TRIAL_QC_TRAITS, collapse = ", "), call. = FALSE)
}
if (length(have) < length(TRIAL_QC_TRAITS)) {
  message("not in the phenotype file, skipped: ",
          paste(setdiff(TRIAL_QC_TRAITS, have), collapse = ", "))
}

message(dplyr::n_distinct(pheno$studyName), " trial(s), ",
        length(have), " trait(s): ", paste(have, collapse = ", "))

# ------------------------------------------------------------
# One trial, one trait
# ------------------------------------------------------------

diagnose <- function(d, trait) {
  y <- d[[trait]]
  n_plots <- nrow(d)
  n_obs   <- sum(!is.na(y))
  base <- tibble::tibble(
    studyName = d$studyName[1], trait = trait,
    n_plots = n_plots, n_obs = n_obs,
    pct_missing = 1 - n_obs / n_plots,
    mean = mean(y, na.rm = TRUE), sd = stats::sd(y, na.rm = TRUE))

  if (n_obs < 10 || !is.finite(base$mean) || base$mean <= 0) {
    return(dplyr::mutate(base, cv = NA_real_, var_geno = NA_real_,
                         var_partner = NA_real_, var_resid = NA_real_,
                         repeat_geno = NA_real_, repeat_partner = NA_real_,
                         singular = NA, fitted = FALSE))
  }

  dd <- d |>
    dplyr::filter(!is.na(.data[[trait]])) |>
    dplyr::mutate(g  = factor(germplasmName),
                  ic = factor(intercropGermplasmName),
                  bl = factor(blockNumber))

  # The conditional block term; see the header.
  rhs <- c("(1 | g)", "(1 | ic)")
  if (nlevels(dd$bl) > 1) rhs <- c(rhs, "(1 | bl)")
  form <- stats::as.formula(paste(trait, "~", paste(rhs, collapse = " + ")))

  fit <- suppressWarnings(suppressMessages(try(
    lme4::lmer(form, data = dd,
               control = lme4::lmerControl(calc.derivs = FALSE)),
    silent = TRUE)))
  if (inherits(fit, "try-error")) {
    return(dplyr::mutate(base, cv = base$sd / base$mean, var_geno = NA_real_,
                         var_partner = NA_real_, var_resid = NA_real_,
                         repeat_geno = NA_real_, repeat_partner = NA_real_,
                         singular = NA, fitted = FALSE))
  }

  vc <- as.data.frame(lme4::VarCorr(fit))
  v  <- stats::setNames(vc$vcov, vc$grp)
  pick <- function(nm) if (nm %in% names(v)) unname(v[[nm]]) else 0

  vg <- pick("g"); vic <- pick("ic"); vr <- pick("Residual")
  dplyr::mutate(
    base,
    cv = base$sd / base$mean,
    var_geno = vg, var_partner = vic, var_resid = vr,
    repeat_geno    = vg / (vg + vr),
    repeat_partner = vic / (vic + vr),
    singular = lme4::isSingular(fit), fitted = TRUE)
}

diag_tbl <- tidyr::expand_grid(
    studyName = sort(unique(pheno$studyName)), trait = have) |>
  purrr::pmap(function(studyName, trait)
    diagnose(dplyr::filter(pheno, studyName == !!studyName), trait)) |>
  purrr::list_rbind()

# rel_mean compares a trial against its peers, because an absolute floor cannot
# be written for a trait whose units are not known in advance.
diag_tbl <- diag_tbl |>
  dplyr::group_by(trait) |>
  dplyr::mutate(rel_mean = mean / stats::median(mean, na.rm = TRUE)) |>
  dplyr::ungroup()

# ------------------------------------------------------------
# Flags, and the verdict
# ------------------------------------------------------------

flag_of <- function(r) {
  f <- character(0)
  if (!isTRUE(r$fitted))                                f <- c(f, "model_failed")
  if (isTRUE(r$n_plots < TRIAL_QC$min_plots))           f <- c(f, "too_few_plots")
  if (isTRUE(r$pct_missing > TRIAL_QC$max_missing))     f <- c(f, "missing_data")
  if (isTRUE(r$cv > TRIAL_QC$max_cv))                   f <- c(f, "high_cv")
  if (isTRUE(r$rel_mean < TRIAL_QC$min_rel_mean))       f <- c(f, "low_mean")
  if (isTRUE(r$repeat_geno < TRIAL_QC$min_repeat))      f <- c(f, "low_repeatability")
  f
}

diag_tbl <- diag_tbl |>
  dplyr::mutate(
    flags = purrr::map_chr(seq_len(dplyr::n()),
                           \(i) paste(flag_of(diag_tbl[i, ]), collapse = "; ")),
    n_flags = purrr::map_int(seq_len(dplyr::n()),
                             \(i) length(flag_of(diag_tbl[i, ]))),
    # Zero genotypic variance is on its own: it is not one problem among
    # several, it is the statement that the trial cannot distinguish any two
    # accessions. AL's oat yield is exactly this.
    no_signal = isTRUE(fitted) & (var_geno <= 0 | !is.finite(var_geno)),
    trait_status = dplyr::case_when(
      no_signal | n_flags >= 2 ~ "fail",
      n_flags == 1             ~ "review",
      TRUE                     ~ "pass"))

qc <- diag_tbl |>
  dplyr::group_by(studyName) |>
  dplyr::summarise(
    n_fail   = sum(trait_status == "fail"),
    n_review = sum(trait_status == "review"),
    status = dplyr::case_when(n_fail > 0 ~ "fail",
                              n_review > 0 ~ "review",
                              TRUE ~ "pass"),
    reason = paste(unique(unlist(strsplit(
      paste(flags[flags != ""], collapse = "; "), "; "))), collapse = "; "),
    .groups = "drop") |>
  dplyr::mutate(keep = status != "fail", decided_by = "diagnostic")

# ------------------------------------------------------------
# The manual override, which wins
# ------------------------------------------------------------

if (file.exists(manual_file)) {
  man <- readr::read_csv(manual_file, show_col_types = FALSE)
  if (!all(c("studyName", "keep") %in% names(man))) {
    stop(manual_file, " must have columns studyName and keep", call. = FALSE)
  }
  unknown <- setdiff(man$studyName, qc$studyName)
  if (length(unknown)) {
    warning("trial(s) in ", basename(manual_file),
            " that are not in the data, ignored: ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }
  man <- dplyr::filter(man, studyName %in% qc$studyName)
  if (nrow(man)) {
    qc <- qc |>
      dplyr::left_join(dplyr::select(man, studyName, keep_manual = keep,
                                     note = dplyr::any_of("note")),
                       by = "studyName") |>
      dplyr::mutate(
        overridden = !is.na(keep_manual) & keep_manual != keep,
        keep = dplyr::coalesce(keep_manual, keep),
        decided_by = dplyr::if_else(!is.na(keep_manual), "manual", decided_by)) |>
      dplyr::select(-keep_manual)

    for (i in which(qc$overridden)) {
      message("OVERRIDE: ", qc$studyName[i], " diagnostic says ",
              qc$status[i], ", ", basename(manual_file), " says keep = ",
              qc$keep[i],
              if ("note" %in% names(qc) && !is.na(qc$note[i]))
                paste0(" (", qc$note[i], ")") else "")
    }
  }
}

readr::write_csv(diag_tbl, file.path(out_dir, "trial_diagnostics.csv"))
readr::write_csv(qc, file.path(out_dir, "trial_qc.csv"))

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)

cat("\n", strrep("=", 86), "\nPer trial and trait\n", strrep("=", 86), "\n",
    sep = "")
diag_tbl |>
  dplyr::transmute(studyName, trait, n = n_obs,
                   mean = fmt(mean, 1), cv = fmt(cv),
                   rel_mean = fmt(rel_mean),
                   rep_geno = fmt(repeat_geno), rep_ptnr = fmt(repeat_partner),
                   status = trait_status, flags) |>
  print(n = 60, width = 200)

cat("\n", strrep("=", 86), "\nVerdict by trial\n", strrep("=", 86), "\n",
    sep = "")
print(dplyr::select(qc, studyName, status, keep, decided_by, reason),
      n = 40, width = 200)

dropped <- dplyr::filter(qc, !keep)
if (nrow(dropped)) {
  cat("\nDROPPED from everything downstream:\n")
  for (i in seq_len(nrow(dropped))) {
    cat("  ", dropped$studyName[i], " -- ", dropped$reason[i], "\n", sep = "")
  }
} else {
  cat("\nNo trial was dropped.\n")
}
if (any(qc$status == "review")) {
  cat("\nFLAGGED but kept -- look at these before trusting the fit:\n")
  for (s in qc$studyName[qc$status == "review"]) {
    cat("  ", s, " -- ", qc$reason[qc$studyName == s], "\n", sep = "")
  }
}

cat("\nTo overrule any of this, write data/trial_qc_manual.csv with columns\n",
    "studyName,keep[,note] and re-run:  Rscript code/validate_refresh.R --from qc\n",
    sep = "")

message("\nwrote:\n  output/trial_diagnostics.csv\n  output/trial_qc.csv")
