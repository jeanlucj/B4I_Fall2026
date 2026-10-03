# ============================================================
# WHAT THE EXPERIMENT ACTUALLY CONTAINS
#
#   Rscript code/summarise_experiment.R
#   Rscript code/summarise_experiment.R --no-qc     # ignore the trial screen
#
# The closing inventory of the intercrop trials: how many accessions of each
# species, how many of them are genotyped, how much of the oat x pea grid was
# evaluated, how often each accession appears, and how often a specific
# combination was repeated -- and when it was, whether the repeat was inside one
# trial or across trials.
#
# Everything is reported three times: both years together, 2025 alone, 2026
# alone. The split matters because the two years are not the same experiment --
# the accession list moved between them, so a combination "repeated across
# trials" can mean repeated across locations within a year or repeated across
# years, and those are different kinds of information.
#
# WHICH TRIALS. Whatever output/trial_qc.csv keeps, which is the trial screen
# from code/curate_trials.R as overridden by data/trial_qc_manual.csv. At the
# time of writing that is 8 trials: IA, IL, ND and NY in both years, with
# B4I_2025_AL dropped as a crop failure. `--no-qc` reports on every trial
# instead, which is the right flag if you want to see what the screen removed.
#
# EVERYTHING HERE IS ON ANALYSIS NAMES, which is what the rest of the pipeline
# uses and what makes the counts mean anything.
#
# assemble_B4I_phenotypes.R:121-138 rewrites germplasmName to the analysis name
# when curation found two entries to be one genotype, so the plot table this
# script reads is already keyed that way. The GRM is NOT: data/GRM_Avena.rds
# carries the original names, and `collapse_grm()` averages the rows and columns
# of lines sharing an analysis name into one -- which for a non-segregating
# family is exactly right, since averaging the relationships of identical lines
# is the GRM of their averaged marker profiles.
#
# EVERY production script already does this -- BGLR_multi_trait_model.R:252,
# validate_crossval.R:86, cross_validate_combinations.R:78,
# megalmm_build_inputs.R:253 -- so nothing in the validate_refresh pipeline
# drops a plot for want of a genotype. An earlier version of THIS script read
# the GRM without collapsing, found 8 oat analysis names absent, and reported
# 301 plots as unusable. That was wrong: after the collapse all 462 oat
# accessions are present and none are dropped. tests/test_grm.R pins it so the
# bare read cannot come back.
#
# Outputs: output/summary_accessions.csv      counts and genotyping, per scope
#          output/summary_ungenotyped.csv     the accessions to chase
#          output/summary_sparsity.csv        grid occupancy, per scope
#          output/summary_partners.csv        partners per accession, per scope
#          output/summary_combination_reps.csv  how often each combination ran
#          output/summary_replication.csv     within-trial vs across-trial
#          output/summary_partners.png
#          output/summary_combination_reps.png
# ============================================================

library(tidyverse)

here::i_am("code/summarise_experiment.R")

source(here::here("code", "dge_ige_functions.R"))

args   <- commandArgs(trailingOnly = TRUE)
use_qc <- !("--no-qc" %in% args)

out_dir <- here::here("output")

# Same labels b4i_plot_table() drops: a plot where one component was not sown is
# a monoculture, not an intercrop, and says nothing about mixing.
MONOCULTURE <- c("oat monoculture", "pea monoculture", "no intercrop", "none",
                 "NO_OATS_PLANTED", "NO_PEAS_PLANTED")

pheno <- readRDS(file.path(out_dir, "B4I_intercrop_pheno.rds"))
if (use_qc) pheno <- apply_trial_qc(pheno)

plots <- pheno |>
  dplyr::filter(!germplasmName %in% MONOCULTURE,
                !intercropGermplasmName %in% MONOCULTURE,
                !is.na(germplasmName), !is.na(intercropGermplasmName)) |>
  dplyr::transmute(
    year  = as.integer(studyYear),
    trial = studyName,
    oat   = as.character(germplasmName),
    pea   = as.character(intercropGermplasmName),
    combo = paste(oat, pea, sep = " :: "))

message(nrow(plots), " intercrop plot(s) from ",
        dplyr::n_distinct(plots$trial), " trial(s): ",
        paste(sort(unique(plots$trial)), collapse = ", "))

name_files <- list(oat = here::here("output", "oat_analysis_names.csv"),
                   pea = here::here("output", "pea_analysis_names.csv"))

# collapse_grm(), not bare read_grm(). See the header.
grms <- list(
  oat = collapse_grm(read_grm(here::here("data", "GRM_Avena.rds")),
                     name_files$oat),
  pea = collapse_grm(read_grm(here::here("data", "GRM_Pisum.rds")),
                     name_files$pea))

# analysis_name -> the original entries it stands for, from the curation step
name_maps <- purrr::map(name_files, \(f) if (file.exists(f))
  readr::read_csv(f, show_col_types = FALSE) else NULL)

#' Why is this accession not in the GRM?
#'
#' Returns one row per missing accession with the diagnosis and the evidence.
diagnose_missing <- function(missing, sp) {
  gr <- rownames(grms[[sp]])
  map <- name_maps[[sp]]
  purrr::map(missing, function(a) {
    members <- if (is.null(map)) character(0) else
      map$germplasmName[map$analysis_name == a]
    n_mem <- length(members)
    n_gt  <- sum(members %in% gr)
    tibble::tibble(
      accession = a, n_members = n_mem, n_members_genotyped = n_gt,
      diagnosis = dplyr::case_when(
        n_mem > 0 && n_gt > 0 ~ "collapsed name present in the map but absent after collapse -- a bug, not a genotyping gap",
        n_mem > 0             ~ "collapsed name; no member genotyped -- needs genotyping",
        TRUE                  ~ "not in the GRM and not a collapsed name -- needs genotyping"),
      members = paste(utils::head(members, 20), collapse = "; "))
  }) |> purrr::list_rbind()
}

# The three scopes, as a list of plot tables.
scopes <- list(
  `both years` = plots,
  `2025`       = dplyr::filter(plots, year == 2025L),
  `2026`       = dplyr::filter(plots, year == 2026L))

# ------------------------------------------------------------
# 1. Accessions, and whether they are genotyped
# ------------------------------------------------------------

accessions_of <- function(d, sp) {
  col <- if (sp == "oat") d$oat else d$pea
  sort(unique(col))
}

acc_tbl <- purrr::imap(scopes, function(d, scope) {
  purrr::map(c("oat", "pea"), function(sp) {
    a  <- accessions_of(d, sp)
    gt <- a %in% rownames(grms[[sp]])
    tibble::tibble(scope = scope, species = sp,
                   n_accessions = length(a),
                   n_genotyped = sum(gt), n_ungenotyped = sum(!gt),
                   pct_genotyped = 100 * mean(gt))
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()

# The list to act on. Reported with how much data each one carries, because an
# ungenotyped accession in forty plots is a more urgent request than one in two.
ungenotyped <- purrr::imap(scopes, function(d, scope) {
  purrr::map(c("oat", "pea"), function(sp) {
    a  <- accessions_of(d, sp)
    missing <- setdiff(a, rownames(grms[[sp]]))
    if (length(missing) == 0) return(NULL)
    self <- if (sp == "oat") "oat" else "pea"
    other <- if (sp == "oat") "pea" else "oat"
    d |>
      dplyr::filter(.data[[self]] %in% missing) |>
      dplyr::group_by(accession = .data[[self]]) |>
      dplyr::summarise(n_plots = dplyr::n(),
                       n_partners = dplyr::n_distinct(.data[[other]]),
                       trials = paste(sort(unique(trial)), collapse = "; "),
                       .groups = "drop") |>
      dplyr::left_join(diagnose_missing(missing, sp), by = "accession") |>
      dplyr::mutate(scope = scope, species = sp, .before = 1)
  }) |> purrr::compact() |> purrr::list_rbind()
}) |> purrr::list_rbind()

# ------------------------------------------------------------
# 2. Sparsity: how much of the grid was evaluated
# ------------------------------------------------------------

sparsity_tbl <- purrr::imap(scopes, function(d, scope) {
  n_oat <- dplyr::n_distinct(d$oat); n_pea <- dplyr::n_distinct(d$pea)
  n_combo <- dplyr::n_distinct(d$combo)
  tibble::tibble(
    scope = scope, n_oat = n_oat, n_pea = n_pea,
    possible_combinations = n_oat * n_pea,
    observed_combinations = n_combo,
    pct_observed = 100 * n_combo / (n_oat * n_pea),
    n_plots = nrow(d),
    plots_per_combination = nrow(d) / n_combo)
}) |> purrr::list_rbind()

# ------------------------------------------------------------
# 3. Partners per accession
#
# "the number of combinations each accession was evaluated in" = the number of
# DISTINCT partners it was grown with. Counting plots instead would mix in
# replication, which is the next question rather than this one.
# ------------------------------------------------------------

partners <- purrr::imap(scopes, function(d, scope) {
  purrr::map(c("oat", "pea"), function(sp) {
    self <- sp; other <- if (sp == "oat") "pea" else "oat"
    d |>
      dplyr::group_by(accession = .data[[self]]) |>
      dplyr::summarise(n_partners = dplyr::n_distinct(.data[[other]]),
                       n_plots = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(scope = scope, species = sp, .before = 1)
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()

partner_summary <- partners |>
  dplyr::group_by(scope, species) |>
  dplyr::summarise(n_accessions = dplyr::n(),
                   min = min(n_partners),
                   q25 = stats::quantile(n_partners, 0.25),
                   median = stats::median(n_partners),
                   mean = mean(n_partners),
                   q75 = stats::quantile(n_partners, 0.75),
                   max = max(n_partners), .groups = "drop")

# ------------------------------------------------------------
# 3b. Analysis names that stand for more than one entry
#
# These behave as one genotype in every model, which is correct -- curation
# found them not to segregate. But they are not one SEED LOT. If one of them is
# selected into a validation pool, somebody has to decide which of its member
# entries to actually sow, source that seed, and then name the line definitively.
# This is that worklist, produced now rather than discovered at planting.
# ------------------------------------------------------------

collapsed <- purrr::imap(name_files, function(f, sp) {
  if (!file.exists(f)) return(NULL)
  map <- readr::read_csv(f, show_col_types = FALSE)
  used <- accessions_of(scopes[["both years"]], sp)
  map |>
    dplyr::filter(analysis_name %in% used) |>
    dplyr::group_by(analysis_name) |>
    dplyr::summarise(n_entries = dplyr::n(),
                     entries = paste(sort(germplasmName), collapse = "; "),
                     reasons = paste(sort(unique(reason)), collapse = "; "),
                     .groups = "drop") |>
    dplyr::filter(n_entries > 1) |>
    dplyr::mutate(species = sp, .before = 1)
}) |> purrr::compact() |> purrr::list_rbind()

if (nrow(collapsed) > 0) {
  use <- dplyr::bind_rows(
    partners |>
      dplyr::filter(scope == "both years") |>
      dplyr::select(species, analysis_name = accession, n_plots, n_partners))
  collapsed <- dplyr::left_join(collapsed, use, by = c("species", "analysis_name")) |>
    dplyr::arrange(species, dplyr::desc(n_plots))
}

# ------------------------------------------------------------
# 4. How often was a specific combination evaluated, and where
#
# A combination repeated inside ONE trial is replication: the same pairing in
# the same environment, which estimates error. A combination repeated ACROSS
# trials is something else -- the same pairing in different environments, which
# is what the interaction term needs. They are counted apart because they buy
# different things.
# ------------------------------------------------------------

combos <- purrr::imap(scopes, function(d, scope) {
  d |>
    dplyr::group_by(combo) |>
    dplyr::summarise(
      n_plots  = dplyr::n(),
      n_trials = dplyr::n_distinct(trial),
      max_in_one_trial = max(table(trial)),
      .groups = "drop") |>
    dplyr::mutate(scope = scope, .before = 1)
}) |> purrr::list_rbind()

combo_reps <- combos |>
  dplyr::count(scope, n_plots, name = "n_combinations") |>
  dplyr::arrange(scope, n_plots)

replication <- combos |>
  dplyr::mutate(kind = dplyr::case_when(
    n_plots == 1                             ~ "evaluated once",
    n_trials == 1                            ~ "repeated within one trial",
    max_in_one_trial == 1                    ~ "repeated across trials only",
    TRUE                                     ~ "repeated both ways")) |>
  dplyr::count(scope, kind, name = "n_combinations") |>
  dplyr::group_by(scope) |>
  dplyr::mutate(pct = 100 * n_combinations / sum(n_combinations)) |>
  dplyr::ungroup()

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

readr::write_csv(acc_tbl,        file.path(out_dir, "summary_accessions.csv"))
readr::write_csv(ungenotyped,    file.path(out_dir, "summary_ungenotyped.csv"))
readr::write_csv(sparsity_tbl,   file.path(out_dir, "summary_sparsity.csv"))
readr::write_csv(partners,       file.path(out_dir, "summary_partners.csv"))
readr::write_csv(combo_reps,     file.path(out_dir, "summary_combination_reps.csv"))
readr::write_csv(replication,    file.path(out_dir, "summary_replication.csv"))
if (nrow(collapsed) > 0) {
  readr::write_csv(collapsed, file.path(out_dir, "summary_collapsed_names.csv"))
}

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

rule <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n",
                        strrep("=", 78), "\n", sep = "")

rule("1. Accessions and genotyping")
print(as.data.frame(acc_tbl), row.names = FALSE, digits = 4)

if (nrow(ungenotyped) > 0) {
  cat("\n--- accessions NOT IN THE GRM, both years ---\n")
  ug <- ungenotyped |>
    dplyr::filter(scope == "both years") |>
    dplyr::arrange(species, dplyr::desc(n_plots))
  print(as.data.frame(dplyr::select(ug, species, accession, n_plots,
                                    n_partners, n_members,
                                    n_members_genotyped, diagnosis)),
        row.names = FALSE)

  need_lab <- dplyr::filter(ug, grepl("needs genotyping", diagnosis))
  need_grm <- dplyr::filter(ug, grepl("rebuild the GRM", diagnosis))
  cat("\n", nrow(ug), " accession(s) cannot enter a kinship model, in ",
      sum(ug$n_plots), " plot(s) carrying ", sum(ug$n_partners),
      " combination(s).\n", sep = "")
  cat("  needing genotyping (a request to the lab): ", nrow(need_lab), "\n",
      sep = "")
  if (nrow(need_lab)) {
    cat("    ", paste(need_lab$accession, collapse = ", "), "\n", sep = "")
  }
  cat("  collapsed names whose members are already genotyped: ", nrow(need_grm),
      "\n", sep = "")
  if (nrow(need_grm)) {
    cat("    these need code/create_GRMs_T3.R re-run so the GRM carries the\n",
        "    analysis names, NOT new marker data:\n      ",
        paste(need_grm$accession, collapse = ", "), "\n", sep = "")
  }
} else {
  cat("\nEvery accession in the kept trials is in the GRM.\n")
}

rule("1b. Analysis names standing for more than one entry")
if (nrow(collapsed) > 0) {
  cat("These are ONE genotype in every model -- curation found them not to\n",
      "segregate -- but they are not one seed lot. If one is selected into a\n",
      "validation pool, its member entries have to be resolved, seed sourced,\n",
      "and the line named definitively.\n\n", sep = "")
  # `reasons` is a paragraph per row; it stays in the CSV and out of the console
  print(as.data.frame(dplyr::select(collapsed, species, analysis_name,
                                    n_entries, n_plots, n_partners)),
        row.names = FALSE)
  cat("\n", nrow(collapsed), " collapsed name(s) covering ",
      sum(collapsed$n_entries), " original entries, ",
      sum(collapsed$n_plots, na.rm = TRUE), " plots.\n", sep = "")
  cat("Member entries are in output/summary_collapsed_names.csv.\n")
} else {
  cat("No analysis name stands for more than one entry.\n")
}

rule("2. Grid occupancy")
print(as.data.frame(sparsity_tbl), row.names = FALSE, digits = 4)

rule("3. Partners per accession")
print(as.data.frame(partner_summary), row.names = FALSE, digits = 3)

rule("4. How often each specific combination was evaluated")
print(combo_reps |>
        tidyr::pivot_wider(names_from = scope, values_from = n_combinations,
                           values_fill = 0) |>
        as.data.frame(), row.names = FALSE)

rule("5. Replication: within a trial, or across trials?")
print(replication |>
        dplyr::mutate(pct = round(pct, 1)) |>
        tidyr::pivot_wider(names_from = scope,
                           values_from = c(n_combinations, pct),
                           values_fill = 0) |>
        as.data.frame(), row.names = FALSE)

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

SURFACE <- "#fcfcfb"; GRID <- "#e3e2de"
FILL    <- "#9ec5f4"; EDGE <- "#2a78d6"
INK     <- "#0b0b0b"; INK_SOFT <- "#52514e"

base_theme <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = SURFACE, colour = NA),
    panel.background = ggplot2::element_rect(fill = SURFACE, colour = NA),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
    plot.title       = ggplot2::element_text(colour = INK, face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(colour = INK_SOFT, size = 9,
                                             margin = ggplot2::margin(b = 10)),
    axis.title       = ggplot2::element_text(colour = INK_SOFT, size = 9),
    axis.text        = ggplot2::element_text(colour = INK_SOFT, size = 8),
    strip.text       = ggplot2::element_text(colour = INK, size = 9, face = "bold"),
    legend.position  = "none")

# The x axis is clipped. A handful of accessions sit far out in the tail -- all
# of them curation collapses, where one analysis name stands for a whole
# non-segregating family and therefore inherits every partner its members had --
# and on a shared axis they stretch the scale until the bulk of the distribution
# occupies a fifth of the panel. The clipped count is stated rather than hidden.
X_MAX <- 25L
beyond <- partners |>
  dplyr::filter(n_partners > X_MAX) |>
  dplyr::count(scope, species, name = "n")
beyond_txt <- if (nrow(beyond) == 0) "none" else
  paste(sprintf("%s %s: %d", beyond$scope, beyond$species, beyond$n),
        collapse = "; ")

p_partners <- partners |>
  dplyr::mutate(scope = factor(scope, c("both years", "2025", "2026"))) |>
  ggplot2::ggplot(ggplot2::aes(n_partners)) +
  ggplot2::geom_histogram(binwidth = 1, fill = FILL, colour = EDGE,
                          linewidth = 0.2) +
  ggplot2::facet_grid(species ~ scope, scales = "free_y") +
  ggplot2::coord_cartesian(xlim = c(0, X_MAX)) +
  ggplot2::labs(
    title = "Partners per accession",
    subtitle = paste0(
      "How many distinct partners of the other species each accession was grown with.\n",
      "Axis clipped at ", X_MAX, "; beyond it (", beyond_txt,
      ") are curation collapses, where one analysis\nname stands for a whole ",
      "non-segregating family and inherits every partner its members had."),
    x = "distinct partners", y = "accessions") +
  base_theme

ggplot2::ggsave(file.path(out_dir, "summary_partners.png"), p_partners,
                width = 8.4, height = 5.0, dpi = 300, bg = SURFACE)

# The once-only bar is removed rather than the axis rescaled: it is an order of
# magnitude taller than everything else and flattens the rest into the baseline.
# Its count is in the subtitle so nothing is hidden, only moved.
once <- combos |> dplyr::filter(n_plots == 1) |> dplyr::count(scope)
once_txt <- paste(sprintf("%s: %s", once$scope, format(once$n, big.mark = ",")),
                  collapse = "   ")

p_reps <- combos |>
  dplyr::filter(n_plots > 1) |>
  dplyr::mutate(scope = factor(scope, c("both years", "2025", "2026"))) |>
  ggplot2::ggplot(ggplot2::aes(factor(n_plots))) +
  ggplot2::geom_bar(fill = FILL, colour = EDGE, linewidth = 0.2, width = 0.8) +
  ggplot2::facet_wrap(~ scope, nrow = 1, scales = "free_y") +
  ggplot2::labs(
    title = "Replication of specific oat × pea combinations",
    subtitle = paste0(
      "Combinations evaluated more than once. The once-only bar is excluded ",
      "because it dwarfs the rest —\nevaluated exactly once: ", once_txt, "."),
    x = "plots in which the combination was evaluated", y = "combinations") +
  base_theme

ggplot2::ggsave(file.path(out_dir, "summary_combination_reps.png"), p_reps,
                width = 8.4, height = 3.8, dpi = 300, bg = SURFACE)

message("\nwrote:\n  output/summary_accessions.csv\n  output/summary_ungenotyped.csv\n",
        "  output/summary_sparsity.csv\n  output/summary_partners.csv\n",
        "  output/summary_combination_reps.csv\n  output/summary_replication.csv\n",
        "  output/summary_collapsed_names.csv\n",
        "  output/summary_partners.png\n  output/summary_combination_reps.png")
