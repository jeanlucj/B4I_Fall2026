# ============================================================
# ASSEMBLE THE PLOT-LEVEL B4I INTERCROP PHENOTYPE TABLE
#
# Builds the table that code/BGLR_multi_trait_model.R fits: one row per
# intercrop plot, carrying the oat accession, the pea accession, and both
# grain yields.
#
# The pea half of a plot is not a second germplasm record -- Breedbase
# stores it on the plot's observation unit as
# additionalInfo$intercropGermplasm -- so the observation units have to be
# fetched separately from the observations and joined by plot.  Rep and
# block come from the unit's observationLevelRelationships.
#
# Accession names are replaced by the analysis names from
# code/curate_oat_accessions.R and code/curate_pea_accessions.R, so that
# lines which are genetically one genotype enter the model once, with all
# of their plots.  Set `apply_analysis_names` to FALSE to keep the names
# as recorded.
#
# Outputs: output/B4I_intercrop_pheno.rds   the model's input table
#          output/B4I_intercrop_pheno.csv   the same, as text
# ============================================================

library(tidyverse)

here::i_am("code/assemble_B4I_phenotypes.R")

library(httr)

source(here::here("code", "t3_functions.R"))

source(here::here("code", "curation_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

db_name <- "T3/Oat"
out_dir <- here::here("output")

trials_file <- here::here("output", "B4I_trials_selected.csv")
obs_file    <- here::here("data", "B4I_observations.rds")

oat_names_file <- here::here("output", "oat_analysis_names.csv")
pea_names_file <- here::here("output", "pea_analysis_names.csv")

cache_dir <- here::here("output", "trial_cache")

oat_yield_trait <- "Grain yield - g/m2|CO_350:0000260"
pea_yield_trait <- "Pea Grain Yield - g/m2|CO_xxx:0003008"

# Only intercrop trials carry both yields; this keeps the selection honest
# rather than hard-coding a list of trial names.
trial_name_pattern <- "^B4I_"

apply_analysis_names <- TRUE

monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

# ------------------------------------------------------------
# Observation units
# ------------------------------------------------------------

# ------------------------------------------------------------
# Driver
# ------------------------------------------------------------

conn <- connect_t3(db_name)

trials <- readr::read_csv(trials_file, show_col_types = FALSE) |>
  dplyr::mutate(trialDbId = as.character(trialDbId)) |>
  dplyr::filter(
    stringr::str_detect(trialName, trial_name_pattern),
    grain_yield, pea_grain_yield
  )

message("intercrop trials with both yields: ", nrow(trials), " (",
        paste(trials$trialName, collapse = ", "), ")")
stopifnot(nrow(trials) > 0)

units <- trials$trialDbId |>
  purrr::map(\(id) fetch_observation_units(conn, id, cache_dir),
             .progress = "Observation units") |>
  purrr::list_rbind()

message("plot-level observation units: ", nrow(units))

# --- the two yields, from the observations already downloaded ---
observations <- readRDS(obs_file)

yields <- observations |>
  dplyr::filter(
    observationVariableName %in% c(oat_yield_trait, pea_yield_trait),
    !is.na(value), value != ""
  ) |>
  dplyr::mutate(
    trait = dplyr::if_else(observationVariableName == oat_yield_trait,
                           "oat_yield", "pea_yield"),
    value = suppressWarnings(as.numeric(value))
  ) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::select(observationUnitDbId, trait, value) |>
  dplyr::distinct(observationUnitDbId, trait, .keep_all = TRUE) |>
  tidyr::pivot_wider(names_from = trait, values_from = value)

pheno <- units |>
  dplyr::inner_join(yields, by = "observationUnitDbId") |>
  dplyr::mutate(
    studyYear = as.integer(stringr::str_extract(studyName, "(?<=_)\\d{4}(?=_)")),
    .before = studyName
  )

message("plots with at least one yield: ", nrow(pheno))

# --- curation: drop monoculture plots, apply the analysis names ---
pheno <- pheno |>
  dplyr::filter(
    !germplasmName %in% monoculture_labels,
    !intercropGermplasmName %in% monoculture_labels,
    !is.na(germplasmName), !is.na(intercropGermplasmName)
  )

if (apply_analysis_names) {
  read_names <- function(path) {
    if (!file.exists(path)) return(tibble::tibble(germplasmName = character(0),
                                                  analysis_name = character(0)))
    readr::read_csv(path, show_col_types = FALSE) |>
      dplyr::select(germplasmName, analysis_name)
  }

  oat_map <- read_names(oat_names_file)
  pea_map <- read_names(pea_names_file)

  pheno <- pheno |>
    dplyr::left_join(oat_map, by = c("germplasmName" = "germplasmName")) |>
    dplyr::mutate(
      oat_renamed   = !is.na(analysis_name),
      germplasmName = dplyr::coalesce(analysis_name, germplasmName)
    ) |>
    dplyr::select(-analysis_name) |>
    dplyr::left_join(pea_map, by = c("intercropGermplasmName" = "germplasmName")) |>
    dplyr::mutate(
      pea_renamed            = !is.na(analysis_name),
      intercropGermplasmName = dplyr::coalesce(analysis_name, intercropGermplasmName)
    ) |>
    dplyr::select(-analysis_name)

  message("plots whose oat accession was renamed: ", sum(pheno$oat_renamed),
          "; pea: ", sum(pheno$pea_renamed))
}

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

cat("\n=== B4I intercrop phenotype table ===\n")
cat("plots: ", nrow(pheno), "\n", sep = "")
print(dplyr::count(pheno, studyYear, studyName, name = "plots"))

cat("\nboth yields present   : ",
    sum(!is.na(pheno$oat_yield) & !is.na(pheno$pea_yield)), "\n", sep = "")
cat("oat yield only        : ",
    sum(!is.na(pheno$oat_yield) & is.na(pheno$pea_yield)), "\n", sep = "")
cat("pea yield only        : ",
    sum(is.na(pheno$oat_yield) & !is.na(pheno$pea_yield)), "\n", sep = "")

cat("\ndistinct oat accessions: ", dplyr::n_distinct(pheno$germplasmName), "\n", sep = "")
cat("distinct pea accessions: ", dplyr::n_distinct(pheno$intercropGermplasmName), "\n", sep = "")
cat("distinct combinations  : ",
    dplyr::n_distinct(paste(pheno$germplasmName, pheno$intercropGermplasmName)), "\n", sep = "")

cat("\nyields:\n")
print(summary(dplyr::select(pheno, oat_yield, pea_yield)))

saveRDS(pheno, file.path(out_dir, "B4I_intercrop_pheno.rds"))
readr::write_csv(pheno, file.path(out_dir, "B4I_intercrop_pheno.csv"))

message("\nwrote output/B4I_intercrop_pheno.rds and .csv")
