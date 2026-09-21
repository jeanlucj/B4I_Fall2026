# ============================================================
# ASSEMBLE THE PILOT OAT-PEA EXPERIMENT
#
# The pilot that preceded B4I: twelve trials named in
# data/Trial_PilotOatPea.txt, five from 2023 and seven from 2024.
#
# Everything here is the B4I pipeline pointed at a different trial list --
# code/t3_functions.R does the downloading and the GRMs, and the model in
# pilot_dge_ige.R is the same one code/BGLR_multi_trait_model.R fits. Only
# the settings differ, which is the point of keeping the pilot in its own
# directory rather than adding a flag to the B4I scripts.
#
# Outputs: output/pilot/pilot_trials.csv        the trials, and what they hold
#          output/pilot/pilot_pheno.rds/.csv    plot table for the model
#          output/pilot/GRM_pilot_{oat,pea}.rds relationship matrices
# ============================================================

library(tidyverse)

here::i_am("code/pilot/pilot_assemble.R")

library(httr)

source(here::here("code", "curation_functions.R"))   # connect_t3
source(here::here("code", "t3_functions.R"))         # downloads, GRMs

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

db_name    <- "T3/Oat"
trial_file <- here::here("data", "Trial_PilotOatPea.txt")

out_dir   <- here::here("output", "pilot")
cache_dir <- here::here("output", "pilot", "cache")

oat_yield_trait <- "Grain yield - g/m2|CO_350:0000260"
pea_yield_trait <- "Pea Grain Yield - g/m2|CO_xxx:0003008"

page_size <- 10000
refresh   <- FALSE

# Oat 3K and the GenoPea 13K array, as for B4I
protocols <- c(oat = "66", pea = "67")

monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Trials
# ------------------------------------------------------------

conn <- connect_t3(db_name)

wanted <- readr::read_lines(trial_file) |>
  stringr::str_trim() |>
  (\(x) x[nzchar(x)])()

# The pilot is named, not discovered: we know which trials it was, so the
# wizard is only used to turn names into ids.
w <- conn$wizard("trials")
all_trials <- tibble::tibble(
  trialDbId = as.character(w$data$ids),
  trialName = as.character(w$data$names)
)

trials <- dplyr::filter(all_trials, trialName %in% wanted)

missing <- setdiff(wanted, trials$trialName)
if (length(missing) > 0) {
  warning(length(missing), " trial(s) not found on ", db_name, ": ",
          paste(missing, collapse = ", "), call. = FALSE)
}
message("pilot trials found: ", nrow(trials), " of ", length(wanted))

# ------------------------------------------------------------
# Observations and observation units
# ------------------------------------------------------------

observations <- download_observations(conn, trials$trialDbId,
                                      cache_dir = cache_dir,
                                      page_size = page_size, refresh = refresh)

units <- trials$trialDbId |>
  purrr::map(\(id) fetch_observation_units(conn, id, cache_dir),
             .progress = "Observation units") |>
  purrr::list_rbind()

message("observations: ", nrow(observations),
        " | plot-level units: ", nrow(units))

# Which trials actually carry both yields. Four of the twelve have no
# phenotypes uploaded at all, which is worth recording rather than silently
# dropping.
have <- observations |>
  dplyr::filter(observationVariableName %in% c(oat_yield_trait, pea_yield_trait)) |>
  dplyr::distinct(trialDbId, observationVariableName) |>
  dplyr::mutate(trait = dplyr::if_else(observationVariableName == oat_yield_trait,
                                       "oat_yield", "pea_yield"), present = TRUE) |>
  dplyr::select(-observationVariableName) |>
  tidyr::pivot_wider(names_from = trait, values_from = present, values_fill = FALSE)

trial_summary <- trials |>
  dplyr::left_join(have, by = "trialDbId") |>
  dplyr::left_join(dplyr::count(observations, trialDbId, name = "n_obs"),
                   by = "trialDbId") |>
  dplyr::mutate(dplyr::across(c(oat_yield, pea_yield), \(x) tidyr::replace_na(x, FALSE)),
                n_obs = tidyr::replace_na(n_obs, 0L),
                usable = oat_yield & pea_yield)

readr::write_csv(trial_summary, file.path(out_dir, "pilot_trials.csv"))

cat("\n=== Pilot trials ===\n")
print(as.data.frame(trial_summary), row.names = FALSE)

usable <- dplyr::filter(trial_summary, usable)
message("\ntrials with both yields: ", nrow(usable), " of ", nrow(trials))
stopifnot(nrow(usable) > 0)

# ------------------------------------------------------------
# Plot table
# ------------------------------------------------------------

yields <- observations |>
  dplyr::filter(trialDbId %in% usable$trialDbId,
                observationVariableName %in% c(oat_yield_trait, pea_yield_trait),
                !is.na(value), value != "") |>
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
  dplyr::mutate(studyYear = as.integer(stringr::str_extract(studyName, "20\\d{2}")),
                .before = studyName)

message("plots with at least one yield: ", nrow(pheno))

# Half the plots in some pilot trials are monocultures -- no partner recorded.
# They carry no information about mixing ability and are dropped, as in B4I.
n_mono <- sum(is.na(pheno$intercropGermplasmName) |
                pheno$intercropGermplasmName %in% monoculture_labels |
                pheno$germplasmName %in% monoculture_labels)

pheno <- pheno |>
  dplyr::filter(
    !is.na(germplasmName), !is.na(intercropGermplasmName),
    nzchar(intercropGermplasmName),
    !germplasmName %in% monoculture_labels,
    !intercropGermplasmName %in% monoculture_labels
  )

message("dropped ", n_mono, " monoculture or partnerless plot(s); ",
        nrow(pheno), " intercrop plots remain")

cat("\n=== Pilot plot table ===\n")
print(dplyr::count(pheno, studyYear, studyName, name = "plots"))
cat("\noat accessions: ", dplyr::n_distinct(pheno$germplasmName),
    " | pea accessions: ", dplyr::n_distinct(pheno$intercropGermplasmName),
    " | combinations: ",
    dplyr::n_distinct(paste(pheno$germplasmName, pheno$intercropGermplasmName)),
    "\n", sep = "")
cat("\nboth yields present: ",
    sum(!is.na(pheno$oat_yield) & !is.na(pheno$pea_yield)), "\n", sep = "")
print(summary(dplyr::select(pheno, oat_yield, pea_yield)))

saveRDS(pheno, file.path(out_dir, "pilot_pheno.rds"))
readr::write_csv(pheno, file.path(out_dir, "pilot_pheno.csv"))

# ------------------------------------------------------------
# Relationship matrices
# ------------------------------------------------------------

cfg <- T3GenoTools::geno_config(db_name, impute = "mean", progress = FALSE)

resolve <- function(names) {
  res <- conn$search("germplasm", body = list(germplasmNames = as.list(names)))
  res$combined_data |>
    purrr::map(\(g) tibble::tibble(germplasmDbId = as.character(g$germplasmDbId),
                                   germplasmName = as.character(g$germplasmName))) |>
    purrr::list_rbind() |>
    dplyr::distinct(germplasmDbId, .keep_all = TRUE)
}

for (sp in c("oat", "pea")) {
  nm <- if (sp == "oat") unique(pheno$germplasmName) else unique(pheno$intercropGermplasmName)
  acc <- resolve(nm)
  message("\n", sp, ": ", nrow(acc), " of ", length(nm), " accessions resolved")

  grm <- build_species_grm(conn, sp, list(accessions = acc), cfg,
                           protocol_id = protocols[[sp]])

  f <- file.path(out_dir, paste0("GRM_pilot_", sp, ".rds"))
  saveRDS(grm, f)

  G <- grm$G
  prior_only <- rownames(G)[rowSums(abs(G - diag(diag(G)))) == 0]
  cat("\n--- ", sp, " GRM ---\n", sep = "")
  cat("dim: ", paste(dim(G), collapse = " x "), "\n", sep = "")
  cat("ungenotyped (prior-only) lines: ", length(prior_only), "\n", sep = "")
  cat("diagonal:\n"); print(summary(diag(G)))
  message("wrote ", f)
}
