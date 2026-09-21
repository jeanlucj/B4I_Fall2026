# ============================================================
# FIND T3/OAT FIELD TRIALS THAT EVALUATED THE B4I ACCESSIONS,
# DOWNLOAD THEIR PHENOTYPES, AND SUMMARISE TRAIT AVAILABILITY
#
# Asks T3/Oat which trials evaluated the accessions in
# data/Acc_B4I_Avena.txt and data/Acc_B4I_Pisum.txt, keeps the trials
# carrying at least `min_b4i_accessions` of them, downloads every
# observation from those trials, and builds a trial x trait matrix of
# how many observations each trial holds for each trait.
#
# Selection is in two stages.  A trial must carry at least
# `min_b4i_accessions` B4I accessions, and it must then have
# `selection_trait` recorded on at least that many of them -- a trial
# with no oat grain yield is of no use here however many accessions it
# holds.  Trait availability is only knowable after downloading, so the
# count-qualified trials are all downloaded and then cut down.
#
# Everything a surviving trial measured is kept, not just the traits
# the model needs.
#
# Credentials come from .Renviron (T3_USERNAME / T3_PASSWORD).
# Per-trial downloads are cached under output/trial_cache/, so a re-run
# only fetches trials it has not seen.
#
# Outputs: output/B4I_trial_search.csv        every candidate trial
#          output/B4I_trials_selected.csv     the trials kept, with metadata
#          output/B4I_trait_availability.csv  trial x trait matrix
#          output/B4I_trait_availability.png  the same, as a heatmap
#          data/B4I_observations.rds          all observations, long
#          data/B4I_observations.csv.gz       the same, as text
# ============================================================

library(tidyverse)

here::i_am("code/find_trials_with_B4I_accessions.R")

# BrAPI.R calls httr::timeout() unqualified, so httr must be attached
library(httr)

source(here::here("code", "t3_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

db_name <- "T3/Oat"

acc_files <- c(
  Avena = here::here("data", "Acc_B4I_Avena.txt"),
  Pisum = here::here("data", "Acc_B4I_Pisum.txt")
)

out_dir   <- here::here("output")
cache_dir <- here::here("output", "trial_cache")

# The downloaded observations are an INPUT to everything downstream and are
# expensive to re-fetch, so they live in data/ and are versioned, unlike the
# derived results in output/ which is gitignored.
data_dir  <- here::here("data")

# A trial is kept only if it evaluated at least this many B4I accessions
min_b4i_accessions <- 20

# The traits the intercrop analysis needs.  Presence of each is reported
# per trial.
required_traits <- c(
  grain_yield     = "Grain yield - g/m2|CO_350:0000260",
  pea_grain_yield = "Pea Grain Yield - g/m2|CO_xxx:0003008"
)

# Absolute filter: a trial is kept only if this trait was recorded on at
# least `min_b4i_accessions` of its B4I accessions.  Set to NA to keep
# every trial that passes the accession count.
selection_trait <- required_traits[["grain_yield"]]

# BrAPI's default page size is 10 records, which turns one trial into
# well over a thousand requests.  Asking for large pages takes the same
# download from ~165 s to ~4 s.
page_size <- 10000

# Names that stand for "this half of the intercrop was not sown"
monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

# TRUE re-downloads trials already in the cache
refresh <- FALSE

# ------------------------------------------------------------
# Connection
# ------------------------------------------------------------

connect_t3 <- function(db_name) {
  readRenviron(here::here(".Renviron"))

  if (!nzchar(Sys.getenv("T3_USERNAME")) || !nzchar(Sys.getenv("T3_PASSWORD"))) {
    stop("T3_USERNAME / T3_PASSWORD not found: check .Renviron in the project root",
         call. = FALSE)
  }

  conn <- BrAPI::getBrAPIConnection(db_name)
  conn$login(
    username = Sys.getenv("T3_USERNAME"),
    password = Sys.getenv("T3_PASSWORD")
  )
  conn
}

# ------------------------------------------------------------
# Accessions
# ------------------------------------------------------------

read_b4i_accessions <- function(paths, drop = monoculture_labels) {
  paths |>
    purrr::map(\(p) readr::read_lines(p)) |>
    purrr::list_c() |>
    stringr::str_trim() |>
    unique() |>
    (\(x) x[nzchar(x)])() |>
    setdiff(drop)
}

# ------------------------------------------------------------
# Trial discovery
#
# The Breedbase search wizard maps a set of accessions to every trial
# that evaluated any of them.  It is asked in batches because the filter
# travels in the request.
# ------------------------------------------------------------

candidate_trials <- function(conn, accessions, batch_size = 200L) {
  batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))

  batches |>
    purrr::map(\(b) {
      w <- conn$wizard("trials", filters = list(accessions = b))
      tibble::tibble(
        trialDbId = as.character(w$data$ids),
        trialName = as.character(w$data$names)
      )
    }, .progress = "Trials from accessions") |>
    purrr::list_rbind() |>
    dplyr::distinct(trialDbId, .keep_all = TRUE)
}

# How many of our accessions does each candidate trial actually carry?
# The wizard is asked the other way round, one trial at a time.
count_b4i_accessions <- function(conn, trial_ids, b4i_names) {
  trial_ids |>
    purrr::map(\(id) {
      a <- conn$wizard("accessions", filters = list(trials = as.integer(id)))
      names_in_trial <- as.character(a$data$names)
      tibble::tibble(
        trialDbId      = as.character(id),
        n_accessions   = length(names_in_trial),
        n_b4i          = length(intersect(names_in_trial, b4i_names))
      )
    }, .progress = "B4I accessions per trial") |>
    purrr::list_rbind()
}

select_trials <- function(trial_counts, min_n = min_b4i_accessions) {
  trial_counts |>
    dplyr::filter(n_b4i >= min_n) |>
    dplyr::arrange(dplyr::desc(n_b4i))
}

# Second selection stage: how many B4I accessions does each trial have a
# non-missing value of `trait` for?  Needs the observations, so it runs
# after the download.
count_trait_coverage <- function(observations, trait, accessions) {
  observations |>
    dplyr::filter(
      observationVariableName == trait,
      germplasmName %in% accessions,
      !is.na(value), value != ""
    ) |>
    dplyr::distinct(trialDbId, germplasmName) |>
    dplyr::count(trialDbId, name = "n_b4i_with_trait")
}

trial_metadata <- function(conn, trial_ids) {
  meta <- try(
    T3BrapiHelpers::get_trial_meta_data_from_trial_vec(as.character(trial_ids), conn),
    silent = TRUE
  )
  if (inherits(meta, "try-error")) {
    warning("trial metadata unavailable; continuing without it", call. = FALSE)
    return(tibble::tibble(trialDbId = as.character(trial_ids)))
  }
  meta |>
    dplyr::rename(dplyr::any_of(c(trialDbId = "study_db_id"))) |>
    dplyr::mutate(trialDbId = as.character(trialDbId))
}

# ------------------------------------------------------------
# Observations
# ------------------------------------------------------------

# ------------------------------------------------------------
# Trait availability
#
# Rows are trials, columns are traits, cells are the number of
# observations that trial holds for that trait (0 = not measured).
# ------------------------------------------------------------

# T3 trait names end in the ontology term ("...|CO_350:0000260",
# "...|COMP:0000114").  Drop it for display, which leaves the part that
# distinguishes traits from each other, including "|timepoint 2".
# The full name is kept everywhere the data is written out.
clean_trait_name <- function(x) {
  stringr::str_remove(x, "\\|[A-Za-z0-9_]+:[A-Za-z0-9]+$")
}

trait_availability <- function(observations, trials = NULL,
                               accessions = NULL) {
  obs <- observations
  if (!is.null(accessions)) {
    obs <- dplyr::filter(obs, germplasmName %in% accessions)
  }

  mat <- obs |>
    dplyr::count(trialDbId, observationVariableName, name = "n_obs") |>
    tidyr::pivot_wider(
      names_from  = observationVariableName,
      values_from = n_obs,
      values_fill = 0
    )

  if (!is.null(trials)) {
    mat <- trials |>
      dplyr::select(trialDbId, trialName) |>
      dplyr::left_join(mat, by = "trialDbId") |>
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) tidyr::replace_na(x, 0)))
  }

  # Traits in descending order of how many trials measured them
  trait_cols <- setdiff(names(mat), c("trialDbId", "trialName"))
  ordered <- trait_cols[order(
    colSums(mat[trait_cols] > 0), decreasing = TRUE
  )]

  dplyr::select(mat, dplyr::any_of(c("trialDbId", "trialName")), dplyr::all_of(ordered))
}

# ============================================================
# Driver
# ============================================================

conn <- connect_t3(db_name)

b4i_names <- read_b4i_accessions(acc_files)
message("B4I accessions: ", length(b4i_names))

# --- which trials evaluated any of them ---
candidates <- candidate_trials(conn, b4i_names)
message("candidate trials: ", nrow(candidates))

counts <- count_b4i_accessions(conn, candidates$trialDbId, b4i_names)

trial_search <- candidates |>
  dplyr::left_join(counts, by = "trialDbId") |>
  dplyr::arrange(dplyr::desc(n_b4i))

count_qualified <- select_trials(trial_search, min_b4i_accessions)
message("trials with >= ", min_b4i_accessions, " B4I accessions: ",
        nrow(count_qualified))

stopifnot(nrow(count_qualified) > 0)

# --- download every observation from those trials ---
all_observations <- download_observations(
  conn, count_qualified$trialDbId,
  cache_dir = cache_dir, page_size = page_size, refresh = refresh
)

# --- second stage: require the selection trait ---
if (!is.na(selection_trait)) {
  coverage <- count_trait_coverage(all_observations, selection_trait, b4i_names)

  count_qualified <- count_qualified |>
    dplyr::left_join(coverage, by = "trialDbId") |>
    dplyr::mutate(n_b4i_with_trait = tidyr::replace_na(n_b4i_with_trait, 0L))

  selected <- count_qualified |>
    dplyr::filter(n_b4i_with_trait >= min_b4i_accessions) |>
    dplyr::arrange(dplyr::desc(n_b4i_with_trait))

  dropped <- dplyr::filter(count_qualified, n_b4i_with_trait < min_b4i_accessions)
  message("of those, ", nrow(selected), " have '", selection_trait,
          "' on >= ", min_b4i_accessions, " B4I accessions")
  if (nrow(dropped) > 0) {
    message("dropped for lack of it: ",
            paste(dropped$trialName, collapse = ", "))
  }
} else {
  selected <- count_qualified
}

stopifnot(nrow(selected) > 0)
print(selected, n = Inf)

# Record the outcome of both stages against every candidate
trial_search <- trial_search |>
  dplyr::left_join(
    dplyr::select(count_qualified, trialDbId, dplyr::any_of("n_b4i_with_trait")),
    by = "trialDbId"
  ) |>
  dplyr::mutate(selected = trialDbId %in% selected$trialDbId)

readr::write_csv(trial_search, file.path(out_dir, "B4I_trial_search.csv"))

observations <- selected |>
  dplyr::select(trialDbId, trialName) |>
  dplyr::inner_join(all_observations, by = "trialDbId")

message("observations: ", nrow(observations),
        " on ", dplyr::n_distinct(observations$germplasmName), " accessions, ",
        dplyr::n_distinct(observations$observationVariableName), " traits")

saveRDS(observations, file.path(data_dir, "B4I_observations.rds"))
readr::write_csv(observations, file.path(data_dir, "B4I_observations.csv.gz"))

# --- trait availability ---
availability <- trait_availability(observations, trials = selected)
readr::write_csv(availability, file.path(out_dir, "B4I_trait_availability.csv"))

cat("\n=== Trait availability (observations per trial x trait) ===\n")
print(as.data.frame(availability))

# --- are the traits the intercrop analysis needs actually there? ---
required_present <- purrr::imap(required_traits, \(trait, label) {
  n <- if (trait %in% names(availability)) availability[[trait]] else 0
  tibble::tibble(trialDbId = availability$trialDbId, !!label := n > 0)
}) |>
  purrr::reduce(dplyr::left_join, by = "trialDbId")

trial_obs_counts <- observations |>
  dplyr::count(trialDbId, name = "n_obs") |>
  dplyr::left_join(
    dplyr::count(observations, trialDbId, observationVariableName) |>
      dplyr::count(trialDbId, name = "n_traits"),
    by = "trialDbId"
  )

trials_out <- selected |>
  dplyr::left_join(trial_obs_counts, by = "trialDbId") |>
  dplyr::mutate(
    n_obs    = tidyr::replace_na(n_obs, 0L),
    n_traits = tidyr::replace_na(n_traits, 0L)
  ) |>
  dplyr::left_join(required_present, by = "trialDbId") |>
  dplyr::left_join(trial_metadata(conn, selected$trialDbId), by = "trialDbId")

empty_trials <- dplyr::filter(trials_out, n_obs == 0)
if (nrow(empty_trials) > 0) {
  message("\n", nrow(empty_trials),
          " selected trial(s) have accessions but no observations yet: ",
          paste(empty_trials$trialName, collapse = ", "))
}

readr::write_csv(trials_out, file.path(out_dir, "B4I_trials_selected.csv"))

cat("\n=== Trials carrying the traits the intercrop model needs ===\n")
trials_out |>
  dplyr::select(trialDbId, trialName, n_b4i, n_obs, n_traits,
                dplyr::any_of(names(required_traits))) |>
  as.data.frame() |>
  print()

purrr::iwalk(required_traits, \(trait, label) {
  n <- sum(required_present[[label]], na.rm = TRUE)
  message(label, " (", trait, "): present in ", n, " of ",
          nrow(required_present), " selected trials")
})

# ------------------------------------------------------------
# Heatmap of the availability matrix
# ------------------------------------------------------------

avail_long <- availability |>
  tidyr::pivot_longer(
    -c(trialDbId, trialName),
    names_to = "trait", values_to = "n_obs"
  ) |>
  dplyr::mutate(
    trait     = forcats::fct_inorder(clean_trait_name(trait)),
    trialName = forcats::fct_rev(forcats::fct_inorder(trialName)),
    measured  = n_obs > 0
  )

p_avail <- ggplot2::ggplot(
  avail_long,
  ggplot2::aes(trait, trialName, fill = ifelse(measured, n_obs, NA_integer_))
) +
  ggplot2::geom_tile(colour = "grey90") +
  ggplot2::scale_fill_viridis_c(
    trans = "log10", na.value = "grey95", name = "observations"
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
    panel.grid  = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title    = "Trait availability in T3/Oat trials evaluating B4I accessions",
    subtitle = paste0("trials with at least ", min_b4i_accessions,
                      " B4I accessions carrying ", clean_trait_name(selection_trait),
                      "; grey = not measured"),
    x = NULL, y = NULL
  )

ggplot2::ggsave(
  file.path(out_dir, "B4I_trait_availability.png"), p_avail,
  width = 12, height = 8, dpi = 150
)

message("\nwrote:\n  ",
        paste(c(file.path("output", c(
          "B4I_trial_search.csv", "B4I_trials_selected.csv",
          "B4I_trait_availability.csv", "B4I_trait_availability.png")),
          file.path("data", c("B4I_observations.rds",
                              "B4I_observations.csv.gz"))),
          collapse = "\n  "))
