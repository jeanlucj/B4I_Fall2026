# ============================================================
# SHARED T3 DATA ACCESS
#
# Downloading phenotypes from a Breedbase instance, and the observation units
# that carry a plot's design and its intercrop partner. Used by the B4I
# pipeline and by code/pilot/, which asks the same questions of the earlier
# pilot trials.
#
# This is also the one place the project opens a connection: connect_t3() lives
# here, and every script that talks to T3 sources this file for it.
#
# Two things about Breedbase that these functions exist to encapsulate:
#
#   * BrAPI's default pageSize is 10 records, which turns one trial into more
#     than a thousand requests. Asking for large pages takes a trial from
#     ~165 s to ~4 s.
#   * A plot's intercrop partner is not a second germplasm record. It lives on
#     the observation unit as additionalInfo$intercropGermplasm, with rep and
#     block in observationLevelRelationships, so units have to be fetched
#     separately from observations and joined by plot.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  # BrAPI.R calls httr::timeout() unqualified, so httr must be attached
  library(httr)
})

# ------------------------------------------------------------
# Connection
#
# readRenviron() on the project's own .Renviron, so the credentials are found
# whatever directory R started in -- including inside wflow_build().
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

# A BrAPI record holds NULLs and nested lists; pull one scalar safely
pluck_chr <- function(rec, field) {
  v <- rec[[field]]
  if (is.null(v) || length(v) == 0) NA_character_ else as.character(v)[1]
}

download_trial_observations <- function(conn, trial_id,
                                        cache_dir = NULL,
                                        page_size = 10000,
                                        refresh = FALSE) {
  cache_file <- if (!is.null(cache_dir)) {
    file.path(cache_dir, paste0("obs_", trial_id, ".rds"))
  } else {
    NULL
  }

  if (!is.null(cache_file) && file.exists(cache_file) && !refresh) {
    return(readRDS(cache_file))
  }

  res <- conn$search(
    "observations",
    body     = list(studyDbIds = list(as.character(trial_id))),
    pageSize = page_size
  )

  records <- res$combined_data

  obs <- if (length(records) == 0) {
    tibble::tibble(
      trialDbId = character(0), observationUnitDbId = character(0),
      observationUnitName = character(0), germplasmDbId = character(0),
      germplasmName = character(0), observationVariableDbId = character(0),
      observationVariableName = character(0), value = character(0),
      season = character(0), observationTimeStamp = character(0)
    )
  } else {
    tibble::tibble(
      trialDbId               = as.character(trial_id),
      observationUnitDbId     = purrr::map_chr(records, pluck_chr, "observationUnitDbId"),
      observationUnitName     = purrr::map_chr(records, pluck_chr, "observationUnitName"),
      germplasmDbId           = purrr::map_chr(records, pluck_chr, "germplasmDbId"),
      germplasmName           = purrr::map_chr(records, pluck_chr, "germplasmName"),
      observationVariableDbId = purrr::map_chr(records, pluck_chr, "observationVariableDbId"),
      observationVariableName = purrr::map_chr(records, pluck_chr, "observationVariableName"),
      value                   = purrr::map_chr(records, pluck_chr, "value"),
      season                  = purrr::map_chr(records, \(r) {
                                  s <- r$season
                                  if (is.null(s)) NA_character_
                                  else as.character(s$year %||% s$season %||% NA)
                                }),
      observationTimeStamp    = purrr::map_chr(records, pluck_chr, "observationTimeStamp")
    )
  }

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    saveRDS(obs, cache_file)
  }

  obs
}

download_observations <- function(conn, trial_ids,
                                  cache_dir = NULL,
                                  page_size = 10000,
                                  refresh = FALSE) {
  trial_ids |>
    purrr::map(\(id) download_trial_observations(conn, id, cache_dir, page_size, refresh),
               .progress = "Downloading observations") |>
    purrr::list_rbind()
}

# A plot's pea partner and its rep/block live on the observation unit, not
# on the observation.  Plot-level units only: the subplot units carry the
# repeated-measure traits and have no yield.
fetch_observation_units <- function(conn, trial_id, cache_dir = NULL,
                                    page_size = 5000, refresh = FALSE) {
  cache_file <- if (!is.null(cache_dir)) {
    file.path(cache_dir, paste0("units_", trial_id, ".rds"))
  } else NULL

  if (!is.null(cache_file) && file.exists(cache_file) && !refresh) {
    return(readRDS(cache_file))
  }

  res <- conn$search(
    "observationunits",
    body     = list(studyDbIds = list(as.character(trial_id))),
    pageSize = page_size
  )

  level_code <- function(u, want) {
    rel <- u$observationUnitPosition$observationLevelRelationships
    if (is.null(rel)) return(NA_character_)
    hit <- purrr::keep(rel, \(r) identical(r$levelName, want))
    if (length(hit) == 0) NA_character_ else as.character(hit[[1]]$levelCode)
  }

  units <- res$combined_data |>
    purrr::keep(\(u) identical(u$observationUnitPosition$observationLevel$levelName,
                               "plot")) |>
    purrr::map(\(u) {
      inter <- u$additionalInfo$intercropGermplasm
      tibble::tibble(
        trialDbId              = as.character(u$studyDbId),
        studyName              = as.character(u$studyName %||% NA),
        observationUnitDbId    = as.character(u$observationUnitDbId),
        observationUnitName    = as.character(u$observationUnitName %||% NA),
        germplasmName          = as.character(u$germplasmName %||% NA),
        intercropGermplasmName = if (is.null(inter) || length(inter) == 0) {
          NA_character_
        } else {
          as.character(inter[[1]]$germplasmName)
        },
        repNumber   = level_code(u, "rep"),
        blockNumber = level_code(u, "block")
      )
    }) |>
    purrr::list_rbind() |>
    # Some trials return a plot-level record once per subplot, so the same
    # observationUnitDbId comes back several times
    dplyr::distinct(observationUnitDbId, .keep_all = TRUE)

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    saveRDS(units, cache_file)
  }
  units
}

# ------------------------------------------------------------
# Relationship matrices
#
# Protocols covering a negligible share of the accessions are skipped: on
# T3/Oat three multi-GB GBS archives each cover a single B4I oat, and
# downloading them costs hours to add one line to the GRM.
# ------------------------------------------------------------

build_species_grm <- function(conn, species, accessions, cfg,
                              protocol_id = NULL, min_coverage = 0,
                              refresh = FALSE) {
  message("\n=== ", species, ": building GRM from ", cfg$db_name, " ===")

  train <- dplyr::select(accessions$accessions, germplasmDbId, germplasmName)

  # Which protocols cover these accessions, and how well
  protocols <- T3GenoTools::find_geno_sources(
    conn, unique(train$germplasmDbId), cfg, "protocol", refresh
  )
  print(as.data.frame(protocols))

  use_ids <- if (!is.null(protocol_id)) {
    as.character(protocol_id)
  } else {
    protocols |>
      dplyr::filter(n_covered >= min_coverage * nrow(train)) |>
      dplyr::pull(dbId)
  }

  if (length(use_ids) == 0) {
    stop(species, ": no protocol covers at least ", min_coverage * 100,
         "% of the accessions; lower min_protocol_coverage", call. = FALSE)
  }
  message(species, ": using protocol(s) ", paste(use_ids, collapse = ", "),
          " of ", nrow(protocols), " covering protocol(s)")

  grm <- T3GenoTools::build_grm(
    conn             = conn,
    train_accessions = train,
    cfg              = cfg,
    protocol_id      = use_ids,
    refresh          = refresh
  )

  message(species, ": G is ", nrow(grm$G), " x ", ncol(grm$G),
          " from ", length(grm$protocol_ids), " protocol(s)")

  grm
}
