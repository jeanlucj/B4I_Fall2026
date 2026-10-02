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

# ------------------------------------------------------------
# Cache freshness
#
# The per-trial download cache had no invalidation of any kind: once a trial was
# fetched, it was never fetched again. That is wrong whenever a collaborator adds
# data to a trial that has already been downloaded, and the symptom is
# indistinguishable from the trial having no data at all. Measured on
# 2026-10-02: B4I_2026_ND and B4I_2026_NY were reported as lacking oat and pea
# yield, when in fact T3 had both and the September cache predated the upload --
# obs_7002.rds held zero rows and obs_6954.rds held 6,901 rows with no yield
# trait among them.
#
# So a cached file is now used only if it is NEWER than the trial's last
# modification on T3.
# ------------------------------------------------------------

#' Is a cached download still good?
#'
#' Pure, so the rule can be tested without a T3 connection.
#'
#' @param cache_file Path; may not exist.
#' @param modified_at When the trial last changed on T3; `NA`/`NULL` when T3 did
#'   not tell us.
#' @param refresh TRUE forces a re-download.
#' @return TRUE to use the cache.
#'
#' AN UNKNOWN MODIFICATION TIME KEEPS THE CACHE. The alternative -- re-downloading
#' whenever T3 is silent -- would re-fetch all 42 cached trials on every run of a
#' script whose whole point is not to. The caller is expected to say so out loud
#' instead; see the `unknown` branch in find_trials_with_B4I_accessions.R.
cache_is_fresh <- function(cache_file, modified_at = NULL, refresh = FALSE) {
  if (isTRUE(refresh)) return(FALSE)
  if (is.null(cache_file) || !file.exists(cache_file)) return(FALSE)
  if (is.null(modified_at) || length(modified_at) == 0 || is.na(modified_at)) {
    return(TRUE)
  }
  file.mtime(cache_file) > as.POSIXct(modified_at, tz = "UTC")
}

#' When did each trial last change on T3?
#'
#' BrAPI's `GET studies/{id}` carries a `lastUpdate` block in v2, but T3 does not
#' always populate it, so several fields are tried in turn and the one that
#' answered is reported alongside the time. `createDate` is last and is a
#' fallback rather than a synonym: on T3 it does appear to move when a study
#' record is rewritten -- B4I_2025_IL reads 2025-10-07 for a trial sown in
#' spring -- but that is an observation about this database, not a guarantee.
#'
#' `endDate` is deliberately NOT used. It is when the field season ended, which
#' has nothing to do with when the data were uploaded.
#'
#' @return tibble(trialDbId, last_modified, modified_source). `last_modified` is
#'   `NA` when nothing usable came back, which is a reportable state, not an error.
trial_last_modified <- function(conn, trial_ids) {
  parse_stamp <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA)
    x <- as.character(x)[1]
    if (is.na(x) || !nzchar(x)) return(NA)
    for (f in c("%Y-%m-%dT%H:%M:%OSZ", "%Y-%m-%dT%H:%M:%OS%z",
                "%Y-%m-%dT%H:%M:%OS", "%Y-%m-%d %H:%M:%OS", "%Y-%m-%d")) {
      t <- suppressWarnings(as.POSIXct(x, format = f, tz = "UTC"))
      if (!is.na(t)) return(t)
    }
    NA
  }

  one <- function(id) {
    res <- try(conn$get(paste0("studies/", id))$content$result, silent = TRUE)
    if (inherits(res, "try-error") || is.null(res)) {
      return(tibble::tibble(trialDbId = as.character(id),
                            last_modified = as.POSIXct(NA, tz = "UTC"),
                            modified_source = "unavailable"))
    }
    candidates <- list(
      lastUpdate_timestamp = res$lastUpdate$timestamp,
      lastUpdate_date      = res$lastUpdate$date,
      additionalInfo       = res$additionalInfo$lastUpdate,
      createDate           = res$createDate %||% res$create_date
    )
    for (nm in names(candidates)) {
      t <- parse_stamp(candidates[[nm]])
      if (!is.na(t)) {
        return(tibble::tibble(trialDbId = as.character(id),
                              last_modified = t, modified_source = nm))
      }
    }
    tibble::tibble(trialDbId = as.character(id),
                   last_modified = as.POSIXct(NA, tz = "UTC"),
                   modified_source = "none")
  }

  purrr::map(trial_ids, one, .progress = "Trial modification times") |>
    purrr::list_rbind()
}

download_trial_observations <- function(conn, trial_id,
                                        cache_dir = NULL,
                                        page_size = 10000,
                                        refresh = FALSE,
                                        modified_at = NULL) {
  cache_file <- if (!is.null(cache_dir)) {
    file.path(cache_dir, paste0("obs_", trial_id, ".rds"))
  } else {
    NULL
  }

  if (cache_is_fresh(cache_file, modified_at, refresh)) {
    return(readRDS(cache_file))
  }
  if (!is.null(cache_file) && file.exists(cache_file) && !isTRUE(refresh)) {
    message("  trial ", trial_id, ": cached copy predates the trial's last ",
            "change on T3 -- re-downloading")
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

#' @param modified Optional tibble(trialDbId, last_modified) from
#'   `trial_last_modified()`. A trial whose cache is older than its last change
#'   on T3 is re-downloaded; without this argument the cache is trusted, which is
#'   the behaviour that went wrong.
#' @param refresh_ids Trial ids to re-download whatever the timestamps say -- the
#'   manual escape hatch for when T3 reports no modification time.
download_observations <- function(conn, trial_ids,
                                  cache_dir = NULL,
                                  page_size = 10000,
                                  refresh = FALSE,
                                  modified = NULL,
                                  refresh_ids = character(0)) {
  mod_of <- function(id) {
    if (is.null(modified)) return(NULL)
    i <- match(as.character(id), as.character(modified$trialDbId))
    if (is.na(i)) NULL else modified$last_modified[i]
  }
  trial_ids |>
    purrr::map(\(id) download_trial_observations(
                 conn, id, cache_dir, page_size,
                 refresh = refresh || as.character(id) %in% as.character(refresh_ids),
                 modified_at = mod_of(id)),
               .progress = "Downloading observations") |>
    purrr::list_rbind()
}

# A plot's pea partner and its rep/block live on the observation unit, not
# on the observation.  Plot-level units only: the subplot units carry the
# repeated-measure traits and have no yield.
fetch_observation_units <- function(conn, trial_id, cache_dir = NULL,
                                    page_size = 5000, refresh = FALSE,
                                    modified_at = NULL) {
  cache_file <- if (!is.null(cache_dir)) {
    file.path(cache_dir, paste0("units_", trial_id, ".rds"))
  } else NULL

  # Observation UNITS carry the pea partner and the rep/block, so a trial whose
  # plots were re-laid out needs these re-fetched for the same reason the
  # observations do.
  if (cache_is_fresh(cache_file, modified_at, refresh)) {
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
