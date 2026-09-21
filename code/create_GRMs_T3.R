# ============================================================
# Build genomic relationship matrices (GRMs) for the B4I oat and
# pea accession lists, using T3GenoTools against the T3/Oat
# Breedbase instance.
#
#   data/Acc_B4I_Avena.txt  -> data/GRM_Avena.rds
#   data/Acc_B4I_Pisum.txt  -> data/GRM_Pisum.rds
#
# Both species live on the T3/Oat database (the B4I intercrop
# trials are hosted there), so one connection and one geno_config()
# serve both; the GRMs themselves are built separately, since no
# genotyping protocol covers both species.
#
# Credentials come from .Renviron (T3_USERNAME / T3_PASSWORD).
# The first run downloads and parses archived VCFs and can take
# hours; everything expensive is cached in the shared T3GenoTools
# cache (see T3GenoTools::geno_cache_root()), so later runs are fast.
# ============================================================

library(tidyverse)

here::i_am("code/create_GRMs_T3.R")

# BrAPI.R calls httr::timeout() unqualified, so httr must be attached
library(httr)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

db_name <- "T3/Oat"

acc_files <- c(
  Avena = here::here("data", "Acc_B4I_Avena.txt"),
  Pisum = here::here("data", "Acc_B4I_Pisum.txt")
)

# The GRMs are an INPUT to every model downstream and take hours to rebuild
# from scratch, so they live in data/ and are versioned, unlike the derived
# results in output/ which is gitignored.
out_files <- c(
  Avena = here::here("data", "GRM_Avena.rds"),
  Pisum = here::here("data", "GRM_Pisum.rds")
)

# Set to TRUE to ignore every cached tier and rebuild from scratch
refresh <- FALSE

# Placeholder germplasm names used in T3 to record that one component of
# the intercrop was not sown.  They are not accessions, so they are
# dropped before anything is asked of the server.
monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

# Ignore genotyping protocols covering less than this fraction of a
# species' accessions.  On T3/Oat the B4I oats are essentially all on
# the Oat 3K array, but three multi-GB GBS archives each cover a single
# accession: downloading them costs hours and adds one line to the GRM.
# Set to 0 to use every covering protocol.
min_protocol_coverage <- 0.05

# ------------------------------------------------------------
# Connect to T3/Oat
# ------------------------------------------------------------

# readRenviron() so the credentials are found even when R was started
# elsewhere (e.g. by workflowr::wflow_build())
readRenviron(here::here(".Renviron"))

if (!nzchar(Sys.getenv("T3_USERNAME")) || !nzchar(Sys.getenv("T3_PASSWORD"))) {
  stop("T3_USERNAME / T3_PASSWORD not found: check .Renviron in the project root")
}

conn <- BrAPI::getBrAPIConnection(db_name)
conn$login(
  username = Sys.getenv("T3_USERNAME"),
  password = Sys.getenv("T3_PASSWORD")
)

# ------------------------------------------------------------
# Resolve accession names to germplasmDbIds
#
# build_grm() wants a tibble of (germplasmDbId, germplasmName).
# The /search/germplasm endpoint also matches synonyms, so the name
# that comes back is the primary name -- which is the one the VCF
# samples and the GRM dimnames will carry.
# ------------------------------------------------------------

resolve_accessions <- function(conn, names, batch_size = 250L) {
  names <- names |>
    stringr::str_trim() |>
    unique() |>
    (\(x) x[nzchar(x)])()

  batches <- split(names, ceiling(seq_along(names) / batch_size))

  found <- batches |>
    purrr::map(\(b) {
      res <- conn$search("germplasm", body = list(germplasmNames = as.list(b)))
      res$combined_data
    }) |>
    purrr::list_flatten() |>
    purrr::map(\(g) tibble::tibble(
      germplasmDbId = as.character(g$germplasmDbId),
      germplasmName = as.character(g$germplasmName),
      genus         = as.character(g$genus %||% NA_character_)
    )) |>
    purrr::list_rbind() |>
    dplyr::distinct(germplasmDbId, .keep_all = TRUE)

  # A requested name may have matched as a synonym, so check both the
  # primary names and the requested strings before calling one missing
  missing <- setdiff(names, found$germplasmName)
  if (length(missing) > 0) {
    warning(
      length(missing), " of ", length(names),
      " accession name(s) did not resolve on ", db_name, ": ",
      paste(utils::head(missing, 10), collapse = ", "),
      if (length(missing) > 10) ", ..." else "",
      call. = FALSE
    )
  }

  list(accessions = found, missing = missing)
}

read_accession_list <- function(path) {
  nm <- path |>
    readr::read_lines() |>
    stringr::str_trim()
  nm <- nm[nzchar(nm)]

  dropped <- intersect(nm, monoculture_labels)
  if (length(dropped) > 0) {
    message("  dropping non-accession label(s): ",
            paste(dropped, collapse = ", "))
  }
  setdiff(nm, monoculture_labels)
}

accessions <- acc_files |>
  purrr::map(read_accession_list) |>
  purrr::map(\(nm) resolve_accessions(conn, nm))

purrr::iwalk(accessions, \(a, sp) {
  message(
    sp, ": ", nrow(a$accessions), " accessions resolved, ",
    length(a$missing), " missing"
  )
  print(dplyr::count(a$accessions, genus))
})

# ------------------------------------------------------------
# Genotyping configuration
#
# Identity (db_name) picks the cache partition; everything else is
# tuning and goes into the cache keys.
# ------------------------------------------------------------

cfg <- T3GenoTools::geno_config(
  db_name        = db_name,
  target_density = 10000,  # thin each VCF to ~this many markers
  max_missing    = 0.5,
  min_maf        = 0.01,
  impute         = "mean", # "glmnet" is better but much slower
  panel_min      = 1000,   # estimate each protocol's GRM on >= this many
  progress       = TRUE
)

print(cfg)

# ------------------------------------------------------------
# Build one GRM per species
#
# build_grm() finds every protocol covering the accessions,
# downloads the archived VCFs, builds a standardized GRM per
# protocol, and EM-combines them into one matrix named by
# germplasmName.
# ------------------------------------------------------------

build_species_grm <- function(species, accessions, cfg,
                              min_coverage = 0, refresh = FALSE) {
  message("\n=== ", species, ": building GRM from ", cfg$db_name, " ===")

  train <- dplyr::select(accessions$accessions, germplasmDbId, germplasmName)

  # Which protocols cover these accessions, and how well
  protocols <- T3GenoTools::find_geno_sources(
    conn, unique(train$germplasmDbId), cfg, "protocol", refresh
  )
  print(as.data.frame(protocols))

  use_ids <- protocols |>
    dplyr::filter(n_covered >= min_coverage * nrow(train)) |>
    dplyr::pull(dbId)

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

grms <- accessions |>
  purrr::imap(\(a, sp) build_species_grm(
    sp, a, cfg,
    min_coverage = min_protocol_coverage,
    refresh      = refresh
  ))

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

purrr::iwalk(grms, \(g, sp) {
  G <- g$G
  cat("\n---", sp, "---\n")
  cat("dim:", paste(dim(G), collapse = " x "), "\n")
  cat("requested but absent from G:",
      length(setdiff(accessions[[sp]]$accessions$germplasmName, rownames(G))),
      "\n")
  # build_grm() injects a training accession with neither genotypes nor
  # pedigree as a prior-only line: diagonal 1, zero covariances
  prior_only <- rownames(G)[rowSums(abs(G - diag(diag(G)))) == 0]
  cat("prior-only (ungenotyped) lines:", length(prior_only), "\n")
  if (length(prior_only) > 0) print(utils::head(prior_only, 10))
  cat("diagonal:\n"); print(summary(diag(G)))
  cat("off-diagonal:\n"); print(summary(G[upper.tri(G)]))
  cat("symmetric:", isSymmetric(unname(G)), "\n")
  cat("min eigenvalue:", min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), "\n")
})

# ------------------------------------------------------------
# Save
#
# The whole build_grm() object is saved (G plus the protocols and
# projects it came from) so a downstream analysis can record its
# provenance; G itself is grms[[sp]]$G.
# ------------------------------------------------------------

purrr::iwalk(grms, \(g, sp) {
  saveRDS(g, out_files[[sp]])
  message("wrote ", out_files[[sp]])
})
