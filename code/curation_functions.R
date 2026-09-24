# ============================================================
# SHARED FUNCTIONS FOR ACCESSION CURATION
#
# Species-agnostic helpers used by code/curate_oat_accessions.R and
# code/curate_pea_accessions.R: read an accession list, confirm
# that one genotyping protocol supplies every marker, fetch dosages,
# correlate marker profiles, group near-identical accessions, and -- where
# pedigrees exist -- summarise full-sib families and resolve the analysis
# names of clonal lines.
#
# connect_t3() is NOT here: it lives in code/t3_functions.R, which every
# script that talks to T3 sources.  One connection function for the project.
#
# Sourced, not run.  Each caller supplies its own settings.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  # BrAPI.R calls httr::timeout() unqualified, so httr must be attached
  library(httr)
})

read_accessions <- function(path,
                            drop = c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")) {
  path |>
    readr::read_lines() |>
    stringr::str_trim() |>
    unique() |>
    (\(x) x[nzchar(x)])() |>
    setdiff(drop)
}

# ------------------------------------------------------------
# Question 1: pedigrees
#
# Breedbase carries the pedigree as a "female/male" string on the
# germplasm record, and the same parents in structured form, with an
# explicit parentType, at /germplasm/{id}/pedigree.  The string is what
# can be fetched in bulk; the structured endpoint is one request per
# accession and is used here only to confirm the string's convention.
# ------------------------------------------------------------

fetch_pedigrees <- function(conn, names, batch_size = 250L) {
  batches <- split(names, ceiling(seq_along(names) / batch_size))

  batches |>
    purrr::map(\(b) {
      res <- conn$search("germplasm", body = list(germplasmNames = as.list(b)))
      res$combined_data
    }, .progress = "Pedigrees") |>
    purrr::list_flatten() |>
    purrr::map(\(g) tibble::tibble(
      germplasmDbId = as.character(g$germplasmDbId),
      germplasmName = as.character(g$germplasmName),
      pedigree      = as.character(g$pedigree %||% NA_character_)
    )) |>
    purrr::list_rbind() |>
    dplyr::distinct(germplasmDbId, .keep_all = TRUE) |>
    dplyr::mutate(
      seed_parent   = dplyr::na_if(
        stringr::str_trim(stringr::str_split_i(pedigree, stringr::fixed("/"), 1)), "NA"),
      pollen_parent = dplyr::na_if(
        stringr::str_trim(stringr::str_split_i(pedigree, stringr::fixed("/"), 2)), "NA"),
      # "NA/NA" is Breedbase for "no pedigree recorded"
      pedigree      = dplyr::if_else(is.na(seed_parent) & is.na(pollen_parent),
                                     NA_character_, pedigree),
      both_parents  = !is.na(seed_parent) & !is.na(pollen_parent)
    )
}

# Confirm that the first element of the pedigree string is the FEMALE
# (seed) parent and the second the MALE (pollen) parent
verify_parent_order <- function(conn, pedigrees, n = 25) {
  if (n <= 0) return(invisible(NULL))

  check <- pedigrees |>
    dplyr::filter(both_parents) |>
    dplyr::slice_sample(n = min(n, sum(pedigrees$both_parents)))

  agree <- check |>
    dplyr::mutate(ok = purrr::map2_lgl(germplasmDbId, seed_parent, \(id, seed) {
      p <- try(conn$get(paste0("/germplasm/", id, "/pedigree"))$content$result$parents,
               silent = TRUE)
      if (inherits(p, "try-error") || length(p) == 0) return(NA)
      female <- purrr::keep(p, \(x) identical(x$parentType, "FEMALE"))
      if (length(female) == 0) return(NA)
      identical(as.character(female[[1]]$germplasmName), seed)
    }, .progress = "Verifying parent order"))

  n_ok <- sum(agree$ok, na.rm = TRUE)
  n_chk <- sum(!is.na(agree$ok))
  message("parent order check: ", n_ok, " of ", n_chk,
          " agree that pedigree = seed_parent/pollen_parent")
  if (n_chk > 0 && n_ok < n_chk) {
    warning("pedigree string order does not always match parentType; ",
            "treat seed/pollen assignment with care", call. = FALSE)
  }
  invisible(agree)
}

# ------------------------------------------------------------
# Question 2: near-identical marker profiles
# ------------------------------------------------------------

marker_dosage <- function(conn, train, cfg, protocol_id) {
  parts <- T3GenoTools::build_partials(
    conn,
    train_accessions = dplyr::select(train, germplasmDbId, germplasmName),
    cfg              = cfg,
    protocol_id      = protocol_id
  )

  dosage <- parts$proto_dosage[[protocol_id]]
  if (is.null(dosage)) {
    stop("no dosage matrix returned for protocol ", protocol_id, call. = FALSE)
  }
  dosage
}

# Correlation of marker dosages between every pair of accessions
identity_matrix <- function(dosage) {
  stats::cor(t(dosage))
}

# Single-linkage groups: accessions joined whenever any pair exceeds the
# threshold.  Single linkage is deliberate -- a chain of near-identical
# accessions is one line however it is labelled.
identity_groups <- function(similarity, threshold) {
  d <- stats::as.dist(1 - similarity)
  cl <- stats::hclust(d, method = "single")
  membership <- stats::cutree(cl, h = 1 - threshold)

  tibble::tibble(
    germplasmName = rownames(similarity),
    group         = membership
  ) |>
    dplyr::add_count(group, name = "group_size") |>
    dplyr::filter(group_size > 1) |>
    # renumber so groups run 1..n in decreasing size
    dplyr::mutate(group = dplyr::dense_rank(dplyr::desc(group_size * 1e6 - group))) |>
    dplyr::arrange(group, germplasmName)
}

# For each group: do the members share a seed parent while naming
# different pollen parents?  That is the self-pollination signature.
classify_groups <- function(groups, pedigrees, min_modal_fraction = 0.5) {
  mode_of <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }

  groups |>
    dplyr::left_join(
      dplyr::select(pedigrees, germplasmName, pedigree, seed_parent, pollen_parent),
      by = "germplasmName"
    ) |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      n_members        = dplyr::n(),
      n_pedigreed      = sum(!is.na(seed_parent)),
      n_seed_parents   = dplyr::n_distinct(seed_parent,   na.rm = TRUE),
      n_pollen_parents = dplyr::n_distinct(pollen_parent, na.rm = TRUE),
      n_pedigrees      = dplyr::n_distinct(pedigree,      na.rm = TRUE),
      modal_seed       = mode_of(seed_parent),
      n_modal_seed     = sum(seed_parent == modal_seed, na.rm = TRUE),
      frac_modal_seed  = dplyr::if_else(n_pedigreed > 0, n_modal_seed / n_pedigreed,
                                        NA_real_),
      # distinct pollen parents recorded against that one female
      n_pollen_on_modal = dplyr::n_distinct(
        pollen_parent[!is.na(seed_parent) & seed_parent == modal_seed], na.rm = TRUE
      ),
      # one female, several recorded males, all genetically identical
      selfing_suspect  = !is.na(frac_modal_seed) &
                         frac_modal_seed >= min_modal_fraction &
                         n_pollen_on_modal > 1,
      same_cross       = n_pedigrees == 1 & n_pedigreed == n_members,
      interpretation   = dplyr::case_when(
        selfing_suspect  ~ paste0("selfed female: ", modal_seed, " with ",
                                  n_pollen_on_modal, " recorded pollen parents"),
        same_cross       ~ "all one cross: indistinguishable sibs or a duplicated line",
        n_pedigreed == 0 ~ "no pedigree recorded",
        TRUE             ~ "mixed pedigrees"
      )
    ) |>
    dplyr::ungroup()
}

# ------------------------------------------------------------
# Full-sib families: accessions sharing BOTH parents
#
# A cross between two inbred oats segregates, so its progeny correlate
# around 0.7-0.8 with each other.  A family whose members all correlate
# above 0.98 never segregated: the "cross" produced one genotype.
# ------------------------------------------------------------

pedigree_families <- function(pedigrees, genotyped, min_size = 2) {
  pedigrees |>
    dplyr::filter(
      germplasmName %in% genotyped,
      !is.na(seed_parent), !is.na(pollen_parent)
    ) |>
    dplyr::mutate(family = paste(seed_parent, pollen_parent, sep = "/")) |>
    dplyr::add_count(family, name = "n_members") |>
    dplyr::filter(n_members >= min_size)
}

# Every within-family pair and its marker correlation
family_pair_correlations <- function(families, similarity) {
  families |>
    dplyr::group_by(family, seed_parent, pollen_parent, n_members) |>
    dplyr::group_modify(\(g, key) {
      m   <- g$germplasmName
      sub <- similarity[m, m, drop = FALSE]
      keep <- upper.tri(sub)
      tibble::tibble(
        a = m[row(sub)[keep]],
        b = m[col(sub)[keep]],
        r = sub[keep]
      )
    }) |>
    dplyr::ungroup()
}

summarise_families <- function(family_pairs, flag) {
  family_pairs |>
    dplyr::group_by(family, seed_parent, pollen_parent, n_members) |>
    dplyr::summarise(
      n_pairs  = dplyr::n(),
      mean_r   = mean(r),
      median_r = stats::median(r),
      min_r    = min(r),
      max_r    = max(r),
      .groups  = "drop"
    ) |>
    dplyr::mutate(
      clonal_family = mean_r > flag,
      verdict = dplyr::if_else(
        clonal_family,
        "did not segregate: one genotype under several names",
        "segregating cross"
      )
    ) |>
    dplyr::arrange(dplyr::desc(mean_r))
}

# ------------------------------------------------------------
# Protocol check and parent genotypes
# ------------------------------------------------------------

# Every marker used here should come from one array.  Stop if a second
# protocol covers enough accessions to be silently mixed in.
confirm_single_protocol <- function(conn, germplasm_db_ids, cfg,
                                    expected_id, max_other_fraction = 0.05) {
  protocols <- T3GenoTools::find_geno_sources(conn, germplasm_db_ids, cfg, "protocol")

  cat("\n=== Genotyping protocols covering these accessions ===\n")
  print(as.data.frame(protocols))

  expected <- dplyr::filter(protocols, dbId == expected_id)
  if (nrow(expected) == 0) {
    stop("protocol ", expected_id, " does not cover these accessions", call. = FALSE)
  }

  others <- dplyr::filter(protocols, dbId != expected_id)
  n_total <- expected$n_total[1]
  if (nrow(others) > 0 && max(others$n_covered) > max_other_fraction * n_total) {
    stop("a second protocol covers more than ", max_other_fraction * 100,
         "% of the accessions; the marker set is not from ", expected_id,
         " alone", call. = FALSE)
  }

  message("confirmed: ", expected$name[1], " (id ", expected_id, ") covers ",
          expected$n_covered[1], " of ", n_total,
          "; no other protocol covers more than ",
          if (nrow(others)) max(others$n_covered) else 0)

  protocols
}

# Every accession named as a parent in the pedigrees
parent_names <- function(pedigrees) {
  c(pedigrees$seed_parent, pedigrees$pollen_parent) |>
    unique() |>
    (\(x) x[!is.na(x)])() |>
    sort()
}

# Resolve the parents to germplasm records so they can be genotyped
fetch_parent_records <- function(conn, parents) {
  res <- conn$search("germplasm", body = list(germplasmNames = as.list(parents)))

  found <- res$combined_data |>
    purrr::map(\(g) tibble::tibble(
      germplasmDbId = as.character(g$germplasmDbId),
      germplasmName = as.character(g$germplasmName)
    )) |>
    purrr::list_rbind() |>
    dplyr::distinct(germplasmDbId, .keep_all = TRUE)

  missing <- setdiff(parents, found$germplasmName)
  if (length(missing) > 0) {
    message(length(missing), " parent(s) not found: ",
            paste(utils::head(missing, 10), collapse = ", "))
  }
  found
}

# Why a parent has no usable genotype.  T3's coverage table can list an
# accession against a protocol while the protocol's archived VCF holds no
# sample of that name, so "covered" and "genotyped" are not the same
# thing and the difference is worth recording.
parent_genotype_status <- function(conn, parent_records, cfg, protocol_id,
                                   genotyped_names) {
  covered <- T3GenoTools::find_geno_sources(
    conn, unique(parent_records$germplasmDbId), cfg, "protocol"
  )
  covered_ids <- if (protocol_id %in% as.character(covered$dbId)) {
    purrr::map(parent_records$germplasmDbId, \(id) {
      p <- T3GenoTools::find_geno_sources(conn, id, cfg, "protocol")
      protocol_id %in% as.character(p$dbId)
    }) |> purrr::list_c()
  } else {
    rep(FALSE, nrow(parent_records))
  }

  vcf_samples <- T3GenoTools::download_protocol_vcfs(conn, protocol_id, cfg) |>
    purrr::map(T3GenoTools::read_vcf_samples) |>
    purrr::list_c() |>
    unique()

  parent_records |>
    dplyr::mutate(
      covered_by_protocol = covered_ids,
      in_vcf_archive      = germplasmName %in% vcf_samples,
      has_dosage          = germplasmName %in% genotyped_names,
      status = dplyr::case_when(
        has_dosage                              ~ "genotyped",
        covered_by_protocol & !in_vcf_archive   ~ "listed against the protocol but absent from its archived VCF",
        in_vcf_archive                          ~ "in the VCF but dropped by QC",
        TRUE                                    ~ "not on this protocol"
      )
    )
}

# Correlation of every accession with its own seed and pollen parent
accession_vs_parents <- function(pedigrees, similarity) {
  genotyped <- rownames(similarity)

  pedigrees |>
    dplyr::filter(germplasmName %in% genotyped) |>
    dplyr::mutate(
      seed_genotyped   = !is.na(seed_parent)   & seed_parent   %in% genotyped,
      pollen_genotyped = !is.na(pollen_parent) & pollen_parent %in% genotyped,
      r_seed_parent = purrr::map2_dbl(germplasmName, seed_parent, \(a, p) {
        if (is.na(p) || !(p %in% genotyped)) NA_real_ else similarity[a, p]
      }),
      r_pollen_parent = purrr::map2_dbl(germplasmName, pollen_parent, \(a, p) {
        if (is.na(p) || !(p %in% genotyped)) NA_real_ else similarity[a, p]
      })
    ) |>
    dplyr::select(germplasmName, pedigree, seed_parent, pollen_parent,
                  seed_genotyped, pollen_genotyped,
                  r_seed_parent, r_pollen_parent)
}

# ------------------------------------------------------------
# Analysis names
#
# Clonal accessions are kept and relabelled, not dropped, so that their
# phenotype records stay in the analysis under one genotype.
# ------------------------------------------------------------

# Two clonal families that share a seed parent may not be two things.  If
# their members correlate with each other as strongly as they do within
# their own family, the pollen parent recorded against them is making no
# genetic difference and they are one selfed line under several cross
# names.  Compare the between-family mean with the weaker of the two
# within-family means, allowing `tolerance` for noise.
pool_clonal_families <- function(clonal, members_by_family, similarity,
                                 tolerance = 0.01) {
  seeds <- clonal |>
    dplyr::count(seed_parent, name = "n_families") |>
    dplyr::filter(n_families > 1, !is.na(seed_parent))

  if (nrow(seeds) == 0) {
    return(tibble::tibble(family = character(0), pool = character(0),
                          between_r = numeric(0), within_r = numeric(0)))
  }

  seeds$seed_parent |>
    purrr::map(\(sp) {
      fams <- dplyr::filter(clonal, seed_parent == sp)
      pairs <- utils::combn(fams$family, 2, simplify = FALSE)

      comparisons <- pairs |>
        purrr::map(\(pr) {
          a <- members_by_family[[pr[1]]]
          b <- members_by_family[[pr[2]]]
          between <- mean(similarity[a, b, drop = FALSE])
          within  <- min(fams$mean_r[fams$family %in% pr])
          tibble::tibble(f1 = pr[1], f2 = pr[2],
                         between_r = between, within_r = within,
                         indistinguishable = between >= within - tolerance)
        }) |>
        purrr::list_rbind()

      merged <- comparisons |> dplyr::filter(indistinguishable)
      if (nrow(merged) == 0) return(NULL)

      tibble::tibble(
        family    = unique(c(merged$f1, merged$f2)),
        pool      = paste0(sp, "_self"),
        between_r = mean(merged$between_r),
        within_r  = min(merged$within_r)
      )
    }) |>
    purrr::list_rbind()
}

resolve_analysis_names <- function(pedigrees, family_summary, similarity,
                                   identity_threshold, family_flag,
                                   pool_tolerance = 0.01) {
  genotyped <- rownames(similarity)

  clonal <- dplyr::filter(family_summary, mean_r > family_flag)

  ped <- pedigrees |>
    dplyr::filter(germplasmName %in% genotyped) |>
    dplyr::mutate(family = paste(seed_parent, pollen_parent, sep = "/"))

  # --- rule 1: members of a family that did not segregate ---
  no_cross <- ped |>
    dplyr::filter(family %in% clonal$family) |>
    dplyr::mutate(
      analysis_name = paste0(stringr::str_replace_all(family, "/", "_"), "_no_cross"),
      reason = paste0("full-sib family ", family, " has mean pairwise r > ",
                      family_flag, ": it did not segregate")
    )

  # --- rule 1b: pool clonal families that share a seed parent and are not
  # distinguishable from each other; they are selfs of that female ---
  members_by_family <- split(no_cross$germplasmName, no_cross$family)

  pooled <- pool_clonal_families(clonal, members_by_family, similarity,
                                 pool_tolerance)

  if (nrow(pooled) > 0) {
    no_cross <- no_cross |>
      dplyr::left_join(pooled, by = "family") |>
      dplyr::mutate(
        reason = dplyr::if_else(
          !is.na(pool),
          paste0("clonal family ", family, " is not distinguishable from the ",
                 "other clonal families of ", seed_parent,
                 " (between-family r = ", round(between_r, 4),
                 " vs within-family r = ", round(within_r, 4),
                 "): one selfed line under several cross names"),
          reason
        ),
        analysis_name = dplyr::coalesce(pool, analysis_name)
      ) |>
      dplyr::select(-pool, -between_r, -within_r)
  }

  # --- rule 2: matches a no_cross line and shares its seed parent ---
  candidates <- setdiff(ped$germplasmName, no_cross$germplasmName)
  seed_of <- stats::setNames(ped$seed_parent, ped$germplasmName)

  selfed <- candidates |>
    purrr::map(\(a) {
      if (is.na(seed_of[[a]])) return(NULL)
      hits <- no_cross$germplasmName[
        similarity[a, no_cross$germplasmName] > identity_threshold &
          !is.na(seed_of[no_cross$germplasmName]) &
          seed_of[no_cross$germplasmName] == seed_of[[a]]
      ]
      if (length(hits) == 0) return(NULL)
      best <- hits[which.max(similarity[a, hits])]
      tibble::tibble(
        germplasmName = a,
        analysis_name = paste0(seed_of[[a]], "_self"),
        reason = paste0("r = ", round(similarity[a, best], 4), " with ", best,
                        " (> ", identity_threshold, ") and shares seed parent ",
                        seed_of[[a]])
      )
    }) |>
    purrr::list_rbind()

  renamed <- dplyr::bind_rows(
    dplyr::select(no_cross, germplasmName, analysis_name, reason),
    selfed
  ) |>
    dplyr::left_join(
      dplyr::select(ped, germplasmName, pedigree, seed_parent, pollen_parent),
      by = "germplasmName"
    )

  # --- rule 3: confirmed identical to the genotyped seed parent ---
  # Overrides both labels above: if the line simply IS the female, call it
  # the female.
  renamed |>
    dplyr::mutate(
      r_seed_parent = purrr::map2_dbl(germplasmName, seed_parent, \(a, p) {
        if (is.na(p) || !(p %in% genotyped)) NA_real_ else similarity[a, p]
      }),
      seed_parent_confirmed = !is.na(r_seed_parent) &
                              r_seed_parent > identity_threshold,
      analysis_name = dplyr::if_else(seed_parent_confirmed, seed_parent, analysis_name),
      reason = dplyr::if_else(
        seed_parent_confirmed,
        paste0("r = ", round(r_seed_parent, 4), " with genotyped seed parent ",
               seed_parent, " (> ", identity_threshold, "): it is that accession"),
        reason
      )
    ) |>
    dplyr::select(germplasmName, analysis_name, reason, pedigree,
                  seed_parent, pollen_parent, r_seed_parent,
                  seed_parent_confirmed) |>
    dplyr::arrange(analysis_name, germplasmName)
}

# ------------------------------------------------------------
# Curation: one representative per group
# ------------------------------------------------------------

choose_representatives <- function(classified, observations_file = NULL) {
  ranked <- classified

  if (!is.null(observations_file) && file.exists(observations_file)) {
    obs_counts <- readRDS(observations_file) |>
      dplyr::filter(!is.na(value), value != "") |>
      dplyr::count(germplasmName, name = "n_observations")
    ranked <- ranked |>
      dplyr::left_join(obs_counts, by = "germplasmName") |>
      dplyr::mutate(n_observations = tidyr::replace_na(n_observations, 0L))
  } else {
    ranked <- dplyr::mutate(ranked, n_observations = 0L)
  }

  ranked |>
    dplyr::group_by(group) |>
    # keep the best-phenotyped member, alphabetical order to break ties
    dplyr::arrange(dplyr::desc(n_observations), germplasmName, .by_group = TRUE) |>
    dplyr::mutate(
      action = dplyr::if_else(dplyr::row_number() == 1, "keep", "drop"),
      representative = dplyr::first(germplasmName)
    ) |>
    dplyr::ungroup()
}
