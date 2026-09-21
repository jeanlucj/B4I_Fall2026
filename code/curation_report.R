# ============================================================
# WRITE THE CURATION REPORT
#
# Turns the outputs of code/curate_oat_accessions.R and
# code/curate_pea_accessions.R into one readable account of what was
# collapsed and on what grounds: output/CURATION.md.
#
# Generated rather than written by hand, and generated from the settings each
# curation run recorded for itself rather than from the current contents of
# the scripts, so the report cannot drift away from the results it describes.
# If a threshold changes, re-run that species' curation and then this.
#
# Outputs: output/CURATION.md
# ============================================================

library(tidyverse)

here::i_am("code/curation_report.R")

out_dir  <- here::here("output")
out_file <- file.path(out_dir, "CURATION.md")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

read_if <- function(path, reader = readr::read_csv) {
  if (!file.exists(path)) return(NULL)
  if (identical(reader, readr::read_csv)) reader(path, show_col_types = FALSE) else reader(path)
}

# Markdown table from a data frame, with no dependency on knitr
md_table <- function(df, digits = 4) {
  if (is.null(df) || nrow(df) == 0) return("_none_\n")
  df <- df |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, digits))) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) tidyr::replace_na(x, "")))
  paste0(
    "| ", paste(names(df), collapse = " | "), " |\n",
    "|", paste(rep("---", ncol(df)), collapse = "|"), "|\n",
    paste(apply(df, 1, \(r) paste0("| ", paste(r, collapse = " | "), " |")),
          collapse = "\n"),
    "\n"
  )
}

species_block <- function(sp, label) {
  settings <- read_if(file.path(out_dir, paste0(sp, "_curation_settings.rds")), readRDS)
  names_tbl <- read_if(file.path(out_dir, paste0(sp, "_analysis_names.csv")))
  groups    <- read_if(file.path(out_dir, paste0(sp, "_identity_groups.csv")))
  families  <- read_if(file.path(out_dir, paste0(sp, "_family_correlations.csv")))

  if (is.null(settings)) {
    return(paste0("## ", label, "\n\n_No curation settings found; run ",
                  "`code/curate_", sp, "_accessions.R` first._\n"))
  }

  n_renamed <- if (is.null(names_tbl)) 0 else nrow(names_tbl)
  n_names   <- if (is.null(names_tbl)) 0 else dplyr::n_distinct(names_tbl$analysis_name)

  # Count what is left by applying the mapping, not by arithmetic on the
  # totals. A pea collapsed onto ND VICTORY joins a name that already exists,
  # so it removes a line; an oat collapsed onto IL18-8562_self joins a name
  # that does not, so it removes several and adds one. Only the mapping knows
  # which, and getting this wrong overstates the panel.
  n_net <- if (is.null(settings$genotyped_names)) {
    NA_integer_
  } else {
    mapped <- settings$genotyped_names
    if (n_renamed > 0) {
      hit <- match(mapped, names_tbl$germplasmName)
      mapped[!is.na(hit)] <- names_tbl$analysis_name[hit[!is.na(hit)]]
    }
    dplyr::n_distinct(mapped)
  }

  out <- paste0(
    "## ", label, "\n\n",
    "Curated ", format(settings$run_at, "%Y-%m-%d %H:%M"),
    " against genotyping protocol ", settings$protocol_id, ".\n\n",
    "- accessions requested: **", settings$n_requested, "**\n",
    "- with usable marker data: **", settings$n_genotyped, "**\n",
    "- collapsed into another name: **", n_renamed, "**, becoming **",
    n_names, "** distinct ", if (n_names == 1) "line" else "lines", "\n",
    "- distinct lines entering an analysis: **", n_net, "**\n\n"
  )

  if (length(settings$ungenotyped) > 0) {
    out <- paste0(out,
      "**Not checkable.** ", length(settings$ungenotyped), " accession(s) have ",
      "no marker data and could not be compared with anything: ",
      paste0("`", settings$ungenotyped, "`", collapse = ", "), ". They enter ",
      "analyses under their own names with a prior-only row in the GRM.\n\n")
  }

  out <- paste0(out, "### Thresholds used\n\n",
    md_table(tibble::tibble(
      setting = c("identity_threshold", "family_mean_r_flag", "pool_tolerance",
                  "min_family_size", "min_modal_seed_fraction"),
      value = as.character(c(
        settings$identity_threshold, settings$family_mean_r_flag,
        settings$pool_tolerance %||% "n/a - no pedigrees, so no families to pool",
        settings$min_family_size, settings$min_modal_seed_fraction)),
      meaning = c(
        "marker-profile correlation above which two accessions are one line",
        "mean within-family correlation above which a full-sib family did not segregate",
        "how close between-family and within-family correlation must be before two clonal families of one female are pooled",
        "smallest full-sib family summarised",
        "share of a group's pedigreed members that must share a seed parent before it is called a selfed female"
      )
    ), digits = 4))

  if (!is.null(families) && nrow(families) > 0) {
    clonal <- dplyr::filter(families, mean_r > settings$family_mean_r_flag)
    out <- paste0(out, "\n### Full-sib families that did not segregate\n\n",
      "Of ", nrow(families), " families with at least ", settings$min_family_size,
      " genotyped members, **", nrow(clonal), "** exceed ",
      settings$family_mean_r_flag, ". A segregating cross averages about ",
      round(mean(families$mean_r[families$mean_r <= settings$family_mean_r_flag]), 3),
      "; these average ", round(mean(clonal$mean_r), 3), ".\n\n",
      md_table(dplyr::select(clonal, family, n_members, n_pairs,
                             mean_r, min_r, max_r)))
  }

  if (n_renamed > 0) {
    tally <- names_tbl |>
      dplyr::count(analysis_name, name = "accessions") |>
      dplyr::arrange(dplyr::desc(accessions))

    out <- paste0(out, "\n### What was collapsed\n\n", md_table(tally),
      "\nEvery accession, with the reason recorded at the time:\n\n",
      md_table(dplyr::select(names_tbl, dplyr::any_of(
        c("germplasmName", "analysis_name", "pedigree", "reason")))))
  } else {
    out <- paste0(out, "\n### What was collapsed\n\n_Nothing._\n")
  }

  if (!is.null(groups) && nrow(groups) > 0) {
    n_groups <- dplyr::n_distinct(groups$group)
    out <- paste0(out, "\n### Near-identical groups\n\n",
      "At r > ", settings$identity_threshold, ": **", n_groups,
      "** group(s) covering **", nrow(groups), "** accessions. ",
      "Full membership in `", sp, "_identity_groups.csv`.\n")
  }

  out
}

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

header <- paste0(
  "# CURATION\n\n",
  "What was collapsed, and on what grounds.\n\n",
  "**Generated by `code/curation_report.R`** from the settings each curation ",
  "run recorded for itself — not written by hand, and not read out of the ",
  "scripts as they stand now, so it describes the run that produced the ",
  "files beside it. Re-run a species' curation and then this script after ",
  "changing a threshold.\n\n",
  "Written ", format(Sys.time(), "%Y-%m-%d %H:%M"), ".\n\n",
  "Accessions that markers show to be one genotype are **renamed, not ",
  "dropped**, so their phenotype records stay in the analysis under a single ",
  "name. The GRM is collapsed the same way, by averaging the rows and ",
  "columns of a group. The machine-readable mapping is ",
  "`oat_analysis_names.csv` and `pea_analysis_names.csv`; the reasoning ",
  "behind the thresholds is in [BACKGROUND.md](../BACKGROUND.md).\n\n",
  "Naming convention:\n\n",
  "| analysis name | means |\n|---|---|\n",
  "| `<seed>_<pollen>_no_cross` | a full-sib family that did not segregate |\n",
  "| `<seed_parent>_self` | matches such a family and shares its female, or a set of that female's clonal families that cannot be told apart |\n",
  "| `<seed_parent>` | confirmed identical to its genotyped parent |\n",
  "| an existing accession name | same line under two names, collapsed onto the better-phenotyped one |\n\n",
  "---\n\n"
)

writeLines(
  paste0(header,
         species_block("oat", "Oat"), "\n---\n\n",
         species_block("pea", "Pea")),
  out_file
)

message("wrote ", out_file)
