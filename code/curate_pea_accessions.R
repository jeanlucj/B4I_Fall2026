# ============================================================
# CURATE THE B4I PEA ACCESSIONS: NEAR-IDENTICAL MARKER PROFILES
#
# The pea counterpart of code/curate_oat_accessions.R, using the same
# functions in code/curation_functions.R and the same thresholds, so the
# two species are judged on the same scale.
#
# One difference is forced by the data.  No B4I pea accession carries a
# pedigree on T3/Oat -- every one reads "NA/NA" -- so the half of the oat
# analysis built on pedigrees cannot be run here: there are no full-sib
# families to average within, no parents to compare against, and no basis
# for a "_no_cross" or "_self" label.  What remains, and what this script
# does, is the marker half: correlate every pair of pea accessions and
# find any that are too alike to be different lines.
#
# The script checks the pedigree situation each time rather than assuming
# it, so it will start reporting families if pedigrees are ever loaded.
#
# Outputs: output/pea_pedigrees.csv          what pedigree data exists
#          output/pea_identity_groups.csv    near-identical groups
#          output/pea_curation_keep_drop.csv one representative per group
#          output/pea_analysis_names.csv     original -> analysis name
#          output/pea_identity_similarity.png distribution of correlations
#          output/pea_identity_largest_group.png
# ============================================================

library(tidyverse)

here::i_am("code/curate_pea_accessions.R")

source(here::here("code", "curation_functions.R"))

# ------------------------------------------------------------
# Settings -- deliberately the same thresholds as the oat script
# ------------------------------------------------------------

db_name  <- "T3/Oat"          # the B4I peas live on the oat instance
acc_file <- here::here("data", "Acc_B4I_Pisum.txt")
out_dir  <- here::here("output")

identity_threshold <- 0.99
family_mean_r_flag <- 0.985
min_family_size    <- 2
min_modal_seed_fraction <- 0.5

# GenoPea 13K Array
protocol_id <- "67"

monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

observations_file <- here::here("output", "B4I_observations.rds")

# ============================================================
# Driver
# ============================================================

conn <- connect_t3(db_name)
cfg  <- T3GenoTools::geno_config(db_name, impute = "mean", progress = FALSE)

pea_names <- read_accessions(acc_file, monoculture_labels)
message("pea accessions requested: ", length(pea_names))

# --- pedigrees, such as they are ---
pedigrees <- fetch_pedigrees(conn, pea_names)
readr::write_csv(pedigrees, file.path(out_dir, "pea_pedigrees.csv"))

cat("\n=== Pedigrees ===\n")
cat("accessions found on ", db_name, ": ", nrow(pedigrees), "\n", sep = "")
cat("with a named seed parent  : ", sum(!is.na(pedigrees$seed_parent)), "\n", sep = "")
cat("with a named pollen parent: ", sum(!is.na(pedigrees$pollen_parent)), "\n", sep = "")

have_pedigrees <- any(pedigrees$both_parents)

if (!have_pedigrees) {
  message("No pea accession has a recorded pedigree, so the full-sib family ",
          "and parent-identity analyses are skipped. Only the marker-profile ",
          "comparison below applies.")
}

# --- one protocol? ---
protocols <- confirm_single_protocol(conn, pedigrees$germplasmDbId, cfg, protocol_id)

# --- marker dosages ---
dosage <- marker_dosage(conn, pedigrees, cfg, protocol_id)
message("marker dosages: ", nrow(dosage), " lines x ", ncol(dosage), " markers")

ungenotyped <- setdiff(pedigrees$germplasmName, rownames(dosage))
if (length(ungenotyped) > 0) {
  message(length(ungenotyped), " accession(s) have no markers and cannot be checked: ",
          paste(utils::head(ungenotyped, 10), collapse = ", "))
}

similarity <- identity_matrix(dosage)
acc_genotyped <- intersect(pedigrees$germplasmName, rownames(similarity))
similarity_acc <- similarity[acc_genotyped, acc_genotyped, drop = FALSE]

# --- full-sib families, if there is anything to work with ---
if (have_pedigrees) {
  families       <- pedigree_families(pedigrees, acc_genotyped, min_family_size)
  family_pairs   <- family_pair_correlations(families, similarity_acc)
  family_summary <- summarise_families(family_pairs, family_mean_r_flag)

  readr::write_csv(family_summary, file.path(out_dir, "pea_family_correlations.csv"))
  readr::write_csv(family_pairs,   file.path(out_dir, "pea_family_pairs.csv"))

  cat("\n=== Full-sib families ===\n")
  print(as.data.frame(family_summary), row.names = FALSE, digits = 4)
}

# --- near-identical accessions ---
groups     <- identity_groups(similarity_acc, identity_threshold)
classified <- classify_groups(groups, pedigrees, min_modal_seed_fraction)
curated    <- choose_representatives(classified, observations_file)

readr::write_csv(
  dplyr::select(curated, group, n_members, germplasmName, pedigree,
                seed_parent, pollen_parent, n_observations, interpretation),
  file.path(out_dir, "pea_identity_groups.csv")
)

readr::write_csv(
  dplyr::select(curated, germplasmName, group, action, representative,
                n_observations, interpretation),
  file.path(out_dir, "pea_curation_keep_drop.csv")
)

# With no pedigrees there is no "_no_cross" or "_self" to assign, so the
# only rename available is to collapse a near-identical group onto its
# representative -- typically a released cultivar and its experimental
# designation.
analysis_names <- curated |>
  dplyr::filter(action == "drop") |>
  dplyr::transmute(
    germplasmName,
    analysis_name = representative,
    reason = paste0("marker profile identical (r > ", identity_threshold,
                    ") to ", representative, ": same line under two names")
  )

readr::write_csv(analysis_names, file.path(out_dir, "pea_analysis_names.csv"))

cat("\n=== Near-identical marker profiles ===\n")
pairwise <- similarity_acc[upper.tri(similarity_acc)]
cat("genotyped pea accessions: ", nrow(similarity_acc), "\n", sep = "")
cat("pairwise marker correlation, quantiles:\n")
print(round(stats::quantile(pairwise, c(0.5, 0.9, 0.99, 0.999, 1)), 3))

for (thr in c(0.9, 0.95, identity_threshold)) {
  cat("pairs with r > ", thr, ": ", sum(pairwise > thr), "\n", sep = "")
}

n_groups <- dplyr::n_distinct(curated$group)
cat("\ngroups at r > ", identity_threshold, ": ", n_groups,
    ", covering ", nrow(curated), " accessions\n", sep = "")

if (n_groups > 0) {
  curated |>
    dplyr::distinct(group, n_members, interpretation) |>
    dplyr::arrange(dplyr::desc(n_members)) |>
    as.data.frame() |>
    print(row.names = FALSE)

  cat("\naccessions proposed for removal: ", sum(curated$action == "drop"), "\n", sep = "")
} else {
  message("No pair of pea accessions reaches r > ", identity_threshold,
          ": nothing to curate on the pea side.")
}

# ------------------------------------------------------------
# Figures -- same form as the oat script, for direct comparison
# ------------------------------------------------------------

p_hist <- tibble::tibble(r = pairwise) |>
  ggplot2::ggplot(ggplot2::aes(r)) +
  ggplot2::geom_histogram(bins = 200, fill = "grey35") +
  ggplot2::geom_vline(xintercept = identity_threshold,
                      colour = "red", linetype = 2) +
  ggplot2::scale_y_continuous(trans = "log1p",
                              breaks = c(0, 10, 100, 1000, 10000, 1e5)) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Marker-profile correlation between B4I pea accessions",
    subtitle = paste0("red line = identity threshold (", identity_threshold, ")"),
    x = "correlation of marker dosages", y = "pairs (log scale)"
  )

ggplot2::ggsave(file.path(out_dir, "pea_identity_similarity.png"),
                p_hist, width = 9, height = 5, dpi = 150)

if (n_groups > 0) {
  largest <- curated |>
    dplyr::count(group, name = "n") |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
    dplyr::pull(group)
  members <- dplyr::filter(curated, group == largest)$germplasmName

  p_heat <- similarity_acc[members, members] |>
    as.data.frame() |>
    tibble::rownames_to_column("a") |>
    tidyr::pivot_longer(-a, names_to = "b", values_to = "r") |>
    ggplot2::ggplot(ggplot2::aes(a, b, fill = r)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      limits = c(min(similarity_acc[members, members]), 1)
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggplot2::labs(
      title = paste0("Largest near-identical pea group (", length(members),
                     " accessions)"),
      x = NULL, y = NULL
    )

  ggplot2::ggsave(file.path(out_dir, "pea_identity_largest_group.png"),
                  p_heat, width = 9, height = 8, dpi = 150)
}

message("\ndone")
