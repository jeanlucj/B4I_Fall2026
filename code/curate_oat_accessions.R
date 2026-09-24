# ============================================================
# CURATE THE B4I OAT ACCESSIONS:
# PEDIGREES AND NEAR-IDENTICAL MARKER PROFILES
#
# Answers two questions about the oat accessions in
# data/Acc_B4I_Avena.txt:
#
#   1. Do we have seed-parent and pollen-parent pedigrees?
#   2. Is there a set of accessions with near-identical marker scores?
#
# and then joins them, because the interesting case is a group of
# accessions with identical markers whose recorded pedigrees share a
# SEED parent but name DIFFERENT pollen parents.  That is the signature
# of a plant in the crossing nursery that self-pollinated: the progeny
# were recorded as several different crosses but are all selfs of the
# female, so they are genetically one line under several names.
#
# Keeping such a group in a genomic prediction training set inflates the
# apparent replication of one genotype and biases the GRM, so the script
# also proposes one representative per group.
#
# It also groups the accessions by pedigree -- same seed parent AND same
# pollen parent, i.e. a full-sib family -- and reports the average
# pairwise marker correlation within each family.  A real cross segregates,
# so its family averages around 0.75; a family that averages above 0.98 did
# not segregate at all and is a clonal set under several names.  That
# contrast is what sets the identity threshold.
#
# Every marker in play comes from one protocol, the Oat 3K array; the
# script checks that and stops if a second protocol has crept in.  The
# pedigree parents are genotyped on the same array where they exist, and
# carried in the same dosage matrix, so an accession can be compared
# directly with its own seed and pollen parent.
#
# Accessions in a clonal family are renamed for analysis rather than
# discarded, so their phenotypes are kept:
#   <seed>_<pollen>_no_cross  a family that did not segregate
#   <seed_parent>_self        matches a no_cross line and shares its female,
#                             or belongs to a set of clonal families of that
#                             female that are not distinguishable from each
#                             other
#   <seed_parent>             confirmed identical to the genotyped parent
#
# Outputs: output/oat_pedigrees.csv          seed/pollen parent per accession
#          output/oat_parent_genotypes.csv   which parents have Oat 3K data
#          output/oat_accession_vs_parents.csv  r to seed and pollen parent
#          output/oat_analysis_names.csv     original -> analysis name + reason
#          output/oat_top_marker_pairs.csv   the most similar pairs
#          output/oat_family_correlations.csv  per full-sib family
#          output/oat_family_pairs.csv       every within-family pair
#          output/oat_identity_groups.csv    the near-identical groups
#          output/oat_curation_keep_drop.csv one representative per group
#          output/oat_family_correlation.png the family histograms
#          output/oat_identity_similarity.png distribution + largest group
# ============================================================

library(tidyverse)

here::i_am("code/curate_oat_accessions.R")

source(here::here("code", "t3_functions.R"))        # connect_t3
source(here::here("code", "curation_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

db_name  <- "T3/Oat"
acc_file <- here::here("data", "Acc_B4I_Avena.txt")
out_dir  <- here::here("output")

# Accessions whose marker profiles correlate above this are treated as
# the same line.  The full-sib family analysis below is what calibrates
# it: within genuine crosses individual sib pairs reach 0.97-0.99, so
# 0.95 is not tight enough to separate a real sib from a duplicate.
# Every family that averages above 0.98 sits well clear of this.
identity_threshold <- 0.99

# A full-sib family whose mean pairwise correlation exceeds this did not
# segregate and is flagged as clonal.
family_mean_r_flag <- 0.985

# Two clonal families of the same female are treated as one selfed line
# when their between-family correlation is within this much of the weaker
# within-family correlation.
pool_tolerance <- 0.01

# Families need at least this many genotyped members to be summarised
min_family_size <- 2

# The genotyping protocol to read marker dosages from.  66 is the Oat 3K
# array, which covers essentially all of the B4I oats.
protocol_id <- "66"

# A group is called a selfing suspect when one female accounts for at
# least this fraction of its pedigreed members and those members name
# more than one different pollen parent.  A single stray member from
# another cross should not mask the pattern, so this is a majority rule
# rather than "exactly one seed parent".
min_modal_seed_fraction <- 0.5

# Check this many pedigree strings against the structured BrAPI endpoint
# to confirm that the string really is seed_parent/pollen_parent.
# 0 skips the check.
n_verify <- 25

monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

# Phenotype counts break ties when choosing which accession in a group
# to keep.  Optional: alphabetical order is used if it is absent.
observations_file <- here::here("data", "B4I_observations.rds")

# How many of the most-similar pairs to report. Longer than the pea list
# because far more oat accessions are collapsed, so the interesting cases --
# pairs that are close but were NOT merged -- are further down.
n_top_pairs <- 100L

# ============================================================
# Driver
# ============================================================

conn <- connect_t3(db_name)
cfg  <- T3GenoTools::geno_config(db_name, impute = "mean", progress = FALSE)

oat_names <- read_accessions(acc_file, monoculture_labels)
message("oat accessions requested: ", length(oat_names))

# --- Question 1 ---
pedigrees <- fetch_pedigrees(conn, oat_names)
readr::write_csv(pedigrees, file.path(out_dir, "oat_pedigrees.csv"))

cat("\n=== Q1: seed and pollen parent pedigrees ===\n")
cat("accessions found on ", db_name, ": ", nrow(pedigrees), "\n", sep = "")
cat("with a pedigree string    : ", sum(!is.na(pedigrees$pedigree)), "\n", sep = "")
cat("with a named seed parent  : ", sum(!is.na(pedigrees$seed_parent)), "\n", sep = "")
cat("with a named pollen parent: ", sum(!is.na(pedigrees$pollen_parent)), "\n", sep = "")
cat("with both parents named   : ", sum(pedigrees$both_parents), "\n", sep = "")

verify_parent_order(conn, pedigrees, n_verify)

# --- All markers from one protocol? ---
protocols <- confirm_single_protocol(conn, pedigrees$germplasmDbId, cfg, protocol_id)

# --- The pedigree parents, genotyped on the same array ---
parents <- parent_names(pedigrees)
message("distinct parents named in the pedigrees: ", length(parents))

parent_records <- fetch_parent_records(conn, parents)

cat("\n=== Genotyping protocols covering the PARENTS ===\n")
parent_protocols <- T3GenoTools::find_geno_sources(
  conn, unique(parent_records$germplasmDbId), cfg, "protocol"
)
print(as.data.frame(parent_protocols))

# --- Question 2: one dosage matrix holding accessions AND parents ---
train <- dplyr::bind_rows(
  dplyr::select(pedigrees, germplasmDbId, germplasmName),
  dplyr::select(parent_records, germplasmDbId, germplasmName)
) |>
  dplyr::distinct(germplasmDbId, .keep_all = TRUE)

dosage <- marker_dosage(conn, train, cfg, protocol_id)
message("marker dosages: ", nrow(dosage), " lines x ", ncol(dosage), " markers")

ungenotyped <- setdiff(pedigrees$germplasmName, rownames(dosage))
if (length(ungenotyped) > 0) {
  message(length(ungenotyped), " accession(s) have no markers and cannot be checked: ",
          paste(utils::head(ungenotyped, 10), collapse = ", "))
}

parents_genotyped <- intersect(parent_records$germplasmName, rownames(dosage))

parent_status <- parent_genotype_status(
  conn, parent_records, cfg, protocol_id, parents_genotyped
) |>
  dplyr::mutate(also_a_b4i_accession = germplasmName %in% pedigrees$germplasmName)

readr::write_csv(parent_status, file.path(out_dir, "oat_parent_genotypes.csv"))

cat("\n=== Parents with usable protocol-", protocol_id, " marker data ===\n", sep = "")
cat(length(parents_genotyped), " of ", nrow(parent_records),
    " parents carry dosages\n\n", sep = "")
print(as.data.frame(dplyr::count(parent_status, status, name = "n_parents")),
      row.names = FALSE)

not_genotyped <- dplyr::filter(parent_status, !has_dosage)
if (nrow(not_genotyped) > 0) {
  cat("\nparents without marker data:\n")
  print(as.data.frame(dplyr::select(not_genotyped, germplasmName, status)),
        row.names = FALSE)
}

# `similarity` spans accessions AND parents; the group analyses below are
# about the B4I accessions, so they use the accession-only block.
similarity <- identity_matrix(dosage)

acc_genotyped  <- intersect(pedigrees$germplasmName, rownames(similarity))
similarity_acc <- similarity[acc_genotyped, acc_genotyped, drop = FALSE]

# --- each accession against its own parents ---
vs_parents <- accession_vs_parents(pedigrees, similarity)
readr::write_csv(vs_parents, file.path(out_dir, "oat_accession_vs_parents.csv"))

cat("\n=== Accessions vs their own parents ===\n")
cat("with a genotyped seed parent  : ", sum(vs_parents$seed_genotyped), "\n", sep = "")
cat("with a genotyped pollen parent: ", sum(vs_parents$pollen_genotyped), "\n", sep = "")
cat("r with seed parent, quantiles:\n")
print(round(stats::quantile(vs_parents$r_seed_parent, c(0, .25, .5, .75, 1),
                            na.rm = TRUE), 3))
cat("r with pollen parent, quantiles:\n")
print(round(stats::quantile(vs_parents$r_pollen_parent, c(0, .25, .5, .75, 1),
                            na.rm = TRUE), 3))
cat("accessions identical (r > ", identity_threshold, ") to their seed parent: ",
    sum(vs_parents$r_seed_parent > identity_threshold, na.rm = TRUE), "\n", sep = "")

# --- full-sib families ---
families     <- pedigree_families(pedigrees, acc_genotyped, min_family_size)
family_pairs <- family_pair_correlations(families, similarity_acc)
family_summary <- summarise_families(family_pairs, family_mean_r_flag)

readr::write_csv(family_summary, file.path(out_dir, "oat_family_correlations.csv"))
readr::write_csv(family_pairs,   file.path(out_dir, "oat_family_pairs.csv"))

cat("\n=== Full-sib families (same seed AND pollen parent) ===\n")
cat("genotyped accessions with both parents: ",
    sum(pedigrees$both_parents & pedigrees$germplasmName %in% acc_genotyped),
    "\n", sep = "")
cat("families with >= ", min_family_size, " members: ", nrow(family_summary),
    ", covering ", dplyr::n_distinct(families$germplasmName), " accessions\n", sep = "")

cat("\nfamily MEAN pairwise correlation, quantiles:\n")
print(round(stats::quantile(family_summary$mean_r, c(0, .25, .5, .75, .9, 1)), 3))
cat("\nall within-family pairwise correlations, quantiles:\n")
print(round(stats::quantile(family_pairs$r, c(0, .25, .5, .75, .9, .99, 1)), 3))

cat("\nfamilies that did not segregate (mean r > ", family_mean_r_flag, "):\n", sep = "")
family_summary |>
  dplyr::filter(clonal_family) |>
  dplyr::select(family, n_members, n_pairs, mean_r, min_r, max_r) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 4)

cat("\nsegregating families, for contrast (lowest 8 mean r):\n")
family_summary |>
  dplyr::filter(!clonal_family) |>
  dplyr::slice_min(mean_r, n = 8) |>
  dplyr::select(family, n_members, n_pairs, mean_r, min_r, max_r) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 4)

# --- near-identical groups ---
groups     <- identity_groups(similarity_acc, identity_threshold)
classified <- classify_groups(groups, pedigrees, min_modal_seed_fraction)

curated <- choose_representatives(classified, observations_file)

readr::write_csv(
  dplyr::select(curated, group, n_members, germplasmName, pedigree,
                seed_parent, pollen_parent, n_observations,
                modal_seed, frac_modal_seed, n_pollen_on_modal,
                selfing_suspect, same_cross, interpretation),
  file.path(out_dir, "oat_identity_groups.csv")
)

readr::write_csv(
  dplyr::select(curated, germplasmName, group, action, representative,
                n_observations, pedigree, interpretation),
  file.path(out_dir, "oat_curation_keep_drop.csv")
)

cat("\n=== Q2: near-identical marker profiles ===\n")
pairwise <- similarity_acc[upper.tri(similarity_acc)]
cat("pairwise marker correlation, quantiles:\n")
print(round(stats::quantile(pairwise, c(0.5, 0.9, 0.99, 0.999, 1)), 3))
cat("\ngroups at r > ", identity_threshold, ": ", dplyr::n_distinct(curated$group),
    ", covering ", nrow(curated), " accessions\n", sep = "")

group_summary <- curated |>
  dplyr::distinct(group, n_members, n_pedigreed, n_seed_parents,
                  n_pollen_parents, modal_seed, n_pollen_on_modal,
                  selfing_suspect, interpretation) |>
  dplyr::arrange(dplyr::desc(n_members))

print(as.data.frame(group_summary), row.names = FALSE)

cat("\nselfing-suspect groups: ", sum(group_summary$selfing_suspect),
    " (", sum(curated$selfing_suspect), " accessions)\n", sep = "")
cat("accessions proposed for removal: ", sum(curated$action == "drop"), "\n", sep = "")

cat("\n=== Selfing-suspect groups in detail ===\n")
curated |>
  dplyr::filter(selfing_suspect) |>
  dplyr::select(group, germplasmName, seed_parent, pollen_parent,
                n_observations, action) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ------------------------------------------------------------
# Analysis names
# ------------------------------------------------------------

analysis_names <- resolve_analysis_names(
  pedigrees, family_summary, similarity,
  identity_threshold = identity_threshold,
  family_flag        = family_mean_r_flag,
  pool_tolerance     = pool_tolerance
)

readr::write_csv(analysis_names, file.path(out_dir, "oat_analysis_names.csv"))

# ------------------------------------------------------------
# The most similar pairs, merged or not
#
# The counts above say how many pairs cleared a threshold; they do not show
# which pairs came close and were left alone. Those are the ones worth a
# human look, so each pair is reported with whether curation actually merged
# it and what the pedigrees say about it -- a pair that is nearly identical
# and shares no parent is a different problem from one that is a pair of
# full sibs.
# ------------------------------------------------------------

final_name <- stats::setNames(rownames(similarity_acc), rownames(similarity_acc))
hit <- match(names(final_name), analysis_names$germplasmName)
final_name[!is.na(hit)] <- analysis_names$analysis_name[hit[!is.na(hit)]]

seed_of   <- stats::setNames(pedigrees$seed_parent,   pedigrees$germplasmName)
pollen_of <- stats::setNames(pedigrees$pollen_parent, pedigrees$germplasmName)

pedigree_link <- function(a, b) {
  sa <- seed_of[a]; sb <- seed_of[b]
  pa <- pollen_of[a]; pb <- pollen_of[b]
  dplyr::case_when(
    is.na(sa) & is.na(pa) | is.na(sb) & is.na(pb) ~ "no pedigree",
    !is.na(sa) & !is.na(sb) & sa == sb & !is.na(pa) & !is.na(pb) & pa == pb ~ "same cross",
    !is.na(sa) & !is.na(sb) & sa == sb ~ "shared seed parent",
    !is.na(pa) & !is.na(pb) & pa == pb ~ "shared pollen parent",
    TRUE ~ "no parent in common"
  )
}

top_pairs <- {
  keep <- upper.tri(similarity_acc)
  a <- rownames(similarity_acc)[row(similarity_acc)[keep]]
  b <- colnames(similarity_acc)[col(similarity_acc)[keep]]
  tibble::tibble(accession_1 = a, accession_2 = b, r = similarity_acc[keep]) |>
    dplyr::arrange(dplyr::desc(r)) |>
    utils::head(n_top_pairs) |>
    dplyr::mutate(
      rank = dplyr::row_number(),
      above_threshold = r > identity_threshold,
      # Curation merges whole clonal families, so a pair can end up with one
      # name without its own correlation clearing the threshold
      merged = final_name[accession_1] == final_name[accession_2],
      analysis_name = dplyr::if_else(merged, final_name[accession_1], NA_character_),
      pedigree_link = pedigree_link(accession_1, accession_2),
      .before = 1
    )
}

readr::write_csv(top_pairs, file.path(out_dir, "oat_top_marker_pairs.csv"))

cat("\n=== Most similar pairs of oat accessions (top ", n_top_pairs,
    ") ===\n", sep = "")
print(as.data.frame(utils::head(top_pairs, 15)), row.names = FALSE, digits = 4)

cat("\nof the top ", nrow(top_pairs), ": ", sum(top_pairs$merged),
    " merged by curation, ", sum(!top_pairs$merged), " left separate\n", sep = "")
cat("pedigree relationship among them:\n")
print(as.data.frame(dplyr::count(top_pairs, pedigree_link, merged,
                                 name = "pairs")), row.names = FALSE)

cat("\nclosest pairs curation did NOT merge:\n")
top_pairs |>
  dplyr::filter(!merged) |>
  utils::head(10) |>
  dplyr::select(rank, accession_1, accession_2, r, pedigree_link) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 4)

saveRDS(
  list(
    species = "oat", protocol_id = protocol_id,
    identity_threshold = identity_threshold,
    family_mean_r_flag = family_mean_r_flag,
    pool_tolerance = pool_tolerance,
    min_family_size = min_family_size,
    min_modal_seed_fraction = min_modal_seed_fraction,
    n_verify = n_verify,
    n_requested = length(oat_names),
    n_genotyped = length(acc_genotyped),
    genotyped_names = acc_genotyped,
    ungenotyped = ungenotyped,
    run_at = Sys.time()
  ),
  file.path(out_dir, "oat_curation_settings.rds")
)

cat("\n=== Accessions renamed for analysis ===\n")
cat("renamed: ", nrow(analysis_names), " of ", length(acc_genotyped),
    " genotyped accessions\n", sep = "")

analysis_names |>
  dplyr::count(analysis_name, name = "n_accessions") |>
  dplyr::arrange(dplyr::desc(n_accessions)) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nconfirmed identical to a genotyped seed parent: ",
    sum(analysis_names$seed_parent_confirmed), "\n", sep = "")

if (sum(analysis_names$seed_parent_confirmed) == 0) {
  message("no clonal accession reaches r > ", identity_threshold,
          " with its own seed parent, so none was renamed to a parent. ",
          "Where the seed parent is genotyped the clonal lines sit near it ",
          "but not on it (see oat_accession_vs_parents.csv).")
}

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

# Histogram of the within-family correlations: one panel for the family
# averages, one for every individual pair
fam_plot_data <- dplyr::bind_rows(
  dplyr::transmute(family_summary, r = mean_r,
                   panel = "Family mean pairwise correlation"),
  dplyr::transmute(family_pairs, r = r,
                   panel = "Every within-family pair")
) |>
  dplyr::mutate(panel = forcats::fct_inorder(panel))

p_family <- ggplot2::ggplot(fam_plot_data, ggplot2::aes(r)) +
  ggplot2::geom_histogram(binwidth = 0.01, fill = "grey35",
                          colour = "white", linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = family_mean_r_flag,
                      colour = "red", linetype = 2) +
  ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Marker correlation within oat full-sib families",
    subtitle = paste0("families share both parents; red line = ",
                      family_mean_r_flag,
                      ", above which a family has not segregated"),
    x = "correlation of marker dosages", y = "count"
  )

ggplot2::ggsave(file.path(out_dir, "oat_family_correlation.png"),
                p_family, width = 9, height = 7, dpi = 150)

p_hist <- tibble::tibble(r = pairwise) |>
  ggplot2::ggplot(ggplot2::aes(r)) +
  ggplot2::geom_histogram(bins = 200, fill = "grey35") +
  ggplot2::geom_vline(xintercept = identity_threshold,
                      colour = "red", linetype = 2) +
  ggplot2::scale_y_continuous(trans = "log1p",
                              breaks = c(0, 10, 100, 1000, 10000, 1e5)) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Marker-profile correlation between B4I oat accessions",
    subtitle = paste0("red line = identity threshold (", identity_threshold, ")"),
    x = "correlation of marker dosages", y = "pairs (log scale)"
  )

largest <- group_summary$group[1]
members <- dplyr::filter(curated, group == largest)$germplasmName

p_heat <- similarity_acc[members, members] |>
  as.data.frame() |>
  tibble::rownames_to_column("a") |>
  tidyr::pivot_longer(-a, names_to = "b", values_to = "r") |>
  ggplot2::ggplot(ggplot2::aes(a, b, fill = r)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(limits = c(min(similarity_acc[members, members]), 1)) +
  ggplot2::theme_minimal(base_size = 8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggplot2::labs(
    title = paste0("Largest near-identical group (group ", largest, ", ",
                   length(members), " accessions)"),
    x = NULL, y = NULL
  )

ggplot2::ggsave(file.path(out_dir, "oat_identity_similarity.png"),
                p_hist, width = 9, height = 5, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "oat_identity_largest_group.png"),
                p_heat, width = 9, height = 8, dpi = 150)

message("\nwrote:\n  ",
        paste(file.path("output", c(
          "oat_pedigrees.csv", "oat_parent_genotypes.csv",
          "oat_accession_vs_parents.csv", "oat_family_correlations.csv",
          "oat_family_pairs.csv", "oat_identity_groups.csv",
          "oat_curation_keep_drop.csv", "oat_analysis_names.csv",
          "oat_family_correlation.png", "oat_identity_similarity.png",
          "oat_identity_largest_group.png"
        )), collapse = "\n  "))
