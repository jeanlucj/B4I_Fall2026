# ============================================================
# BIVARIATE DGE-IGE (PRODUCER-ASSOCIATE) MODEL FOR OAT-PEA INTERCROPS
#
# Fits the joint model of docs/B4I_Proposal_Models.docx:
#
#   (y_oat, y_pea)' = (I2 (x) Z_oat)(Pr_oat, As_oat->pea)'
#                   + (I2 (x) Z_pea)(As_pea->oat, Pr_pea)'
#                   + (I2 (x) Z_mix)(S_oatxpea->oat, S_peaxoat->pea)'
#                   + trial + block + (e_oat, e_pea)'
#
#   (Pr_oat, As_oat->pea)' ~ N(0, Sigma_oat (x) G_oat)
#   (As_pea->oat, Pr_pea)' ~ N(0, Sigma_pea (x) G_pea)
#   (S_...->oat, S_...->pea)' ~ N(0, Sigma_mix (x) G_oat (x) G_pea)
#   (e_oat, e_pea)' ~ N(0, R), R unstructured
#
# The off-diagonal of Sigma_oat is sigma_PrAs,oat, the producer-associate
# covariance that motivates the joint analysis; Sigma_pea likewise.
#
# Each genetic term is fitted as a BRR on X = Z L, where L L' = G.
# That is the same model as an RKHS term with kernel Z G Z' -- the
# implied covariance is Z L L' Z' (x) Sigma = Z G Z' (x) Sigma -- but
# the coefficients stay on the accession scale, so the producer and
# associate effects come back named by accession.  (An RKHS term's
# $beta lives in the eigenvector basis of the plot-level kernel and
# carries no accession identity at all.)
#
# Inputs : a plot-level phenotype table (see `pheno_file` below)
#          output/GRM_Avena.rds, output/GRM_Pisum.rds from
#          code/create_GRMs_T3.R
# Outputs: output/BGLR_multitrait_*.rds  (fitted model)
#          output/BGLR_{oat,pea}_effects_all_seeds.csv
#          output/BGLR_{oat,pea}_rank_stability.csv
#          output/BGLR_variance_components.csv
# ============================================================

library(tidyverse)

here::i_am("code/BGLR_multi_trait_model.R")

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

# Plot-level phenotypes.  Required columns:
#   studyYear, studyName, blockNumber,
#   germplasmName           (oat accession),
#   intercropGermplasmName  (pea accession),
#   oat_yield, pea_yield
pheno_file <- here::here("data", "B4I_2025_intercrop_pheno.rds")

grm_files <- c(
  oat = here::here("output", "GRM_Avena.rds"),
  pea = here::here("output", "GRM_Pisum.rds")
)

out_dir <- here::here("output")

study_year <- 2025

trials <- c(
  "B4I_2025_IA",
  "B4I_2025_IL",
  "B4I_2025_ND",
  "B4I_2025_NY"
)

# The first seed gives the fit that is saved; the rest are there to
# check that the accession rankings are stable across chains.
seeds <- c(12567, 129, 456, 789)

nIter  <- 20000
burnIn <- 3000
thin   <- 10

# Fit the specific-combination (SMA / direct x associate) term?  It is
# estimable only if combinations are replicated -- see the design
# diagnostics below.
fit_mix_term <- TRUE

# Placeholder germplasm names used in T3 to record that one component of
# the intercrop was not sown.  A plot carrying one of these is a
# monoculture, not an intercrop, so it says nothing about mixing ability
# and is dropped.
monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

# ------------------------------------------------------------
# Phenotypes
# ------------------------------------------------------------

if (!file.exists(pheno_file)) {
  stop(
    "Phenotype file not found:\n  ", pheno_file,
    "\nPoint `pheno_file` at the plot-level B4I table. It needs the columns ",
    "studyYear, studyName, blockNumber, germplasmName, ",
    "intercropGermplasmName, oat_yield, pea_yield.",
    call. = FALSE
  )
}

pheno <- if (grepl("\\.rds$", pheno_file, ignore.case = TRUE)) {
  readRDS(pheno_file)
} else {
  readr::read_csv(pheno_file, show_col_types = FALSE)
}

needed <- c("studyYear", "studyName", "blockNumber", "germplasmName",
            "intercropGermplasmName", "oat_yield", "pea_yield")
if (!all(needed %in% names(pheno))) {
  stop("Phenotype file is missing column(s): ",
       paste(setdiff(needed, names(pheno)), collapse = ", "), call. = FALSE)
}

# ------------------------------------------------------------
# Curation: drop monoculture plots
#
# The model is about what an oat and a pea do to each other, so a plot
# where only one of the two was sown carries no information about any
# producer, associate or specific-combination effect.  T3 records those
# plots with a placeholder germplasm name on the missing side.
# ------------------------------------------------------------

monoculture <- pheno |>
  dplyr::filter(
    germplasmName %in% monoculture_labels |
      intercropGermplasmName %in% monoculture_labels
  )

if (nrow(monoculture) > 0) {
  message("dropping ", nrow(monoculture), " monoculture plot(s)")
  monoculture |>
    dplyr::count(germplasmName, intercropGermplasmName, name = "plots") |>
    as.data.frame() |>
    print()
}

pheno <- pheno |>
  dplyr::filter(
    !germplasmName %in% monoculture_labels,
    !intercropGermplasmName %in% monoculture_labels
  )

# One tidy set of names is used from here on: oatAcc / peaAcc / mixID.
# mixID is built from the two accession names rather than taken from the
# data, so that the combination kernel below is exactly G_oat (x) G_pea
# on the observed combinations.
grainWgt <- pheno |>
  dplyr::mutate(
    oatYield     = oat_yield,
    peaYield     = pea_yield,
    oatAcc       = as.character(germplasmName),
    peaAcc       = as.character(intercropGermplasmName),
    mixID        = paste(oatAcc, peaAcc, sep = "::"),
    trialF       = factor(studyName),
    blockNumberF = factor(paste(studyYear, studyName, blockNumber))
  ) |>
  dplyr::filter(
    as.character(studyYear) == as.character(study_year),
    studyName %in% trials,
    !is.na(oatYield), !is.na(peaYield),
    !is.na(oatAcc), !is.na(peaAcc)
  ) |>
  dplyr::mutate(
    trialF       = droplevels(trialF),
    blockNumberF = droplevels(blockNumberF)
  )

message("plots: ", nrow(grainWgt),
        " | trials: ", nlevels(grainWgt$trialF),
        " | blocks: ", nlevels(grainWgt$blockNumberF),
        " | oat: ", dplyr::n_distinct(grainWgt$oatAcc),
        " | pea: ", dplyr::n_distinct(grainWgt$peaAcc),
        " | combinations: ", dplyr::n_distinct(grainWgt$mixID))

print(table(grainWgt$studyName))
print(summary(dplyr::select(grainWgt, oatYield, peaYield)))

stopifnot(nrow(grainWgt) > 0)

# ------------------------------------------------------------
# Design diagnostics
#
# Producer and associate effects within a species are separable only if
# each accession meets several partners; the specific-combination term
# is separable from the residual only if combinations are replicated.
# ------------------------------------------------------------

partners_oat <- grainWgt |>
  dplyr::distinct(oatAcc, peaAcc) |>
  dplyr::count(oatAcc, name = "n_partners")

partners_pea <- grainWgt |>
  dplyr::distinct(peaAcc, oatAcc) |>
  dplyr::count(peaAcc, name = "n_partners")

message("distinct pea partners per oat accession: median ",
        stats::median(partners_oat$n_partners),
        ", range ", paste(range(partners_oat$n_partners), collapse = "-"))
message("distinct oat partners per pea accession: median ",
        stats::median(partners_pea$n_partners),
        ", range ", paste(range(partners_pea$n_partners), collapse = "-"))

if (max(partners_oat$n_partners) == 1 || max(partners_pea$n_partners) == 1) {
  warning("Some accessions appear with a single partner: producer and ",
          "associate effects are then aliased for those accessions.",
          call. = FALSE)
}

combo_reps <- as.integer(table(grainWgt$mixID))
message("plots per combination: median ", stats::median(combo_reps),
        ", range ", paste(range(combo_reps), collapse = "-"))

if (fit_mix_term && mean(combo_reps == 1) > 0.9) {
  warning("More than 90% of combinations occur in a single plot: the ",
          "specific-combination term is nearly confounded with the residual. ",
          "Consider fit_mix_term <- FALSE.", call. = FALSE)
}

# ------------------------------------------------------------
# GRMs
#
# create_GRMs_T3.R saves the whole build_grm() object; G is $G.
# ------------------------------------------------------------

read_grm <- function(path) {
  if (!file.exists(path)) {
    stop("GRM not found:\n  ", path,
         "\nRun code/create_GRMs_T3.R first.", call. = FALSE)
  }
  g <- readRDS(path)
  if (is.list(g) && !is.null(g$G)) g$G else g
}

G_oat_all <- read_grm(grm_files[["oat"]])
G_pea_all <- read_grm(grm_files[["pea"]])

missing_oat <- setdiff(unique(grainWgt$oatAcc), rownames(G_oat_all))
missing_pea <- setdiff(unique(grainWgt$peaAcc), rownames(G_pea_all))

if (length(missing_oat) > 0 || length(missing_pea) > 0) {
  stop("Accessions phenotyped but absent from the GRMs -- oat: ",
       length(missing_oat), ", pea: ", length(missing_pea), "\n  ",
       paste(utils::head(c(missing_oat, missing_pea), 10), collapse = ", "),
       call. = FALSE)
}

# Subset and order the GRMs by the accessions actually in the trial
oatAccs <- sort(unique(grainWgt$oatAcc))
peaAccs <- sort(unique(grainWgt$peaAcc))
mixIDs  <- sort(unique(grainWgt$mixID))

G_oat <- G_oat_all[oatAccs, oatAccs, drop = FALSE]
G_pea <- G_pea_all[peaAccs, peaAccs, drop = FALSE]

# Combination kernel: G_oat (x) G_pea restricted to the observed
# combinations, i.e. K[k, l] = G_oat[oat_k, oat_l] * G_pea[pea_k, pea_l]
mix_parts <- stringr::str_split_fixed(mixIDs, stringr::fixed("::"), 2)
G_mix <- G_oat[mix_parts[, 1], mix_parts[, 1], drop = FALSE] *
         G_pea[mix_parts[, 2], mix_parts[, 2], drop = FALSE]
dimnames(G_mix) <- list(mixIDs, mixIDs)

for (nm in c("G_oat", "G_pea", "G_mix")) {
  G <- get(nm)
  cat("\n", nm, ": ", nrow(G), " x ", ncol(G), "\n", sep = "")
  cat("  diagonal      : "); print(summary(diag(G)))
  cat("  symmetric     : ", isSymmetric(unname(G)), "\n", sep = "")
  cat("  min eigenvalue: ",
      min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), "\n", sep = "")
}

# ------------------------------------------------------------
# Response and design matrices
# ------------------------------------------------------------

Y <- as.matrix(grainWgt[, c("peaYield", "oatYield")])
colnames(Y) <- c("peaYield", "oatYield")
stopifnot(!anyNA(Y))

# Full set of trial dummies WITHOUT a separate intercept.  A full dummy
# set plus an intercept is rank-deficient: BGLR does not error, it
# samples along the ridge, so neither the intercept nor the trial
# effects mean anything on their own and the chain mixes badly.
incTrials <- stats::model.matrix(~ 0 + trialF, grainWgt)
colnames(incTrials) <- levels(grainWgt$trialF)

incBlocks <- stats::model.matrix(~ 0 + blockNumberF, grainWgt)
colnames(incBlocks) <- levels(grainWgt$blockNumberF)

# Incidence matrices, with factor levels forced to the GRM row order
incidence <- function(x, levels) {
  Z <- stats::model.matrix(~ 0 + factor(x, levels = levels))
  colnames(Z) <- levels
  Z
}

Z_oat <- incidence(grainWgt$oatAcc, rownames(G_oat))
Z_pea <- incidence(grainWgt$peaAcc, rownames(G_pea))
Z_mix <- incidence(grainWgt$mixID,  rownames(G_mix))

# ------------------------------------------------------------
# G = L L', so that a BRR on Z L is the RKHS model with kernel Z G Z'
# while keeping the coefficients on the accession scale
# ------------------------------------------------------------

grm_factor <- function(G, tol = 1e-8) {
  e <- eigen(G, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  if (sum(keep) < ncol(G)) {
    message("  ", ncol(G) - sum(keep), " of ", ncol(G),
            " dimensions dropped as numerically null")
  }
  L <- sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), "*")
  rownames(L) <- rownames(G)
  L
}

L_oat <- grm_factor(G_oat)
L_pea <- grm_factor(G_pea)
L_mix <- grm_factor(G_mix)

stopifnot(
  max(abs(tcrossprod(L_oat) - G_oat)) < 1e-8,
  max(abs(tcrossprod(L_pea) - G_pea)) < 1e-8,
  max(abs(tcrossprod(L_mix) - G_mix)) < 1e-8
)

# ------------------------------------------------------------
# ETA
# ------------------------------------------------------------

ETA <- list(
  trial = list(X = incTrials,      model = "FIXED"),
  block = list(X = incBlocks,      model = "BRR"),
  G_pea = list(X = Z_pea %*% L_pea, model = "BRR"),
  G_oat = list(X = Z_oat %*% L_oat, model = "BRR")
)

if (fit_mix_term) {
  ETA$G_mix <- list(X = Z_mix %*% L_mix, model = "BRR")
}

# ------------------------------------------------------------
# Fit
# ------------------------------------------------------------

fit_one <- function(seed_value) {
  message("\n--- fitting, seed ", seed_value, " ---")
  set.seed(seed_value)
  BGLR::Multitrait(
    y         = Y,
    ETA       = ETA,
    intercept = FALSE,
    resCov    = list(df0 = 4, S0 = NULL, type = "UN"),
    nIter     = nIter,
    burnIn    = burnIn,
    thin      = thin,
    saveAt    = file.path(out_dir, paste0("BGLR_seed_", seed_value, "_")),
    verbose   = FALSE
  )
}

fits <- rlang::set_names(purrr::map(seeds, fit_one), seeds)

saveRDS(
  fits[[1]],
  file.path(out_dir, "BGLR_multitrait_pea_oat_yield_Gpea_Goat_Gmix.rds")
)

# ------------------------------------------------------------
# Variance components
#
# This is what the bivariate model is for: Sigma's off-diagonal is the
# producer-associate covariance, which the two separate single-species
# analyses cannot reach.
# ------------------------------------------------------------

# Trait 1 is peaYield, trait 2 is oatYield, so for the oat kernel the
# effect on oat yield is the producer effect and the effect on pea yield
# is the associate effect; for the pea kernel it is the other way round.
effect_roles <- list(
  G_oat = c(Pr = "oatYield", As = "peaYield"),
  G_pea = c(Pr = "peaYield", As = "oatYield"),
  G_mix = c(Pr = "oatYield", As = "peaYield")  # S->oat, S->pea
)

covariance_components <- function(fit, term) {
  Omega <- fit$ETA[[term]]$Cov$Omega
  dimnames(Omega) <- list(colnames(Y), colnames(Y))
  roles <- effect_roles[[term]]
  tibble::tibble(
    term        = term,
    var_Pr      = Omega[roles[["Pr"]], roles[["Pr"]]],
    var_As      = Omega[roles[["As"]], roles[["As"]]],
    cov_PrAs    = Omega[roles[["Pr"]], roles[["As"]]],
    cor_PrAs    = cov_PrAs / sqrt(var_Pr * var_As)
  )
}

genetic_terms <- intersect(c("G_oat", "G_pea", "G_mix"), names(ETA))

varcomp <- purrr::imap(fits, \(fit, seed) {
  purrr::map(genetic_terms, \(tm) covariance_components(fit, tm)) |>
    purrr::list_rbind() |>
    dplyr::mutate(seed = as.integer(seed), .before = 1)
}) |>
  purrr::list_rbind()

cat("\n=== Variance components (Pr = producer, As = associate) ===\n")
print(as.data.frame(varcomp), digits = 4)

resid_cov <- purrr::imap(fits, \(fit, seed) {
  R <- fit$resCov$R
  dimnames(R) <- list(colnames(Y), colnames(Y))
  tibble::tibble(
    seed       = as.integer(seed),
    var_pea    = R["peaYield", "peaYield"],
    var_oat    = R["oatYield", "oatYield"],
    cov_pea_oat = R["peaYield", "oatYield"],
    cor_pea_oat = cov_pea_oat / sqrt(var_pea * var_oat)
  )
}) |>
  purrr::list_rbind()

cat("\n=== Residual covariance between oat and pea yield on a plot ===\n")
print(as.data.frame(resid_cov), digits = 4)

readr::write_csv(
  dplyr::bind_rows(
    dplyr::mutate(varcomp, component = "genetic"),
    dplyr::mutate(resid_cov, term = "residual", component = "residual")
  ),
  file.path(out_dir, "BGLR_variance_components.csv")
)

# ------------------------------------------------------------
# Accession-level producer and associate effects
#
# beta is on the L basis; L %*% beta puts the effects back on the
# accession scale, where they carry the accession names.
# ------------------------------------------------------------

accession_effects <- function(fit, term, L) {
  g <- L %*% fit$ETA[[term]]$beta
  dimnames(g) <- list(rownames(L), colnames(Y))
  roles <- effect_roles[[term]]

  tibble::as_tibble(g, rownames = "accession") |>
    dplyr::transmute(
      accession,
      PrEff = .data[[roles[["Pr"]]]],
      AsEff = .data[[roles[["As"]]]],
      GMA   = PrEff + AsEff   # GMA = Pr + As, docs/B4I_Proposal_Models.docx
    )
}

species_effects <- function(term, L, id_name) {
  purrr::imap(fits, \(fit, seed) {
    accession_effects(fit, term, L) |>
      dplyr::mutate(seed = as.integer(seed), .before = 1)
  }) |>
    purrr::list_rbind() |>
    dplyr::rename({{ id_name }} := accession)
}

oat_all_seeds <- species_effects("G_oat", L_oat, "oatAcc")
pea_all_seeds <- species_effects("G_pea", L_pea, "peaAcc")

rank_summary <- function(effects, id) {
  effects |>
    dplyr::group_by(seed) |>
    dplyr::mutate(rank_GMA = dplyr::min_rank(dplyr::desc(GMA))) |>
    dplyr::group_by(.data[[id]]) |>
    dplyr::summarise(
      mean_rank = mean(rank_GMA),
      sd_rank   = stats::sd(rank_GMA),
      mean_GMA  = mean(GMA),
      sd_GMA    = stats::sd(GMA),
      mean_Pr   = mean(PrEff),
      mean_As   = mean(AsEff),
      .groups   = "drop"
    ) |>
    dplyr::arrange(mean_rank)
}

oat_rank_summary <- rank_summary(oat_all_seeds, "oatAcc")
pea_rank_summary <- rank_summary(pea_all_seeds, "peaAcc")

cat("\n=== Top oat accessions by mixing ability ===\n")
print(head(oat_rank_summary, 20))
cat("\n=== Top pea accessions by mixing ability ===\n")
print(head(pea_rank_summary, 20))

# Between-chain agreement of the accession rankings
chain_agreement <- function(effects, id) {
  wide <- effects |>
    dplyr::select(seed, dplyr::all_of(id), GMA) |>
    tidyr::pivot_wider(names_from = seed, values_from = GMA) |>
    dplyr::select(-dplyr::all_of(id))
  stats::cor(wide, method = "spearman")
}

cat("\n=== Rank correlation of oat GMA between chains ===\n")
print(round(chain_agreement(oat_all_seeds, "oatAcc"), 3))
cat("\n=== Rank correlation of pea GMA between chains ===\n")
print(round(chain_agreement(pea_all_seeds, "peaAcc"), 3))

readr::write_csv(oat_all_seeds,    file.path(out_dir, "BGLR_oat_effects_all_seeds.csv"))
readr::write_csv(pea_all_seeds,    file.path(out_dir, "BGLR_pea_effects_all_seeds.csv"))
readr::write_csv(oat_rank_summary, file.path(out_dir, "BGLR_oat_rank_stability.csv"))
readr::write_csv(pea_rank_summary, file.path(out_dir, "BGLR_pea_rank_stability.csv"))

# ------------------------------------------------------------
# Producer vs associate effects
# ------------------------------------------------------------

pr_as_plot <- function(effects, species) {
  labels <- effects |>
    dplyr::group_by(seed) |>
    dplyr::summarise(
      r = stats::cor(PrEff, AsEff, use = "complete.obs"),
      .groups = "drop"
    ) |>
    dplyr::mutate(label = paste0("r = ", round(r, 3)))

  ggplot2::ggplot(effects, ggplot2::aes(PrEff, AsEff)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "red") +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 4
    ) +
    ggplot2::facet_wrap(~ seed) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste(species, "producer vs associate effects across chains"),
      x = "Producer effect", y = "Associate effect"
    )
}

p_oat <- pr_as_plot(oat_all_seeds, "Oat")
p_pea <- pr_as_plot(pea_all_seeds, "Pea")

print(p_oat)
print(p_pea)

ggplot2::ggsave(file.path(out_dir, "BGLR_oat_Pr_vs_As.png"), p_oat,
                width = 8, height = 6, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "BGLR_pea_Pr_vs_As.png"), p_pea,
                width = 8, height = 6, dpi = 150)
