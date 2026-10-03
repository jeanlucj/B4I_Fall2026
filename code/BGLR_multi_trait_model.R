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
#          data/GRM_Avena.rds, data/GRM_Pisum.rds from
#          code/create_GRMs_T3.R
# Outputs: output/BGLR_multitrait_*.rds  (fitted model)
#          output/BGLR_{oat,pea}_effects_all_seeds.csv
#          output/BGLR_{oat,pea}_rank_stability.csv
#          output/BGLR_variance_components.csv
# ============================================================

library(tidyverse)

here::i_am("code/BGLR_multi_trait_model.R")

source(here::here("code", "dge_ige_functions.R"))
# read_beta_draws(): the streamed-draw reader, which already knows that BGLR
# writes every thinned iteration INCLUDING burn-in. Reused rather than
# rewritten, so there is one place that knows about that trap.
source(here::here("code", "interaction_decomp.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

# Plot-level phenotypes.  Required columns:
#   studyYear, studyName, blockNumber,
#   germplasmName           (oat accession),
#   intercropGermplasmName  (pea accession),
#   oat_yield, pea_yield
pheno_file <- here::here("output", "B4I_intercrop_pheno.rds")

grm_files <- c(
  oat = here::here("data", "GRM_Avena.rds"),
  pea = here::here("data", "GRM_Pisum.rds")
)

# The curation scripts rename accessions that are genetically one line.
# The phenotype table already carries those names, so the GRMs have to be
# collapsed the same way or the two will not join.
analysis_name_files <- c(
  oat = here::here("output", "oat_analysis_names.csv"),
  pea = here::here("output", "pea_analysis_names.csv")
)

out_dir <- here::here("output")

# WHICH TRIALS ARE FITTED IS NOT SET HERE. It comes from output/trial_qc.csv
# via apply_trial_qc(), inside b4i_fit_frame(). There used to be a hard-coded
# `trials` whitelist at this spot, and because it was never updated when the
# 2026 trials landed it silently held the fit to five trials and 1,985 plots
# while every table in VALIDATION_DESIGN.md reported the phenotype file's nine
# and 3,567. Nothing compared the two, so nothing noticed for two vintages.
# See SELF_CRITIQUE.md finding A. To include or exclude a trial, edit
# data/trial_qc_manual.csv and re-run code/curate_trials.R.
#
# study_years remains only as a guard against a stray year arriving from T3.
study_years <- c(2025L, 2026L)

# The first seed gives the fit that is saved; the rest are there to
# check that the accession rankings are stable across chains.
seeds <- c(12567, 129, 456, 789)

nIter  <- 20000
burnIn <- 3000
thin   <- 10

# Precision of the streamed per-draw coefficients. "single" halves the files
# (~13 MB per seed for both genetic terms rather than ~26) and is far more
# precision than a posterior variance needs. The tight colMeans(draws) == beta
# identity in tests/test_fits.R uses "double", on a 400-iteration chain where
# the file is tiny -- so do not "fix" the two to agree.
draw_storage_mode <- "single"

# Fit the specific-combination (SMA / direct x associate) term?  It is
# estimable only if combinations are replicated.  In the B4I trials about
# 91% of oat-pea combinations occur in a single plot, so the term would be
# confounded with the residual and it is left out of this first fit.  The
# design diagnostics below report the replication either way.
fit_mix_term <- FALSE

# Placeholder germplasm names used in T3 to record that one component of
# the intercrop was not sown.  A plot carrying one of these is a
# monoculture, not an intercrop, so it says nothing about mixing ability
# and is dropped.
monoculture_labels <- c("NO_OATS_PLANTED", "NO_PEAS_PLANTED")

# ------------------------------------------------------------
# Phenotypes
# ------------------------------------------------------------

# ------------------------------------------------------------
# Phenotypes
#
# One filter chain, shared with validation_inputs() through b4i_fit_frame():
# trial QC, monoculture plots, missing yields, and the tidy modelling names.
# Sharing it is what makes the plot set the fit saw assertable downstream
# rather than re-derived and hoped equal.
# ------------------------------------------------------------

grainWgt <- b4i_fit_frame(pheno_file = pheno_file,
                          study_years = study_years,
                          monoculture_labels = monoculture_labels)

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

n_single_oat <- sum(partners_oat$n_partners == 1)
n_single_pea <- sum(partners_pea$n_partners == 1)

if (n_single_oat > 0 || n_single_pea > 0) {
  warning(n_single_oat, " oat and ", n_single_pea,
          " pea accession(s) appear with a single partner. Their producer and ",
          "associate effects are aliased and rest on the genomic covariance ",
          "with better-connected relatives rather than on their own plots.",
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

G_oat_all <- collapse_grm(read_grm(grm_files[["oat"]]), analysis_name_files[["oat"]])
G_pea_all <- collapse_grm(read_grm(grm_files[["pea"]]), analysis_name_files[["pea"]])

# The GRMs must cover every phenotyped accession. fit_producer_associate()
# checks this too, but checking here first means the message names the file to
# fix rather than arriving from inside a fitting function.
missing_oat <- setdiff(unique(grainWgt$oatAcc), rownames(G_oat_all))
missing_pea <- setdiff(unique(grainWgt$peaAcc), rownames(G_pea_all))

if (length(missing_oat) > 0 || length(missing_pea) > 0) {
  stop("Accessions phenotyped but absent from the GRMs -- oat: ",
       length(missing_oat), ", pea: ", length(missing_pea), "\n  ",
       paste(utils::head(c(missing_oat, missing_pea), 10), collapse = ", "),
       "\n  Re-run code/create_GRMs_T3.R, or check the analysis names.",
       call. = FALSE)
}

for (nm in c("oat", "pea")) {
  accs <- sort(unique(grainWgt[[paste0(nm, "Acc")]]))
  G <- if (nm == "oat") G_oat_all[accs, accs] else G_pea_all[accs, accs]
  cat("\nG_", nm, ": ", nrow(G), " x ", ncol(G), "\n", sep = "")
  cat("  diagonal      : "); print(summary(diag(G)))
  cat("  symmetric     : ", isSymmetric(unname(G)), "\n", sep = "")
  cat("  min eigenvalue: ",
      min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), "\n", sep = "")
}

# ------------------------------------------------------------
# Fit
#
# Through fit_producer_associate(), which the cross-validation folds in
# code/validate_crossval.R also call. This script used to build its own ETA and
# call BGLR::Multitrait() directly, so VALIDATION_DESIGN.md's claim that "the
# fitting itself is shared" was not true of the production fit. It is now, and
# BGLR_fit_provenance.csv records which function ran so the claim stays
# checkable.
#
# save_effects streams the per-draw coefficients of the two genetic terms, which
# is the only correct route to a PER-ACCESSION prediction error variance:
# SD.beta is the elementwise posterior SD on the L basis, and Var(L b) needs the
# full p x p posterior covariance of b, which BGLR does not accumulate.
# ------------------------------------------------------------

seed_prefix <- function(s) file.path(out_dir, paste0("BGLR_seed_", s, "_"))

fit_one <- function(seed_value) {
  message("\n--- fitting, seed ", seed_value, " ---")
  fit_producer_associate(
    grainWgt, G_oat_all, G_pea_all,
    seed         = seed_value,
    nIter        = nIter,
    burnIn       = burnIn,
    thin         = thin,
    fit_mix_term = fit_mix_term,
    saveAt       = seed_prefix(seed_value),
    save_effects = c("G_oat", "G_pea"),
    storage_mode = draw_storage_mode,
    verbose      = FALSE
  )
}

fp <- rlang::set_names(purrr::map(seeds, fit_one), seeds)

# The raw BGLR objects, so the reporting below is unchanged by the refactor.
fits  <- purrr::map(fp, "fit")
# L is a deterministic function of the GRM, so it is identical across seeds.
L_oat <- fp[[1]]$L$oat
L_pea <- fp[[1]]$L$pea
Y     <- as.matrix(grainWgt[, DGE_IGE_TRAITS])
ETA   <- fits[[1]]$ETA

stopifnot(
  "the seeds disagree about which trials were fitted" =
    length(unique(purrr::map(fp, "trials"))) == 1L,
  "the fitted trial set is not the QC-kept set" =
    setequal(fp[[1]]$trials, levels(grainWgt$trialF))
)

# ------------------------------------------------------------
# What the fit actually saw
#
# Written out because nothing recorded it before, which is why a five-trial fit
# could be reported as nine for two vintages. validation_inputs() reads
# BGLR_fit_trials.csv and refuses to proceed if it disagrees with the
# QC-filtered phenotype table.
# ------------------------------------------------------------

fit_trials <- grainWgt |>
  dplyr::group_by(studyName, studyYear) |>
  dplyr::summarise(
    n_plots        = dplyr::n(),
    n_oat          = dplyr::n_distinct(oatAcc),
    n_pea          = dplyr::n_distinct(peaAcc),
    n_combinations = dplyr::n_distinct(mixID),
    .groups        = "drop"
  ) |>
  dplyr::arrange(studyName)

readr::write_csv(fit_trials, file.path(out_dir, "BGLR_fit_trials.csv"))

# Timestamps are taken BEFORE the tibble, because tibble() evaluates its
# arguments in order and `pheno_file = basename(pheno_file)` would rebind the
# name before file.mtime() saw it -- which silently wrote NA.
qc_path     <- file.path(out_dir, "trial_qc.csv")
pheno_mtime <- format(file.mtime(pheno_file))
qc_mtime    <- if (file.exists(qc_path)) format(file.mtime(qc_path)) else NA_character_
stopifnot("the phenotype file has no modification time" = !is.na(pheno_mtime))

readr::write_csv(
  tibble::tibble(
    fit_date   = format(Sys.Date()),
    fit_fn     = "fit_producer_associate",
    pheno_file = basename(pheno_file),
    pheno_mtime = pheno_mtime,
    qc_file    = "trial_qc.csv",
    qc_mtime   = qc_mtime,
    n_trials   = nrow(fit_trials),
    n_plots    = nrow(grainWgt),
    n_oat      = dplyr::n_distinct(grainWgt$oatAcc),
    n_pea      = dplyr::n_distinct(grainWgt$peaAcc),
    nIter = nIter, burnIn = burnIn, thin = thin,
    seeds      = paste(seeds, collapse = ","),
    fit_mix_term = fit_mix_term,
    draw_storage_mode = draw_storage_mode
  ),
  file.path(out_dir, "BGLR_fit_provenance.csv")
)

cat("\n=== Trials in the fit ===\n")
print(as.data.frame(fit_trials), row.names = FALSE)

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

genetic_terms <- intersect(c("G_oat", "G_pea", "G_mix"), names(ETA))

varcomp <- purrr::imap(fits, \(fit, seed) {
  purrr::map(genetic_terms, \(tm) covariance_components(fit, tm, colnames(Y), effect_roles[[tm]])) |>
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

species_effects <- function(term, L, id_name) {
  purrr::imap(fits, \(fit, seed) {
    accession_effects(fit, term, L, colnames(Y), effect_roles[[term]]) |>
      dplyr::mutate(seed = as.integer(seed), .before = 1)
  }) |>
    purrr::list_rbind() |>
    dplyr::rename({{ id_name }} := accession)
}

oat_all_seeds <- species_effects("G_oat", L_oat, "oatAcc")
pea_all_seeds <- species_effects("G_pea", L_pea, "peaAcc")

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
# Per-accession prediction error variance
#
# WHY THE STREAMED DRAWS ARE THE ONLY ROUTE. The returned fit carries
# `SD.beta`, the elementwise posterior SD of the coefficients on the L basis.
# That is not enough: an accession effect is `L %*% beta`, so its variance is
# `L Var(beta) L'` and needs the full p x p posterior covariance of beta, which
# BGLR does not accumulate. Streaming each draw and rotating it to the
# accession scale is exact.
#
# WHY IT MATTERS. Everything downstream used ONE reliability per species,
# derived from var(BLUP)/sigma2. That gives an accession seen with twenty
# partners the same prediction error as one seen with a single partner, and the
# single-partner accessions are exactly the ones whose producer and associate
# effects are aliased. A per-accession PEV turns the validation pools'
# eligibility rule from a partner count into a reliability threshold.
# ------------------------------------------------------------

acc_pev <- function(seed_value, term, L) {
  roles <- DGE_IGE_ROLES[[term]]
  draws <- read_beta_draws(seed_prefix(seed_value), term = paste0("ETA_", term),
                           nIter = nIter, burnIn = burnIn, thin = thin,
                           p = ncol(L), traits = length(DGE_IGE_TRAITS),
                           storage_mode = draw_storage_mode)

  # Rotate every draw to the accession scale at once: [n_draws x n_acc] per
  # trait. tcrossprod(draws[, , k], L) is draws %*% t(L), and L's rows are the
  # accessions, so column j of the result is accession j across draws.
  g <- purrr::map(seq_along(DGE_IGE_TRAITS),
                  \(k) tcrossprod(draws[, , k], L)) |>
    rlang::set_names(DGE_IGE_TRAITS)

  Pr <- g[[roles[["Pr"]]]]
  As <- g[[roles[["As"]]]]

  # A cheap consistency check against BGLR's own posterior mean. This is the
  # same identity tests/test_fits.R asserts at 1e-10 with double storage; here
  # the draws are single-precision, so the tolerance is relative and loose. It
  # catches a wrong burn-in offset, a wrong p, or a byte-order mistake, all of
  # which would quietly corrupt every PEV.
  beta_mean <- L %*% fits[[as.character(seed_value)]]$ETA[[term]]$beta
  recovered <- cbind(colMeans(g[[DGE_IGE_TRAITS[1]]]),
                     colMeans(g[[DGE_IGE_TRAITS[2]]]))
  rel <- max(abs(recovered - beta_mean)) / max(abs(beta_mean))
  if (!is.finite(rel) || rel > 1e-3) {
    stop("seed ", seed_value, ", ", term, ": the streamed draws average to ",
         "something other than BGLR's posterior mean (relative error ",
         signif(rel, 3), "). Check burn-in, p and storage mode before ",
         "trusting any PEV.", call. = FALSE)
  }

  n   <- nrow(Pr)
  Prc <- sweep(Pr, 2, colMeans(Pr))
  Asc <- sweep(As, 2, colMeans(As))

  tibble::tibble(
    seed     = as.integer(seed_value),
    acc      = rownames(L),
    PEV_Pr   = colSums(Prc^2) / (n - 1),
    PEV_As   = colSums(Asc^2) / (n - 1),
    cov_PrAs = colSums(Prc * Asc) / (n - 1),
    # GMA = Pr + As, so its PEV carries the covariance. This is what the
    # total-yield estimand needs and it is free here.
    PEV_GMA  = PEV_Pr + PEV_As + 2 * cov_PrAs,
    n_draws  = n
  )
}

pev <- purrr::map(c("oat", "pea"), \(sp) {
  term <- paste0("G_", sp)
  L    <- if (sp == "oat") L_oat else L_pea
  out  <- purrr::map(seeds, \(s) acc_pev(s, term, L)) |> purrr::list_rbind()
  readr::write_csv(out, file.path(out_dir, paste0("BGLR_", sp, "_pev.csv")))
  out
}) |> rlang::set_names(c("oat", "pea"))

cat("\n=== Per-accession PEV, averaged over chains ===\n")
purrr::imap(pev, \(x, sp) {
  m <- x |>
    dplyr::group_by(acc) |>
    dplyr::summarise(PEV_As = mean(PEV_As), PEV_GMA = mean(PEV_GMA),
                     .groups = "drop")
  tibble::tibble(
    species = sp, n_acc = nrow(m),
    PEV_As_min = min(m$PEV_As), PEV_As_median = stats::median(m$PEV_As),
    PEV_As_max = max(m$PEV_As),
    fold_range = max(m$PEV_As) / min(m$PEV_As)
  )
}) |>
  purrr::list_rbind() |>
  as.data.frame() |>
  print(digits = 4, row.names = FALSE)

cat("\n  fold_range is the point: a single global PEV would give every one of\n",
    "  these accessions the middle value.\n", sep = "")

# The accession-scale associate-effect draws for the first seed, kept so the
# validation chain can compute the exact posterior variance of a pool contrast
# -- including the posterior COVARIANCE between pool members, which a
# per-accession PEV alone cannot give.
for (sp in c("oat", "pea")) {
  term <- paste0("G_", sp)
  L    <- if (sp == "oat") L_oat else L_pea
  d    <- read_beta_draws(seed_prefix(seeds[1]), term = paste0("ETA_", term),
                          nIter = nIter, burnIn = burnIn, thin = thin,
                          p = ncol(L), traits = length(DGE_IGE_TRAITS),
                          storage_mode = draw_storage_mode)
  As <- tcrossprod(d[, , match(DGE_IGE_ROLES[[term]][["As"]], DGE_IGE_TRAITS)], L)
  colnames(As) <- rownames(L)
  saveRDS(As, file.path(out_dir, paste0("BGLR_", sp, "_As_draws_seed",
                                        seeds[1], ".rds")))
}


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
