# ============================================================
# BIVARIATE DGE-IGE MODEL FOR THE PILOT OAT-PEA EXPERIMENT
#
# The same model code/BGLR_multi_trait_model.R fits to B4I, pointed at the
# pilot: the shared machinery is in code/dge_ige_functions.R and only the
# settings below differ.
#
# One substantive difference, forced by the data rather than chosen. B4I has
# 91% of its oat x pea combinations in a single plot, so its
# specific-combination term is confounded with the residual and is left out.
# The pilot has 48 oats x 12 peas with 86% of the matrix filled and half the
# combinations replicated, so that term is estimable here and is fitted. This
# is the first fit in the project that includes it.
#
# Outputs: output/pilot/pilot_variance_components.csv
#          output/pilot/pilot_{oat,pea}_effects.csv
#          output/pilot/pilot_{oat,pea}_Pr_vs_As.png
#          output/pilot/pilot_fit.rds
# ============================================================

library(tidyverse)

here::i_am("code/pilot/pilot_dge_ige.R")

source(here::here("code", "dge_ige_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

out_dir    <- here::here("output", "pilot")
pheno_file <- file.path(out_dir, "pilot_pheno.rds")

grm_files <- c(oat = file.path(out_dir, "GRM_pilot_oat.rds"),
               pea = file.path(out_dir, "GRM_pilot_pea.rds"))

# The pilot has not been curated for clonal accessions; its entries are
# largely released cultivars rather than the breeding lines where the
# crossing-nursery problem showed up. collapse_grm() no-ops on a missing file,
# so pointing these at real mappings is all it would take.
analysis_name_files <- c(oat = file.path(out_dir, "pilot_oat_analysis_names.csv"),
                         pea = file.path(out_dir, "pilot_pea_analysis_names.csv"))

seeds  <- c(12567, 129, 456, 789)
nIter  <- 20000
burnIn <- 3000
thin   <- 10

# Estimable here; see the header
fit_mix_term <- TRUE

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------

pheno <- readRDS(pheno_file)

grainWgt <- pheno |>
  dplyr::mutate(
    oatYield     = oat_yield,
    peaYield     = pea_yield,
    oatAcc       = as.character(germplasmName),
    peaAcc       = as.character(intercropGermplasmName),
    mixID        = paste(oatAcc, peaAcc, sep = "::"),
    trialF       = factor(studyName),
    blockNumberF = factor(paste(studyName, blockNumber))
  ) |>
  dplyr::filter(!is.na(oatYield), !is.na(peaYield)) |>
  dplyr::mutate(trialF = droplevels(trialF),
                blockNumberF = droplevels(blockNumberF))

message("plots: ", nrow(grainWgt),
        " | trials: ", nlevels(grainWgt$trialF),
        " | blocks: ", nlevels(grainWgt$blockNumberF),
        " | oat: ", dplyr::n_distinct(grainWgt$oatAcc),
        " | pea: ", dplyr::n_distinct(grainWgt$peaAcc),
        " | combinations: ", dplyr::n_distinct(grainWgt$mixID))

# --- design diagnostics, as for B4I ---
partners_oat <- dplyr::count(dplyr::distinct(grainWgt, oatAcc, peaAcc), oatAcc, name = "n")
partners_pea <- dplyr::count(dplyr::distinct(grainWgt, peaAcc, oatAcc), peaAcc, name = "n")
combo_reps   <- as.integer(table(grainWgt$mixID))

message("pea partners per oat: median ", stats::median(partners_oat$n),
        ", range ", paste(range(partners_oat$n), collapse = "-"))
message("oat partners per pea: median ", stats::median(partners_pea$n),
        ", range ", paste(range(partners_pea$n), collapse = "-"))
message("plots per combination: median ", stats::median(combo_reps),
        "; singletons ", round(100 * mean(combo_reps == 1)), "%")

if (mean(combo_reps == 1) > 0.9) {
  warning("more than 90% of combinations are unreplicated; the ",
          "specific-combination term will be confounded with the residual",
          call. = FALSE)
}

# ------------------------------------------------------------
# Relationship matrices
# ------------------------------------------------------------

G_oat_all <- collapse_grm(read_grm(grm_files[["oat"]]), analysis_name_files[["oat"]])
G_pea_all <- collapse_grm(read_grm(grm_files[["pea"]]), analysis_name_files[["pea"]])

missing_oat <- setdiff(unique(grainWgt$oatAcc), rownames(G_oat_all))
missing_pea <- setdiff(unique(grainWgt$peaAcc), rownames(G_pea_all))
if (length(missing_oat) + length(missing_pea) > 0) {
  message("dropping ", length(missing_oat), " oat and ", length(missing_pea),
          " pea accession(s) with no GRM row")
  grainWgt <- dplyr::filter(grainWgt, !oatAcc %in% missing_oat,
                            !peaAcc %in% missing_pea)
}

oatAccs <- sort(unique(grainWgt$oatAcc))
peaAccs <- sort(unique(grainWgt$peaAcc))
mixIDs  <- sort(unique(grainWgt$mixID))

G_oat <- G_oat_all[oatAccs, oatAccs, drop = FALSE]
G_pea <- G_pea_all[peaAccs, peaAccs, drop = FALSE]

# The combination kernel is G_oat (x) G_pea on the observed combinations
mix_parts <- stringr::str_split_fixed(mixIDs, stringr::fixed("::"), 2)
G_mix <- G_oat[mix_parts[, 1], mix_parts[, 1], drop = FALSE] *
         G_pea[mix_parts[, 2], mix_parts[, 2], drop = FALSE]
dimnames(G_mix) <- list(mixIDs, mixIDs)

for (nm in c("G_oat", "G_pea", "G_mix")) {
  G <- get(nm)
  cat("\n", nm, ": ", nrow(G), " x ", ncol(G),
      "  min eigenvalue ",
      signif(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), 3),
      "\n", sep = "")
}

# ------------------------------------------------------------
# Model matrices
# ------------------------------------------------------------

Y <- as.matrix(grainWgt[, c("peaYield", "oatYield")])
colnames(Y) <- c("peaYield", "oatYield")
stopifnot(!anyNA(Y))

# Full trial dummies with no separate intercept: a full dummy set plus an
# intercept is rank-deficient and BGLR samples along the ridge rather than
# erroring.
incTrials <- stats::model.matrix(~ 0 + trialF, grainWgt)
colnames(incTrials) <- levels(grainWgt$trialF)

# Every pilot trial has four blocks, so unlike B4I none of them is aliased
# with its own trial effect and all are fitted.
incBlocks <- stats::model.matrix(~ 0 + blockNumberF, grainWgt)
colnames(incBlocks) <- levels(grainWgt$blockNumberF)

Z_oat <- incidence(grainWgt$oatAcc, rownames(G_oat))
Z_pea <- incidence(grainWgt$peaAcc, rownames(G_pea))
Z_mix <- incidence(grainWgt$mixID,  rownames(G_mix))

L_oat <- grm_factor(G_oat)
L_pea <- grm_factor(G_pea)
L_mix <- grm_factor(G_mix)

ETA <- list(
  trial = list(X = incTrials,       model = "FIXED"),
  block = list(X = incBlocks,       model = "BRR"),
  G_pea = list(X = Z_pea %*% L_pea, model = "BRR"),
  G_oat = list(X = Z_oat %*% L_oat, model = "BRR")
)
if (fit_mix_term) ETA$G_mix <- list(X = Z_mix %*% L_mix, model = "BRR")

# ------------------------------------------------------------
# Fit
# ------------------------------------------------------------

fit_one <- function(seed_value) {
  message("\n--- fitting, seed ", seed_value, " ---")
  set.seed(seed_value)
  BGLR::Multitrait(
    y = Y, ETA = ETA, intercept = FALSE,
    resCov = list(df0 = 4, S0 = NULL, type = "UN"),
    nIter = nIter, burnIn = burnIn, thin = thin,
    saveAt = file.path(out_dir, paste0("pilot_seed_", seed_value, "_")),
    verbose = FALSE
  )
}

fits <- rlang::set_names(purrr::map(seeds, fit_one), seeds)
saveRDS(fits[[1]], file.path(out_dir, "pilot_fit.rds"))

# ------------------------------------------------------------
# Variance components
# ------------------------------------------------------------

effect_roles <- list(
  G_oat = c(Pr = "oatYield", As = "peaYield"),
  G_pea = c(Pr = "peaYield", As = "oatYield"),
  G_mix = c(Pr = "oatYield", As = "peaYield")
)

genetic_terms <- intersect(c("G_oat", "G_pea", "G_mix"), names(ETA))

varcomp <- purrr::imap(fits, \(fit, seed) {
  purrr::map(genetic_terms,
             \(tm) covariance_components(fit, tm, colnames(Y), effect_roles[[tm]])) |>
    purrr::list_rbind() |>
    dplyr::mutate(seed = as.integer(seed), .before = 1)
}) |>
  purrr::list_rbind()

resid_cov <- purrr::imap(fits, \(fit, seed) {
  R <- fit$resCov$R
  dimnames(R) <- list(colnames(Y), colnames(Y))
  tibble::tibble(seed = as.integer(seed),
                 var_pea = R["peaYield", "peaYield"],
                 var_oat = R["oatYield", "oatYield"],
                 cov_pea_oat = R["peaYield", "oatYield"],
                 cor_pea_oat = cov_pea_oat / sqrt(var_pea * var_oat))
}) |>
  purrr::list_rbind()

cat("\n=== Variance components (Pr = producer, As = associate) ===\n")
print(as.data.frame(varcomp), row.names = FALSE, digits = 4)
cat("\n=== Residual covariance between the two yields on a plot ===\n")
print(as.data.frame(resid_cov), row.names = FALSE, digits = 4)

readr::write_csv(
  dplyr::bind_rows(dplyr::mutate(varcomp, component = "genetic"),
                   dplyr::mutate(resid_cov, term = "residual", component = "residual")),
  file.path(out_dir, "pilot_variance_components.csv")
)

# --- credible intervals from the saved samples ---
n_burn_rows <- burnIn / thin
omega_file <- function(seed, j) file.path(out_dir, sprintf("pilot_seed_%d_Omega_%d.dat", seed, j))

cat("\n=== 95% credible intervals, pooled over chains ===\n")
term_index <- stats::setNames(seq_along(ETA), names(ETA))
purrr::iwalk(effect_roles[genetic_terms], \(roles, tm) {
  files <- purrr::map_chr(seeds, \(s) omega_file(s, term_index[[tm]]))
  if (!all(file.exists(files))) return(invisible(NULL))
  S <- purrr::map(files, \(f) {
    m <- as.matrix(utils::read.table(f))
    m[-seq_len(n_burn_rows), , drop = FALSE]
  }) |> purrr::reduce(rbind)

  # A 2x2 unstructured covariance is stored as its three unique elements
  pr_first <- roles[["Pr"]] == colnames(Y)[1]
  v_pr <- if (pr_first) S[, 1] else S[, ncol(S)]
  v_as <- if (pr_first) S[, ncol(S)] else S[, 1]
  cv   <- S[, 2]
  q <- \(x) round(stats::quantile(x, c(.025, .5, .975)), 3)

  cat("\n", tm, "\n", sep = "")
  cat("  var Pr   : ", paste(q(v_pr), collapse = "  "), "\n", sep = "")
  cat("  var As   : ", paste(q(v_as), collapse = "  "), "\n", sep = "")
  cat("  cov Pr,As: ", paste(q(cv), collapse = "  "),
      "   P(cov < 0) = ", round(mean(cv < 0), 3),
      "   excludes zero: ",
      ifelse(stats::quantile(cv, .025) > 0 | stats::quantile(cv, .975) < 0, "YES", "no"),
      "\n", sep = "")
})

R_files <- purrr::map_chr(seeds, \(s) file.path(out_dir, sprintf("pilot_seed_%d_R.dat", s)))
if (all(file.exists(R_files))) {
  R <- purrr::map(R_files, \(f) {
    m <- as.matrix(utils::read.table(f)); m[-seq_len(n_burn_rows), , drop = FALSE]
  }) |> purrr::reduce(rbind)
  q <- round(stats::quantile(R[, 2], c(.025, .5, .975)), 3)
  cat("\nresidual cov oat,pea: ", paste(q, collapse = "  "),
      "   excludes zero: ", ifelse(q[1] > 0 | q[3] < 0, "YES", "no"), "\n", sep = "")
}

# ------------------------------------------------------------
# Accession effects
# ------------------------------------------------------------

species_effects <- function(term, L, id_name) {
  purrr::imap(fits, \(fit, seed) {
    accession_effects(fit, term, L, colnames(Y), effect_roles[[term]]) |>
      dplyr::mutate(seed = as.integer(seed), .before = 1)
  }) |>
    purrr::list_rbind() |>
    dplyr::rename({{ id_name }} := accession)
}

oat_eff <- species_effects("G_oat", L_oat, "oatAcc")
pea_eff <- species_effects("G_pea", L_pea, "peaAcc")

oat_rank <- rank_summary(oat_eff, "oatAcc")
pea_rank <- rank_summary(pea_eff, "peaAcc")

cat("\n=== Oat accessions by mixing ability ===\n")
print(head(oat_rank, 15))
cat("\n=== Pea accessions by mixing ability (all 12) ===\n")
print(pea_rank, n = Inf)

readr::write_csv(oat_rank, file.path(out_dir, "pilot_oat_effects.csv"))
readr::write_csv(pea_rank, file.path(out_dir, "pilot_pea_effects.csv"))

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

pr_as_plot <- function(effects, species) {
  labels <- effects |>
    dplyr::group_by(seed) |>
    dplyr::summarise(r = stats::cor(PrEff, AsEff), .groups = "drop") |>
    dplyr::mutate(label = paste0("r = ", round(r, 3)))

  ggplot2::ggplot(effects, ggplot2::aes(PrEff, AsEff)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "red") +
    ggplot2::geom_text(data = labels,
                       ggplot2::aes(x = Inf, y = Inf, label = label),
                       inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 4) +
    ggplot2::facet_wrap(~ seed) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(title = paste(species, "producer vs associate effects, pilot"),
                  x = "Producer effect", y = "Associate effect")
}

ggplot2::ggsave(file.path(out_dir, "pilot_oat_Pr_vs_As.png"),
                pr_as_plot(oat_eff, "Oat"), width = 8, height = 6, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "pilot_pea_Pr_vs_As.png"),
                pr_as_plot(pea_eff, "Pea"), width = 8, height = 6, dpi = 150)

message("\nwrote pilot outputs to output/pilot/")
