# ============================================================
# CAN THE DGE-IGE INTERACTION BE READ AS SCORES AND LOADINGS?
#
#   Rscript code/sim_decomp_run.R --check          # positive control; run first
#   Rscript code/sim_decomp_run.R --pilot          # two cheap cells, end to end
#   Rscript code/sim_decomp_run.R                  # the design
#   Rscript code/sim_decomp_run.R --task 3 --ntasks 20
#   Rscript code/sim_decomp_run.R --combine        # rebuild the CSVs from cache
#
# THE QUESTION. `dge_ige` predicts the interaction better than MegaLMM when
# combinations are sparse and -- measured at full rank in
# output/simulation_kron_study.csv -- also when the interaction is higher-rank.
# But it returns a number per combination, not a score for an oat and a loading
# for a pea. Its fitted surface is nevertheless an exact bilinear form, so one
# SVD per MCMC draw recovers the interpretable object. This design asks whether
# what comes back is the TRUTH or an artefact of the kinship basis.
#
# WHAT IS SCORED, AND AGAINST WHAT.
#
#   share1_post       component 1's share of the interaction, with a 95%
#                     interval. Taken from the per-DRAW spectra, because each
#                     draw is a posterior sample of the true surface. The
#                     corresponding figure for the posterior-MEAN surface is
#                     carried as share1_meansurf and runs systematically higher
#                     -- singular values are convex, so averaging draws cancels
#                     their idiosyncratic directions and leaves something that
#                     looks more low-rank than anything real. `share1` is the
#                     mean-surface value for both models, and is the only share
#                     comparable with MegaLMM, which has no streamed draws.
#   participation     (sum d^2)^2 / sum d^4 -- how many components there really
#                     are, against a known n_factors of 1, 3 or 5
#   sub_cor           canonical correlations between the recovered top-n_factors
#                     score subspace and sim$truth$U_oat. Subspaces, not
#                     columns: for n_factors > 1 the truth is identified only up
#                     to rotation, so a column-wise correlation would understate
#                     recovery for reasons that have nothing to do with the
#                     model.
#   frac_of_ceiling   sub_cor divided by what the TRUNCATED BASIS allows. At
#                     kron_rank = q the recovered scores are confined to
#                     span(kron_A) while the truth has mass on every direction
#                     of G, so there is a hard analytic ceiling -- computable
#                     with nothing fitted. Recovery that is not expressed as a
#                     fraction of it cannot be read, and the apparent n_acc
#                     effect is partly 30-of-200 versus 30-of-400 directions.
#
# THE NULL CELLS ARE NOT OPTIONAL. BGLR puts one scalar prior variance over all
# kron_rank^2 coefficients, which is misspecified for a low-rank truth and
# should FLATTEN the recovered spectrum. So `participation` has no meaning in
# the abstract, only against the value it takes when there is no interaction at
# all. `--null-cells` (on by default) adds interaction_pct = 0 scenarios to
# measure that floor. If the n_factors = 1 cells come back at the same
# participation as the null, the spectrum carries no rank information and the
# honest deliverable shrinks to the leading direction plus a caveat.
#
# MegaLMM IS SCORED WITH THE SAME OPERATOR, not with its own factors. Its
# factors are not orthogonal -- per-factor shares of its fitted interaction do
# not sum to one -- so "MegaLMM factor 1's share" and "dge_ige component 1's
# share" are not comparable quantities until both surfaces go through
# decompose_surface(). tests/test_decomp.R asserts the non-orthogonality, so the
# reason for this is on the record rather than a matter of taste.
#
# Outputs: output/simulation_decomp_results.csv    one row per cell x trait x model
#          output/simulation_decomp_spectrum.csv   one row per component
#          output/simulation_decomp/               per-cell cache, and the
#                                                  streamed draws
# ============================================================

library(tidyverse)

here::i_am("code/sim_decomp_run.R")

source(here::here("code", "dge_ige_functions.R"))
source(here::here("code", "sim_config.R"))
source(here::here("code", "sim_int_config.R"))
source(here::here("code", "sim_generate.R"))
source(here::here("code", "sim_fit.R"))
source(here::here("code", "megalmm_setup.R"))
source(here::here("code", "interaction_decomp.R"))

out_dir    <- here::here("output")
cache_dir  <- here::here("output", "simulation_decomp")
draws_dir  <- file.path(cache_dir, "draws")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
has_flag <- function(flag) flag %in% args

n_reps       <- as.integer(arg_value("--reps", as.character(SIM_INT_REPS)))
kron_rank    <- as.integer(arg_value("--rank", as.character(SIM_KRON_RANK)))
n_iter       <- as.integer(arg_value("--n-iter", "3000"))
burn_in      <- as.integer(arg_value("--burn-in", "600"))
thin         <- as.integer(arg_value("--thin", "10"))
check        <- has_flag("--check")
pilot        <- has_flag("--pilot")
combine_only <- has_flag("--combine")
refresh      <- has_flag("--refresh")
keep_draws   <- !has_flag("--drop-draws")
null_cells   <- !has_flag("--no-null-cells")
with_megalmm <- !has_flag("--no-megalmm")
task         <- as.integer(arg_value("--task", NA))
ntasks       <- as.integer(arg_value("--ntasks", NA))
filter_expr  <- arg_value("--filter", NULL)

# The cache scheme marker, following the rest of the project: a change to what
# is stored has to move it, or bind_rows() over a mixed cache produces both
# schemas with half the rows NA in each.
DECOMP_SCHEME <- "d1"

# The null cells are new scenarios, so they need their own seed block. Clear of
# SIM_BASE_SEED (20260921), SIM_INT_BASE_SEED (20270000) and sim_kron_study's
# 20280000 -- a collision would make two different scenarios share a simulated
# panel and look more alike than they are.
DECOMP_NULL_BASE_SEED <- 20290000L

# ------------------------------------------------------------
# One cell
# ------------------------------------------------------------

#' Recovery of a recovered factor matrix against the simulated truth.
#'
#' Reported three ways because they answer different questions: the leading
#' canonical correlation (is the main direction there), the mean over the first
#' `n_factors` (is the whole subspace there), and both divided by what the
#' truncated basis allows.
score_recovery <- function(recovered, truth, basis, n_factors) {
  if (is.null(truth) || is.null(recovered)) {
    return(tibble::tibble(cor1 = NA_real_, cor_mean = NA_real_,
                          ceil1 = NA_real_, ceil_mean = NA_real_,
                          frac1 = NA_real_, frac_mean = NA_real_))
  }
  k <- min(n_factors, ncol(recovered), ncol(truth))
  got  <- subspace_cors(recovered[, seq_len(k), drop = FALSE],
                        truth[, seq_len(k), drop = FALSE])
  ceil <- recovery_ceiling(truth[, seq_len(k), drop = FALSE], basis)
  ceil <- ceil[seq_len(min(k, length(ceil)))]
  tibble::tibble(
    cor1 = got[1], cor_mean = mean(got),
    ceil1 = ceil[1], ceil_mean = mean(ceil),
    frac1 = got[1] / ceil[1], frac_mean = mean(got) / mean(ceil))
}

#' Decompose one model's fitted interaction, and score what came back.
summarise_model <- function(dec, model, trait, sim, basis, n_factors,
                            draws_summary = NULL, seconds = NA_real_) {
  truth_U   <- if (trait == "oat") sim$truth$U_oat      else sim$truth$U_pea
  truth_Lam <- if (trait == "oat") sim$truth$Lambda_oat else sim$truth$Lambda_pea

  sc <- score_recovery(dec$scores, truth_U, basis$A, n_factors)
  ld <- score_recovery(dec$loadings, truth_Lam, basis$B, n_factors)

  share_at <- function(i) if (i <= length(dec$share)) dec$share[i] else NA_real_
  out <- tibble::tibble(
    model = model, trait = trait, seconds = seconds,
    rank_kept = dec$rank_kept, total_ss = dec$total_ss,
    share1 = share_at(1), share2 = share_at(2), share3 = share_at(3),
    cum_share3 = dec$cum_share[min(3L, length(dec$cum_share))],
    participation = dec$participation, n90 = dec$n90,
    score_cor1 = sc$cor1, score_cor_mean = sc$cor_mean,
    score_ceil1 = sc$ceil1, score_ceil_mean = sc$ceil_mean,
    score_frac1 = sc$frac1, score_frac_mean = sc$frac_mean,
    load_cor1 = ld$cor1, load_cor_mean = ld$cor_mean,
    load_frac1 = ld$frac1, load_frac_mean = ld$frac_mean)

  if (!is.null(draws_summary)) {
    out <- dplyr::bind_cols(out, dplyr::select(
      draws_summary, n_draw, participation_draw_mean, order_stable,
      order_stable_all, rot_diag, gap_ratio,
      share1_post = share1, share1_lo, share1_hi,
      share1_meansurf, cum_share3_meansurf,
      proj1_mean, proj1_lo, proj1_hi))
  }
  out
}

run_cell <- function(row) {
  key <- paste0(row$scenario, "_rep", row$rep, "_s", row$seed, "_",
                DECOMP_SCHEME)
  f <- file.path(cache_dir, paste0(key, ".rds"))
  if (!refresh && file.exists(f)) {
    message("  cached: ", row$scenario, " rep ", row$rep)
    return(readRDS(f))
  }
  message("\n### ", row$scenario, " rep ", row$rep, " ###")

  grms <- sim_grms(row$n_acc, seed = row$n_acc)
  set.seed(row$seed)
  sim <- simulate_experiment(grms$G_oat, grms$G_pea, sparsity = row$sparsity,
                             n_factors = row$n_factors,
                             interaction_pct = row$interaction_pct,
                             n_envs = row$n_envs, gxe_cor = row$gxe_cor)
  sp <- prepare_scenario(sim)

  prefix <- file.path(draws_dir, paste0(key, "_"))
  t0 <- Sys.time()
  fd <- fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, with_interaction = TRUE,
                    kron_rank = kron_rank, nIter = n_iter, burnIn = burn_in,
                    seed = row$seed, save_effects = TRUE, saveAt = prefix,
                    storage_mode = "single")
  secs_bglr <- as.numeric(difftime(Sys.time(), t0, units = "s"))

  ka <- ncol(fd$kron$A); kb <- ncol(fd$kron$B)
  draws <- read_beta_draws(prefix, nIter = n_iter, burnIn = burn_in,
                           thin = thin, p = ka * kb, traits = 2L,
                           storage_mode = "single")

  rows <- list(); spec <- list()
  for (trait in c("oat", "pea")) {
    ti <- match(if (trait == "oat") "oatYield" else "peaYield",
                fd$kron$traits)

    # Guard, not decoration: the additive part of the dge_ige surface is an
    # outer sum, whose interaction_part is exactly zero. So decomposing the
    # bilinear term must give the same spectrum as decomposing the whole
    # surface. If it does not, the padding or the reshape is wrong.
    dec_b <- decompose_bilinear(fd$kron$A,
                                beta_from_vector(fd$beta[, ti], ka, kb),
                                fd$kron$B)
    dec_s <- decompose_surface(fd$surface[[trait]])
    k_cmp <- min(5L, length(dec_b$d), length(dec_s$d))
    if (max(abs(dec_b$d[seq_len(k_cmp)] - dec_s$d[seq_len(k_cmp)])) > 1e-6) {
      stop("the bilinear decomposition disagrees with the surface ",
           "decomposition for ", trait, ": the reshape or the padding is wrong",
           call. = FALSE)
    }

    dd <- decompose_draws(fd$kron$A, draws, fd$kron$B, trait = ti,
                          mean_beta = fd$beta[, ti])
    cs <- component_summary(dd)
    rows[[length(rows) + 1L]] <- summarise_model(
      dd$reference, "dge_ige", trait, sim, fd$kron, row$n_factors,
      draws_summary = cs$summary, seconds = secs_bglr)
    spec[[length(spec) + 1L]] <- dplyr::mutate(
      cs$components, model = "dge_ige", trait = trait, .before = 1)
  }

  if (with_megalmm) {
    t0 <- Sys.time()
    rd <- file.path(tempdir(), paste0("decomp_", key))
    dir.create(rd, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(rd, recursive = TRUE), add = TRUE)
    mm <- fit_megalmm_both(sp$train, sim$G_oat, sim$G_pea, runID = rd,
                           K = SIM_INT_MEGALMM$K,
                           eigen_variance = SIM_INT_MEGALMM$eigen_variance,
                           fixed_main_effect = FALSE)
    secs_mm <- as.numeric(difftime(Sys.time(), t0, units = "s"))
    for (trait in c("oat", "pea")) {
      dec <- decompose_surface(mm$surface[[trait]])
      rows[[length(rows) + 1L]] <- summarise_model(
        dec, "megalmm", trait, sim,
        # MegaLMM is not basis-truncated, so its ceiling is the full span: the
        # identity of the right dimension, which recovery_ceiling() reads as 1.
        list(A = diag(nrow(sim$G_oat)), B = diag(nrow(sim$G_pea))),
        row$n_factors, seconds = secs_mm)
      spec[[length(spec) + 1L]] <- tibble::tibble(
        model = "megalmm", trait = trait,
        component = seq_along(dec$share),
        share_mean = NA_real_, share_sd = NA_real_,
        share_lo = NA_real_, share_hi = NA_real_, d_mean = dec$d,
        share_ref = dec$share, cum_share_ref = dec$cum_share)
    }
  }

  design_cols <- dplyr::select(row, scenario, rep, seed, n_acc, sparsity,
                               n_factors, interaction_pct, n_envs, gxe_cor)
  out <- list(
    results = dplyr::bind_cols(
      design_cols[rep(1, length(rows)), ],
      purrr::list_rbind(rows)) |>
      dplyr::mutate(kron_rank = kron_rank,
                    n_combos = dplyr::n_distinct(paste(sim$obs$oat,
                                                       sim$obs$pea))),
    spectrum = dplyr::bind_cols(
      design_cols[rep(1, nrow(purrr::list_rbind(spec))), ],
      purrr::list_rbind(spec)))

  saveRDS(out, f)
  if (!keep_draws) unlink(paste0(prefix, "ETA_G_mix_beta.bin"))
  out
}

# ------------------------------------------------------------
# Positive control
#
# A dense, rank-1 scenario is the easiest case there is: one true component, half
# the combinations observed. If the leading component does not dominate there,
# and does not look like the simulated scores, the decomposition is mis-wired
# and every number the design reports is an artefact. A negative result is only
# worth reporting if the pipeline can produce a positive one.
# ------------------------------------------------------------

if (check) {
  dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)
  message("dense rank-1 scenario: component 1 should dominate and match truth")
  grms <- sim_grms(120L, seed = 120L)
  set.seed(1L)
  sim <- simulate_experiment(grms$G_oat, grms$G_pea, sparsity = 0.48,
                             n_factors = 1L, interaction_pct = 0.20,
                             n_envs = 1L)
  sp <- prepare_scenario(sim)
  prefix <- file.path(draws_dir, "check_")
  fd <- fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, TRUE,
                    kron_rank = kron_rank, nIter = 1500L, burnIn = 300L,
                    seed = 1L, save_effects = TRUE, saveAt = prefix,
                    storage_mode = "single")
  ka <- ncol(fd$kron$A); kb <- ncol(fd$kron$B)
  draws <- read_beta_draws(prefix, nIter = 1500L, burnIn = 300L, thin = 10L,
                           p = ka * kb, traits = 2L, storage_mode = "single")
  ti <- match("oatYield", fd$kron$traits)

  check_mean <- max(abs(colMeans(draws[, , ti]) - fd$beta[, ti]))
  dd <- decompose_draws(fd$kron$A, draws, fd$kron$B, trait = ti,
                        mean_beta = fd$beta[, ti])
  cs <- component_summary(dd)
  ceil <- recovery_ceiling(sim$truth$U_oat, fd$kron$A)[1]
  got  <- subspace_cors(dd$reference$scores[, 1, drop = FALSE],
                        sim$truth$U_oat)[1]

  mmshare <- NA_real_
  mm_ok <- TRUE
  if (with_megalmm) {
    rd <- file.path(tempdir(), "decomp_check")
    dir.create(rd, showWarnings = FALSE, recursive = TRUE)
    mm <- fit_megalmm_both(sp$train, sim$G_oat, sim$G_pea, runID = rd,
                           K = SIM_INT_MEGALMM$K,
                           eigen_variance = SIM_INT_MEGALMM$eigen_variance,
                           fixed_main_effect = FALSE)
    mmshare <- sum(decompose_surface(mm$surface$oat)$share)
    mm_ok <- abs(mmshare - 1) < 1e-8
  }

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf("draw mean vs posterior mean   %.3g  (must be ~0)\n", check_mean))
  cat(sprintf("share of component 1          %.3f  [95%% %.3f - %.3f]\n",
              cs$summary$share1, cs$summary$share1_lo, cs$summary$share1_hi))
  cat(sprintf("  mean-surface summary        %.3f  (overstates concentration)\n",
              cs$summary$share1_meansurf))
  cat(sprintf("  draw variance in that dir.  %.3f  (identifiability)\n",
              cs$summary$proj1_mean))
  cat(sprintf("participation ratio           %.2f (ref) / %.2f (draws)\n",
              cs$summary$participation_ref, cs$summary$participation_draw_mean))
  cat(sprintf("order stability (top 3)       %.2f   rot_diag %.2f\n",
              cs$summary$order_stable, cs$summary$rot_diag))
  cat(sprintf("score vs truth                %.3f of a %.3f ceiling = %.2f\n",
              got, ceil, got / ceil))
  if (with_megalmm) {
    cat(sprintf("MegaLMM shares sum to         %.6f  (same operator)\n", mmshare))
  }
  cat(strrep("=", 70), "\n", sep = "")

  # Threshold set for SINGLE-precision draws, which is what the sweep writes:
  # the round trip is then exact to ~1e-9, not to 1e-16. tests/test_fits.R
  # asserts the double-precision version at 1e-10.
  if (check_mean > 1e-7) {
    stop("the streamed draws do not average to BGLR's posterior mean: the ",
         "burn-in offset or the storage mode is wrong", call. = FALSE)
  }
  if (cs$summary$share1 < 0.3) {
    stop(sprintf("component 1 carries only %.3f of a rank-1 interaction. ",
                 cs$summary$share1),
         "The decomposition is mis-wired.", call. = FALSE)
  }
  if (!is.finite(got / ceil) || got / ceil < 0.5) {
    stop(sprintf("recovered score reaches %.2f of its %.2f ceiling. ",
                 got, ceil),
         "The scores do not match the simulated truth.", call. = FALSE)
  }
  if (!mm_ok) {
    stop("MegaLMM's shares do not sum to 1 under the shared operator",
         call. = FALSE)
  }
  message("\nOK: the decomposition recovers a rank-1 truth.")
  quit(save = "no")
}

# ------------------------------------------------------------
# The design
# ------------------------------------------------------------

scenarios <- sim_int_grid(n_reps = n_reps)

# The null cells: no interaction at all, so the fitted spectrum is pure noise
# and prior. One per n_acc x sparsity, which is all that is needed -- n_factors
# and the GxE levels cannot matter when the interaction variance is zero.
if (null_cells) {
  nulls <- tidyr::expand_grid(
    n_acc = SIM_INT_LEVELS$n_acc, sparsity = SIM_INT_LEVELS$sparsity,
    n_factors = 1L, interaction_pct = 0, n_envs = 1L) |>
    dplyr::mutate(
      gxe_cor = unname(SIM_INT_GXE[as.character(n_envs)]),
      scenario = sprintf("null_n%d_sp%03d", n_acc, round(sparsity * 1000)))
  n_null <- nrow(nulls)
  nulls <- tidyr::expand_grid(rep = seq_len(n_reps), nulls) |>
    dplyr::mutate(seed = DECOMP_NULL_BASE_SEED + (rep - 1L) * n_null +
                    rep(seq_len(n_null), times = n_reps)) |>
    dplyr::relocate(scenario, rep, seed)
  scenarios <- dplyr::bind_rows(scenarios, nulls)
}

message(nrow(scenarios), " cell(s) at kron_rank ", kron_rank,
        if (null_cells) " (including the interaction_pct = 0 nulls)" else "")

if (pilot) {
  scenarios <- scenarios |>
    dplyr::filter(rep == 1, n_acc == 200, sparsity == 0.48) |>
    dplyr::slice_head(n = 2)
}
if (!is.null(filter_expr)) {
  scenarios <- dplyr::filter(scenarios, !!rlang::parse_expr(filter_expr))
}
if (nrow(scenarios) == 0) { message("nothing to run"); quit(save = "no") }

dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

if (!is.na(task) && !is.na(ntasks)) {
  keep <- which((seq_len(nrow(scenarios)) - 1L) %% ntasks == (task - 1L))
  scenarios <- scenarios[keep, ]
  message("task ", task, " of ", ntasks, ": ", nrow(scenarios), " cell(s)")
  if (nrow(scenarios) == 0) quit(save = "no")
}

if (!combine_only) {
  for (i in seq_len(nrow(scenarios))) {
    tryCatch(run_cell(scenarios[i, ]),
             error = function(e) message("  FAILED: ", conditionMessage(e)))
  }
}

# ------------------------------------------------------------
# Combine and report
# ------------------------------------------------------------

files <- list.files(cache_dir, pattern = paste0("_", DECOMP_SCHEME, "\\.rds$"),
                    full.names = TRUE)
if (length(files) == 0) { message("no cached cells"); quit(save = "no") }

cached   <- purrr::map(files, readRDS)
results  <- purrr::map(cached, "results")  |> purrr::list_rbind()
spectrum <- purrr::map(cached, "spectrum") |> purrr::list_rbind()

readr::write_csv(results,  file.path(out_dir, "simulation_decomp_results.csv"))
readr::write_csv(spectrum, file.path(out_dir, "simulation_decomp_spectrum.csv"))

fmt <- function(x, k = 3) formatC(x, format = "f", digits = k)

cat("\n", strrep("=", 78),
    "\nIs the spectrum low-rank, and does it track the truth?",
    "\n(participation against a known n_factors; the nulls are the floor)\n",
    strrep("=", 78), "\n", sep = "")
results |>
  dplyr::filter(model == "dge_ige", trait == "oat") |>
  dplyr::mutate(cell = dplyr::if_else(interaction_pct == 0, "null",
                                      as.character(n_factors))) |>
  dplyr::group_by(n_acc, sparsity, cell) |>
  dplyr::summarise(n = dplyr::n(),
                   share1 = mean(share1_post),
                   part = mean(participation_draw_mean),
                   .groups = "drop") |>
  tidyr::pivot_wider(names_from = cell,
                     values_from = c(share1, part, n)) |>
  print(n = 30, width = 220)

cat("\n", strrep("=", 78),
    "\nDoes the recovered score match the truth, as a fraction of its ceiling?\n",
    strrep("=", 78), "\n", sep = "")
results |>
  dplyr::filter(interaction_pct > 0, trait == "oat") |>
  dplyr::group_by(model, sparsity, n_factors) |>
  dplyr::summarise(cor1 = mean(score_cor1, na.rm = TRUE),
                   ceil1 = mean(score_ceil1, na.rm = TRUE),
                   frac1 = mean(score_frac1, na.rm = TRUE),
                   .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = c(cor1, ceil1, frac1)) |>
  print(n = 40, width = 220)

cat("\n", strrep("=", 78),
    "\nIs 'component 1' a stable object, or only the leading subspace?\n",
    strrep("=", 78), "\n", sep = "")
results |>
  dplyr::filter(model == "dge_ige", trait == "oat", interaction_pct > 0) |>
  dplyr::group_by(sparsity, n_factors) |>
  dplyr::summarise(order_stable = mean(order_stable),
                   rot_diag = mean(rot_diag),
                   share1_width = mean(share1_hi - share1_lo),
                   proj1 = mean(proj1_mean),
                   .groups = "drop") |>
  print(n = 40)

message("\nwrote:\n  output/simulation_decomp_results.csv\n",
        "  output/simulation_decomp_spectrum.csv")
