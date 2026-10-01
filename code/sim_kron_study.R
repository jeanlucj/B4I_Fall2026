# ============================================================
# DOES SIM_KRON_RANK DECIDE THE MegaLMM-vs-DGE-IGE RESULT?
#
#   Rscript code/sim_kron_study.R --check        # one cell, fast
#   Rscript code/sim_kron_study.R                # the design
#   Rscript code/sim_kron_study.R --n-acc 60,100 --ranks 10,20,30,full
#
# THE CONFOUND. docs/BGLR_vs_MegaLMM.md reports MegaLMM beating dge_ige at
# recovering the oat x pea interaction above a density threshold -- by 0.37 at
# rank 1 and 16% observed. But dge_ige's interaction term is RANK-TRUNCATED and
# MegaLMM's is not. At SIM_KRON_RANK = 30 on a 400-accession panel the basis
# carries 62.7% of the oat GRM's trace and 54.8% of the pea's, so about
#
#     0.627 x 0.548 = 34% of the Kronecker covariance, 900 of 160,000 dimensions
#
# and the directions it keeps are chosen by GRM eigenvalue, not by where the
# interaction actually is. So the finding might read, strictly, "MegaLMM beats a
# rank-30-truncated Kronecker kernel".
#
# WHY A RANK SWEEP ALONE CANNOT SETTLE IT. Raising the rank and watching the gap
# narrow tells you the direction but not the destination: you cannot reach full
# rank at n = 400, because the basis is rank^2 columns -- rank 90 is 8,100
# columns on 61,440 plots, and full rank is 160,000. Any conclusion from a sweep
# is an extrapolation, and an extrapolation needs anchoring.
#
# WHAT "FULL RANK" HAS TO MEAN HERE. Not the exact kernel. The exact kernel
# G_oat[i,i'] * G_pea[j,j'] is defined only over OBSERVED combinations, so it
# returns no value for a combination nobody grew -- and the never-observed cells
# are exactly where r_int is scored. Measured: fitting it and asking for r_int
# gives NA, because the returned surface has no interaction component at all.
# Predicting an unobserved combination from it would need a BLUP extension
# (cross-covariance between the new combination and the training ones), which is
# not implemented and is not cheap.
#
# The low-rank basis gives it for free instead. At `kron_rank = n_acc` the two
# bases span their whole spaces, so their row-wise Kronecker product spans the
# ENTIRE n_acc^2-dimensional Kronecker space: trace retained = 1, nothing
# truncated, and still a full oat x pea surface. So `full` in the rank list means
# rank = n_acc, and at small and middling panels it is affordable -- n = 60 is
# 3,600 columns, n = 100 is 10,000.
#
# THE DESIGN. Two parts, and the second is what makes the first interpretable.
#
# 1. ANCHOR: panels small enough to run at FULL rank, so what truncation costs
#    is measured rather than extrapolated.
#
# 2. SCALING: the same sweep at several panel sizes. The thing to learn is
#    whether the truncation cost tracks the ABSOLUTE rank or the RATIO rank/n --
#    because rank 30 is half the directions at n = 60 and 7.5% of them at
#    n = 400. If it tracks the retained trace (which the script records), the
#    anchor at small n licenses a prediction at large n. If it does not, the
#    honest answer is that the confound cannot be resolved by this route and the
#    interaction comparison has to be restricted to panels where the exact
#    kernel is affordable.
#
# THE DECISION RULE, which is what was actually asked for. For each cell define
#
#     delta(rank) = r_int(megalmm) - r_int(dge_ige at that rank)
#
# The published conclusion is delta > 0. The script reports delta at every rank
# and, where `exact` was run, delta at FULL rank -- which is the ceiling, since
# no kron_rank can do better than the exact kernel. If delta(exact) > 0 with a
# margin, then **no value of SIM_KRON_RANK can flip that cell** and the
# conclusion is safe as stated. If delta(exact) <= 0, the published finding is an
# artefact of the truncation and must be withdrawn for that cell.
#
# It also reports the BREAKEVEN: the r_int that dge_ige would have to reach to
# erase the gap, and the trace retention that would imply on the fitted
# trace-to-accuracy curve. A breakeven above 100% retention means unreachable.
# ============================================================

library(tidyverse)

here::i_am("code/sim_kron_study.R")

source(here::here("code", "dge_ige_functions.R"))
source(here::here("code", "sim_config.R"))
source(here::here("code", "sim_generate.R"))
source(here::here("code", "sim_fit.R"))
source(here::here("code", "megalmm_setup.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
check   <- "--check" %in% args
refresh <- "--refresh" %in% args

# "full" means rank = n_acc, which spans the whole Kronecker space -- the
# ceiling no truncated basis can beat.
ranks   <- strsplit(arg_value("--ranks", "10,20,30,full"), ",")[[1]]
n_accs  <- as.integer(strsplit(arg_value("--n-acc", "60,100,200"), ",")[[1]])
sparsities <- as.numeric(strsplit(arg_value("--sparsity", "0.16,0.48"), ",")[[1]])
n_factors  <- as.integer(strsplit(arg_value("--n-factors", "1,5"), ",")[[1]])
n_reps  <- as.integer(arg_value("--reps", "3"))
n_iter  <- as.integer(arg_value("--n-iter", "6000"))
burn_in <- as.integer(arg_value("--burn-in", "1000"))

if (check) {
  ranks <- c("10", "full"); n_accs <- 60L; sparsities <- 0.48
  n_factors <- 1L; n_reps <- 1L; n_iter <- 2000L; burn_in <- 400L
}

out_dir   <- here::here("output")
cache_dir <- here::here("output", "simulation_kron")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

BASE_SEED <- 20280000L   # clear of both other designs

design <- tidyr::expand_grid(rep = seq_len(n_reps), n_acc = n_accs,
                             sparsity = sparsities, n_factors = n_factors) |>
  dplyr::mutate(
    scenario = sprintf("kron_n%d_sp%03d_f%d", n_acc, round(sparsity * 1000),
                       n_factors),
    seed = BASE_SEED + dplyr::row_number())

message(nrow(design), " scenario-replicate(s) x ", length(ranks),
        " rank setting(s): ", paste(ranks, collapse = ", "))

# How much of the Kronecker trace a given rank retains -- the x axis the
# extrapolation runs on, and computable without fitting anything.
trace_retained <- function(G_oat, G_pea, rank) {
  share <- function(G) {
    v <- pmax(eigen(G, symmetric = TRUE, only.values = TRUE)$values, 0)
    k <- min(rank, length(v))
    sum(v[seq_len(k)]) / sum(v)
  }
  share(G_oat) * share(G_pea)
}

run_cell <- function(row) {
  f <- file.path(cache_dir, paste0(row$scenario, "_rep", row$rep, "_s",
                                   row$seed, ".rds"))
  if (!refresh && file.exists(f)) {
    message("  cached: ", row$scenario, " rep ", row$rep)
    return(readRDS(f))
  }
  message("\n### ", row$scenario, " rep ", row$rep, " ###")
  grms <- sim_grms(row$n_acc, seed = row$n_acc)
  set.seed(row$seed)
  sim <- simulate_experiment(grms$G_oat, grms$G_pea, sparsity = row$sparsity,
                             n_factors = row$n_factors, interaction_pct = 0.20,
                             n_envs = 1)
  sp <- prepare_scenario(sim)
  n_combos <- dplyr::n_distinct(paste(sim$obs$oat, sim$obs$pea))

  # MegaLMM once per scenario: it does not depend on kron_rank
  rd <- file.path(tempdir(), paste0("kron_", row$scenario, "_", row$rep))
  dir.create(rd, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(rd, recursive = TRUE), add = TRUE)
  mm <- fit_megalmm_both(sp$train, sim$G_oat, sim$G_pea, runID = rd,
                         K = 5L, eigen_variance = 0.75,
                         fixed_main_effect = FALSE)
  s_mm <- score_predictions(mm, sim, sp$idx, "megalmm", train = sp$train)

  rows <- purrr::map(ranks, function(rk) {
    is_full <- identical(rk, "full")
    kr <- if (is_full) as.integer(row$n_acc) else as.integer(rk)
    if (kr > row$n_acc) {
      message("  skipping rank ", kr, ": above the panel size")
      return(NULL)
    }
    if (kr^2 > 40000L) {
      message("  skipping rank ", kr, ": ", kr^2,
              " basis columns is beyond what is affordable here")
      return(NULL)
    }
    message("  kron_rank = ", kr, if (is_full) " (FULL)" else "",
            " -> ", kr^2, " columns")
    t0 <- Sys.time()
    fit <- tryCatch(
      fit_dge_ige(sp$train, sim$G_oat, sim$G_pea, with_interaction = TRUE,
                  kron_rank = kr, nIter = n_iter, burnIn = burn_in,
                  seed = row$seed),
      error = function(e) { message("    FAILED: ", conditionMessage(e)); NULL })
    if (is.null(fit)) return(NULL)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "s"))
    s <- score_predictions(fit, sim, sp$idx, "dge_ige", train = sp$train)
    dplyr::mutate(s,
      kron_rank = as.character(kr), is_full = is_full,
      n_basis_cols = kr^2,
      trace_kept = trace_retained(sim$G_oat, sim$G_pea, kr),
      seconds = secs, .after = model)
  }) |> purrr::compact() |> purrr::list_rbind()

  out <- dplyr::bind_rows(
    dplyr::mutate(s_mm, kron_rank = NA_character_, is_full = NA,
                  n_basis_cols = NA_integer_,
                  trace_kept = NA_real_, .after = model),
    rows) |>
    dplyr::bind_cols(dplyr::select(row, scenario, rep, seed, n_acc, sparsity,
                                   n_factors)[rep(1, nrow(rows) + 1), ]) |>
    dplyr::mutate(n_combos = n_combos)
  saveRDS(out, f)
  out
}

res <- purrr::map(seq_len(nrow(design)), \(i) run_cell(design[i, ])) |>
  purrr::list_rbind()
readr::write_csv(res, file.path(out_dir, "simulation_kron_study.csv"))

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

if (nrow(res) == 0) { message("nothing ran"); quit(save = "no") }

cat("\n", strrep("=", 76), "\nr_int_oat by panel size, density, rank and kron_rank\n",
    strrep("=", 76), "\n", sep = "")
print(res |>
        dplyr::group_by(n_acc, sparsity, n_factors, model, kron_rank) |>
        dplyr::summarise(trace = round(mean(trace_kept), 3),
                         r_int_oat = round(mean(r_int_oat), 3),
                         r_int_pea = round(mean(r_int_pea), 3),
                         secs = round(mean(seconds)), .groups = "drop"),
      n = 60, width = 200)

# delta = megalmm - dge_ige, the quantity the published conclusion is about
mm_r <- res |> dplyr::filter(model == "megalmm") |>
  dplyr::select(scenario, rep, mm_oat = r_int_oat, mm_pea = r_int_pea)
gap <- res |>
  dplyr::filter(model == "dge_ige") |>
  dplyr::left_join(mm_r, by = c("scenario", "rep")) |>
  dplyr::mutate(delta_oat = mm_oat - r_int_oat,
                delta_pea = mm_pea - r_int_pea)

cat("\n", strrep("=", 76),
    "\ndelta = MegaLMM - dge_ige on r_int. delta > 0 is the published finding.\n",
    strrep("=", 76), "\n", sep = "")
summ <- gap |>
  dplyr::group_by(n_acc, sparsity, n_factors, kron_rank) |>
  dplyr::summarise(n = dplyr::n(),
                   trace = round(mean(trace_kept), 3),
                   delta_oat = round(mean(delta_oat), 3),
                   se = round(stats::sd(delta_oat) / sqrt(dplyr::n()), 3),
                   .groups = "drop")
print(summ, n = 60, width = 200)

cat("\n", strrep("=", 76), "\nTHE VERDICT, per cell\n", strrep("=", 76), "\n", sep = "")
verdict <- gap |>
  dplyr::group_by(n_acc, sparsity, n_factors) |>
  dplyr::group_modify(function(d, key) {
    ex <- dplyr::filter(d, isTRUE(is_full) | is_full, !is.na(delta_oat))
    at30 <- dplyr::filter(d, kron_rank == "30", !is.na(delta_oat))
    d_ex <- if (nrow(ex)) mean(ex$delta_oat) else NA_real_
    se <- if (nrow(ex) > 1) stats::sd(ex$delta_oat) / sqrt(nrow(ex)) else NA_real_
    tibble::tibble(
      delta_at_30 = if (nrow(at30)) round(mean(at30$delta_oat), 3) else NA_real_,
      delta_full  = round(d_ex, 3),
      se_full     = round(se, 3),
      verdict = dplyr::case_when(
        is.na(d_ex)                      ~ "full rank not run -- cannot certify",
        !is.na(se) & d_ex - 2 * se > 0   ~ "SAFE: MegaLMM ahead at FULL rank",
        d_ex > 0                         ~ "probably safe, but within noise",
        TRUE ~ "WITHDRAW: the gap is an artefact of truncation"))
  }) |> dplyr::ungroup()
print(verdict, n = 40, width = 200)

cat("\nReading this: `delta_full` is the gap at rank = n_acc, where the basis\n",
    "spans the WHOLE Kronecker space. That is the ceiling for dge_ige, so if it\n",
    "is positive, no value of\n",
    "SIM_KRON_RANK can overturn that cell -- the truncation costs dge_ige\n",
    "something, but not enough to matter. If it is negative, the published\n",
    "finding for that cell is an artefact and should be withdrawn.\n")

# the breakeven, for cells where exact could not be run
cat("\n", strrep("=", 76),
    "\nBreakeven where full rank was not affordable\n", strrep("=", 76), "\n", sep = "")
be <- gap |>
  dplyr::filter(!isTRUE(is_full)) |>
  dplyr::group_by(n_acc, sparsity, n_factors) |>
  dplyr::summarise(
    ranks_run = paste(sort(unique(as.integer(kron_rank))), collapse = "/"),
    trace_lo = round(min(trace_kept), 3), trace_hi = round(max(trace_kept), 3),
    delta_lo = round(mean(delta_oat[trace_kept == min(trace_kept)]), 3),
    delta_hi = round(mean(delta_oat[trace_kept == max(trace_kept)]), 3),
    .groups = "drop") |>
  dplyr::mutate(
    slope_per_trace = round((delta_hi - delta_lo) / (trace_hi - trace_lo), 3),
    trace_at_breakeven = round(trace_hi - delta_hi / slope_per_trace, 3),
    reachable = dplyr::case_when(
      is.na(slope_per_trace) ~ "only one rank run",
      slope_per_trace >= 0 ~ "gap does not shrink with rank -- safe",
      trace_at_breakeven > 1 ~ "breakeven needs >100% of the trace: UNREACHABLE",
      TRUE ~ "breakeven is reachable: run `full` on this cell"))
print(be, n = 40, width = 220)
cat("\nThe linear extrapolation in `trace_at_breakeven` is crude and is a\n",
    "SCREEN, not a result: use it to decide which cells need the exact kernel.\n")

message("\nwrote ", file.path(out_dir, "simulation_kron_study.csv"))
