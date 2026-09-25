# Simulation evaluation checklist

Tick these off over time. Each line names the `arm_evaluation()` group; the
walkthrough is in **[EVALUATION_SIMULATION.md](EVALUATION_SIMULATION.md)** at
the matching level (S1…S11). Everything here is offline and most of it is
instant, so there is no reason to skip ahead. `disarm_evaluation()` at the end
of each level.

Bootstrap once:

```r
library(tidyverse); here::i_am("code/evaluation.R")
source(here::here("code", "evaluation.R"))
eval_load("simulation"); eval_conflicts()          # expect: no duplicated names
panel <- sim_grms(100, seed = 100)
set.seed(7); sim <- simulate_experiment(panel$G_oat, panel$G_pea,
                                        sparsity = 0.20, n_factors = 1,
                                        interaction_pct = 0.20, n_envs = 1)
peek(sim)
```

`eval_load("simulation")` sources `dge_ige_functions.R` first: the simulation
uses the same `grm_factor()` and `read_grm()` the DGE-IGE model fits with, so
loading the analysis alongside it is now harmless.

## The generator — does it produce the truth it advertises?

- [ ] **S1 `sim_design`** — **120** scenarios, 120 distinct names; the cells
      where a level scales nothing collapsed (`interaction_pct` at
      `n_factors == 0`, `gxe_cor` at `n_envs == 1`); `mm_grid` 12 rows → 1440
      MegaLMM fits. The **seed is now part of the cache key**, so a changed grid
      recomputes instead of silently reusing files from a different seed.
- [ ] **S2 `sim_panel`** — `tcrossprod(L)` reproduces `G` to ~1e-15;
      `draw_effect()` hits its target variance exactly with mean 0; count the
      prior-only rows that came along from the real GRM.
- [ ] **S3 `sim_cells`** — `nrow(cl)` equals `round(sparsity * n_acc^2)`
      **exactly** at all four levels on both panel sizes, now asserted in
      `sample_combinations()`; every accession has at least `SIM_MIN_PER_ACC`.
      Confirm the assertion fires by asking for 0.015 at `n_acc = 200`, which
      hugs the floor and used to overshoot to 601–607.
- [ ] **S4 `sim_truth`** — `sum(truth$V) == 1` (now six components, two of them
      zero at `gxe_cor = 1`); realised variances equal their targets;
      `qr(truth$interaction)$rank == n_factors`; `sim$obs$y` has sd ≈ 1, so it is
      not a yield.
- [ ] **S4 GxE** — `gxe_cor = 1` reproduces the pre-GxE generator **bit for
      bit** at both `n_envs` values (the check that the axis is an extension,
      not a rewrite); at `gxe_cor = 0.6`, `truth$realised_gxe_cor` ≈ 0.6 and
      `r_total` falls.

## The scorer — does it measure what it names?

- [ ] **S5 `sim_split`** — **two calls give identical held sets** (both
      frameworks must see one split); train and held disjoint; floors hold;
      count environments that took the centre-only fallback at `n_envs = 10`.
- [ ] **S6 `sim_basis`** — the Kronecker reshape identity holds to ~1e-17:
      `kron_basis(A,B,i,j) %*% as.vector(t(Beta))` equals
      `(A %*% Beta %*% t(B))[cbind(i,j)]`. Re-run after touching either
      function.
- [ ] **S9 `sim_score`** — `M == additive_part(M) + interaction_part(M)` to
      machine precision; `interaction_part()` residualises a purely additive
      surface to zero; know which of `r_total` / `r_gma` / `r_interaction` /
      `r_observed` you are quoting, and that `r_pea_assoc` is the **pea's**
      effect on oat yield; count the `NA`s feeding each mean.

## Bivariate and both orientations

- [ ] **Four effects, correlated as asked** — `cor(oat_prod, oat_assoc)` ≈
      −0.065 and `cor(pea_prod, pea_assoc)` ≈ −0.234 to four decimals (the
      whitening in `draw_effect_pair()` makes these exact, not approximate);
      `var()` of each equals its target; both per-trait budgets sum to 1.
- [ ] **Two interaction surfaces, independent** — `qr(I_oat)$rank` and
      `qr(I_pea)$rank` both equal `n_factors`; `cor(I_oat, I_pea)` ≈ 0.
- [ ] **The comparator is the real model** — `fit_dge_ige()` calls
      `fit_producer_associate()`, so the simulation and the production analysis
      fit one implementation. Its mix term uses the low-rank `kron_rank` basis,
      not the exact kernel.
- [ ] **Both orientations run and the pea side is transposed** —
      `dim(surface$oat) == dim(surface$pea)`, both oat-rows × pea-cols;
      `n_dropped == 0`.
- [ ] **The structural prediction** — `r_*_assoc` degrades faster than
      `r_*_prod` as sparsity falls from 48% to 1.6%, because only the producer
      effects get kinship. On `megalmm_U`, which excludes the per-column
      intercept, the associate effects should be near zero at any density
      (measured: 0.13 and 0.16 on a dense panel against 0.94 for the producer).

## The fractional design

- [ ] **The recoding is what makes it estimable** — a main-effects-plus-2FI
      model is singular on the FULL grid, because `interaction_pct` is nested in
      `n_factors` and `gxe_cor` in `n_envs`. Confirm `sim_design()$model_rank`
      equals `$model_terms` (82).
- [ ] **150 runs covers everything** — all 120 data scenarios touched, all 8
      MegaLMM settings present with roughly equal counts.
- [ ] **Read it as a model, not a table** — `simulation_design_effects.csv`,
      not cell means. There is no per-scenario best-of-settings to take any
      more, which is the point.

## The models

- [ ] **S7 `sim_bglr`** — additive model's `interaction` is exactly 0;
      `dimnames` carry accession names; note that both fits share one BGLR
      `saveAt` prefix; BGLR convergence is unchecked.
- [ ] **S8 `sim_megalmm`** — reference run on the dense 100 × 100:
      `r_total` ≈ 0.89, `r_interaction` ≈ 0.78, `fixed_ok = TRUE` in ~6 s.
      `all(rowSums(!is.na(Y)) > 0)` before trusting `r_mainfactor_truePr`,
      and `abs()` it before averaging (a free factor's sign is arbitrary —
      3 of 6 recorded rows are negative). The fixed factor works (mean |r|
      0.067 → 0.674) but **still loses to a plain row mean (0.760)**; say so.

## The grid

- [ ] **S10 driver** — one `_bglr` file per scenario and twelve `_mm_*` files;
      `--task 1 --ntasks 2` partitions the grid; cache deleted or `--refresh`
      used after any change to `sim_generate.R` or `sim_fit.R`.
- [ ] **`--check` before any sweep** — `Rscript code/sim_run.R --check` reaches
      `r_total` ≈ 0.91 against additive 0.73, comfortably over its 0.5 floor.
      A bare pass is itself a finding.
- [ ] **S11 results** — 12 MegaLMM settings per scenario present;
      **the headline MegaLMM number is a max over those 12 chosen on the
      held-out data while DGE-IGE is fitted once with no tuning.** Measure the
      gap (best-of-12 vs the default setting) before reporting the two side by
      side; on the cached scenario it is +0.135, which flips the verdict.
- [ ] **`--trace` once** — accuracy has plateaued by the last chunk, so
      `SIM_MEGALMM_SAMPLE = 250` is not truncating the chains.

------------------------------------------------------------------------

*Notes / anomalies noticed:*

-

-

-
