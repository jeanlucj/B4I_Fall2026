# EVALUATION — stepping through the simulation

The companion to [EVALUATION.md](EVALUATION.md), for the simulation framework rather than the analysis. [SIMULATION.md](SIMULATION.md) says what the simulation is for and what it found; this file is how to satisfy yourself that it is measuring what it claims to. The tick list is [EVALUATION_SIMULATION_CHECKLIST.md](EVALUATION_SIMULATION_CHECKLIST.md).

The simulation is easier to evaluate than the analysis, and more important to. Easier because **it is entirely offline and it knows its own truth**: nothing needs T3, and every fitted number can be compared against the quantity that was simulated rather than against another model. More important because the simulation is what licenses the interpretation of the negative MegaLMM result — "the factor model needs about 15% of cells observed and B4I has 3%" is a claim about the simulation's calibration, not about the real data. If the simulation is mis-wired in the direction that handicaps MegaLMM, the analysis conclusion inherits the error and looks independently confirmed.

Two failure modes to hunt, in this order:

1.  **The generator does not produce the truth it advertises** — the variance budget does not add up, the interaction is not the rank it says, the sparsity is not the sparsity requested. Then every scenario label is wrong.
2.  **The scorer does not measure what it names** — `r_interaction` picking up main effects, the two frameworks seeing different splits, the "best" MegaLMM setting chosen on the test set.

------------------------------------------------------------------------

## 1. The modules, in one screen

Four files, one of which is a driver. The **group** column is what you pass to `arm_evaluation()`.

| File (`code/…`) | Owns | Group |
|------------------------|------------------------|------------------------|
| `sim_config.R` | the factorial, the MegaLMM sweep, and every parameter taken from the B4I data | `sim_design` |
| `sim_generate.R` | GRM panels, effect drawing, which cells are observed, one simulated experiment | `sim_panel`, `sim_cells`, `sim_truth` |
| `sim_fit.R` | preprocessing, masking, design bases, the three models, the scoring | `sim_split`, `sim_basis`, `sim_bglr`, `sim_megalmm`, `sim_score` |
| `sim_run.R` | the grid driver: caching, job-array slicing, `--check`, the summary tables | *(driver; use `eval_defs`)* |
| `megalmm_setup.R` | shared with the analysis — see EVALUATION.md A6 | `megalmm` |

The chain, once: `sim_grid()` → `sim_grms()` → `simulate_experiment()` → `split_observations()` → {`fit_dge_ige()` ×2, `fit_megalmm()`} → `score_predictions()`.

`arm_evaluation("sim_pipeline")` arms the whole middle of that chain at once, so a single `run_scenario_bglr(sim)` breaks at every stage in order.

------------------------------------------------------------------------

## 2. Bootstrap (paste once per session)

``` r
library(tidyverse)
here::i_am("code/evaluation.R")
source(here::here("code", "evaluation.R"))

eval_load("simulation")          # sim_config, sim_generate, sim_fit, megalmm_setup
eval_conflicts()                 # expect: no name defined twice
eval_groups("simulation")
```

`eval_conflicts()` should answer **no FUNCTION is defined by more than one loaded file**, then list the settings names each driver defines in its own block. Note what `eval_load("simulation")` sources: `dge_ige_functions.R` comes first, because the simulation draws its effects with the **same** `grm_factor()` and reads the GRMs with the same `read_grm()` that the DGE-IGE model fits with. Until 2026-09-23 `sim_generate.R` carried its own copies of both; loading the analysis and the simulation together would silently have picked one. Now there is one of each and it does not matter.

`grm_basis()` in `sim_fit.R` is a genuinely different function — it truncates by rank or by cumulative variance and is used to build design matrices, not to factor a GRM whole — so it stays separate.

One reference object for the whole document:

``` r
panel <- sim_grms(100, seed = 100)
set.seed(7)
sim <- simulate_experiment(panel$G_oat, panel$G_pea, sparsity = 0.20,
                           n_factors = 1, interaction_pct = 0.20, n_envs = 1)
peek(sim)
```

```         
[peek] sim: simulation: 100 x 100 panel, 2000 obs (20.0% of cells), 1 env, 1 factor(s)
        V(true): producer=0.169  associate=0.128  interaction=0.200  residual=0.503
        obs per oat: 20.0  per pea: 20.0  (min 9 / 11)
```

------------------------------------------------------------------------

## 3. Level-by-level walkthrough

Same annotations as EVALUATION.md: **→ returns**, **🔍 eyeball**, **🚩 red flag**, *Assumption:*. Every number quoted below was produced by running the line above it on 2026-09-23.

------------------------------------------------------------------------

### S1 — `sim_design` (instant)

``` r
arm_evaluation("sim_design")
g <- sim_grid(SIM_LEVELS, n_reps = 1)
disarm_evaluation()

nrow(g); dplyr::n_distinct(g$scenario)
dplyr::count(g, n_factors, interaction_pct)
mm_grid <- tidyr::expand_grid(!!!SIM_MEGALMM_LEVELS); nrow(mm_grid)
```

- **→ returns** **120** scenarios, 120 distinct names. `mm_grid` has **12** rows (K 2 × eigen_variance 3 × fixed_main_effect 2), so the full grid is 120 × 2 BGLR fits + 120 × 12 = **1440 MegaLMM fits**.

- **🔍 eyeball** the redundant cells really are collapsed: all 12 rows with `n_factors == 0` have `interaction_pct == 0`, not two copies at 0.10 and 0.20.

- **🚩 red flag — the cache key does not include the seed.** `seed` is `SIM_BASE_SEED + row_number()` over the *whole grid*, but the cache filename is built from the scenario name and rep only. Change the grid and the same file name now corresponds to a different seed:

  ``` r
  a <- sim_grid(SIM_LEVELS, 1); b <- sim_grid(SIM_LEVELS_EXTENDED, 1)
  j <- dplyr::inner_join(dplyr::select(a, scenario, s1 = seed),
                         dplyr::select(b, scenario, s2 = seed), by = "scenario")
  sum(j$s1 != j$s2)     # 60 of 60
  ```

  So `--extended` silently reuses all 60 standard-grid cache files, whose `seed` column records a seed that the extended grid would never have used. The *design* is the same, so the results remain valid replicates — but `simulation_results.csv` then contains seeds that do not reproduce their own rows. Same applies to `--reps`. Use `--refresh` when changing the grid, or treat the seed column as decorative.

- **🚩 red flag** `SIM_BASE_SEED` is defined *after* `sim_grid()` in `sim_config.R`. Harmless (R resolves it at call time) but it means `sim_grid()` reads a global rather than taking a seed argument.

- *Assumption:* `SIM_VAR_SHARES` is grounded in the fitted B4I model. It is — verified: 394 (oat producer on oat yield) / 299 (`G_pea`'s `var_As`, the pea's effect **on oat yield**) / 1181 (residual oat variance), normalising to 0.211 / 0.160 / 0.630. Re-derive with `sim_observed_parameters()` if the analysis is re-run.

------------------------------------------------------------------------

### S2 — `sim_panel` (instant)

``` r
arm_evaluation("sim_panel")
panel <- sim_grms(200, seed = 200)
L <- grm_factor(panel$G_oat)
set.seed(1); v <- draw_effect(L, 0.25)
disarm_evaluation()

max(abs(tcrossprod(L) - panel$G_oat)); var(v); mean(v)
```

- **→ returns** two 200 × 200 matrices taken from the real GRMs; `L` is 200 × 200; `tcrossprod(L)` reproduces `G` to 8.5e-15; `var(v)` is **exactly** 0.25 and `mean(v)` is 0.

- **🔍 eyeball** the exactness is the design: `draw_effect()` rescales the realised vector to the target variance, so a single replicate has the variance the scenario asked for rather than having it in expectation. That removes one source of between-replicate noise and is worth knowing when reading replicate spread.

- **🚩 red flag** the panel is drawn from the **real** GRMs, which carry 8 oat and 20 pea prior-only rows (EVALUATION.md A1) — accessions related to nothing. At `n_acc = 200` of 508/435 a few will be in any given panel, and for them `L` contributes an independent column, so kinship-borrowing has nothing to borrow. Count them:

  ``` r
  sum(rowSums(abs(panel$G_pea - diag(diag(panel$G_pea)))) == 0)
  ```

- *Assumption:* `sim_grms()` samples accessions at random so a panel is not systematically the most or least related part of the collection, and uses a fixed `seed = n_acc` in `sim_run.R` so every scenario of a given size shares one panel. That is deliberate: the panel is not part of the factorial.

------------------------------------------------------------------------

### S3 — `sim_cells` (instant; the sparsity claim)

``` r
arm_evaluation("sim_cells")
set.seed(2); cl <- sample_combinations(200, 200, 0.05)
disarm_evaluation()

nrow(cl); round(0.05 * 200 * 200)
min(table(cl$oat)); min(table(cl$pea))
try(sample_combinations(200, 200, 0.002))      # must stop, not silently thin
```

- **→ returns** exactly **2000** cells for a request of 2000; minimum 4 observations per oat and per pea against a floor of `SIM_MIN_PER_ACC = 3`.
- **🔍 eyeball** "5% observed" means 5% observed, exactly. The guaranteed minimum is built first and the remainder filled at random, so the floor does not inflate the count — check this at the *lowest* sparsity you use, where the floor is closest to binding: at 3% and 200 × 200 it is 1200 requested against \~600 forced.
- **🚩 red flag** if the forced cells ever exceeded `n_obs`, `extra` would be `sample(remaining, max(0, negative))` = nothing, and the achieved sparsity would silently exceed the requested one with no warning. The `stop()` guard catches the clearly-impossible case only. Always confirm `nrow(sim$obs)` equals `round(sparsity * n_acc^2)` rather than trusting the label — `peek(sim)` prints the achieved percentage for exactly this reason.
- *Assumption:* guaranteeing every accession `min_per_acc` observations is both more realistic than uniform sampling and necessary, because a pea column holding one oat has no estimable residual variance and takes MegaLMM's ARD sampler to `NaN` rather than to an error.

------------------------------------------------------------------------

### S4 — `sim_truth` (instant; does the budget add up?)

``` r
arm_evaluation("sim_truth")
set.seed(7)
sim <- simulate_experiment(panel$G_oat, panel$G_pea, sparsity = 0.20,
                           n_factors = 1, interaction_pct = 0.20, n_envs = 1)
disarm_evaluation()

sum(sim$truth$V)                                  # 1
var(sim$truth$producer);  sim$truth$V[["producer"]]
var(as.vector(sim$truth$interaction)); sim$truth$V[["interaction"]]
qr(sim$truth$interaction)$rank                    # == n_factors
```

- **→ returns** `V` sums to 1; `var(producer)` equals its target to machine precision (0.1686…); interaction variance over all cells equals `interaction_pct` exactly; the interaction matrix has rank `n_factors`.
- **🔍 eyeball** the rank check is the axis the whole comparison turns on. With `n_factors = 1` the interaction is `u λ'`, rank 1 — MegaLMM's home ground. With `n_factors = 5` it approaches the full Kronecker structure DGE-IGE assumes. Confirm the rank is what the scenario name says before believing any crossover.
- **🚩 red flag — the simulated phenotype is not on a yield scale.** `y = env_mean + env_scale * (genetic + resid)` with the variance budget normalised to 1, so at `n_envs = 1` you get mean 167.3 g/m² and **sd ≈ 0.98**, where the real trials have sd ≈ 34. Harmless — everything is standardised within environment before fitting — but `sim$obs$y` must not be read as a yield, and `SIM_GRAND_MEAN` is doing nothing except sitting in front of a unit-variance signal.
- **🚩 red flag** environments carry **mean and variance heterogeneity only**: `env_scale` multiplies signal and noise alike, so genetic correlations across environments stay exactly 1. There is no rank-changing G×E in this simulation. If a conclusion is about environments, that is the assumption it rests on. Verify: at `n_envs = 10`, `cor` of the genetic value between any two environments is 1 by construction.
- *Assumption:* drawing interaction scores from `G_oat` and loadings from `G_pea` makes one factor's covariance exactly `G_oat[i,i'] * G_pea[j,j']`, so neither framework is favoured by construction. This is the fairness claim of the whole exercise; `qr()$rank` plus the covariance check is how you test it.

------------------------------------------------------------------------

### S5 — `sim_split` (instant; both frameworks must see one split)

``` r
arm_evaluation("sim_split")
a <- split_observations(sim, SIM_CV_FRACTION, seed = 7)
disarm_evaluation()
b <- split_observations(sim, SIM_CV_FRACTION, seed = 7)

identical(a$held$cell, b$held$cell)               # TRUE -- the load-bearing one
nrow(a$train); nrow(a$held)
min(table(a$train$oat)); min(table(a$train$pea))
length(intersect(a$train$cell, a$held$cell))      # 0
```

- **→ returns** identical splits from repeated calls; 1600 train / 400 held; minimum 7 training observations per oat and per pea; train and held disjoint.
- **🔍 eyeball** this is what makes the head-to-head fair. `run_scenario_bglr()` and `run_scenario_megalmm()` are separate functions, cached separately, often run in different processes — and they agree only because both call `split_observations(sim, cv_fraction, seed)` with the same seed and it calls `set.seed()` internally. **If this check ever fails, every head-to-head comparison in `simulation_results.csv` is between models scored on different data.**
- **🚩 red flag** `standardize_within_env()` falls back to centring only when an environment has ≤ 2 observations or zero variance. At `n_envs = 10` and the sparsest cells, check how many environments take the fallback — a centred-but-unscaled environment sits on a different scale from its neighbours, and that is a real cost of the `n_envs` axis rather than a bug.
- *Assumption:* the floor in `mask_observations()` (`SIM_FLOOR_OBS = 2`) matches the real analysis's `mask_cells()` floor, so the simulation and the analysis are masking under the same rule. They are, but they are two separate implementations of it.

------------------------------------------------------------------------

### S6 — `sim_basis` (instant; the reshape that must be exact)

The specific-combination term is fitted in a low-rank Kronecker basis and then reshaped into a full surface. The reshape is index arithmetic, it is not checked anywhere in the code, and getting it wrong transposes the interaction without changing any dimension.

**To understand this more fully, read [docs/specific-combination_kronecker.md](docs/specific-combination_kronecker.md) first** — it derives both from scratch with a worked 3-oat × 4-pea example, and shows what the wrong reshape looks like (no error, right shape, r = 0.983 against the truth). This level is the check; that document is the explanation.

``` r
arm_evaluation("sim_basis")
A <- grm_basis(panel$G_oat, rank = 5)
B <- grm_basis(panel$G_pea, rank = 5)
disarm_evaluation()

set.seed(3); Beta <- matrix(rnorm(25), 5, 5)
oi <- c(1, 4, 9); pj <- c(2, 7, 11)
lhs <- as.vector(kron_basis(A, B, oi, pj) %*% as.vector(t(Beta)))
rhs <- (A %*% Beta %*% t(B))[cbind(oi, pj)]
max(abs(lhs - rhs))                               # 6.9e-18
```

- **→ returns** agreement to 7e-18. This is the oracle: the prediction you get by multiplying the fitted basis by the coefficient vector must equal the prediction you get from the reshaped surface, or `fit_dge_ige()`'s `interaction` matrix is not the model that was fitted.
- **🔍 eyeball** `kron_basis` puts `A[,a] * B[,b]` in column `(a-1)*rank + b`, and `matrix(..., byrow = TRUE)` reads the coefficients back in that same order. Swap `each =` and `times =` in `kron_basis`, or drop `byrow`, and the check above fails immediately — which is why it is worth running once after any edit to either function.
- **🚩 red flag** `SIM_KRON_RANK = 30` per species means the interaction basis has 900 columns and spans a *subspace* of the exact Kronecker kernel. The DGE-IGE interaction is therefore fitted at reduced rank while MegaLMM's is not, at every cell above `SIM_KRON_EXACT_MAX = 3000` observations. When the comparison is close, check whether raising `SIM_KRON_RANK` moves it.
- *Assumption:* `grm_basis(G, variance = 0.999)` for the main effects keeps essentially the whole space, while `rank = 30` for the interaction does not. Two different truncation rules in one model.

------------------------------------------------------------------------

### S7 — `sim_bglr` (seconds to a minute)

``` r
arm_evaluation("sim_bglr")
add <- fit_dge_ige(a$train, sim$G_oat, sim$G_pea, with_interaction = FALSE)
dge <- fit_dge_ige(a$train, sim$G_oat, sim$G_pea, with_interaction = TRUE)
disarm_evaluation()

dim(add$total); max(abs(add$interaction))         # 100 x 100; exactly 0
peek(as.vector(dge$interaction))
```

- **→ returns** full 100 × 100 prediction surfaces. The additive model's `interaction` is an exact zero matrix, which is the correct answer for it and makes `r_interaction` for `additive` `NA` rather than misleading.
- **🔍 eyeball** `total = outer(pr, as_, "+")` — an additive surface from the two main effects, plus the interaction when fitted. Confirm `dimnames(add$total)` carry the accession names, because `score_predictions()` indexes by integer position and a name/position mismatch would be invisible.
- **🚩 red flag** BGLR writes MCMC sample files to `saveAt = file.path(tempdir(), "sim_")`. Two `fit_dge_ige()` calls in one session share that prefix and overwrite each other's `.dat` files. It does not affect the returned fit (the coefficients are in memory), but any attempt to go back and inspect the chains afterwards reads the *second* model's samples under the first model's expectations.
- *Assumption:* `nIter = 6000, burnIn = 1000` is enough. Nothing in the simulation checks convergence of the BGLR half — the `--trace` machinery covers MegaLMM only. If in doubt, refit one scenario at 20,000 and compare.

------------------------------------------------------------------------

### S8 — `sim_megalmm` (5–40 s per fit)

``` r
arm_evaluation("sim_megalmm")
d <- file.path(tempdir(), "mmcheck"); dir.create(d, showWarnings = FALSE)
r <- run_scenario_megalmm(sim, K = SIM_MEGALMM_K,
                          eigen_variance = SIM_EIGEN_VARIANCE, seed = 7,
                          run_dir = d, fixed_main_effect = TRUE)
disarm_evaluation()
as.data.frame(r$scores)
```

On the 100 × 100 panel at 20% observed this takes about 6 s and returns:

| model | r_total | r_interaction | r_observed | r_mainfactor_truePr | r_rowmean_truePr | fixed_ok |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| megalmm | 0.886 | 0.780 | 0.597 | 0.739 | 0.939 | TRUE |
| megalmm_U | 0.689 | 0.808 | 0.455 | 0.739 | 0.939 | TRUE |

- **🔍 eyeball** `fixed_ok` is `TRUE`: `sd(Lambda[1, ]) < 1e-8`. The loadings are constant across columns but **not equal to 1** — `remove_nuisance_parameters` rescales `Lambda` by `sqrt(var(F))`, so row 1 comes back at the main-effect SD. Testing `all(Lambda[1,] == 1)` would fail on a working model.

- **🔍 eyeball — the fixed factor does work, and still loses to a row average.** Over the recorded results the main-effect recovery is:

  |   | mean `|r_mainfactor_truePr|` | mean `r_rowmean_truePr` |
  |----|----|----|
  | `fixed_main_effect = FALSE` | 0.067 | 0.760 |
  | `fixed_main_effect = TRUE` | 0.674 | 0.760 |

  So pinning the first factor's loadings is doing exactly what it was added for — a free factor recovers essentially nothing of the oat main effect, a fixed one recovers most of it — and it still does not reach a plain `tapply(y, oat, mean)`. Both halves of that are worth saying; `SIMULATION.md` currently says neither.

- **🚩 red flag — the sign of `r_mainfactor_truePr` is only meaningful when the factor is fixed.** A free factor's sign is arbitrary, so with `fixed_main_effect = FALSE` the column comes back negative about half the time (3 of 6 recorded rows; `--check` returns −0.939). Averaging the signed column then cancels: mean signed 0.021 against mean absolute 0.067. Take `abs()` before summarising, or restrict to the fixed-factor rows.

- **🚩 red flag — an alignment that works by luck.** `r_mainfactor_truePr` is `cor(mm$U_F[, 1], sim$truth$producer[seq_len(nrow(mm$U_F))])`. `U_F`'s rows are the oats that survived `keep_row` (those with at least one training observation), in `rownames(G_oat)` order; `truth$producer` is indexed 1…n_oat. Taking the **first** `nrow(U_F)` producers is only correct when `keep_row` is all `TRUE`. It is, today, because `SIM_MIN_PER_ACC = 3` and `SIM_FLOOR_OBS = 2` guarantee it — but the correlation would silently misalign, not error, if that ever stopped holding. Check `all(rowSums(!is.na(Y)) > 0)` before trusting this column.

- **🚩 red flag** `place()` writes the fitted submatrix into a full `n_oat × n_pea` matrix of **zeros**. Any oat or pea dropped by `keep_row` / `keep_col` gets a prediction of exactly 0 for every held-out cell, dragging the correlation down with no diagnostic.

- **🔍 eyeball** expect one `Warning: In cor(F) : the standard deviation is zero` when `fixed_main_effect = TRUE`. It comes from `reorder_factors()` seeing the pinned factor and is expected; it should *not* appear with `fixed_main_effect = FALSE`.

- *Assumption:* `Eta_mean` (the predicted phenotype) is the right target, not `U`. The table above is the awkward case: `megalmm` (Eta) wins on `r_total` but `megalmm_U` wins on `r_interaction`. Both are reported, which is the honest choice — but the headline picks one.

------------------------------------------------------------------------

### S9 — `sim_score` (instant; the most load-bearing arithmetic)

Every model returns a full surface, and the surface splits **exactly** into an
additive and an interaction part:

``` r
M == additive_part(M) + interaction_part(M)          # to machine precision
rowMeans(additive_part(M)) == rowMeans(M)
```

That is what makes the metrics comparable across frameworks, and it is where
`r_gma` (the additive part against `Pr + As`), `r_oat_prod` (`rowMeans` against
the true producer effect) and `r_pea_assoc` (`colMeans` against the true
associate effect) come from. Note which effect `r_pea_assoc` is: the response
is **oat yield**, so the column margin is the **pea's** effect on the oat. The
oat's own associate effect is not simulated — see SIMULATION.md's caveats.

**🚩 red flag — the factor route to a main effect does not work, and it is worth
knowing why.** One might try to read the column margin out of the factor
structure — mean factor score times loading, summed over factors. Measured on an
80 × 80 panel at 50% observed, K = 8: the scores are *not* centred
(`colMeans(U_F)` up to −0.54) so the route is non-zero, but it explains **7.5%**
of the column-margin variance and correlates **0.051** with the true associate
effect against **0.938** for the plain `colMeans`. The effect sits in the
per-column intercept, which MegaLMM fits as a **fixed, unshrunk** term with no
kinship — whereas the row margin has no intercept at all and must come from the
kinship-shrunk latent structure. That asymmetry predicts `r_pea_assoc` degrading
faster than `r_oat_prod` as sparsity falls, which the 1.6% level should show.

``` r
arm_evaluation("sim_score")
M  <- matrix(rnorm(20), 4, 5) + outer(1:4, 1:5, "+")
IP <- interaction_part(M)
disarm_evaluation()

max(abs(rowMeans(IP))); max(abs(colMeans(IP)))          # ~1e-16 both
max(abs(interaction_part(outer(rnorm(4), rnorm(5), "+"))))   # ~3e-16
```

- **→ returns** row and column means of the residualised matrix are zero, and a purely additive surface residualises to exactly zero.

- **🔍 eyeball** that second line is the point: `r_interaction` for the additive model must be `NA` (zero variance), not a small positive number. If an additive surface leaves a non-zero residual, `interaction_part()` has its recycling backwards — `M - rowMeans(M)` recycles down columns, which is what is wanted, but `- rep(colMeans(M), each = nrow(M))` has to be written exactly that way to match.

- **🚩 red flag** `score_predictions()` compares against `truth$producer + truth$associate + truth$interaction` at the held cells — the **noise-free genetic value**, without the residual and without the environment scaling. So `r_total` has no ceiling imposed by heritability and is not comparable to an accuracy measured against observed phenotypes. `r_observed` is the one to quote against a phenotype. Three correlations with three different targets live in the same table; label them carefully when reporting.

- **🚩 red flag** `safe_cor()` returns `NA` when either vector has sd \< 1e-10. `NA`s then flow into `mean(r_interaction)` in the summaries — some with `na.rm = TRUE`, some without. Check which scenarios contribute `NA` before reading a mean:

  ``` r
  res <- readr::read_csv(here::here("output", "simulation_results.csv"))
  res |> dplyr::group_by(model) |>
    dplyr::summarise(na_total = sum(is.na(r_total)),
                     na_int = sum(is.na(r_interaction)))
  ```

- *Assumption:* residualising both prediction and truth puts the two frameworks on the same footing — DGE-IGE has an explicit interaction term, MegaLMM's factors carry main effects and interaction together, and neither hands back something directly comparable.

------------------------------------------------------------------------

### S10 — the driver (`sim_run.R`)

``` r
eval_defs("code/sim_run.R")       # arg_value, has_flag, cache_path, run_one
cache_path("n200_sp016_f1_i20_e01_g100", 1, 20260922, "bglr")
list.files(here::here("output", "simulation")) |> head()
```

- **→ returns** one `_bglr.rds` and twelve `_mm_K*_ev*_fx*.rds` per scenario × rep.
- **🔍 eyeball** the two halves are cached apart *because* the MegaLMM settings change nothing about the simulated data — so sweeping them must not refit BGLR. Confirm the BGLR half is written once per scenario, not twelve times.
- **🚩 red flag** the combined CSV is always rebuilt by globbing the cache, not from the current run, so **a stale cache file silently joins the results**. The glob is `_(bglr|mm_K[0-9]+_ev[0-9]+_fx[01])\.rds$`. After changing a generator or a scorer, delete `output/simulation/` or use `--refresh`; nothing versions the cache against the code that wrote it.
- **🚩 red flag** `--check` must be run before trusting any sweep. It fits a dense, strongly structured scenario (100 × 100 at 50% observed) where MegaLMM should clearly win, and stops if `r_total < 0.5`. It currently reaches **0.912** against an additive 0.731, so the margin is comfortable; a `--check` that lands just over 0.5 is itself a finding.
- *Assumption:* job-array slicing by scenario (not by row) keeps both halves of a scenario in one task, so a simulation is never generated twice with different results. Confirm with `--task 1 --ntasks 2` that the slices partition the grid.

------------------------------------------------------------------------

### S11 — reading the results

``` r
res <- readr::read_csv(here::here("output", "simulation_results.csv"))
dplyr::count(res, model)
res |> dplyr::filter(model == "megalmm") |>
  dplyr::count(scenario, rep) |> dplyr::count(n)      # 12 settings per scenario
```

- **🚩 red flag — the headline MegaLMM number is a maximum over 12 settings chosen on the held-out data, and DGE-IGE gets no such selection.** `sim_run.R` builds `best_mm` with `slice_max(r_interaction)` per scenario and compares *that* to a `dge_ige` fitted once at fixed settings. The two are not being asked the same question: one is "how well does this model do", the other is "how well does the best of twelve tunings of this model do, judged on the same cells it is scored on". Measure the gap rather than assuming it is small — on the one scenario currently in the cache (`n200_sp016_f1_i20_e01_g100`):

  ``` r
  mm <- dplyr::filter(res, model == "megalmm")
  mean(dplyr::slice_max(mm, r_interaction, n = 1)$r_interaction)   #  0.109  best of 12
  mean(dplyr::filter(mm, K == 10, eigen_variance == 0.8,
                     !fixed_main_effect)$r_interaction)            # -0.026  the default
  mean(mm$r_interaction)                                           #  0.028  averaged
  mean(res$r_interaction[res$model == "dge_ige"])                  #  0.044
  ```

  The selection moves MegaLMM from −0.026 to +0.109, a swing larger than the DGE-IGE value it is being compared with — and it flips the verdict for that scenario from "MegaLMM loses to DGE-IGE" to "MegaLMM wins". One scenario is an illustration, not an estimate, but the mechanism is real and the magnitude is not negligible.

  The direction matters for how much of the published conclusion survives. The bias favours MegaLMM, and the conclusion is that MegaLMM *loses* below \~15% observed, so the headline is conservative: the true crossover is at least as high as reported. What does not survive is any individual cell where MegaLMM narrowly wins.

  Two ways to make the comparison symmetric, if it is ever reported as a like-for-like: fix one MegaLMM setting in advance and report that (the cheapest fix), or give DGE-IGE the same treatment by sweeping `SIM_KRON_RANK` and taking its best. Reporting best-of-12 against fitted-once is the one thing to avoid.

- **🔍 eyeball** the `--check` scenario is not in the grid, so a grid that produces all-zero MegaLMM numbers and a `--check` that passes are consistent with each other and with "the method does not work here". Only the pair of them supports that reading.

- **🔍 eyeball** `SIM_TRACE_CHUNKS` with `--trace` answers the chain-length question without a sweep: if accuracy is still climbing at the last chunk, `SIM_MEGALMM_SAMPLE = 250` is too short and every MegaLMM number is an underestimate.

------------------------------------------------------------------------

## 4. What each module should generate

| Script | Writes | Expected |
|------------------------|------------------------|------------------------|
| `sim_config.R` | nothing | `sim_grid()` → 60 rows; `SIM_MEGALMM_LEVELS` → 12 combinations |
| `sim_generate.R` | nothing | `simulate_experiment()` → `list(obs, truth, G_oat, G_pea, settings)`; `truth$V` sums to 1; `obs` has exactly `round(sparsity * n_acc^2)` rows |
| `sim_fit.R` | nothing | `run_scenario_bglr()` → 4 rows (additive, dge_ige, oat_mean, oat_plus_pea); `run_scenario_megalmm()` → `$scores` 2 rows (megalmm, megalmm_U) + `$trace` |
| `sim_run.R` | `output/simulation/<scenario>_rep<k>_bglr.rds` | 4 rows each |
|  | `output/simulation/<scenario>_rep<k>_mm_K*_ev*_fx*.rds` | 12 files per scenario, 2 rows each |
|  | `output/simulation/*_trace.rds` | only with `--trace` |
|  | `output/simulation_results.csv` | `4 + 24` rows per scenario × rep |
|  | `output/simulation_summary.png` | r_interaction vs sparsity, faceted panel × factors |
|  | `output/simulation_runs/` | MegaLMM run state, deleted after each fit |

`--check` writes nothing and exits non-zero if MegaLMM fails to clear 0.5.

------------------------------------------------------------------------

## 5. Subtle-bug catalogue

Ranked by consequence.

1.  **The "best" MegaLMM setting is picked on the test score, DGE-IGE is fitted once** (S11). MegaLMM gets a max over 12 tunings judged on the held-out cells; DGE-IGE gets no tuning at all. On the one cached scenario the selection is worth +0.135 in `r_interaction` — more than the whole DGE-IGE value. Biases the comparison toward MegaLMM, so the "MegaLMM loses when sparse" conclusion is conservative and survives; individual narrow MegaLMM wins do not.
2.  **Cache keys omit the seed** (S1). `--extended` and `--reps` reuse files whose recorded seed is not the one the current grid assigns.
3.  **The results CSV is rebuilt by globbing the cache** (S10), so results from before a code change join silently.
4.  **`r_mainfactor_truePr` aligns by position** (S8), correct only while `keep_row` is all `TRUE`.
5.  **`place()` fills dropped rows and columns with zero** (S8) rather than `NA`, so a dropped accession degrades the correlation invisibly.
6.  **Three correlations, three targets** (S9): `r_total` against noise-free genetic value, `r_interaction` against the residualised interaction, `r_observed` against the standardised phenotype. Easy to quote the wrong one.
7.  **The DGE-IGE interaction is rank-limited to 30² and MegaLMM's is not** (S6).
8.  **Both `fit_dge_ige()` calls write to the same `saveAt` prefix** (S7).
9.  **The simulated phenotype has a realistic mean and a unit variance** (S4). Harmless, but `y` is not g/m².
10. **BGLR convergence is never checked** (S7); only MegaLMM has `--trace`.

------------------------------------------------------------------------

## 6. Where things live

- `output/simulation/` — one cache file per scenario × rep × half. Regenerable; delete after changing the generator or the scorer.
- `output/simulation_runs/` — MegaLMM run state, created and deleted per fit.
- `output/simulation_results.csv`, `simulation_summary.png` — the combined results, always rebuilt from the cache.
- `data/GRM_Avena.rds`, `data/GRM_Pisum.rds` — the real relationship matrices the panels are drawn from. The simulation needs no T3 login and no phenotypes.
- `code/scinet/` — the SLURM job array for running the grid on SciNet; see [code/scinet/README.md](code/scinet/README.md).
