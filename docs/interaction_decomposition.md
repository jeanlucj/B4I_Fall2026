# Reading the DGE-IGE interaction as scores and loadings

`dge_ige` predicts the oat × pea interaction better than MegaLMM when
combinations are sparse, and — measured at full rank in
`output/simulation_kron_study.csv` — also when the interaction is higher-rank.
But what it returns is a number per combination, `interaction$oat[i, j]`, not a
score for an oat and a loading for a pea. MegaLMM returns the interpretable
object and loses the accuracy.

This note is about recovering the interpretable object from the better model.
The code is `code/interaction_decomp.R`, driven by `code/sim_decomp_run.R`, with
the algebra tested in `tests/test_decomp.R` and the BGLR plumbing in
`tests/test_fits.R`.

## Why there is no deregression step

The obvious worry about a second analysis of model output is that the inputs are
shrunken BLUPs: you would need to deregress them and weight each by its
expected error before fitting anything to them.

That worry applies to *fitting a new model*. Nothing is fitted here.
`code/dge_ige_functions.R:392` builds the fitted surface as

```
I_hat = kron_A %*% Beta %*% t(kron_B)
```

where `kron_A = grm_basis(G_oat, rank = q)` is `U_a Λ_a^½` and `Beta` is the
`q × q` reshape of the fitted coefficients. **The surface is already an exact
bilinear form.** So one SVD rewrites it in its own singular basis:

```
double_centre(I_hat) = (C_a kron_A) Beta (C_b kron_B)'      centring is linear
                     = scores %*% diag(d) %*% t(loadings)
```

with `scores` orthonormal over the oats and `loadings` orthonormal over the
peas. That is an algebraic rotation of a term the model already estimated — no
information is created, none is double-counted, and there is no deregression
step to get wrong. Uncertainty comes from doing the same rotation inside each
MCMC draw rather than from correcting a point estimate.

Because both outer factors are orthonormal, `sum(d_k²) = ||interaction||_F²`
**exactly**, so `d_k² / Σd_j²` is an exact variance share rather than an
approximation. `tests/test_decomp.R` checks that, the reconstruction, the
orthonormality, the zero column means, and scale invariance, all as algebra at
`tol = 1e-10`.

Deregression does become necessary one step later — if you regress the recovered
scores on biological covariates to ask what a component *means*. The reliability
and PEV machinery at `code/validation_functions.R:152-172` is what to reuse
then.

## The three numbers, and why the obvious one is wrong

`component_summary()` returns three quantities per component. Choosing the wrong
one as the headline is easy and was caught only by noticing that a point
estimate sat outside its own interval.

| quantity | what it is |
|---|---|
| **`share1`** (+ `share1_lo/hi`) | **The estimand.** Each MCMC draw is a posterior sample of the *true* surface, so the posterior of a draw's own component 1 is the posterior of component 1's share of the real interaction. The only one of the three with an honest interval. |
| `share1_meansurf` | Component 1's share of the posterior-**mean** surface. Convenient, and the only share comparable with MegaLMM (which has no streamed draws) — but **not** an estimate of the above. |
| `proj1_mean` (+ interval) | How much of a typical draw's interaction lies in the direction the mean surface picked out. An identifiability diagnostic, not a share of anything real. |

`share1_meansurf` runs systematically high because singular values are convex in
the matrix: averaging the draws cancels their idiosyncratic directions and
leaves a surface that looks more low-rank than anything real. Measured on the
dense rank-1 positive control:

```
share of component 1          0.448  [95% 0.408 - 0.489]
  mean-surface summary        0.691  (overstates concentration)
  draw variance in that dir.  0.417  (identifiability)
participation ratio           2.06 (mean surface) / 4.46 (draws)
```

A rank-1 truth, and the mean surface says 0.69 while the posterior says 0.45.
Quote the wrong one and the interaction looks half again as concentrated as the
model actually believes.

## Rotation: what can and cannot be interpreted

When two singular values are close, "component 2" is not a well-defined object —
the pair spans a plane and any rotation within it fits equally well. No
alignment scheme fixes that, so:

- **rotation-invariant quantities are primary**: cumulative share over the top
  *k*, and `effective_rank()`'s participation ratio `(Σd²)²/Σd⁴`, which is 1 for
  a rank-1 surface and *r* for *r* equal singular values and needs no threshold.
- per-component scores are **secondary**, to be read with `order_stable` (the
  fraction of draws whose top-3 order matched the reference) and `rot_diag`
  (mean `|diag|` of the block Procrustes rotation; 1 = individually stable,
  0 = only the block as a whole is meaningful) beside them.

`order_stable` deliberately watches only the leading components. Over all of
them it is uninformative: the trailing ones are noise, permute freely, and with
enough of them at least one always swaps, so the statistic reads 0 whatever the
leading structure does.

## The truncation ceiling, which has to be quoted

At `kron_rank = q` the recovered scores are confined **by construction** to
`span(kron_A)`, while `sim_generate.R`'s truth has mass on every direction of
`G`. So there is a hard analytic limit on any recovery number, and
`recovery_ceiling()` computes it from the two spans with nothing fitted.

Recovery must be reported as a fraction of it. From the pilot — two cells at
n = 200, 48% observed, rank-1 truth:

| model | leading score vs truth | ceiling | fraction of ceiling |
|---|---|---|---|
| `dge_ige` | 0.718 | 0.763 | **0.94** |
| `megalmm` | 0.970 | 1.000 | 0.97 |

Read without the ceiling, `dge_ige` looks much worse than MegaLMM at recovering
the scores. Read with it, the two are close and the gap is almost entirely the
truncated basis — which is a property of `SIM_KRON_RANK`, not of the
decomposition. The same confound is why `code/sim_kron_study.R` exists; it now
reaches the decomposition and not just the accuracy.

One consequence worth stating in any write-up: a recovered component is
necessarily a smooth function of kinship, so "biologically interpretable" is
bounded here. That is also why the scores extend to accessions that were never
grown, which is the property the validation experiment needs.

## Why the spectrum needs a null, and MegaLMM needs the same operator

**The nulls.** BGLR puts one scalar prior variance over all `q²` coefficients,
which is misspecified for a low-rank truth and should flatten the recovered
spectrum. The control above shows it: a participation ratio of 4.46 on a
genuinely rank-1 interaction. So the participation ratio has no meaning in the
abstract, only against what it reads when there is no interaction at all.
`sim_decomp_run.R` adds `interaction_pct = 0` scenarios by default to measure
that floor. If the `n_factors = 1` cells come back at the same participation as
the null, the spectrum carries no rank information and the honest deliverable
shrinks to the leading direction plus a caveat.

**MegaLMM.** Its factors are *not* orthogonal — per-factor shares of
`interaction_part(U_F %*% Lambda)` do not sum to one, asserted in
`tests/test_decomp.R` — so "MegaLMM factor 1's share" and "dge_ige component 1's
share" are not comparable quantities. Both models' surfaces therefore go through
the same `decompose_surface()`, which also means the MegaLMM side needs no
`U_F`/`Lambda` at all (and `U_F` is `NULL` when no `envCov` was supplied).

## The BGLR trap, in case it is ever re-plumbed

`BGLR::Multitrait` honours `ETA[[j]]$saveEffects = TRUE` and streams every draw
of that term's coefficients to `<saveAt><term>_beta.bin`. Three things about it:

1. **`saveAt` is a filename prefix, not a directory.** Pointing two array tasks
   at the same prefix makes them overwrite each other.
2. **The file contains the burn-in.** The write sits inside
   `if (iter %% thin == 0)` while the posterior mean accumulates under the
   separate gate `(iter > burnIn) & (iter %% thin == 0)`, and the header records
   `nRow = nIter/thin`. `read_beta_draws()` drops the first
   `floor(burnIn/thin)` rows. Keeping them moved the mean by 0.0994 against a
   3e-16 match when dropped — and an inflated posterior SD with a flattened
   spectrum is *exactly* what the prior-bias analysis is looking for, so this
   failure would have been read as a result.
3. **A truncated file reads back as NA without erroring**, which is what a
   killed job leaves. Hence the file-size assertion.

`tests/test_fits.R` pins all three with one assertion:
`colMeans(draws) == fit$ETA$G_mix$beta`. BGLR's returned `beta` *is* the running
posterior mean over exactly the post-burn-in thinned draws, so that equality can
only hold if the offset, the byte order, the trait-major layout and the storage
mode are all right.

## Running it

```
Rscript code/sim_decomp_run.R --check     # positive control; run first
Rscript code/sim_decomp_run.R --pilot     # two cheap cells, end to end
Rscript code/sim_decomp_run.R             # the design
sbatch -A <account> code/scinet/sim_decomp_array.sbatch
Rscript code/sim_decomp_run.R --combine   # after the array: one complete CSV
```

600 cells plus the nulls, about 120 core-hours, roughly 6 h per task over 20
tasks with a 24 h wall clock. Cells are cached per scenario-replicate and the
job is resumable. Outputs are
`output/simulation_decomp_results.csv` (one row per cell × trait × model) and
`output/simulation_decomp_spectrum.csv` (one row per component).

**`--combine` is the step it is easiest to forget.** Every array task writes
both CSVs from whatever was in the cache when *it* finished, so the files left
behind by a 20-task job are each a partial view. Running `--combine` on the
login node fits nothing, globs the whole cache and rewrites the two CSVs
complete, reprinting the three summary tables. Do it before reading anything.

### Every flag

| flag | default | what it does |
|---|---|---|
| `--check` | — | The positive control: one dense rank-1 cell at `n_acc = 120`, 48% observed, short chain. Prints the diagnostics and `stop()`s if the draws do not average to BGLR's posterior mean, if the leading component carries under 0.3, if the recovered score misses half its ceiling, or if MegaLMM's shares do not sum to 1. Run it first; a negative result from the sweep is only worth having if the pipeline can produce a positive one. |
| `--pilot` | — | Two cheap cells (`rep 1`, `n_acc = 200`, 48% observed), end to end. They are genuine design cells with design seeds, so the full run reuses their cache. |
| `--combine` | — | Rebuild both CSVs from the cache and reprint the tables, fitting nothing. See above. |
| `--refresh` | off | Ignore the cache and refit. Needed after any change to what is scored or stored — otherwise move `DECOMP_SCHEME`. |
| `--reps N` | `SIM_INT_REPS` (5) | Replicates per scenario. Replicate is the slow seed index, so raising it is additive: existing cells stay cached rather than being renumbered. |
| `--rank N` | `SIM_KRON_RANK` (30) | `kron_rank`, the truncation axis — and the one flag most worth varying. It sets the ceiling on recovery, costs `(N/30)^2` in basis columns, and `recovery_ceiling()` reports what it allows. |
| `--n-iter N` | 3000 | MCMC iterations. |
| `--burn-in N` | 600 | Burn-in. `read_beta_draws()` uses it to strip the leading rows BGLR streams but does not average. |
| `--thin N` | 10 | Thinning. Controls the fit and the reader together — they must agree, or the file-size assertion fires. Kept draws are `nIter/N − burnIn/N`. |
| `--filter EXPR` | — | An R expression over the design, e.g. `--filter "sparsity <= 0.09"`. |
| `--task N --ntasks M` | — | Array slicing: task `N` of `M`, taking every `M`th cell. Supplied by the sbatch wrapper. |
| `--no-null-cells` | nulls on | Drop the `interaction_pct = 0` scenarios. **Do not, unless they are already run.** They are the floor the participation ratio is read against, and without them a flattened spectrum cannot be told from a genuinely high-rank one. |
| `--no-megalmm` | on | Skip the MegaLMM pair. Saves about 32 of the 120 core-hours and is a reasonable first stage, but it drops the head-to-head. |
| `--drop-draws` | keep | Delete each cell's `.bin` once summarised. Saves roughly 1 GB over the sweep; the cost is that redoing any decomposition choice — centred or not, a different rank, another alignment scheme — then needs a refit rather than seconds. |

## Scope: simulation only

`code/BGLR_multi_trait_model.R:94` sets `fit_mix_term <- FALSE`, and when TRUE
it builds the **exact** combination kernel inline and calls `BGLR::Multitrait`
directly rather than going through `fit_producer_associate`. So there is no
fitted interaction surface for the real experiment, the exact kernel is not
bilinear and this decomposition does not apply to it, and `save_effects` does
not reach that script — it refuses outright when `kron_rank` is `NA`, rather
than silently producing coefficients that do not reshape. The guard at
`BGLR_multi_trait_model.R:233` notes that over 90% of real combinations are
unreplicated, so the term is near-inestimable there anyway.

Making this reach the real data means switching the production script to
`fit_producer_associate(kron_rank = ...)`, which should wait until the
simulation says the decomposition recovers anything at the densities the real
trials achieve.
