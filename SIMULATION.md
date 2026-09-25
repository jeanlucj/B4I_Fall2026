# SIMULATION

A framework for asking when the MegaLMM factor framework beats the DGE-IGE
framework, and when it does not. For usage see [README.md](README.md); for
structure see [DESIGN.md](DESIGN.md); for the models themselves see
[BACKGROUND.md](BACKGROUND.md); for how accuracy is measured on real data see
[CROSS_VALIDATION.md](CROSS_VALIDATION.md).

## The question

On the B4I data neither framework finds oat × pea interaction, and the two
agree it is not there to find. That leaves the more useful question open:
**under what conditions would either work?** A simulation can answer it,
because the truth is known and the design can be moved.

The two frameworks differ in exactly one thing — how they model the
interaction:

- **DGE-IGE** gives it covariance `G_oat ⊗ G_pea`: full rank, every pair of
  combinations related through both parents' relationships.
- **MegaLMM** gives it a small number of factors: low rank, a few axes along
  which peas differ in how they rank oats.

So the axis that should decide between them is the **true rank** of the
interaction, and the axis that should decide whether either works at all is
**how much of the matrix is observed**.

## The design

Five factors, in `SIM_LEVELS` in `code/sim_config.R`:

| factor | levels | why |
|---|---|---|
| panel size | 200×200, 400×400 | more accessions means a bigger matrix but not more data per cell |
| sparsity | 5%, 15%, 45% observed | see below — these are not the levels first proposed |
| interaction rank | 0, 1, 5 factors | 0 = no interaction; 1 = MegaLMM's home ground; 5 = heading towards Kronecker |
| interaction variance | 10%, 20% of total | is there enough signal to find |
| environments | 1, 10 | one site, or ten each holding a tenth of the combinations |

`2 × 3 × 3 × 2 × 2 = 72`, less the twelve cells where `n_factors = 0` would be
run twice at two interaction variances that scale nothing: **60 scenarios**.

### Why 5 / 15 / 45% and not 3 / 6 / 12%

The first pass used 3 / 6 / 12%, bracketing where the B4I experiment actually
sits (1.1% raw, 3.2% after trimming). A pilot showed all three are deep inside
the region where the factor model has nothing to work with. At 100×100:

| observed | MegaLMM, total genetic value | additive DGE-IGE |
|---|---|---|
| 50% | 0.79 | 0.76 |
| 20% | 0.51 | 0.74 |

The crossover is somewhere between. Three levels that straddle it say more than
three that agree it has not happened yet, so the default is now 5 / 15 / 45%,
with 5% still bracketing the real experiment. `--extended` resolves the axis
more finely (3, 5, 10, 15, 25, 45, 60%) for locating the crossover precisely.

## What is simulated

For oat *i* with pea *j* in environment *k*:

```
y_ijk = mu_k + s_k * ( Pr_i + As_j + I_ij + e_ijk )

Pr  ~ N(0, G_oat * V_Pr)          producer: what the oat does for itself
As  ~ N(0, G_pea * V_As)          associate: what the pea does to the oat
I   = sum_f u_if * lam_jf         interaction, rank n_factors
        u_.f ~ N(0, G_oat),  lam_.f ~ N(0, G_pea)
e   ~ N(0, V_e)
```

### The interaction is built to be fair to both frameworks

Drawing the oat scores from `G_oat` and the pea loadings from `G_pea` makes a
single factor's covariance exactly `G_oat[i,i'] * G_pea[j,j']`. So as
`n_factors` grows the interaction converges on the **full Kronecker structure
the DGE-IGE model assumes**, while at `n_factors = 1` it is as low-rank as
MegaLMM could wish for.

This matters. Simulating the interaction as a factor model and then reporting
that the factor model wins would be circular. Here the rank axis slides between
the two frameworks' home ground, and `n_factors = 0` gives a case where the
correct answer is that neither should find anything.

### Everything else comes from the B4I data

- **Variance shares** — producer 0.211, associate 0.160, residual 0.630, from
  the fitted bivariate model (395.1 / 299.3 / 1181.5 for oat yield). The
  interaction takes its share off the top and the rest is split in these
  proportions.
- **Relationship matrices** — the real `GRM_Avena.rds` and `GRM_Pisum.rds`,
  subsampled to the panel size, so the relatedness structure and its
  unevenness are real rather than idealised.
- **Environment heterogeneity** — environments differ in mean and in spread
  only, so genetic correlations across them stay 1. Scale factors are
  log-normal with SD 0.669, the observed SD of log within-trial SD across the
  six B4I trials. That figure is driven by `B4I_2025_AL`, a near-total crop
  failure (mean 8.3 g/m² against 127–383 elsewhere); excluding it gives 0.167.
  The default is the honest "as observed" value. Set `SIM_ENV_LOG_SD <- 0.167`
  to ask what happens in a season where nothing fails.

### Design constraints that are not cosmetic

Every accession is guaranteed `SIM_MIN_PER_ACC` (3) observations. A pea
environment holding a single oat has no estimable residual variance, and
MegaLMM's ARD sampler answers that with `NaN` several frames deep in
`sample_Lambda_prec_ARD` rather than with an error. Guaranteeing three is also
the more realistic design — nobody grows an entry once — and it keeps
`sparsity` meaning what it says, since the panel is never trimmed afterwards.
Masking for cross-validation carries the same floor.

## What is fitted

Three models, all shown identical data and scored on identical held-out cells:

| model | interaction |
|---|---|
| `additive` | none — producer + associate only. The model used on the real B4I data, and the floor the others must clear. |
| `dge_ige` | the specific-combination term, covariance `G_oat ⊗ G_pea`. The proposal's Eqn 2 in full. |
| `megalmm` | factors, with pea GRM eigenvectors as environmental covariates. Scored on `Eta_mean`, the predicted phenotype. |
| `megalmm_U` | the same fit scored on `U` instead, for comparison. |

plus `oat_mean` and `oat_plus_pea` as margin-only baselines.

### Score `Eta_mean`, not `U`

`U = U_F %*% Lambda + U_R` is a genetic value and **excludes the per-column
intercept**, which is where MegaLMM keeps the pea main effect. Scoring `U`
against a truth containing that effect asks it to predict a component it
structurally cannot hold — worth about 16% of the total variance here. An
earlier version of this framework did exactly that, and MegaLMM's `r_total`
was correspondingly understated; `megalmm_U` is retained as a row so the size
of the difference stays visible.

`r_interaction` is unaffected either way, because row and column means are
stripped from prediction and truth alike.

### MegaLMM settings are swept, and cached apart

`K` (5, 10) and `eigen_variance` (0.20, 0.50, 0.80) are crossed with the data
design, giving 360 fits. Both are claims MegaLMM makes about itself — that the
ARD prior makes surplus factors and irrelevant covariates cheap — so the sweep
tests those claims as much as it tunes anything.

Neither changes the simulated experiment, so refitting the BGLR half for each
would waste five sixths of the compute. The two halves are cached separately:
`<scenario>_rep<k>_bglr.rds` and `<scenario>_rep<k>_mm_K<K>_ev<pct>.rds`.

K matters more than the shrinkage argument suggests, at least at low density.
The binding constraint is **observations per pea column relative to K**: a
column's K loadings are estimated from that column's data, and at 5% observed
with 200 oats there are 10 observations per column. The crossover reported
below at 15% observed is 30 per column, which is 3 x K at K = 10 — so it is
better read as a statement about that ratio than about density as such.

### The Kronecker term is low-rank by necessity

The exact specific-combination kernel is `n_obs × n_obs`, which at 19,200
observations is a 2.9 GB matrix to eigen-decompose. The same space is spanned
by products of the two species' own eigenvectors, so the term is built from the
leading `SIM_KRON_RANK` (30) of each: column `(a-1)*rank + b` is
`A[,a] * B[,b]`, and the fitted coefficients reshape to a `rank × rank` matrix
whose full interaction surface is `A %*% Beta %*% t(B)`. That reshaping is also
what makes predicting every cell cheap.

This is an approximation, and it is the DGE-IGE model's handicap in the
comparison — worth remembering when reading a result where it loses.

The basis and the reshape are derived from scratch, with a worked example, in
[docs/specific-combination_kronecker.md](docs/specific-combination_kronecker.md).
Read it before editing either `grm_basis()` or `kron_basis()`: the index
arithmetic fails silently, and the check that catches it is level S6 of
[EVALUATION_SIMULATION.md](EVALUATION_SIMULATION.md).

## How accuracy is measured

20% of observations held out, with the floor above. Three correlations:

- **`r_total`** — against the true genetic value `Pr + As + I`. What a breeder
  ranking on predicted performance would care about.
- **`r_interaction`** — against the true interaction alone. **The metric the
  comparison turns on.** Neither framework returns an "interaction" on the same
  terms — DGE-IGE has an explicit term, MegaLMM's factors carry main effects
  and interaction together — so row and column means are stripped from the
  whole predicted surface and from the truth alike. An additive model
  residualises to exactly zero, which is the correct answer for it.
- **`r_observed`** — against the held-out observation, noise included. The
  ceiling any model faces in practice.

## Running it

```bash
Rscript code/sim_run.R --check                          # sanity check, always first
Rscript code/sim_run.R                                  # the 60-scenario grid
Rscript code/sim_run.R --reps 5                         # replicated
Rscript code/sim_run.R --extended                       # finer sparsity axis
Rscript code/sim_run.R --filter "n_acc == 200 & n_envs == 1"
Rscript code/sim_run.R --refresh                        # ignore the cache
Rscript code/sim_run.R --trace                          # record a convergence trace
Rscript code/sim_run.R --task 3 --ntasks 20             # one slice, for a job array
Rscript code/sim_run.R --combine                        # rebuild the CSV from the cache
```

On a cluster the grid runs as a SLURM job array; see
[code/scinet/README.md](code/scinet/README.md).

### What it writes

| path | what it is |
|---|---|
| `output/simulation/<scenario>_rep<k>.rds` | one scenario's scores: a tibble with one row per model, the same shape as a slice of the combined CSV. This is the cache — the presence of this file is what lets a scenario be skipped. |
| `output/simulation_results.csv` | every scenario run so far, combined. Columns as described under [Reading the results](#reading-the-results). |
| `output/simulation_summary.png` | interaction accuracy against sparsity, faceted by panel size and interaction rank. |
| `output/simulation_runs/` | scratch. MegaLMM needs its run state on disk, so each scenario gets a subdirectory here and it is deleted as soon as the scenario's scores are cached. The parent directory stays behind; nothing in it is worth keeping. |

Everything under `output/` is gitignored, and everything here regenerates.

Because the cache is per scenario, the grid can be run in pieces, interrupted
and resumed. Re-running picks up where it stopped; `--refresh` ignores the
cache and refits.

### How long it takes

Measured at 200 × 200, mean seconds per scenario over all three fits:

| observed | observations | additive | dge_ige | megalmm | per scenario |
|---|---|---|---|---|---|
| 5% | 2,000 | 2.8 | 10.0 | 10.2 | 23 s |
| 15% | 6,000 | 7.3 | 28.6 | 10.0 | 46 s |
| 45% | 18,000 | 20.4 | 70.6 | 9.8 | 101 s |

The thirty 200 × 200 scenarios take about half an hour; `n_envs` does not
change the amount of data, only how the standardisation is grouped, so it
costs nothing. The thirty at 400 × 400 carry four times the observations at a
given sparsity and should take around two and a half hours, putting **one
replicate of the full grid at roughly three hours**. That figure is
extrapolated from the 200 × 200 timings rather than measured, so treat it as
give or take half.

Note that `megalmm` is flat across sparsity while the two BGLR models are not:
MegaLMM's cost is driven by the size of the matrix and `K`, not by how much of
it is filled.

The cell most likely to give trouble is 400 × 400 at 45%, where the Kronecker
design matrix is 72,000 × 900 — 518 MB, with two temporaries of that size
built before they are multiplied. If the full grid fails anywhere it will be
there, on memory rather than time. Halving `SIM_KRON_RANK` to 20 cuts that
term to 400 columns.

### Convergence: checked, not swept

`--trace` splits the sampling into `SIM_TRACE_CHUNKS` pieces and scores after
each, so a single run reports whether accuracy was still moving when the chain
stopped — no chain-length axis needed.

**Caveat: the trace is not yet trustworthy.** Its values do not reconcile with
the final posterior mean from the same run (one case: −0.116 at the last chunk
against 0.087 scored at the end). That points at how `save_posterior_chunk()`
and `load_posterior_param()` accumulate a `posteriorMean` parameter across
chunks — plausibly each chunk's mean rather than the running mean — which
would make the trace a sequence of chunk estimates rather than a convergence
curve. Useful for spotting drift, not for reading off a final number, and it
needs verifying before either use.

### `--check` exists for a reason

It fits a dense, strongly structured scenario in which MegaLMM should clearly
recover the signal, and stops if it does not. The failure mode it guards
against — a wiring mistake that makes every MegaLMM number approximately
zero — is indistinguishable by eye from "the method does not work here", and
this project has already spent time on a MegaLMM result that looked like the
latter. Run it before trusting any sweep.

The check is also a result in miniature. On a 100×100 matrix at 50% observed
with a single interaction factor:

| model | `r_total` | `r_interaction` | seconds |
|---|---|---|---|
| additive | 0.731 | — | 3 |
| dge_ige | **0.787** | 0.444 | 19 |
| megalmm | 0.705 | **0.869** | 5 |
| oat_plus_pea | 0.715 | — | — |

Read that carefully, because it is the shape of the whole answer. MegaLMM is
far better at the **interaction** — 0.87 against 0.44 — which is what it is
for. DGE-IGE is better at the **total**, because it models the additive main
effects explicitly while MegaLMM absorbs the pea main effect into a per-column
intercept that never reaches `U`. A framework can win the question it was built
for and still lose the one the breeder asks.

## Reading the results

`simulation_results.csv` has one row per scenario × replicate × model, with
`r_total`, `r_interaction`, `r_observed`, `seconds`, `n_train` and `n_held`.
The runner prints three summaries: accuracy for the total, accuracy for the
interaction, and a head-to-head of MegaLMM minus DGE-IGE.

Things to look for, given what the pilot already shows:

- **Sparsity should dominate everything.** Find the density at which
  `r_interaction` for MegaLMM lifts off zero; that is the number a redesign
  has to hit. Why it dominates, and why MegaLMM loses to a row average rather
  than merely tying with it, is worked through in
  [docs/MegaLMM_sparsity_challenge.md](docs/MegaLMM_sparsity_challenge.md).
- **Rank should decide the winner.** At 1 factor MegaLMM should win the
  interaction comfortably; at 5 it should narrow as the truth approaches the
  Kronecker structure DGE-IGE assumes.
- **`n_factors = 0` is the false-positive check.** `r_interaction` is
  undefined there — there is no true interaction to correlate against, and the
  scorer returns `NA` rather than a number. The readout is instead that
  neither interaction-fitting model should beat `additive` on `r_total`.
- **Panel size is not sample size.** 400×400 at a given sparsity has four
  times the observations but the same number per cell. Whether that helps
  separates "needs more data" from "needs more data per combination".
- **Ten environments should cost something.** Standardising within
  environment removes the heterogeneity but spends degrees of freedom to do
  it, and with a tenth of the combinations each, that estimate is noisy.

## Results so far

**These numbers predate the `Eta_mean` fix and the K sweep**, and were produced
with MegaLMM scored on `U` at a fixed K = 10. They are kept because the
interaction comparison is unaffected, but `r_total` for MegaLMM is understated
throughout and the crossover location is conditional on K = 10. Re-running the
grid will replace them.

The 200 × 200 single-environment slice, one replicate per cell
(`--filter "n_acc == 200 & n_envs == 1"`). The full grid adds the panel-size
and environment axes.

### Recovering the interaction — the crossover is at about 15% observed

| observed | rank | interaction var | dge_ige | megalmm | MegaLMM gain |
|---|---|---|---|---|---|
| 5% | 1 | 10% | **0.341** | 0.112 | −0.229 |
| 15% | 1 | 10% | 0.440 | **0.475** | +0.035 |
| 45% | 1 | 10% | 0.510 | **0.888** | +0.378 |
| 5% | 1 | 20% | **0.044** | −0.080 | −0.125 |
| 15% | 1 | 20% | 0.519 | **0.803** | +0.284 |
| 45% | 1 | 20% | 0.522 | **0.947** | +0.425 |
| 5% | 5 | 10% | **0.113** | −0.026 | −0.139 |
| 15% | 5 | 10% | **0.367** | 0.266 | −0.101 |
| 45% | 5 | 10% | 0.435 | **0.633** | +0.198 |
| 5% | 5 | 20% | **0.254** | 0.027 | −0.227 |
| 15% | 5 | 20% | 0.327 | **0.383** | +0.056 |
| 45% | 5 | 20% | 0.497 | **0.794** | +0.297 |

Three things, all as the design predicted:

- **Sparsity decides whether MegaLMM works at all.** At 5% it is at or below
  zero in every cell; by 45% it is at 0.63–0.95. DGE-IGE degrades far more
  gracefully — 0.04–0.34 at 5%, 0.44–0.52 at 45% — because a Kronecker kernel
  borrows through relatedness rather than needing a column's own data.
- **The crossover sits near 15% observed**, which is where the two are within
  a few hundredths of each other in three of four cells. That is roughly five
  times the density of the B4I experiment.
- **Rank decides the size of the win.** At rank 1 MegaLMM's advantage at 45%
  is +0.38 to +0.43; at rank 5 it falls to +0.20 to +0.30, as the truth moves
  towards the full Kronecker structure DGE-IGE assumes.

More interaction variance helps MegaLMM disproportionately: at 15% observed
and rank 1, going from 10% to 20% interaction variance takes it from 0.475 to
0.803 while DGE-IGE moves 0.440 to 0.519.

### The breeder-facing number tells a different story

`r_total`, against the true genetic value:

| observed | rank | interaction var | additive | dge_ige | megalmm |
|---|---|---|---|---|---|
| 5% | 0 | — | **0.909** | 0.898 | 0.137 |
| 45% | 0 | — | **0.978** | 0.975 | 0.732 |
| 5% | 1 | 20% | 0.670 | 0.654 | −0.011 |
| 45% | 1 | 20% | 0.753 | 0.825 | **0.828** |
| 45% | 5 | 20% | 0.763 | **0.825** | 0.739 |

DGE-IGE or the plain additive model wins almost everywhere. MegaLMM only draws
level in the single most favourable cell — densest, lowest rank, most
interaction variance. The reason is structural and was visible in `--check`:
MegaLMM absorbs the pea main effect into a per-column intercept that never
reaches `U`, so it is giving away a variance component that the DGE-IGE model
estimates explicitly.

**A framework can win the question it was built for and still lose the one the
breeder asks.** If the goal is ranking oats on expected performance, DGE-IGE is
the better choice at every density tested. If the goal is understanding which
specific oat × pea combinations do something unusual, MegaLMM is much better —
but only once the matrix is dense enough.

### The false-positive check passes

At `n_factors = 0` neither interaction-fitting model beats the additive model
on `r_total` (0.898 and 0.136 against 0.909 at 5% observed; 0.975 and 0.732
against 0.978 at 45%). Fitting an interaction that is not there costs a little
and gains nothing, which is the correct behaviour.

### What this says about B4I

The real experiment sits at 1.1% observed, 3.2% after trimming — below the
lowest level simulated, where MegaLMM is already at or below zero and even
DGE-IGE recovers little. The simulated crossover near 15% is roughly five times
the current density. Reaching it means replicating specific oat × pea pairs
rather than spreading singleton combinations thinly, which is the same design
change the bivariate model needs for `σ_PrAs`.

## Caveats

- The interaction is generated as a factor model with genetic loadings. At
  `n_factors = 5` that is close to, but not identical to, the Kronecker
  structure DGE-IGE assumes, and the approach to it is asymptotic in rank. A
  cleaner test of DGE-IGE's home ground would draw the interaction directly
  from `G_oat ⊗ G_pea`; the generator would need one more branch.
- The Kronecker term is rank-limited (see above), so `dge_ige` is not fitting
  the exact model it claims. `SIM_KRON_RANK` controls this and is worth a
  sensitivity check at small sizes where the exact kernel is tractable.
- Environments carry variance heterogeneity only. Genuine genotype ×
  environment rank changes would be a second kind of interaction and would
  confound the axis the simulation is built around.
- The convergence trace does not reconcile with the final posterior (above).
- Only the oat orientation is fitted. MegaLMM gives the oat axis a relationship
  matrix and the pea axis only a per-column intercept and factor loadings, so
  the pea side is structurally under-specified relative to the DGE-IGE model,
  which shrinks both species through their GRMs. Part of MegaLMM's deficit on
  `r_total` is that asymmetry rather than the method.
  [BOTH_ORIENTATIONS.md](BOTH_ORIENTATIONS.md) sets out what to do about it.
- Only oat yield is simulated. The bivariate model's producer–associate
  covariance, which is what the real analysis is ultimately after, is not in
  this framework at all.
