# CROSS-VALIDATION

How predictive accuracy is measured here, what the numbers mean, and what is
not yet validated. For usage see [README.md](README.md); for structure see
[DESIGN.md](DESIGN.md); for the reasoning behind the models see
[BACKGROUND.md](BACKGROUND.md); for the simulated comparison of the two
frameworks see [SIMULATION.md](SIMULATION.md).

## Status, in one table

| model | cross-validated? | where |
|---|---|---|
| MegaLMM factor model (oat × pea matrix) | **yes** | `code/megalmm_oat_pea.R` |
| Bivariate DGE-IGE model (BGLR) | **no** | — |

Only the MegaLMM analysis is cross-validated. The bivariate model reports
posterior variance components and per-accession effects, and it checks that its
chains agree, but it has never been asked to predict something it was not shown.
The distinction matters and is spelled out under
[What is not validated](#what-is-not-validated).

This is worth stating plainly because Sub-objective 1.4 of the proposal asks for
exactly the missing half: *"Estimate the accuracy of genomic prediction of
intercrop-related breeding values of each crop and the ability to predict overall
performance of specific oat-pea combinations using cross-validation."* The first
clause — breeding values per crop — has no cross-validation yet. The second —
specific combinations — is what the MegaLMM cross-validation below measures, and
it is measured as failing.

---

## The MegaLMM cross-validation

### What is being predicted

The matrix is oat accessions × pea "environments", where a cell holds the oat
grain yield of that oat grown with that pea. A cell is either observed or the
pair was never grown together. The question is whether the model can fill in a
cell it was not shown — the intercrop equivalent of the **CV2** scenario in
multi-environment trial work: predict a missing genotype × environment
combination from other cells in the same rows and columns.

This is a test of **specific-combination** prediction. It is not a test of
whether an oat is a good oat, which is why the baselines below matter so much.

### Masking

`mask_cells()` in `code/megalmm_oat_pea.R`.

Held-out cells are drawn from the **filled** cells only, never from the empty
ones, so every held-out cell has a value to score against. The default is 20%,
repeated with 3 different seeds (`cv_fraction`, `cv_reps`, `cv_seed`).

Masking is constrained, not purely random. Cells are visited in random order and
held out only while **both** their row and their column would keep more than
`floor_obs` (default 2) training cells. Without that floor an oat or pea can be
emptied entirely, and that does not produce a graceful failure: MegaLMM's ARD
sampler returns `NaN` and the run dies several frames deep in
`sample_Lambda_prec_ARD` with "missing value where TRUE/FALSE needed". The same
constraint is why the matrix is trimmed before any of this begins (≥ 3 oats per
pea, ≥ 2 peas per oat — see [BACKGROUND.md](BACKGROUND.md)).

The floor costs nothing in practice here: all 18 fits held out exactly 336 of
1,682 cells, the full 20% asked for.

A consequence worth being explicit about: because rows and columns always retain
some data, this measures **interpolation** into a sparse matrix, not
extrapolation to an oat or pea never seen. Predicting a wholly new pea
environment is the CV0 scenario, which the machinery supports — pass `newEnv` to
`setup_megalmm_state()` and the `U_CV0` / `Eta_CV0` posterior functions appear —
but it is not part of this sweep.

### What is scored

Three MegaLMM predictions and three baselines, all against the same held-out
cells, by Pearson correlation and RMSE.

**MegaLMM targets** — all three are reported rather than one being chosen,
so that a weak result cannot be dismissed as the wrong target having been picked:

| target | what it is |
|---|---|
| `MegaLMM_Eta` | `Eta_mean`, the predicted phenotype. The like-for-like comparison with the baselines, since it carries the per-environment intercept. |
| `MegaLMM_U` | `U_F %*% Lambda + U_R` — genetic value, factor part plus environment-specific genetic residual. |
| `MegaLMM_U_noU_R` | `U_F %*% Lambda` — the factor part alone. |

The two `U` matrices are genetic values and omit the per-column intercept, so
across columns they are at a disadvantage. They are scored anyway because they
are the quantities a breeder would rank on.

**Baselines** — `margin_baselines()`. These are the point of the exercise:

| baseline | prediction for cell (oat *i*, pea *j*) |
|---|---|
| `oat_mean` | mean of oat *i*'s training cells |
| `pea_mean` | mean of pea *j*'s training cells |
| `oat_plus_pea` | `oat_mean + pea_mean − grand_mean` — additive margins, no interaction |

`oat_plus_pea` is the honest null for this question. It is the best you can do
knowing only that some oats are better than others and some peas are easier
partners than others, with **no** genotype-by-genotype term at all. A factor
model exists to beat that. If it does not, it has not found interaction
structure, whatever its internal fit looks like.

An oat or pea whose training cells all vanished would fall back to the grand
mean; the floor above means this does not arise.

### The comparison grid

`3 fills × 2 covariate settings × 3 replicates = 18 fits`, each a full burn-in
and sampling run (5 × 60 burn-in iterations, 250 sampling; the final fit uses
8 × 100 and 400). The three factors are crossed deliberately:

- **fill** (`raw` / `centered` / `standardized`) — settles empirically what
  belongs in a cell rather than arguing it.
- **covariates on/off** — isolates what the *extended* MegaLMM contributes. With
  ~7 oats per pea environment, this is the question of whether covariates are
  carrying the model.
- **replicates** — different masking seeds, so the spread across replicates
  shows how much of any difference is noise.

### Running it

```bash
Rscript code/megalmm_build_inputs.R    # writes all three fills
Rscript code/megalmm_oat_pea.R         # 18 CV fits + the final fit, ~20 min
```

Outputs:

- `output/megalmm_cv_results.csv` — one row per fill × covariates × replicate ×
  predictor, with `r`, `rmse` and `n_held`.
- `output/megalmm_cv_accuracy.png` — the same, mean ± se.
- `output/megalmm_fit_<fill>.rds` — the final fit, with the CV table attached.

To change the design, edit the settings block: `cv_fraction`, `cv_reps`,
`cv_seed`, `fills`, and the chain lengths `cv_burn_rounds` / `cv_burn_iter` /
`cv_sample_iter`.

### How to read the result

Mean correlation with held-out cells, 3 replicates at 20% masking:

| cell value | covariates | best MegaLMM | `oat_mean` | `oat_plus_pea` |
|---|---|---|---|---|
| centered | no | 0.031 | 0.259 | **0.279** |
| centered | **yes** | **0.072** | 0.259 | **0.279** |
| raw | no | 0.010 | 0.024 | 0.009 |
| raw | yes | 0.040 | 0.024 | 0.009 |
| standardized | no | 0.021 | 0.201 | 0.210 |
| standardized | yes | 0.002 | 0.201 | 0.210 |

Two readings, both supported:

- **The covariates work.** Turning them on roughly doubles MegaLMM's accuracy in
  every filling. That is the extended part of the framework doing its job.
- **The model does not.** 0.072 against 0.279 means the additive margins beat the
  factor model by a wide margin. There is no recoverable oat × pea interaction
  signal in this design.

The second reading is corroborated from a different direction: 91% of oat × pea
combinations occur in a single plot, which is also why the bivariate model leaves
its specific-combination term out. Two methods, one answer.

**A caveat on the baselines.** `oat_mean` reaching 0.26 does *not* mean oat
main effects explain 26% of anything — the held-out cells are trial-centered
yields, so this is the correlation between an oat's average deviation and a
single new plot's deviation, residual variance included. It is a floor for the
prediction problem, not an estimate of heritability.

---

## What is not validated

### The bivariate DGE-IGE model has no cross-validation

`code/BGLR_multi_trait_model.R` fits the model and reports posterior variance
components, producer and associate effects per accession, and credible intervals.
None of that is out-of-sample.

It does run four chains from different seeds and report the rank correlation of
accession GMA between them (0.999 throughout). **That is a convergence
diagnostic, not cross-validation.** It says the sampler explores the same
posterior from different starts; it says nothing about whether the model predicts
an accession it has not seen. Four chains agreeing perfectly on a wrong answer is
entirely possible.

### What a cross-validation of it would look like

Sub-objective 1.4 asks for accuracy of *intercrop-related breeding values of each
crop*, which is a different hold-out from the one above: mask whole
**accessions**, not cells.

A design that would answer it:

- **Mask accessions, not plots.** Drop every plot of a held-out oat accession,
  fit, then predict its producer and associate effects from the GRM through its
  relatives — the genomic-prediction question proper. Repeat for pea.
- **Score producer and associate separately.** Correlate predicted effects with
  the effects estimated when the accession *is* in the model. They are different
  quantities with different accuracies, and the associate effect is the harder
  and more interesting one.
- **Fold by relatedness, not at random.** Random folds leak: a clonal group or a
  full-sib family split across folds makes prediction look better than it will be
  for a genuinely new line. Folds should keep families together, which the
  curation output (`oat_analysis_names.csv`, `oat_family_correlations.csv`) makes
  straightforward.
- **Baseline against the trial mean and the GMA mean**, for the same reason the
  MegaLMM sweep does: an accuracy figure with nothing to beat is not informative.

Two things to expect from the design as it stands. 99 oat and 84 pea accessions
have a single partner, so their producer and associate effects are aliased and
will be predicted badly no matter what; they should be reported separately rather
than pooled into one accuracy. And with the specific-combination term omitted,
this would validate GMA-type effects only — the ability to predict *specific*
combinations is the thing the MegaLMM sweep already measures, and already finds
absent.

This is not implemented. See [docs/B4I_followups.md](docs/B4I_followups.md).
