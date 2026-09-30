# SIMULATION GLOSSARY — every column in `simulation_results.csv`

What each `model` is and what each `r_*` column measures. For the design the
columns come from see [SIMULATION.md](SIMULATION.md); for stepping through the
code that produces them see
[EVALUATION_SIMULATION.md](EVALUATION_SIMULATION.md).

Everything here is defined by `score_predictions()` in `code/sim_fit.R` and, for
the MegaLMM diagnostics, by `run_scenario_megalmm()` in the same file.

---

## The one thing to understand first

**Every model returns two full oat × pea surfaces**, one per trait, both laid out
oat-rows × pea-cols. A "surface" is a prediction for every cell of the matrix,
including cells never observed. The pea-side MegaLMM fit is transposed on the way
out so that all six models are in the same orientation.

Every number below is then a correlation computed from those two surfaces. That
uniformity is the point: it asks three quite different frameworks the same
question in the same way, rather than reading each model's own idea of what it
estimated.

Which margin carries which effect follows from the layout, and is the single
easiest thing to get backwards:

| surface | row margin (`rowMeans`) | column margin (`colMeans`) |
|---|---|---|
| **oat yield** | the **oat's producer** effect | the **pea's associate** effect |
| **pea yield** | the **oat's associate** effect | the **pea's producer** effect |

Read "producer" as *what this accession does for its own yield* and "associate"
as *what it does to its partner's yield*. So an oat's associate effect is read on
**pea** yield — which is why both yields have to be simulated, and why MegaLMM
has to be run in both orientations.

---

## The models (first column)

Six rows per data scenario, plus two more for each MegaLMM setting the design
assigned to that scenario.

| `model` | what it is |
|---|---|
| `additive` | The bivariate `BGLR::Multitrait` producer–associate model **without** a specific-combination term. Four genetic effects — oat producer, oat associate, pea producer, pea associate — with kinship on all four through `Σ ⊗ G`, and a residual covariance between the two yields. Its surfaces are pure outer sums, so it has no interaction component at all. |
| `dge_ige` | The **same** `Multitrait` fit **with** the low-rank specific-combination term: the Khatri–Rao basis of the two GRMs, truncated at `SIM_KRON_RANK = 30` per species and so 900 columns. This is the production model — `fit_producer_associate()` from `code/dge_ige_functions.R`, the one the real analysis runs, not a stand-in. It is the only BGLR model with a non-`NA` `r_int_*`. |
| `row_mean` | Baseline. Each accession's own mean over its training plots, and **nothing else** — no borrowing across relatives, no partner information. On the oat-yield surface every column is identical, so it has no column effect by construction and its column-margin metrics come back `NA` rather than poor. |
| `both_means` | Baseline. Both margins: `outer(row mean, column mean, "+") − grand mean`. This is the floor a model has to clear, and at high sparsity it is a surprisingly hard one — at 48% observed it ties `dge_ige` on the main effects. |
| `megalmm` | The MegaLMM factor model, fitted in **both** orientations and scored on `Eta_mean`, the predicted **phenotype**. |
| `megalmm_U` | The same fit scored on `U`, the **genetic value**. |

### Why `megalmm` and `megalmm_U` are both reported

They differ by the per-column intercept, and that difference is the most
informative thing in the table.

MegaLMM has a per-column intercept and **no per-row one**. So in either
orientation:

| | row species | column species |
|---|---|---|
| where its effect lives | latent factors + `U_R`, **kinship-shrunk** | the per-column intercept, **fixed and unshrunk** |
| which effect that is | **producer** | **associate** |

`U` excludes the per-column intercept, so scoring `U` against a truth that
contains the column species' main effect asks the model to predict a component it
structurally cannot hold. That is why `megalmm_U` typically shows a near-zero or
negative `r_*_assoc` while `megalmm` does not — it is a property of the
parameterisation, not a failure of the fit. `megalmm` is the row to quote;
`megalmm_U` is there to show where the effect went.

---

## Per-trait accuracy, at the held-out cells

Four columns per trait, computed **only at the 20% of observed cells held out**
(`idx` in the code). Suffix `_oat` means the oat-yield surface, `_pea` the
pea-yield surface.

| column | correlation between | what it answers |
|---|---|---|
| `r_total_oat` | predicted oat-yield surface, and **true producer + associate + interaction** | Everything the model could in principle know about a held-out cell. The headline accuracy. |
| `r_gma_oat` | `additive_part(surface)`, and **true producer + associate** | Accuracy on the additive part alone, with the interaction residualised out of both sides. |
| `r_int_oat` | `interaction_part(surface)`, and `interaction_part(true I_oat)` | Accuracy on the specific-combination part alone. **`NA` for `additive`, `row_mean` and `both_means`**, which have no interaction: an exactly additive surface residualises to exactly zero, so the correlation is undefined rather than small. |
| `r_obs_oat` | predicted surface, and the **standardised observed yield** | The only column comparable to a real-data cross-validation number, because it is scored against a noisy phenotype rather than against a noise-free truth. Bounded above by the square root of heritability, so it is always much lower than `r_total_oat` and the two must never be quoted as if they were the same quantity. |

`r_total_pea`, `r_gma_pea`, `r_int_pea` and `r_obs_pea` are the exact mirror on
the pea-yield surface.

**The truth these are scored against is the stable effects plus the
interaction.** Environment-specific deviations are unpredictable by construction
and so count against every model equally.

`additive_part(M)` is `M` reduced to its two margins; `interaction_part(M)` is
what is left, `M − additive_part(M)`, which has exactly zero row and column
means. The two sum to `M` to machine precision — `tests/test_surface.R` asserts
it at `1e-12`, because every column in this section is a margin or a residual of
that decomposition and a wrong split would corrupt all of them together.

---

## Effect recovery, over the whole panel

These are **not** restricted to held-out cells: they are margins of the full
surface, compared with the simulated effect for every accession. They answer
"would this model have ranked the accessions correctly", which is the question
breeding actually asks.

| column | predicted | true |
|---|---|---|
| `r_oat_prod` | `rowMeans(oat surface)` | oat producer effect |
| `r_pea_assoc` | `colMeans(oat surface)` | pea associate effect |
| `r_oat_assoc` | `rowMeans(pea surface)` | oat associate effect |
| `r_pea_prod` | `colMeans(pea surface)` | pea producer effect |
| `r_oat_gma` | `rowMeans(oat surface) + rowMeans(pea surface)` | oat producer **+** oat associate |
| `r_pea_gma` | `colMeans(pea surface) + colMeans(oat surface)` | pea producer **+** pea associate |

`r_*_gma` is **general mixing ability**: what an accession contributes to the
total, its own yield plus its effect on its partner's. Each species' GMA is
assembled across the two orientations — one component from each — which is the
main reason both orientations are fitted.

### The prediction to test here

Because the row species' effect is kinship-shrunk and the column species' is a
fixed unshrunk intercept, **`r_*_assoc` should degrade faster than `r_*_prod` as
sparsity falls**. The 1.6% level is where it should be most visible. If it does
not, the reasoning in the table above is wrong and worth knowing.

---

## MegaLMM diagnostics

Present on `megalmm` and `megalmm_U` rows only; `NA` elsewhere.

| column | what it is |
|---|---|
| `r_mainfactor_oat` | Correlation between the **first latent factor's scores** on the oat-side fit and the true oat producer effect. When `fixed_main_effect = TRUE` the factor is pinned and the sign is meaningful, so the value is signed; when the factor is free its sign is arbitrary, so the **absolute** value is reported. |
| `r_mainfactor_pea` | The same on the pea-side fit, against the true pea producer effect. |
| `r_rowmean_oat` | Correlation between each oat's **training row mean** of standardised oat yield and its true producer effect. The no-model reference for the same quantity. |
| `r_rowmean_pea` | The same for pea. |

**These are diagnostics of the pinned factor, not estimates of a main effect.**
A main effect is not confined to factor 1 — other factors pick some of it up —
so the estimate to use is `r_oat_prod` / `r_pea_prod` above, which take the
margin of the whole surface and therefore collect the main effect wherever the
model put it. `r_mainfactor_*` exists to answer a narrower question: did pinning
the first factor actually make it carry the main effect?

Comparing `r_mainfactor_*` with `r_rowmean_*` is the useful reading. The row mean
is free and uses no model; if the pinned factor cannot beat it, the factor
machinery is not earning its keep for that quantity.

---

## Everything else in the table

| column | meaning |
|---|---|
| `seconds` | Wall-clock time for that fit. `NA` on rows that share a fit with the row above — the baselines, and `megalmm_U`. |
| `n_train`, `n_held` | Plots in the training and held-out sets, per trait. |
| `n_dropped` | Accessions MegaLMM dropped for having too few observations. |
| `fixed_ok_oat`, `fixed_ok_pea` | Whether the fixed main-effect factor was successfully pinned in that orientation. `FALSE` makes `r_mainfactor_*` uninterpretable for that row. |
| `scenario` | The data scenario's name, `n<acc>_sp<sparsity×1000>_f<n_factors>_i<interaction_pct>_e<n_envs>_g<gxe_cor×100>`. Encodes only the design, never the seed. A name **without** the trailing `_g###` is from the pre-September-2026 single-trait grid and does not belong in the same table. |
| `rep` | Replicate index. |
| `seed` | The RNG seed that generated this scenario's data. See the note below. |
| `n_acc`, `sparsity` | Panel size per species, and the fraction of the oat × pea matrix observed. |
| `interaction` | Composite level: `none`, `f1_i10`, `f1_i20`, `f5_i10`, `f5_i20` — interaction rank and the share of variance it carries. |
| `environment` | Composite level: `one`, `ten_stable`, `ten_gxe`. |
| `n_factors`, `interaction_pct`, `n_envs`, `gxe_cor` | The underlying levels the two composites decode to, as the generator receives them. |
| `K`, `eigen_variance`, `fixed_main_effect` | The three MegaLMM levers. **`NA` on the BGLR and baseline rows**, which do not depend on them and are fitted once per data scenario — so any filter that keys on these columns for every row will delete the entire comparator side of the experiment. |

`seed` identifies which draw produced the row, and is what distinguishes two
cache files that would otherwise both call themselves the same `(scenario, rep)`.
It is `SIM_BASE_SEED + (rep - 1) * n_scenarios + scenario_index`, so it changes
when the grid changes but **not** when more replicates are requested — which is
what makes `--reps` additive.

> **Results generated before 30 September 2026 do not have this column**, and
> were produced under a numbering that varied `rep` fastest. In those tables two
> different draws can both appear as `rep 1` and cannot be told apart; the
> September run has 138 such rows. `code/sim_filter_results.R` detects and counts
> them, and recovers the seed from the cache filename where the column is
> missing, but it cannot resolve a collision in a CSV that never recorded one.

---

## Reading the table

Because the MegaLMM sweep is a **D-optimal fraction**, most scenarios carry only
some of the settings, so a table of cell means compares unlike with unlike — and
taking a per-scenario maximum over settings reintroduces the selection bias the
fraction was meant to remove (S11 in the evaluation checklist). Fit the design
model instead: main effects and two-way interactions in the composite factors,
on each outcome column. `sim_run.R` writes
`output/simulation_design_effects.csv` as a first pass.

Before any of that, run `code/sim_filter_results.R`: the results CSV is rebuilt
by globbing the cache, so rows from an earlier grid join it silently.
