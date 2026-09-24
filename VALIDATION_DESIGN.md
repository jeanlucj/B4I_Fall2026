# VALIDATION TRIAL — design

How we propose to test whether the producer and associate effects estimated
from the B4I intercrop trials are real. For the project as a whole see
[README.md](README.md); for the model these effects come from see
[BACKGROUND.md](BACKGROUND.md); for what has and has not been validated so far
see [CROSS_VALIDATION.md](CROSS_VALIDATION.md).

**Generated from the fit in `output/`** by `code/validate_pool_selection.R`,
`code/validate_power.R` and `code/validate_design.R`. Every number below is a
snapshot of one data vintage, not a fixed parameter — see §7. Numbers here are
from the **6-trial fit, vintage 2026-09-24**.

---

## 1. What we are testing, and why it sets the whole design

We have estimated, for ~440 oat and ~420 pea accessions, a **producer effect**
(`Pr`, what an accession does for its own yield) and an **associate effect**
(`As`, what it does to its partner's yield). Only a tiny fraction of the
oat × pea combination space was ever sown, so these rest heavily on genomic
covariance rather than on direct observation, and nothing has tested them.

The validation asks one question: **in new plots, does a set of accessions we
predict to have high associate effects actually outperform a set we predict to
have low ones?**

Finding accessions with a *low producer* effect is not interesting — poor
lines are easy to find. So both pools are drawn from **high-producer**
accessions and contrasted only on `As`. The claim being tested is therefore
the useful one: *among good lines, can we tell the good neighbours from the
bad ones?*

Three structural facts follow, and they determine everything else.

**An associate effect is read on the partner's yield.** The oat contrast is
measured on **pea yield**; the pea contrast on **oat yield**.

**So one set of plots validates both species.** Every plot carries an oat and
a pea and both yields are recorded. A 2 × 2 factorial of oat-pool × pea-pool
gives each species' test the *full* plot budget instead of half of it.

**A producer imbalance between the pools is not a confound.** `Pr` acts on own
yield and `As` on partner yield, so unequal pool `Pr` would bias the *other*
trait's test, not this one. We still match them, because "these are all good
producers" is part of the claim — but it is a presentational requirement, not
a statistical one.

---

## 2. Building the pools

### The obvious rule does not work

The natural rule — take the top *n* on `Pr + As` for the As+ pool and the top
*n* on `Pr − As` for the As− pool — **does not produce disjoint pools**. A
sufficiently high-`Pr` accession makes both lists whatever its `As`:

| n per pool | oat accessions in *both* pools |
|---|---|
| 10 | 1 |
| 20 | 8 |
| 30 | 13 |
| 50 | 20 |

At n = 30 that is 13 of 30 shared, which guts the contrast.

### What we do instead

A constrained optimisation: **maximise `mean(As+) − mean(As−)` subject to**

1. the pools being disjoint,
2. both pool means of `Pr` above a threshold (the median of eligible
   candidates),
3. `|mean Pr+ − mean Pr−|` within 1 g/m².

Implemented with a Lagrange multiplier on `Pr`, pushed in opposite directions
for the two pools and bisected until the producer difference crosses zero
(`build_pools()` in `code/validation_functions.R`).

**Candidates are filtered first.** 99 oat and 84 pea accessions appear with
only a single partner, and 8 oat / 20 pea have no marker data at all. Their
effects are almost pure shrinkage, so accessions with **fewer than 3 distinct
partners are excluded** — leaving 265 oat and 254 pea, ample for pools of
15–30.

### What it delivers

| species | n per pool | ΔAs (g/m²) | ΔPr | mean Pr, As+ | mean Pr, As− | population mean Pr |
|---|---|---|---|---|---|---|
| oat | 20 | **14.3** | 0.15 | 14.0 | 13.8 | 0.06 |
| pea | 30 | **16.5** | 0.01 | 10.5 | 10.5 | 0.26 |

The producer constraint is essentially free: it costs about 1 g/m² of
associate contrast and buys two pools whose producer means differ by less than
0.2 g/m² while both sit 10–14 g/m² above the population mean.

---

## 3. The field design

**2 × 2 factorial**, oat pool × pea pool, all four cells equally represented
at **every** location. Complete balance at each site means losing one entirely
— as AL was lost to crop failure in 2025 — costs a fifth of the plots and
confounds nothing.

**Pairing within a cell** is by permutation, redrawn at each location, so an
accession meets a different partner every time. Averaging an accession over
many partners is what makes its pool mean precise; it is worth more than
repeating combinations.

**Anchors.** One combination per cell is sown **twice at each location**, the
same combinations everywhere — 8 plots per location, 40 in total. Within-
location duplication is the only thing that separates plot error from
combination × location; repeating a combination across locations alone
confounds the two. This also buys a first estimate of specific-combination
variance, the quantity the project has concluded is missing from the existing
design.

**Blocks.** Resolvable incomplete blocks of ~10 plots within each location,
each carrying the four cells in equal proportion.

At 5 locations × 80 plots, with oat pools of 20 and pea pools of 30, the
generated design gives:

| property | value |
|---|---|
| plots per location | 80 (exactly) |
| cells per location | 20 / 20 / 20 / 20, of which 4 are anchor duplicates |
| plots per oat accession | median 10 (range 8–12) |
| plots per pea accession | median 6 (range 6–10) |
| distinct partners per oat | median 9 (range 6–10) |
| **partner-pool imbalance** | **0 for every accession, both species** |
| replicated combinations | 41 |

That last row is the one the analysis depends on: **every focal accession
meets the two partner pools equally often**, so partner effects are orthogonal
to the focal contrast and drop out of it.

---

## 4. Analysis

Primary, one per species:

```
y_pea ~ location + block(location) + oatPool + peaPool + oatPool:peaPool
        + (1|oat_acc within oatPool) + (1|pea_acc within peaPool)
        + (1|combination) + e
```

Accessions are **random within pool**, so the inference is about the class of
accessions we predicted, not about these particular 40 — which is what makes
the result generalise. The experimental unit for the contrast is therefore the
**accession**, not the plot (§5).

- **Primary tests**: oat-pool effect on pea yield; pea-pool effect on oat
  yield. **One-sided, α = 0.05**, direction pre-registered. Two primary tests;
  the multiplicity position should be stated before the trial goes in.
- **Secondary**: pool × pool interaction; the correlation between predicted
  and realised effects across all accessions (a continuous validation, more
  informative than the binary contrast); specific-combination variance from
  the anchors.

---

## 5. Power, and why the plot budget matters so little

### The result

One-sided α = 0.05, at the design geometry (n = plots-per-location ÷ 4 for
oat, 30 for pea):

| | plots/loc | n per pool | total plots | λ = 1.0 | λ = 0.8 | λ = 0.6 |
|---|---|---|---|---|---|---|
| **oat** | 60 | 15 | 300 | 0.95 | 0.85 | 0.64 |
| **oat** | 80 | 20 | 400 | 0.98 | 0.91 | 0.72 |
| **oat** | 100 | 25 | 500 | 0.98 | 0.92 | 0.74 |
| **pea** | 60 | 30 | 300 | 0.88 | 0.72 | 0.51 |
| **pea** | 80 | 30 | 400 | 0.91 | 0.77 | 0.56 |
| **pea** | 100 | 30 | 500 | 0.93 | 0.80 | 0.59 |

λ is the **attenuation**: how much of the predicted contrast actually
materialises in new environments, combining model calibration with the
stability of associate effects across sites. It matters more than everything
else combined, so it has been measured rather than assumed — **oat 1.28,
pea 0.80** — see §6, which also shows what the fold-to-fold spread costs.

Using the **one-sided** test rather than two-sided is worth 7–11 points and
costs nothing, because the hypothesis is directional and pre-registered:
oat 0.84 → 0.91, pea 0.66 → 0.77 at λ = 0.8.

### Why more plots buy so little

The contrast's variance has two terms and **only one of them contains the plot
count**:

```
SE(Δ)² = 2·σ²_within / n   +   4·σ²_e / P
         └───────────────┘      └────────┘
      variation among the        ordinary
      accessions in a pool       plot noise
```

More plots re-measure the *same n accessions* more precisely. They do not add
new accessions, so the first term is untouched. The experimental unit for a
pool contrast is the accession, and the trial is really a two-sample
comparison with 15–30 per side.

That term dominates, because `σ²_within = PEV + within-pool spread` and the
reliabilities are low — oat `As` 0.31, pea `As` 0.24 — so even the *true*
effects of accessions we picked as "high As" are widely scattered around their
pool mean. At λ = 0.8:

| | SE at P=300 | SE at P=500 | SE at P=∞ | irreducible | power 300 → 500 → ceiling |
|---|---|---|---|---|---|
| oat n=20 | 4.01 | 3.61 | **2.91** | 77% | 0.88 → 0.93 → **0.99** |
| oat n=25 | 3.80 | 3.38 | **2.61** | 74% | 0.85 → 0.92 → **0.99** |
| pea n=20 | 6.53 | 6.03 | **5.19** | 83% | 0.70 → 0.76 → **0.86** |
| pea n=30 | 5.83 | 5.26 | **4.27** | 78% | 0.72 → 0.80 → **0.92** |

**Roughly three quarters of the standard error cannot be bought down with
plots.** With infinitely many plots, pea at n = 20 still caps at 0.86.

What actually moves the needle, in order:

1. **λ** — worth ~40 points across its plausible range.
2. **More accessions per pool** — attacks the dominant term. Pea's *ceiling*
   goes from 0.86 to 0.92 moving n from 20 to 30.
3. **Better estimates** — incoming trials raise reliability, which shrinks PEV,
   which shrinks the dominant term. This is why pools are selected as late as
   possible.
4. **A one-sided test** — 7–11 points, free.
5. **More plots** — 5–10 points across the whole 300 → 500 range.

### Why pool size is not the hard question

Power is nearly flat in *n*: the contrast shrinks as pools grow at almost
exactly the rate the standard error falls. Oat sits at 0.85–0.88 from n = 15
to n = 50. So pool size is set by the **design geometry** — n = plots-per-
location ÷ 4 makes every accession appear exactly twice at every location —
rather than by the power curve. Pea is given larger pools only because its
accession term is the one that dominates.

### Verification

Both an analytic formula and a simulation of the actual design are run
(`code/validate_power.R`). The simulation draws each accession's *true* effect
from its posterior and then selects pools on the *estimates*, which is the
only way to capture selection on noisy predictions:

| | predicted Δ | simulated Δ | analytic SE | empirical SE |
|---|---|---|---|---|
| oat, λ=1.0 | 14.3 | 14.3 | 3.77 | 3.67 |
| pea, λ=1.0 | 16.5 | 16.7 | 5.48 | 5.15 |

The simulation recovers the predicted contrast, and its standard errors come
in slightly *below* the analytic ones — the design's balance removes partner
variance a little better than the formula assumes. Where the two diverge the
analytic number is the conservative one, and it is the one quoted above.

Under a **null** truth (all associate effects set to zero) the rejection rate
returns **0.028 (oat) and 0.022 (pea)** against a nominal α = 0.05. The
experiment can produce a negative result, which is what licenses believing a
positive one. It comes in below nominal because replication is not equal
across accessions — anchors carry twice the plots, and where a pool is larger
than a cell its members appear at only some locations — and the equal-variance
test used in the simulation assumes that heterogeneity away. It errs toward
not rejecting, so the power above is if anything understated; the mixed model
in the real analysis weights by precision and should recover the nominal rate.

---

## 6. λ, measured

λ dominates every other decision, so rather than assume it we estimated it
from the trials already in hand (`code/validate_crossval.R`): hold out one
trial, refit the model on the rest with the *same* fitting code the production
script uses, and ask the held-out trial what the predictions were worth.

**Not every fold asks the same question.** The five 2025 trials share nearly
all their accessions, so holding one out asks *"same lines, new environment"* —
which is what the validation trial will do. B4I_2026_IL is almost disjoint
from the rest (only 13% of its oat accessions appear elsewhere), so holding it
out asks *"new germplasm"* instead. Only the former bears on this design.

### The rehearsal — the validation trial in miniature

The most direct check: build As+/As− pools from the **training trials alone**,
then measure the contrast those pools actually show in the held-out trial.

| held out | oat: predicted → realised | pea: predicted → realised |
|---|---|---|
| B4I_2025_AL | 14.8 → **−1.8** | 15.3 → **0.9** |
| B4I_2025_IA | 14.2 → **−0.9** | 14.6 → **3.4** |
| B4I_2025_IL | 11.7 → **15.4** | 11.7 → **11.9** |
| B4I_2025_ND | 11.7 → **24.0** | 15.6 → **19.1** |
| B4I_2025_NY | 9.9 → **26.1** | 15.2 → **11.1** |

In three of five trials the pools show a large, correctly signed contrast in an
environment the model never saw — **often larger than predicted**. In AL and IA
they show nothing.

AL is explicable: it was a near-total crop failure, mean yield 8 g/m² against
127–383 elsewhere. A trial cannot exhibit a 15 g/m² effect when the whole trial
spans 8. Those folds are excluded from the pooled λ for having too little
variation to show an effect at all, which is recorded rather than quietly done.
IA is not explicable that way and is worth understanding.

### The number

Over the informative "same lines, new environment" folds:

| | λ | 95% CI | across-fold spread | accession-level *r* |
|---|---|---|---|---|
| oat | **1.28** | 0.42 – 2.14 | 0.69 | 0.28 |
| pea | **0.80** | 0.35 – 1.24 | 0.57 | 0.18 |

Two things follow.

**The effects do transfer.** λ at or above 0.8 for both species is at the
optimistic end of the range the design was sized against. λ > 1 for oat means
the realised contrast *exceeds* the predicted one — the BLUPs are over-shrunk,
which is what low reliability plus BGLR's shrinkage would produce.

**But the fold-to-fold spread is larger than we assumed.** The across-fold SD
of the slope is 0.57–0.69 of its mean, against the 0.5 used in the sensitivity.
Read as associate × environment interaction, it costs real power:

| | power ignoring the spread | allowing for it | at the lower CI bound for λ |
|---|---|---|---|
| oat | 1.00 | **0.84** | 0.41 |
| pea | 0.77 | **0.65** | 0.26 |

So: **oat is adequately powered; pea is marginal**, and both depend on λ being
near its point estimate rather than near the bottom of its interval.

### What this rests on

- λ comes from four usable folds, and its own confidence interval is wide.
- The across-fold spread is treated as genuine interaction, but part of it is
  estimation noise in each fold's slope; that split cannot be made with five
  trials, and the three incoming ones will help.
- **The model assumes associate effects are constant in absolute g/m²**, while
  the trials differ several-fold in mean yield. That assumption is what makes
  AL uninformative, and it is questionable generally — note that
  [SIMULATION.md](SIMULATION.md) models environments as *scaling* the whole
  signal, which is the opposite convention. Worth reconciling.

---

## 7. This is re-issued as data arrives

Three further trials are expected. Every input above moves when they land —
the BLUPs shift, reliability rises, PEV falls, fewer accessions fail the
connectivity filter, power rises. **Nothing in the scripts is hard-coded**;
all of it is re-derived at run time, and each run writes a dated vintage under
`output/validation/<date>/`.

Each refresh also reports **which accessions entered and left each pool** since
the previous vintage. That churn is a result in its own right: if a few new
trials reshuffle most of a pool, the effects are not stable enough to be worth
validating — and that is far better known before seed is ordered than after.

Pool membership should be frozen at the last moment compatible with seed
logistics, on the largest dataset available. Everything before that is
provisional.

---

## 8. Open decisions

| # | decision | status |
|---|---|---|
| 1 | Pea pool size: 30–40 (recommended) or p/4 for symmetry with oat | proposed at 30 |
| 2 | Seed availability for 2n accessions per species at 5 locations | unknown — the real constraint on n |
| 3 | Whether to proceed if λ < 0.5 | **decide before seeing λ** |
| 4 | Multiplicity position across the two primary tests | to state pre-trial |
| 5 | Monoculture checks (for land-equivalent ratio) | omitted; would need extra plots |

---

## Reproducing this

```bash
Rscript code/validate_crossval.R         # measures lambda (do this first)
Rscript code/validate_pool_selection.R   # pools + diff vs previous vintage
Rscript code/validate_power.R            # analytic + simulation + null check
Rscript code/validate_design.R           # field book + design checks
```

Outputs land in `output/validation/<date>/`: `crossval_folds.csv`,
`crossval_summary.csv`, `crossval_power.csv`, `pools.csv`, `pool_summary.csv`,
`pool_diff.csv`, `power_grid.csv`, `power_ceiling.csv`,
`power_simulation.csv`, `field_book.csv`, `design_check.txt`, and figures.

The fitting itself is shared: `fit_producer_associate()` in
`code/dge_ige_functions.R` is used both by the production model and by every
cross-validation fold, so the folds refit the same model rather than one that
merely resembles it.
