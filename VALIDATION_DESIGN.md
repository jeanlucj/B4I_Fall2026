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

Implemented with a Lagrange multiplier `θ` on `Pr` and **one selection index
serving both pools**: rank the candidates by `As + θ·Pr`, take the top *n* and
the bottom *n*. Both selections then move along a single axis, so θ buys producer
balance efficiently, and disjointness is structural — the top and bottom of one
ordering cannot collide. θ is chosen by a dense scan that minimises |ΔPr|, not by
bisection: ΔPr(θ) is a step function and bisecting for a sign change can land on
a local crossing (`build_pools()` in `code/validation_functions.R`).

Until 2026-10-03 the two pools were scored by *different* indices — `As + θ·Pr`
and `As − θ·Pr` — which rewarded high Pr in the As⁻ pool and so dragged both
pools upward instead of equalising them. Correcting it raised the contrast by
30% for oat and 42% for pea while *tightening* the producer balance; see
[docs/B4I_followups.md](docs/B4I_followups.md) item 13.

**Candidates are filtered first.** 99 oat and 84 pea accessions appear with
only a single partner, and 8 oat / 20 pea have no marker data at all. Their
effects are almost pure shrinkage, so accessions with **fewer than 3 distinct
partners are excluded** — leaving 265 oat and 254 pea, ample for pools of
15–30.

### What it delivers

<!-- BEGIN GENERATED: pools -->
| species | n per pool | candidates | θ | ΔAs (g/m²) | ΔPr | mean Pr, As+ | mean Pr, As− |
|---|---|---|---|---|---|---|---|
| oat | 20 | 215 | 0.194 | **19.59** | 0.226 | 10.02 | 9.80 |
| pea | 30 | 204 | 0.718 | **25.07** | -0.131 | 7.63 | 7.77 |

*Selection rule: `single-v1`.*
<!-- END GENERATED: pools -->

The producer constraint is close to free: it costs a little associate contrast
and buys two pools whose producer means differ by a fraction of a g/m² while
both sit well above the population mean.

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

**Anchors.** One combination per cell is sown **twice at each location**, and it
is the *same* four combinations everywhere — 8 plots per location, 40 in total.
Both halves of that matter and for different reasons. Within-location duplication
is the only thing that separates plot error from combination × location;
repeating a combination across locations alone confounds the two. Holding the
combinations fixed across locations is what makes combination × location
estimable at all. Getting either half wrong leaves a design that still looks
replicated and answers neither question — so `check_validation_design()` asserts
both, and `tests/test_validation.R` asserts that it does. This also buys a first estimate of specific-combination
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
| plots per oat accession | median 10 (range 8–20) |
| plots per pea accession | median 6 (range 6–20) |
| distinct partners per oat | median 9 (range 2–10) |
| **partner-pool imbalance** | **0 for every accession, both species** |
| replicated combinations | 23: four anchors at 10 plots each, the rest at 2–3 |

The two ranges that look lopsided are the anchors. An anchor combination is sown
twice at each of the five locations, so its four accessions appear in 20 plots
against everyone else's 8–12, and they meet only their four fixed partners rather
than nine. That is the cost of making combination × location estimable, and it is
paid by four accessions per species. It is benign for the contrast — the analysis
averages to one value per accession before comparing pools, so an anchor
contributes a single, slightly more precise mean — and the unevenness is what the
null-rate paragraph below accounts for.

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

<!-- BEGIN GENERATED: vintage -->
*Numbers below are from vintage **2026-10-03**.*

| species | trials | plots | accessions | eligible | reliability of As |
|---|---|---|---|---|---|
| oat | 9 | 3567 | 442 | 430 | 0.35 |
| pea | 9 | 3567 | 423 | 407 | 0.27 |

Trials kept by the QC screen (8): B4I_2025_IA, B4I_2025_IL, B4I_2025_ND, B4I_2025_NY, B4I_2026_IA, B4I_2026_IL, B4I_2026_ND, B4I_2026_NY
Dropped (1): B4I_2025_AL
<!-- END GENERATED: vintage -->

<!-- BEGIN GENERATED: power -->
One-sided α = 0.05, at the configured pool size, with λ and the
pool × location interaction as measured (see §6).

| | plots/loc | n per pool | total plots | λ | interaction | power |
|---|---|---|---|---|---|---|
| **oat** | 60 | 20 | 300 | 0.83 | 0.73 | 0.76 |
| **oat** | 80 | 20 | 400 | 0.83 | 0.73 | 0.78 |
| **oat** | 100 | 20 | 500 | 0.83 | 0.73 | 0.79 |
| **pea** | 60 | 30 | 300 | 0.78 | 0.53 | 0.81 |
| **pea** | 80 | 30 | 400 | 0.78 | 0.53 | 0.84 |
| **pea** | 100 | 30 | 500 | 0.78 | 0.53 | 0.85 |
<!-- END GENERATED: power -->

λ is the **attenuation**: how much of the predicted contrast actually
materialises in new environments, combining model calibration with the
stability of associate effects across sites. It matters more than everything
else combined, so it is measured rather than assumed — see §6.

The grid no longer sweeps λ or the pool × location interaction. Both are
estimated from the trials in hand, and sweeping a measured quantity reports
uncertainty the data have already resolved. What is still swept is what is still
a choice: the plot budget and the pool size.

The test is **one-sided** throughout, because the hypothesis is directional and
pre-registered. The two-sided arm was dropped rather than reported alongside: it
was never this design's test, and having it in the table invited reading the
wrong column.

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

Power is nearly flat in *n*, and the current defaults are within half a point of
the best available. At P = 400 plots:

| n per pool | oat power | pea power |
|---|---|---|
| 10 | 0.750 | 0.804 |
| 15 | 0.771 | 0.834 |
| **20** | 0.779 | **0.840** |
| 25 | 0.782 | 0.836 |
| **30** | **0.784** | 0.839 |
| 40 | 0.780 | 0.826 |
| 50 | 0.770 | 0.807 |

The whole range spans **0.034 for oat and 0.036 for pea**, and the optimum sits
at n = 30 (oat) and n = 20 (pea) — so the configured 20 and 30 are each within
0.005 of the best. Pool size is therefore set by the **design geometry**,
n = plots-per-location ÷ 4, which makes every accession appear exactly twice at
every location, rather than by the power curve.

#### Why it is flat — the two forces cancel

Halving the pool does two opposite things.

**It sharpens selection.** Pools are the extremes of the associate-effect
distribution, so a smaller pool reaches further into the tail. For a normal the
mean of the top *p* fraction is `φ(z_p)/p`, and with 215 oat and 204 pea
candidates after the Pr filter, halving *n* from 40 to 20 buys a contrast about
**1.24×** larger.

**It costs precision.** The standard error of a pool mean goes as 1/√n, so
halving *n* multiplies the SE by **√2 = 1.41**.

Those nearly cancel, and the achieved contrast now behaves as the theory says:
ΔAs falls from 23.2 to 14.9 for oat across n = 10 → 50, almost exactly the 1/√n
the SE follows. The residual curvature is why there is a shallow optimum rather
than a perfectly level line.

> **This section said the opposite for one day, and the reason is worth keeping.**
> Until 2026-10-03 the pool builder used a different selection index for each
> pool, which was inefficient enough that enlarging the pool helped mainly by
> escaping that inefficiency — power then rose from 0.654 to 0.724 for oat and
> 0.414 to 0.695 for pea, and pool size looked like a major lever. Fixing the
> index (§2) raised power at every *n*, and restored the flat curve. The lesson
> is that a power curve shaped by a defect in the selection rule will argue for
> changing the design when what needs changing is the rule.

#### What the Pr constraints cost

**Candidates must have Pr above the median**, which halves the candidate set —
430 eligible oat become 215, 407 pea become 204 — so n = 20 is the top 9% of
candidates rather than the top 4.7% of all eligible.

**The two pools must have nearly equal mean Pr**, and Pr and As are negatively
correlated (−0.51 across all accessions). The As⁺ extreme is systematically
low-Pr and the As⁻ extreme high-Pr, so balancing costs As extremity. With the
corrected index that cost is small: the achieved producer gap is 0.23 g/m² for
oat and −0.13 for pea against a tolerance of 1.0, and both pools sit about 10
(oat) and 8 (pea) g/m² above the population mean Pr.

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

(These two rates were measured before the anchor combinations were fixed across
locations, which made replication slightly less even — four accessions per
species now carry 20 plots rather than 12. The direction of the effect is the
same and the argument is unchanged, but expect the rates to move a little when
`validate_power.R` is next re-run.)

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

<!-- BEGIN GENERATED: lambda -->
| | λ | across-fold spread (interaction_frac) | folds | accession-level *r* |
|---|---|---|---|---|
| oat | **0.83** | 0.73 | 8 | 0.17 |
| pea | **0.78** | 0.53 | 8 | 0.20 |
<!-- END GENERATED: lambda -->

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

The 2026 trials — `B4I_2026_IA`, `B4I_2026_ND`, `B4I_2026_NY` — have now landed,
alongside `B4I_2026_IL`, and every input above moved when they did: the BLUPs
shifted, reliability rose, PEV fell, fewer accessions failed the connectivity
filter. Power did **not** simply rise, because λ is now estimated over nine folds
rather than four and the newer environments reproduce the predicted effects less
well than the 2025 ones did. More trials bought a better-measured λ, not a
larger one. Further trials will move all of it again. **Nothing in the scripts is hard-coded**;
all of it is re-derived at run time, and each run writes a dated vintage under
`output/validation/<date>/`.

Each refresh also reports **which accessions entered and left each pool** since
the previous vintage. That churn is a result in its own right: if a few new
trials reshuffle most of a pool, the effects are not stable enough to be worth
validating — and that is far better known before seed is ordered than after.

Pool membership should be frozen at the last moment compatible with seed
logistics, on the largest dataset available. Everything before that is
provisional.

### How to re-issue it

**One command runs the whole chain:**

```bash
Rscript code/validate_refresh.R --dry-run    # see the plan, run nothing
Rscript code/validate_refresh.R              # the whole chain
```

It drives the existing scripts in dependency order — trial discovery, accession
curation, assembly, the **trial QC screen**, the BGLR fit, then the four
`validate_*` scripts — stops at the first failure, and writes a per-step log
under `output/refresh_logs/`. The thing to read on a failure is the step's own
log, not the driver.

`validate_crossval.R` runs **before** `validate_power.R`, because power now
reads the λ and interaction it measures rather than sweeping assumed values.

#### `--refresh-doc`, and what it does not do

```bash
Rscript code/validate_refresh.R --refresh-doc
```

This rewrites the **numbers** in this document from the vintage the run just
produced — the vintage summary in §5, the power table in §5, the λ table in §6.
Those sit inside `<!-- BEGIN GENERATED: key -->` … `<!-- END GENERATED: key -->`
markers, and `code/validation_report.R` replaces what is between them and
nothing else. `Rscript code/validation_report.R --check` reports whether the
document has drifted from the latest vintage without writing anything.

**It does not touch the prose, and cannot.** Everything outside those markers is
argument, and an argument should not change because a number moved. So a
sentence like "three further trials are expected" goes stale and stays stale
until a person edits it — that exact sentence survived two data vintages before
anyone noticed. If the numbers move enough to change what the document *claims*,
the claim needs rewriting by hand; the flag only keeps the tables honest. Review
the diff before committing it.

Steps 1–3 talk to T3 and need `T3_USERNAME` / `T3_PASSWORD` in `.Renviron`; the
driver checks for them **before** starting a twenty-minute download rather than
after. Everything from `assemble` onward is offline, so if the download has
already happened:

```bash
Rscript code/validate_refresh.R --from assemble
Rscript code/validate_refresh.R --only validate   # just the four validate_* scripts
```

Budget roughly an hour, most of it in trial discovery and the BGLR fit.

**What it deliberately does not touch: the simulation's parameters.**
`SIM_VAR_SHARES`, `SIM_PR_AS_COR` and `SIM_RESID_COR` in `code/sim_config.R` are
frozen at the values the six-trial fit gave on 2026-09-21, so that simulation
results stay comparable across data vintages. New variance components will appear
in `output/BGLR_variance_components.csv` and the power calculation *will* use
them — that is the point of refreshing — but nothing in the chain edits
`sim_config.R`. Refreshing those constants is a separate, deliberate act that
invalidates every cached simulation result.

When it finishes, read in this order:

| file | what it tells you |
|---|---|
| `output/validation/<date>/pool_diff.csv` | **who entered and left each pool.** Read this first: it is the stability check |
| `output/validation/<date>/pool_summary.csv` | the new pools |
| `output/validation/<date>/power_grid.csv` | power at the new variance components |
| `output/validation/<date>/field_book.csv` | the layout to hand to the stations |

---

## 7b. Judging a trial by what it does to power — don't

A tempting test, when the QC screen drops a trial and you want to know whether
it should have: put it back, re-run, and see what happens to the power numbers.
If power improves, the trial had signal; if power drops, it was noise.

**This reasoning runs backwards, and in the dangerous direction.** Power here is
not a measurement of data quality. It is a function of parameters estimated from
the same data, and adding a trial moves several of them at once.

**Noise can raise power.** The predicted contrast ΔAs is the gap between the
mean associate effect of the two pools, and the pools are chosen as the extremes
of the BLUP distribution. A trial that is mostly noise, in a design where most
combinations appear once, inflates the estimated associate variance — noise with
nowhere else to go is absorbed as genotype variance. Wider BLUPs mean a wider
gap between the extremes, a larger ΔAs, and **higher** computed power. The
validation would then fail in the field, because the contrast being predicted is
not there. This is the direction that matters: the test as stated would keep the
bad trial.

**Signal can lower power.** More good data shrinks the BLUPs toward the truth.
The spread narrows, ΔAs falls, and power falls with it — while the prediction
has become more honest, not less. More trials also means λ is estimated over
more environments, which usually lowers it, lowering power again for the same
reason.

**And the comparison is not clean anyway.** Power is computed on the data used
to choose the pools, so it is an in-sample statement throughout. Adding a trial
also moves the across-trial median that `curate_trials.R` compares against, so
other trials' QC verdicts can change in the same run — the two configurations
differ by more than the one trial.

### What to use instead

The honest question is out-of-sample: does the trial's information
**reproduce**? `validate_crossval.R` holds out one trial at a time and answers
it per trial, per species.

**Where to look.** `output/validation/<vintage>/crossval_jackknife.csv`, one row
per trial × species. Find the row whose `held_out` is the trial in question.

**The three columns that decide it**, all about that trial's own fold:

| column | read it as | excludes when |
|---|---|---|
| `lambda` | slope of realised on predicted associate effect in the held-out trial — how much of what this trial's data predict actually shows up | near zero, or negative, with `lambda_p` not significant |
| `r_accession` | correlation of predicted with realised at the accession level in that fold | near zero |
| `narrows_by` | how much the across-fold spread **falls** when this trial's fold is removed | large and positive |

`narrows_by` is the answer to "where do I see whether interaction_frac widens?".
You cannot read it off a column, because `interaction_frac` is a property of the
*set* of folds — `sd(λ)/|mean λ|` — so no single trial has one. It is computed
by removing that trial's fold and recomputing the spread, which is what the
jackknife file now does. **`narrows_by > 0` means this trial is what widens the
interaction.**

The three travel together, and that is not a coincidence: a fold whose λ sits far
from the others both fails to predict itself and is what inflates the spread.

**The rule.** Exclude a trial when its own fold has λ at or near zero, a
non-significant `lambda_p`, `r_accession` near zero, and a large positive
`narrows_by`. Keep it otherwise. Power then follows from whatever that gives,
rather than being the thing consulted.

**Judge it per species.** A trial can fail on one and carry the other, and the
species are separate decisions because the pools are.

### The case in hand, worked

From the 2026-10-02 vintage, nine same-line folds:

| held out | species | λ | *p* | r_accession | mean λ without it | spread narrows by |
|---|---|---|---|---|---|---|
| **B4I_2025_AL** | oat | **−0.04** | 0.26 | 0.00 | 0.73 → 0.83 | **0.139** |
| **B4I_2025_AL** | pea | **0.01** | 0.78 | 0.03 | 0.69 → 0.78 | **0.137** |
| **B4I_2025_IA** | oat | **−0.03** | 0.83 | 0.01 | 0.73 → 0.83 | **0.135** |
| B4I_2025_IA | pea | 0.66 | 0.011 | 0.15 | 0.69 → 0.70 | −0.041 |
| B4I_2025_IL | oat | 1.33 | 3e-7 | 0.28 | — | −0.100 |
| B4I_2025_ND | oat | 1.60 | 5e-6 | 0.29 | — | −0.070 |

**AL fails on both species** — λ indistinguishable from zero, r at zero, and the
two largest `narrows_by` in the table. Dropping it raises mean λ by 0.10 on each
species. That is a clean exclusion, and it agrees with the trial QC screen.

**B4I_2025_IA splits.** On **oat** it looks exactly like AL: λ = −0.03, *p* =
0.83, r = 0.01, and removing it narrows the spread by 0.135. On **pea** it is a
contributing fold: λ = 0.66 at *p* = 0.011, and removing it would *widen* the
spread. So its pea data predict out of sample and its oat data do not.

That is the nuance a power comparison would have hidden, and it is why the trial
is kept in `data/trial_qc_manual.csv` rather than dropped: the cost is a weaker
oat λ, the benefit is a real pea fold. If the oat contrast is the primary
validation, revisit that.

A trial can also be worth keeping while failing both — it contributes plots,
partners and connectivity to accessions that would otherwise be unestimable,
even if its own yields are poor. `curate_trials.R` cannot see that, which is why
its verdict is overridable in `data/trial_qc_manual.csv` and why the override is
recorded rather than silent.

Both are dropped by the QC screen at the current thresholds — AL on high CV, low
mean and low repeatability; IA on high CV and low mean, the latter acquired only
when the 2026 trials raised the across-trial median it is compared against.

## 7c. Collapsed analysis names are one genotype but not one seed lot

Curation pools entries it finds not to segregate under a single **analysis
name** — `<seed>_<pollen>_no_cross` for a full-sib family that did not
segregate, `<line>_self` for a selfed line appearing under several cross names
(CURATION.md). Eight oat names do this, standing for 50 original entries across
301 plots.

Everything downstream is keyed on those names and nothing is lost by it.
`assemble_B4I_phenotypes.R` rewrites `germplasmName` to the analysis name, and
every script that reads a GRM collapses it to match with `collapse_grm()`, which
averages the rows and columns of the members — exactly right for lines that are
genetically identical. All 462 oat and 434 pea accessions are present after the
collapse; no plot is dropped for want of a genotype.

**But a pooled name is not a seed lot.** If one is selected into a validation
pool, somebody has to decide which member entry to sow, source that seed, and
name the line definitively before the trial is planted. The largest,
`IL18-735_IL17-10306_no_cross`, stands for 13 entries.

`output/summary_collapsed_names.csv` is the worklist: every pooled name with its
member entries, the curation reason, and how many plots and partners it carries.
Check the selected pools against it (`output/validation/<vintage>/pools.csv`)
before ordering seed.

A second consequence worth remembering when reading any per-accession figure:
a pooled name inherits every partner its members had, so it looks far better
connected than a real accession. The largest has 59 partners against a median of
5, which is why the partners histogram in
`output/summary_partners.png` clips its axis.

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
