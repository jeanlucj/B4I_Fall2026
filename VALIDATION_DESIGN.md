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

**Candidates are filtered first**, on **how well their own associate effect is
estimated** rather than on a partner count. The threshold is the first quartile
of associate reliability among accessions sitting at exactly three partners —
so it is the bar the old count rule was already accepting, stated as a
reliability instead of a proxy for one — **and** a hard floor of two distinct
partners in the data the model was fitted to.

The floor is not redundant. Reliability measures posterior precision, not
identifiability: an accession grown with a single partner has its producer and
associate effects *perfectly aliased*, and can still score a respectable
reliability by borrowing from well-genotyped relatives. The pools exist to
separate exactly those two effects, so such an accession cannot enter them
whatever its reliability. Nine accessions per species are excluded on that
ground alone.

<!-- BEGIN GENERATED: eligibility -->
| species | rule | reliability threshold | eligible | of |
|---|---|---|---|---|
| oat | reliability | 0.605 | 398 | 462 |
| pea | reliability | 0.614 | 389 | 434 |

Candidates must clear that reliability **and** have met at least two distinct partners in the data the model was fitted to. The partner floor is not redundant: reliability measures posterior precision, not identifiability, and an accession grown with a single partner has its producer and associate effects perfectly aliased while still scoring well by borrowing from genotyped relatives.
<!-- END GENERATED: eligibility -->

Counting partners in the right data matters as much as the rule. Until
2026-10-03 the count was taken over every trial in the phenotype file while the
effects came from a five-trial fit, so 11 of 40 oat and 20 of 60 pea selected
accessions were below the project's own floor in the data the model had
actually seen — and twelve had met exactly one partner. See `SELF_CRITIQUE.md`
finding A1.

### What it delivers

<!-- BEGIN GENERATED: pools -->
| species | n per pool | candidates | θ | ΔAs (g/m²) | ΔPr | mean Pr, As+ | mean Pr, As− |
|---|---|---|---|---|---|---|---|
| oat | 20 | 199 | 0.131 | **21.15** | -0.033 | 13.08 | 13.12 |
| pea | 30 | 195 | 0.414 | **30.17** | -0.127 | 10.93 | 11.05 |

*Selection rule: `single-v2`.*
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

## 4. Analysis, pre-registered

One analysis, written down once, and **every power number in §5 is computed
through it**. Until 2026-10-03 this section specified a mixed model while both
power routes used a two-stage *t*-test, so no reported power described the
model that would actually be fitted. The model now lives as data in
`PREREG_MODEL` (`code/validation_functions.R`) and the block below is generated
from it, so the two cannot drift apart again.

<!-- BEGIN GENERATED: analysis -->
*Analysis `prereg-v1`, fixed 2026-10-03. Generated from `PREREG_MODEL` in
`code/validation_functions.R`, which is what every power number is
computed through.*

One random-effects structure, three instantiations:

```
y ~ location + (1|location:block) + <FIXED> + <AsxE>
      + (1|acc_focal) + (1|acc_partner) + (1|combination) + e
```

| estimand | role | response | `<FIXED>` | `<AsxE>` |
|---|---|---|---|---|
| `slope` | PRIMARY, one per species | the PARTNER's yield | `x_focal + x_partner   (predicted associate effects, centred)` | `(0 + x_focal | location) + (0 + x_partner | location)` |
| `total` | CO-PRIMARY | oat_yield + w_pea * pea_yield   (w_pea = 1, the physical total) | `x_total = GMA_oat + w_pea * GMA_pea, centred` | `(0 + x_total | location)` |
| `pool` | DESCRIPTIVE, not a separate test | the PARTNER's yield | `pool_focal + pool_partner + pool_focal:pool_partner` | `(1|location:pool_focal) + (1|location:pool_partner)` |

**Multiplicity.** HOLM over the THREE primaries: the two per-species slopes and the total-yield slope. Directions are pre-registered, so every test is one-sided at alpha = 0.05. The pool contrasts are descriptive summaries of the first two and are reported unadjusted; counting them as separate tests would penalise reporting two views of the same regression. This settles open decision 4.

**Secondary:**
- pool x pool interaction
- specific-combination variance from the anchors -- estimable almost only from the four anchor combinations, so expect it near the boundary
- the realised associate effect regressed on the focal accession's phenology BLUE from PRIOR trials (not from this trial: in-trial phenology is a MEDIATOR of the associate effect, so adjusting for it would bias the slope toward zero)
- an economically weighted total, w_pea != 1
- accession-level correlation between predicted and realised, for comparison with the cross-validation folds

**Considered and not done:**
- a third pool selected on a phenotypic proxy -- considered and declined; it would cost about 25% more plots
- monoculture checks, so no land-equivalent ratio
- extra plots replicating specific combinations beyond the four anchors -- to be argued at the field-design stage via n_anchor_per_cell
<!-- END GENERATED: analysis -->

**Accessions are random within pool**, so the inference is about the class of
accessions we predicted, not about these particular 40 — which is what makes
the result generalise, and why the experimental unit for the contrast is the
**accession** rather than the plot (§5).

### Why `pool × location` is in it

It was not, and that was the single largest defect in this section. The power
formula says the pool × location term is the *largest* of its three variance
components, and an omitted pool-level term that varies by location does not
vanish — it inflates the test. Simulated with associate effects that are
entirely location-specific and have a population mean of exactly zero, the
model this section used to specify rejects a true null at:

| denominator | oat | pea |
|---|---|---|
| no `pool × location` — the old §4 | **0.170** | **0.120** |
| with it, variance re-estimated from the trial's own locations | 0.062 | 0.046 |
| **with it, variance taken as measured** | **0.027** | **0.034** |
| location-level contrasts, 4 df (narrow inference) | 0.052 | 0.046 |

against a nominal 0.05, over 2,000 replicates each
(`output/validation/<vintage>/power_type1.csv`). That is the concrete reason
for the term, and it is measured rather than argued.

**And the third row is the pre-registered test.** Re-estimating the interaction
variance from the trial's own five locations gives four degrees of freedom to
do it with, and dividing by a variance that noisy makes the statistic
heavy-tailed — which is where the second row's residual inflation comes from.
Taking the cross-validation's eight-fold estimate as *known* instead holds the
size, conservatively. That is already what `contrast_power()`, `slope_power()`
and `total_yield_power()` do, so the power in §5 describes a test that holds
its size rather than one that beats it.

The fourth row is a reference, not a bound. The location-level test holds size
exactly, but it conditions on the particular accessions sown and so answers the
narrow question — "do these 40 lines differ?" — rather than the broad one about
the class they stand for, which is the inference target §4 is built around.

### Why the pool contrast is not a separate test

§4 used to make the pool contrast primary and the predicted-versus-realised
regression secondary, while describing the latter as "more informative" — a
contradiction. The regression is now primary, for three reasons: its slope **is
λ measured in new data**, on the same scale as the `interaction_frac` the design
is conditioned on, so the trial measures the quantity it was sized against; it
does not depend on where the pool boundary fell; and it is directly comparable
with the cross-validation folds in §6.

It is not, however, more *powerful*. Selection has already removed most of the
within-pool spread in the predictor that a regression would exploit: Var(x) over
the 2*n* selected accessions is (Δ/2)² + the within-pool variance, and the
second term is a few percent of the first. The two statistics come out within
half a point of each other (§5), which is why counting both as primaries would
penalise reporting two views of one regression.

---

## 5. Power, and why the plot budget matters so little

### The result

<!-- BEGIN GENERATED: vintage -->
*Numbers below are from vintage **2026-10-03**.*

| species | trials | plots | accessions | eligible | reliability of As |
|---|---|---|---|---|---|
| oat | 8 | 3181 | 462 | 398 | 0.29 |
| pea | 8 | 3181 | 434 | 389 | 0.36 |

Trials kept by the QC screen (8): B4I_2025_IA, B4I_2025_IL, B4I_2025_ND, B4I_2025_NY, B4I_2026_IA, B4I_2026_IL, B4I_2026_ND, B4I_2026_NY
Dropped (1): B4I_2025_AL
<!-- END GENERATED: vintage -->

### The three pre-registered estimands

<!-- BEGIN GENERATED: estimands -->
At **80 plots per location**, one-sided α = 0.05, λ as measured.
Locations are swept because the pool × location term is the largest of the
three variance components and carries only `n_loc − 1` degrees of freedom.

| estimand | role | target | locations | effect | SE | df | power | power at λ's lower 95% |
|---|---|---|---|---|---|---|---|---|
| `pool` | descriptive | oat | 4 | 16.34 | 7.280 | 6.7 | **0.64** | 0.46 |
| `pool` | descriptive | oat | 5 | 16.34 | 6.586 | 9.3 | **0.74** | 0.52 |
| `pool` | descriptive | pea | 4 | 22.68 | 8.622 | 11.5 | **0.80** | 0.62 |
| `pool` | descriptive | pea | 5 | 22.68 | 7.820 | 16.1 | **0.87** | 0.69 |
| `slope` | primary | oat | 4 | 0.77 | 0.341 | 6.5 | **0.65** | 0.47 |
| `slope` | primary | oat | 5 | 0.77 | 0.309 | 9.0 | **0.75** | 0.53 |
| `slope` | primary | pea | 4 | 0.75 | 0.277 | 10.1 | **0.81** | 0.65 |
| `slope` | primary | pea | 5 | 0.75 | 0.250 | 14.0 | **0.89** | 0.73 |
| `total` | co-primary | both | 4 | 0.76 | 0.298 | 6.9 | **0.74** | 0.59 |
| `total` | co-primary | both | 5 | 0.76 | 0.268 | 9.4 | **0.84** | 0.67 |

The pool × location term is 61% of the variance.
<!-- END GENERATED: estimands -->

The three come out within a few points of each other, which is the point: the
continuous slope is adopted for interpretability (§4), not because it is more
powerful, and the total-yield co-primary is comparable rather than better
because the 2 × 2 factorial doubles the predictor's spread and so offsets the
larger plot variance of a sum.

**Locations matter and plots do not.** Losing one site costs about ten points
on every estimand. That is not a quirk — it follows directly from the variance
split below, where the only large term that more plots cannot touch is the one
divided by the number of locations.

### The older per-species grid, for comparison

<!-- BEGIN GENERATED: power -->
One-sided α = 0.05, at the configured pool size, with λ and the
pool × location interaction as measured (see §6).

| | plots/loc | n per pool | total plots | λ | interaction | power |
|---|---|---|---|---|---|---|
| **oat** | 60 | 20 | 300 | 0.77 | 0.73 | 0.72 |
| **oat** | 80 | 20 | 400 | 0.77 | 0.73 | 0.74 |
| **oat** | 100 | 20 | 500 | 0.77 | 0.73 | 0.75 |
| **pea** | 60 | 30 | 300 | 0.75 | 0.54 | 0.84 |
| **pea** | 80 | 30 | 400 | 0.75 | 0.54 | 0.87 |
| **pea** | 100 | 30 | 500 | 0.75 | 0.54 | 0.89 |
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

### Why more plots buy so little, and what does

The contrast's variance has **three** terms, and only the last contains the
plot count:

```
SE(Δ)² = (σ²_within⁺ + σ²_within⁻)/n  +  (int_sd)²/n_loc  +  4·σ²_e/P
         └──────────────────────────┘     └────────────┘     └───────┘
            variation among the            associate ×        ordinary
            accessions in a pool            location         plot noise
```

This section said "two terms" until 2026-10-03, and the missing middle one is
the largest of the three — **65% of the variance for oat and 50% for pea** at
the measured interaction, at 80 plots across 5 locations. That changes the
conclusion's *direction*, not its force:

**More plots re-measure the same *n* accessions.** They do not add new
accessions, so the first term is untouched; the experimental unit for a pool
contrast is the accession, and the trial is really a two-sample comparison with
15–30 per side. `σ²_within` is the per-accession prediction error variance plus
the within-pool spread of the predictions, and it is now computed **per pool**,
because the two extremes of the distribution are not equally well estimated —
oat's As⁺ pool carries 56.9 against As⁻'s 40.4.

**And more plots do nothing at all to the middle term**, which is divided by the
number of *locations*. So **88% of the oat standard error and 80% of the pea
cannot be bought with plots**, and most of what cannot be bought is
associate × location. **The one remaining lever is more sites**, which is why
`n_loc` is the sweep in the table above and the plot budget is simply fixed at
80 per location. It is also why the middle term drives the *degrees of freedom*
as well as the variance: it carries only `n_loc − 1`, which is what takes the
effective df from 38 down to 9.3 for oat.

A note on `σ²_within`, because the number moved a long way. It used to be
backed out of the identity Var(true) = Var(BLUP) + E[PEV] by solving for the
second term. That identity does not hold in this fit: for oat,
Var(BLUP) = 46.9 and the measured E[PEV] = 49.0 sum to 95.9 against a fitted
component of 160.9, with mean(diag(G)) = 1.00 so the GRM scale is not the
explanation. The backed-out PEV was therefore **2.3× the posterior's own value**
(the same ratio appears independently in pea), which inflated this term and
*understated* power. PEV is now measured per accession from the streamed
coefficient draws; both values and their ratio are reported in `vintage.csv`.

At λ = 0.8:

The current numbers are in `power_ceiling.csv`; the shape of the thing is what
matters here. **Roughly 80–90% of the standard error cannot be bought down with
plots**, and most of what cannot be bought is the pool × location term.

What actually moves the needle, in order:

1. **λ** — worth ~40 points across its plausible range.
2. **More locations.** At the same per-site effort, going from four sites to
   five takes oat from 0.643 to 0.743 — ten points. That is more than the
   entire 60-to-100 plots-per-location range buys, because locations are the
   only lever on the dominant term, and because that term carries just
   `n_loc − 1` degrees of freedom.
3. **Better estimates** — more trials raise reliability, which shrinks PEV,
   which shrinks the accession term. This is why pools are frozen as late as
   seed logistics allow.
4. **A one-sided test** — 7–11 points, free.
5. **More plots per location** — about 2–3 points across the whole 60-to-100
   range. The weakest lever by a wide margin, which is why the budget is now
   simply fixed at 80.

**A note on reading the budget.** Until 2026-10-03 this grid crossed a *total*
plot count with a location count, which produces cells that are not the same
design: 400 plots is 100 per site at four locations and 80 at five, and those
give 0.652 and 0.743. The sweep is now over **plots per location**, with the
total derived, so the two location arms are compared at equal per-site effort.
`power_curves.png` saw-toothed for exactly this reason — it was joining points
from two different designs.

### Why pool size is not the hard question

Power is nearly flat in *n*. At 80 plots per location and five locations:

| n per pool | oat power | pea power |
|---|---|---|
| 10 | 0.742 | 0.882 |
| 15 | **0.744** | **0.887** |
| **20** | 0.743 | 0.880 |
| 25 | 0.737 | 0.879 |
| **30** | 0.733 | 0.871 |
| 40 | 0.719 | 0.851 |
| 50 | 0.701 | 0.823 |

The whole range spans **0.041 for oat and 0.064 for pea**, the optimum sits at
n = 15 for both, and the configured 20 and 30 are within 0.001 and 0.016 of it.
Pool size is therefore set by the **design geometry**,
n = plots-per-location ÷ 4, which makes every accession appear exactly twice at
every location, rather than by the power curve.

#### Why it is flat — the two forces cancel

Halving the pool does two opposite things.

**It sharpens selection.** Pools are the extremes of the associate-effect
distribution, so a smaller pool reaches further into the tail. For a normal the
mean of the top *p* fraction is `φ(z_p)/p`, and with 199 oat and 195 pea
candidates after the Pr filter, halving *n* from 40 to 20 buys a contrast about
**1.24×** larger.

**It costs precision.** The standard error of a pool mean goes as 1/√n, so
halving *n* multiplies the SE by **√2 = 1.41**.

Those nearly cancel, and the achieved contrast behaves as the theory says: ΔAs
falls from 24.5 to 15.7 for oat across n = 10 → 50, and 39.7 to 24.2 for pea —
close to the 1/√n the SE follows. The residual curvature is why there is a
shallow optimum rather than a perfectly level line, and why that optimum has
moved to n = 15: the cancellation is not exact, and with the pool × location
term added the SE now falls off a little more slowly in *n* than the contrast
does.

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
398 eligible oat become 199, 389 pea become 195 — so n = 20 is the top 10% of
candidates rather than the top 5% of all eligible.

**The two pools must have nearly equal mean Pr**, and Pr and As are negatively
correlated (−0.47 for oat and −0.56 for pea across all accessions). The As⁺
extreme is systematically low-Pr and the As⁻ extreme high-Pr, so balancing
costs As extremity. With the corrected index that cost is small: the achieved
producer gap is **−0.03 g/m² for oat and −0.13 for pea** against a tolerance of
1.0, and both pools sit about 10 (oat) and 11 (pea) g/m² above the population
mean Pr.

### Verification

Two routes are run (`code/validate_power.R`): the analytic formula above, and a
simulation of the **actual** generated design.

**What the simulation checks, and what it cannot.** It checks the standard-error
formula and the design's balance — that partner effects really do orthogonalise,
that the anchors' unequal replication behaves, and that a denominator omitting
`pool × location` over-rejects when that term is real (§4). It is **not**
evidence that the predicted contrast materialises in the field, and this section
used to claim otherwise: the truth is drawn as `λ · (BLUP + error)` with the
error centred on zero, so `E[true | BLUP] = λ · BLUP` *by construction* and the
simulated contrast is pinned to `λ · ΔAs` whatever the BLUPs are actually worth.
A table showing "predicted 14.3, simulated 14.3" was reporting an identity, not
an agreement. The only evidence about whether the contrast transfers is λ, and
λ comes from §6.

Two things the simulation does now that it did not:

- **The error is drawn per accession.** Each accession's own posterior SD,
  measured from the streamed coefficient draws, rather than one global value
  that gave a line seen with twenty partners the same prediction error as one
  seen with two. Across accessions the PEV spans a 6.9-fold range for oat.
- **Associate × location is generated**, and *not* centred to sum zero across
  the sites sown. Centring looks careful and silently destroys the test: the
  estimate averages over locations, so deviations summing to zero leave the
  contrast untouched and there is no inflation left to find. The locations in a
  trial are a *sample* — their deviations have mean zero in the population, not
  across five particular sites — and that realised, non-zero mean is exactly
  what the `(int_sd)²/n_loc` term describes.

The type I rates that follow are in §4. The headline: the model this document
specified until 2026-10-03 rejects a true null at 0.170 for oat, and the
pre-registered one at 0.027.

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

Eight folds, one per trial the QC screen keeps. The crop-failure trial at AL no
longer becomes a fold at all, because the screen now runs upstream of the
cross-validation (§7).

| held out | oat: predicted → realised | pea: predicted → realised |
|---|---|---|
| B4I_2025_IA | 12.4 → **1.5** | 17.3 → **0.7** |
| B4I_2025_IL | 13.1 → **12.8** | 17.1 → **21.7** |
| B4I_2025_ND | 12.0 → **16.9** | 18.0 → **30.2** |
| B4I_2025_NY | 11.2 → **20.7** | 18.2 → **12.9** |
| B4I_2026_IA | 13.4 → **18.0** | 23.7 → **48.7** |
| B4I_2026_IL | 16.4 → **4.2** | 33.5 → **23.4** |
| B4I_2026_ND | 17.3 → **−0.9** | 33.7 → **16.0** |
| B4I_2026_NY | 17.0 → **4.3** | 27.3 → **22.2** |

Read the two columns differently. **Pea delivers in seven of eight** folds, and
in three of them by more than predicted. **Oat is split**: four folds deliver
and four come in near zero, and the four that fail are not the four you would
guess — they are IA 2025 and the three 2026 sites other than IA. That is the
spread §6 quantifies as `interaction_frac`, seen one fold at a time, and it is
why oat's λ interval is the wider of the two.

One thing the table no longer shows, and it matters: the predicted contrast
*grows* down the oat column, from about 12 to about 17, while the realised one
does not. §7b sets out the mechanism — more data in a design where most
combinations appear once can widen the BLUPs faster than it sharpens them, which
raises the predicted contrast and lowers λ. The pooled λ absorbs both, which is
the argument for freezing pools as late as seed logistics allow.

### The number

Over the informative "same lines, new environment" folds:

<!-- BEGIN GENERATED: lambda -->
| | λ | across-fold spread (interaction_frac) | folds | accession-level *r* |
|---|---|---|---|---|
| oat | **0.77** | 0.73 | 8 | 0.18 |
| oat | **0.02** | 0.81 | 8 | 0.18 |
| pea | **0.75** | 0.54 | 8 | 0.21 |
| pea | **0.01** | 0.37 | 8 | 0.21 |
<!-- END GENERATED: lambda -->

### Absolute or proportional? The question §6 used to leave open

This section used to end by noting that the model assumes associate effects are
constant in absolute g/m² while the trials differ several-fold in spread, and
that `SIMULATION.md` uses the opposite convention — "worth reconciling". It is
now reconciled, by estimating λ both ways. Because `interaction_frac` is
`sd(λ)/|mean(λ)|` it is scale-free, so the two scales are directly comparable.

<!-- BEGIN GENERATED: lambda_scales -->
λ is estimated twice: in absolute g/m², which is what the trial is sized
in, and with the held-out trial's response divided by its own SD.
`interaction_frac` is `sd(λ)/|mean(λ)|`, so it is scale-free and the two
rows per species are directly comparable — which is what settles whether
the across-fold spread is interaction or just the trials differing in
spread. `interaction_frac_corrected` additionally removes fold-level
estimation noise, using the accession-clustered standard errors.

| species | scale | folds | λ | λ in g/m² | SE of λ | interaction_frac | corrected |
|---|---|---|---|---|---|---|---|
| oat | raw | 8 | 0.773 | 0.77 | 0.199 | 0.728 | 0.624 |
| oat | z | 8 | 0.022 | 0.77 | 0.006 | 0.811 | 0.731 |
| pea | raw | 8 | 0.752 | 0.75 | 0.144 | 0.543 | 0.417 |
| pea | z | 8 | 0.015 | 0.75 | 0.002 | 0.371 | 0.130 |

By year, as a **diagnostic** — the global estimate over all folds is what
the power table uses, and nothing in the chain re-sizes on a subset of folds.

| species | year | folds | λ | SD |
|---|---|---|---|---|
| oat | 2025 | 4 | 0.940 | 0.685 |
| oat | 2026 | 4 | 0.605 | 0.440 |
| pea | 2025 | 4 | 0.709 | 0.361 |
| pea | 2026 | 4 | 0.794 | 0.503 |
<!-- END GENERATED: lambda_scales -->

The answer differs by species, which is why it was worth asking. For **oat**,
standardising makes the spread slightly *worse* (0.73 → 0.81), so oat's
across-fold variation is genuine associate × environment interaction and not an
artefact of trials differing in spread. For **pea**, standardising removes about
a third of it (0.54 → 0.37), so a real part of pea's apparent interaction is
scale heterogeneity. Removing fold-level estimation noise as well takes the raw
figures to 0.62 and 0.42; neither floors at zero, so the interaction is real in
both species, just smaller than the uncorrected numbers suggest.

The trial is sized on the **raw** scale, because that is the scale it is sown
and harvested in.

Two things follow.

**The effects do transfer.** λ near 0.77 for both species, over eight folds,
with the oat 95% interval 0.38–1.16 and pea's 0.47–1.03. Individual folds do
exceed 1, which would mean the realised contrast beats the predicted one and
the BLUPs are over-shrunk; the *mean* does not, and it is the mean the trial is
sized on.

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

- λ comes from **eight** folds, one per trial the QC screen keeps, and its own
  confidence interval is still wide — oat 0.38 to 1.16.
- The across-fold spread is **no longer** treated as entirely genuine
  interaction. Each fold's slope now carries an accession-clustered standard
  error, and subtracting the mean squared error leaves the part that is real:
  0.73 → 0.62 for oat, 0.54 → 0.42 for pea. Neither floors at zero.
- The **absolute-versus-proportional** question is answered above rather than
  deferred.
- λ is no longer conditional on a fold-inclusion rule chosen by looking at λ.
  The trial QC screen now runs upstream of the folds, so the crop-failure trial
  never becomes a fold at all, and the "all folds" and "informative folds" sets
  are **identical** — which is what removes the circularity `SELF_CRITIQUE.md`
  finding G describes.

---

## 7. This is re-issued as data arrives

The 2026 trials — `B4I_2026_IA`, `B4I_2026_ND`, `B4I_2026_NY` — landed alongside
`B4I_2026_IL`, and for two vintages **they did not reach the fit at all**.
`BGLR_multi_trait_model.R` carried a hard-coded trial whitelist that was never
updated, so the production fit ran on five trials and 1,985 plots while every
table here reported the phenotype file's nine and 3,567. The reliability rise
that §7 used to attribute to the new trials was caused by *dropping* the
crop-failure trial at AL, not by adding anything. `SELF_CRITIQUE.md` finding A
has the full account.

That whitelist is gone. The trial set is now derived from `output/trial_qc.csv`,
the fit writes `BGLR_fit_trials.csv` and `BGLR_fit_provenance.csv` recording
what it actually saw, and `validation_inputs()` **refuses to run** if those
disagree with the QC-filtered plot table. Both now go through one shared filter
chain, `b4i_fit_frame()`, so they agree by construction rather than by
discipline.

The current fit: **8 trials, 3,181 plots, 462 oat and 434 pea accessions**.
Further trials will move all of it again; each run writes a dated vintage under
`output/validation/<date>/`.

Each refresh also reports **which accessions entered and left each pool** since
the previous vintage. It is reported because somebody ordering seed will want to
know what moved — **not as a criterion**. This section used to say that a pool
retaining under 70% of its members was a sign the effects were not stable enough
to be worth validating. That is the wrong standard for this project, and the
guidance has been removed.

Two years have gone into estimation trials. If the honest conclusion were
"estimate more", that is itself the result — a statement that the approach is
not cost-effective — not a reason to postpone. The year-3 trial is a validation
trial regardless. And while the analysis method is still being settled, churn
between vintages mostly measures changes to the *method*: the `pool_index` stamp
in `pool_summary.csv` records which rule built a vintage, so a diff across a
rule change can be recognised as incomparable rather than read as instability.

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
| `lambda` | slope of realised on predicted associate effect in the held-out trial — how much of what this trial's data predict actually shows up | near zero, or negative, with `lambda_p` not significant. **Use the clustered `lambda_p`, not `lambda_p_plot`** — see below |
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
non-significant clustered `lambda_p`, `r_accession` near zero, and a large positive
`narrows_by`. Keep it otherwise. Power then follows from whatever that gives,
rather than being the thing consulted.

**Judge it per species.** A trial can fail on one and carry the other, and the
species are separate decisions because the pools are.

**Which *p*-value.** The slope is fitted at the plot level, where plots are not
independent: accessions recur across plots and plots sit in blocks. The
plot-level *p* is therefore too small, and this section leans on it to *exclude
a trial* — a decision error, not a cosmetic one. The fold table now carries both:
`lambda_p` is clustered by focal accession and is the one to use;
`lambda_p_plot` is retained only so the worked example below stays readable.
The inflation is modest in practice — the clustered standard error is about 6%
larger — but the rule should not depend on that having been checked.

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
| 3 | Whether to proceed if λ < 0.5 | **moot**: λ is 0.77 (oat) and 0.75 (pea) over eight folds. The lower 95% bound, 0.38, puts oat at 0.53 — that is the scenario to decide about |
| 4 | Multiplicity position | **settled**: Holm over three primaries, §4 |
| 5 | Monoculture checks (for land-equivalent ratio) | omitted; would need extra plots. Note total yield itself needs **none** — both yields are already on every plot, and it is a co-primary (§4) |
| 6 | Whether the primary test estimates the interaction variance or takes it as measured | **recommend: as measured.** Re-estimating it from five locations gives four degrees of freedom and a test that rejects a true null at 0.062 (oat); taking the cross-validation's eight-fold estimate as known holds size at 0.032 and 0.044. The reported power already assumes the latter |
| 7 | Extra plots replicating specific combinations, within and across locations | to argue at the field-design stage, via `n_anchor_per_cell`. Power here assumes the current one anchor per cell |

---

## Reproducing this

```bash
Rscript code/validate_crossval.R         # measures lambda (do this first)
Rscript code/validate_pool_selection.R   # pools + diff vs previous vintage
Rscript code/validate_power.R            # analytic + simulation + null check
Rscript code/validate_design.R           # field book + design checks
```

Outputs land in `output/validation/<date>/`: `crossval_folds.csv`,
`crossval_summary.csv`, `crossval_lambda_scales.csv`,
`crossval_lambda_year.csv`, `crossval_power.csv`, `pools.csv`,
`pool_summary.csv`, `pool_diff.csv`, `estimand_power.csv`, `power_grid.csv`,
`power_ceiling.csv`, `power_simulation.csv`, `power_type1.csv`,
`field_book.csv`, `design_check.txt`, and figures. Beside them in `output/`:
`BGLR_fit_trials.csv` and `BGLR_fit_provenance.csv`, which record what the fit
saw, and `BGLR_{oat,pea}_pev.csv`, the per-accession prediction error variances.

**The fitting itself is shared** — and as of 2026-10-03 that is actually true.
`fit_producer_associate()` in `code/dge_ige_functions.R` is called by the
production model and by every cross-validation fold, so the folds refit the same
model rather than one that merely resembles it. Until then this sentence was
wrong: `BGLR_multi_trait_model.R` built its own design matrices and called
`BGLR::Multitrait()` directly. The two constructions were verified identical
(`max|diff| = 0` on every term) before the switch, and
`BGLR_fit_provenance.csv` records `fit_fn` so the claim stays checkable rather
than asserted.
