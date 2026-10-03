# SELF-CRITIQUE — the pathway from data to validation design

A step back from the code. Not "do the functions do what the comments say" —
the test suite asks that, and mostly gets a yes. This asks the question the
test suite cannot: **is the logic of the pathway right, and is each task on it
done the right way?**

Written 2026-10-03, against the `2026-10-03` vintage, before seed is ordered.

For the pathway itself see [README.md](README.md) §3 and
[DESIGN.md](DESIGN.md); for the reasoning see [BACKGROUND.md](BACKGROUND.md);
for the design being critiqued see
[VALIDATION_DESIGN.md](VALIDATION_DESIGN.md).

---

## 0. How to read this

**A note on trial counts.** This document was written against a vintage in
which the assembled phenotype table held **9** trials and the QC screen kept
**8**, and it uses 9 in places where it should say 8. The analysis set is 8 —
`B4I_2025_{IA,IL,ND,NY}` and `B4I_2026_{IA,IL,ND,NY}`, with `B4I_2025_AL`
dropped by `data/trial_qc_manual.csv`. The number 9 is only ever the *raw*
table, and the gap between the two is finding A itself. Where a count below
reads 9, check whether the raw or the filtered set is meant.

Findings are ordered by **what they would change**, not by how interesting
they are. Three categories:

| | meaning |
|---|---|
| **BLOCKING** | fix before seed is ordered; the current numbers are wrong or mean something other than what they say |
| **MATERIAL** | the design still stands but a stated conclusion needs weakening or a cheap addition is being left on the table |
| **NOTED** | real, worth recording, does not change a decision now |

Every quantitative claim below was recomputed from `output/` during this
review. Where I could not check something I say so (§8).

**The headline.** The pathway's *logic* is sound and in places better than
standard practice — §7 says what I think should not be touched. But the
**numbers currently attached to it describe an experiment other than the one
that was run**: the BLUPs the pools are built from come from a fit on five
trials and 1,985 plots, while every table in VALIDATION_DESIGN.md advertises
nine trials and 3,567 plots. The consequences reach the pool membership
itself. That is finding A, and nothing else matters until it is fixed.

---

## 1. The pathway, stated plainly

Stripped of implementation, the chain is:

```
T3/Oat (BrAPI)
  │
  ├─1─ trial discovery          find_trials_with_B4I_accessions.R
  ├─2─ accession curation       curate_{oat,pea}_accessions.R
  ├─3─ phenotype assembly       assemble_B4I_phenotypes.R
  ├─4─ trial QC                 curate_trials.R
  ├─5─ the model                BGLR_multi_trait_model.R  →  Pr, As BLUPs
  │
  ├─6─ does it transfer?        validate_crossval.R       →  λ, interaction
  ├─7─ who to sow               validate_pool_selection.R →  As+ / As− pools
  ├─8─ is it worth doing?       validate_power.R          →  power
  └─9─ how to sow it            validate_design.R         →  field book
```

**Does that order make sense? Yes, and one thing about it is unusually
good:** step 6 comes before step 8. Measuring transferability from data in
hand, and then conditioning the power claim on the measurement, is the right
way round and is not how most validation proposals are written. Most assume a
λ and sweep it.

**Is anything missing from the chain? Two things, and both are real:**

- **Nothing checks that step 5's fit used the data step 3 assembled.** That
  gap is finding A. It is a structural hole in the pathway, not a bug in a
  function: no script compares the fit's trial set against the phenotype
  table's, and `validation_vintage()` reports the *phenotype table's* trial
  count as if it described the fit.
- **Nothing connects the validation back to the breeding objective.** The
  chain ends at "is the As contrast real", which is a model question. The
  breeding question is "does selecting on this index produce a better
  intercrop". Those are not the same, and the second is nearly free to add
  (finding J).

---

## 2. BLOCKING findings

### A. The fit never saw three of the nine trials, and the pools inherit it

`code/BGLR_multi_trait_model.R:72` carries a hard-coded trial whitelist:

```r
trials <- c("B4I_2025_AL", "B4I_2025_IA", "B4I_2025_IL",
            "B4I_2025_ND", "B4I_2025_NY", "B4I_2026_IL")
```

`B4I_2026_IA`, `B4I_2026_ND` and `B4I_2026_NY` are **not in it**. The filter
at line 181 (`studyName %in% trials`) drops them, and `apply_trial_qc()` at
line 127 then drops `B4I_2025_AL`. So the production fit ran on **five
trials and 1,985 plots**.

The evidence, recomputed:

| | oat accessions | pea accessions | plots |
|---|---|---|---|
| assembled phenotype table (9 trials, **before** QC) | 462 | 434 | 3,567 |
| the hard-coded 6 | 442 | 423 | 2,371 |
| **the 6 minus AL = what was fitted (5)** | **442** | **423** | **1,985** |
| `BGLR_{oat,pea}_effects_all_seeds.csv` | **442** | **423** | — |

The effects files carry 442 and 423 — the five-trial figure. The fit was run
at 17:40 on 2026-10-02, four minutes after the nine-trial phenotype table was
written at 17:36, and still produced five trials' worth of accessions. The
whitelist silently absorbed the new data.

**Why nothing caught it.** `validation_inputs()` computes
`n_trials = n_distinct(pheno$studyName)` from the phenotype table
(`code/validation_functions.R:170`) and never from the fit. So
`vintage.csv` says `n_trials = 9, n_plots = 3567` while describing BLUPs that
saw five trials and 1,985 plots, and §5 of VALIDATION_DESIGN.md prints that
row verbatim. The one number that would have exposed the problem is
structurally incapable of doing so. `validate_refresh.R` runs the fit step,
gets exit status 0, and reports success.

**Four consequences, in order of severity.**

**A1 — the connectivity filter admits the accessions it exists to exclude.**
`min_partners = 3` is justified in VALIDATION_DESIGN.md §2 because below it
"their effects are almost pure shrinkage, so there is nothing in it to
validate". But partner counts are taken from the nine-trial table while the
effects come from five trials. Counting partners in the data the fit actually
used:

| | eligible by the 9-trial count | eligible in the fitted data |
|---|---|---|
| oat | 430 | **248** |
| pea | 407 | **244** |

And in the selected pools themselves:

| partners in the fitted data | 1 | 2 | 3 | 4 | ≥5 |
|---|---|---|---|---|---|
| oat pool members (of 40) | **4** | 7 | 5 | 7 | 17 |
| pea pool members (of 60) | **8** | 12 | 6 | 5 | 29 |

**11 of 40 oat and 20 of 60 pea pool members fall below the project's own
eligibility threshold** in the data the model saw. Twelve of them met exactly
one partner — and BACKGROUND.md states the consequence itself: "a genotype
grown with one partner has the two aliased." The design's entire premise is
pools *matched on Pr and contrasted on As*. Twelve selected accessions have no
Pr/As separation at all. Their position in the pool is an artefact of which
single partner they happened to draw.

**A2 — λ and ΔAs come from different training sets.** `validate_crossval.R`
reads `B4I_intercrop_pheno.rds` directly, with no `apply_trial_qc()` and no
whitelist, so it runs **nine folds on all nine trials**, each fold training on
eight. λ = 0.83 therefore describes eight-trial BLUPs; ΔAs = 19.59 and
PEV = 81.1 describe five-trial BLUPs. `contrast_power()` multiplies them
together. More training data means less shrinkage and less attenuation, so a λ
measured on eight trials is the optimistic one for BLUPs built from five.
Every fold's training set also includes AL, which production discards — so the
claim that "the folds refit the same model rather than one that merely
resembles it" is true of the *model* and false of the *data policy*.

**A3 — §7 attributes a change to the wrong cause.** VALIDATION_DESIGN.md §7
states: "The 2026 trials … have now landed … and every input above moved when
they did: the BLUPs shifted, reliability rose, PEV fell." The BLUPs did not
shift. The vintage history shows what actually happened:

| vintage | ΔAs oat | reliability_As oat | what changed |
|---|---|---|---|
| 2026-09-24 | 14.3 | 0.307 | 6-trial fit incl. AL |
| 2026-09-27 | 14.3 | 0.307 | nothing |
| 2026-10-02 | 15.2 | 0.349 | **AL dropped by the QC screen** |
| 2026-10-03 | 19.6 | 0.349 | **the pool index was fixed** |

Reliability rose because a crop-failure trial was *removed*, not because three
trials were *added*. The sentence "**Nothing in the scripts is hard-coded**;
all of it is re-derived at run time" is the opposite of the truth about the one
script that matters most.

**A4 — the stability check has not been exercised.** §7 makes pool churn the
go/no-go: "if a few new trials reshuffle most of a pool, the effects are not
stable enough to be worth validating." Under a fixed selection rule, the only
data change so far has been dropping AL, and it retained **27 of 40 oat
(68%)** and **43 of 60 pea (72%)** — already at the 70% line, for removing one
trial. The three new trials have never been put through this test. My
expectation is that they will breach it.

**Fix.** Replace the whitelist with the QC verdict as the single source of
truth (`trials` should be derived from `output/trial_qc.csv`, not typed), and
make `validation_vintage()` read the trial set *from the fit* — write it out
alongside the effects and assert it matches. Then re-run the chain and read
`pool_diff.csv` before anything else.

---

### B. λ has a year structure, and the mean sits in a gap between two regimes

This is the finding I would most want someone to argue with me about, because
it is a judgement about what λ means rather than an error.

The per-fold oat slopes, from `crossval_folds.csv`:

| held out | λ | *p* | r_accession | predicted → realised |
|---|---|---|---|---|
| B4I_2025_AL | −0.04 | 0.26 | 0.00 | 14.1 → −1.8 |
| B4I_2025_IA | −0.03 | 0.83 | 0.00 | 13.3 → −4.3 |
| B4I_2025_IL | **1.33** | 4e-7 | 0.28 | 11.8 → 11.6 |
| B4I_2025_ND | **1.60** | 5e-6 | 0.29 | 11.2 → 17.0 |
| B4I_2025_NY | **1.22** | 6e-5 | 0.34 | 10.4 → 21.6 |
| B4I_2026_IA | 1.32 | 0.002 | 0.17 | 13.4 → 20.6 |
| B4I_2026_IL | 0.53 | 0.014 | 0.10 | 17.1 → 9.0 |
| B4I_2026_ND | 0.26 | 0.22 | 0.12 | 18.2 → −4.2 |
| B4I_2026_NY | 0.40 | 0.20 | 0.05 | 17.2 → 1.5 |

The 2025 folds that have any signal average **λ = 1.38**. The 2026 folds
average **λ = 0.63**, and three of the four are at 0.26–0.53 with two not
significant. The pooled 0.83 is an average across two regimes and describes
neither.

Two readings, and they have different consequences:

1. **Training-set composition.** Four of the five fitted trials are 2025, so a
   2025 fold is predicted largely by its own year's siblings and a 2026 fold is
   not. Then λ ≈ 0.63 is the honest number for a *new* year — which is what the
   validation trial is.
2. **The predicted contrast is inflating.** The `predicted` column rises from
   10–14 in the 2025 folds to 17–34 in the 2026 folds while `realised` does
   not. §7b of VALIDATION_DESIGN.md predicts exactly this mechanism — "noise
   with nowhere else to go is absorbed as genotype variance … Wider BLUPs mean
   a wider gap between the extremes, a larger ΔAs, and **higher** computed
   power. The validation would then fail in the field." The project wrote the
   warning and is now showing the symptom.

**Resolved 2026-10-03, and the second reading was wrong.** "Forward-looking" was
muddled reasoning on my part: holding out 2026 while training mostly on 2025 is
not "predicting the future", and predicting 2025 from 2026 is just as valid a
question. The split was reading 1 — training-set composition. After the refit on
8 trials, with 4 from each year and every fold training on 7, it has collapsed:
oat 0.94 for 2025 folds against 0.60 for 2026, a 0.34 gap against a within-year
SD of 0.56, and pea shows no split at all (0.71 against 0.79). **λ is estimated
globally over all 8 folds — 0.77 for oat, 0.75 for pea — and by-year is reported
as a homogeneity diagnostic, not as an estimate to re-size on.**

The table below is kept because it shows how much the λ choice is worth, but its
"2026 folds only" row should be read as a defunct hypothesis rather than a
recommendation. Power at the design geometry (n = 20/30, P = 400, one-sided
α = 0.05):

| λ used | source | oat power | pea power |
|---|---|---|---|
| 0.83 / 0.78 | informative folds, as quoted in §5–6 | **0.78** | **0.84** |
| 0.73 / 0.69 | *all* same-line folds — in `crossval_summary.csv`, never quoted | 0.66 | 0.72 |
| 0.63 | oat, 2026 folds only | **0.67** | — |
| 0.41 / 0.49 | lower 95% bound, as §6 reports | 0.51 | 0.60 |

Note the second row: the unconditional estimate is computed, written to disk,
and never shown in the document. The one that is shown is conditional on a
fold-inclusion rule (finding I).

**Fix, as implemented.** Report λ by year as a diagnostic and size on the global
estimate over all 8 folds. Do not re-size on a subset of folds. §8 item 3's
question is now answerable: λ is 0.77 and 0.75, and the scenario worth deciding
about is the lower 95% bound of 0.38, which puts oat at 0.53.

---

## 3. MATERIAL findings

### C. The analysis model, the power formula and the simulation are three different experiments

Three places describe how the trial will be analysed, and no two agree:

| | model | denominator for the pool test |
|---|---|---|
| VALIDATION_DESIGN.md §4 | mixed model: `location + block + oatPool + peaPool + oatPool:peaPool + (1\|oat_acc) + (1\|pea_acc) + (1\|combination)` | accession within pool |
| `contrast_se()` | `2σ²within/n + 4σ²e/P + (int·Δ)²/n_loc` | accession + plot + pool×location |
| `simulate_validation()` | centre on location mean, average to accession, equal-variance *t* | accession + plot |

**No power number describes the model that will actually be fitted.** Two
specific problems follow.

**C1 — §4 has no `pool × location` term, and that is the largest variance
component in the power formula.** For oat at the design geometry:

```
2σ²within/n  =  8.91      accessions
4σ²e/P       =  6.55      plots
(int·Δ)²/5   = 28.17      pool x location   <- 65% of the total
```

The power calculation says this term dominates. The analysis model omits it
entirely. An omitted pool-level term that varies by location does not vanish;
it inflates the test. The worst case is concrete and plausible: associate
effects that are strongly location-specific with a near-zero mean across
locations would produce a large `pool:location` variance and no real pool main
effect, and the §4 model would reject it as a main effect. The null simulation
cannot see this, because setting all As to zero removes the interaction along
with the effect — so the reported false-positive rates (0.028 oat, 0.025 pea)
are measured under a null that omits the one thing worth worrying about.

**C2 — the dominant term has 4 df, and the test uses 38.** `contrast_power()`
adds the interaction variance and keeps `df = 2n − 2`. Satterthwaite on the
actual variance split:

| | naive df | effective df | power as quoted | power corrected |
|---|---|---|---|---|
| oat (n = 20) | 38 | **9.3** | 0.779 | **0.735** |
| pea (n = 30) | 58 | 21.4 | 0.839 | 0.825 |

4.4 points for oat, 1.4 for pea. Smaller than the structure suggests, and I
checked rather than assuming — but it compounds with B, and it bites harder if
a location is lost: at four locations and 320 plots, oat power falls from 0.78
to **0.64**.

**Fix.** Pick one analysis, pre-register it, and compute power for that one.
I would add `(1|location:oatPool)` and `(1|location:peaPool)` to §4, drop
`(1|combination)` to a secondary model (it is estimable almost entirely from
four anchor combinations and will sit near the boundary), and extend the null
simulation to generate As × location variance with zero mean.

### D. The simulation does not check what it is said to check

`VALIDATION_DESIGN.md` §5 and the header of `validate_power.R` both present
the simulation as the test of whether selection on noisy predictions delivers
the predicted contrast — "the only way to capture selection on noisy
predictions."

It cannot be. `draw_truth()` sets

```r
true = lambda * (BLUP + rnorm(0, sqrt(PEV)))
```

so `E[true | BLUP] = λ · BLUP` **by construction**. The simulated contrast
must come back at λ·ΔAs, and the table reporting "predicted 14.3, simulated
14.3" is reporting an identity, not an agreement. Mean-zero error around the
BLUP adds variance; it does not attenuate the mean. The simulation *assumes*
the BLUPs are calibrated, which is precisely the assumption at issue.

This matters because it makes the evidence base thinner than it looks. The
simulation is a genuine and valuable check on **the SE formula and the
design's balance** — that partner effects really do orthogonalise, that the
layout behaves. It is no evidence at all about whether the contrast
materialises. The only evidence for that is λ, from eight folds with a 95%
interval of 0.41–1.25 and the year structure of finding B.

A secondary symptom: §5 explains the simulation's lower SE as "the design's
balance removes partner variance a little better than the formula assumes."
That is backwards — the analytic formula contains no partner term at all
(finding E), so it should be the *smaller* of the two. The gap is the λ
convention, which `contrast_power()`'s own comment explains correctly.

**Fix.** Rewrite the claim to what the simulation does. If regression to the
mean is to be probed, it needs a truth generated independently of the
BLUPs — e.g. simulate a population, sample a training design like B4I's, fit,
select, and score. That is a bigger job and probably not worth it before field
season; saying plainly that λ is the only evidence is worth doing today.

### E. Two variance components are set to zero without being named

`contrast_se()` has three terms. Two real contributors are absent.

**Specific combination (oat × pea).** Each accession meets about nine partners,
so σ²_comb enters each accession's mean at roughly σ²_comb/9 — not negligible.
It is set to zero because the current design cannot estimate it, which is the
project's central finding. Sizing the experiment that is meant to *measure*
σ²_comb on the assumption that σ²_comb = 0 deserves at least a sensitivity
band (say σ²_comb at 10%, 20%, 30% of σ²_e) so the exposure is bounded.

**Partner identity.** The design balances partner *pool*, not partner
*identity*, so each oat's mean pea yield carries the mean of its partners'
producer effects, contributing about Var(pea Pr within pool)/9. The simulation
includes this; the analytic formula, which is the one quoted, does not. It is
small here — because the pools are Pr-selected and therefore Pr-restricted —
but it should be in the formula rather than discovered by the simulation and
then misexplained.

### F. Absolute versus proportional scale is flagged and left open, and it decides λ

§6 ends with: "The model assumes associate effects are constant in absolute
g/m², while the trials differ several-fold in mean yield … note that
SIMULATION.md models environments as *scaling* the whole signal, which is the
opposite convention. Worth reconciling."

This is not a footnote. It determines:

- whether AL is uninformative or merely low-signal (and therefore whether the
  SD floor is a scale correction or an outcome filter);
- what λ means, since a slope in raw g/m² pooled across trials spanning
  8–383 g/m² is dominated by the high-yielding trials;
- **whether `interaction_frac` is interaction at all.** `response_sd` across
  the oat folds runs 5 to 51 — a tenfold range. A slope estimated on an
  absolute scale across environments that differ tenfold in spread will vary
  across folds for reasons that have nothing to do with genotype ×
  environment. Since this term is 65% of the oat variance in the power
  formula, the question is load-bearing.

**Fix, and it is cheap and decisive.** Re-estimate λ with the response
standardised within trial, and compare `interaction_frac`. If it collapses,
much of what is being charged as interaction is a scale artefact, power rises,
and "pea is marginal" may not survive. If it does not collapse, the
interaction is real and finding C1 becomes more urgent, not less. Either
result is worth having, and the code already records `response_sd` per fold —
it is used only for the SD floor and then discarded.

### G. Trial inclusion is partly selected on the outcome

§7b is the best-argued section in the repository. Its case against judging a
trial by what it does to power is correct, well reasoned, and the right call.

But the rule it puts in place is still a rule on the outcome variable:
"Exclude a trial when its own fold has λ at or near zero, a non-significant
`lambda_p`, `r_accession` near zero, and a large positive `narrows_by`."
`narrows_by` is defined as how much the across-fold spread *falls* when the
trial is removed. A rule that removes folds which disagree with the others
will, applied as stated, raise mean λ and lower `interaction_frac` — which is
what it did (AL: "Dropping it raises mean λ by 0.10 on each species"). That is
not a criticism of excluding AL, which is amply justified on agronomic grounds
and by the QC screen independently. It is a point about what the resulting
number means: **λ = 0.83 is conditional on a trial-inclusion rule that was
chosen partly by looking at λ.**

`lambda_p`, which the rule uses as a criterion, is also too small. `calib()`
fits `lm(peaYield ~ oat_As + pea_Pr)` at the **plot** level, treating plots as
independent when accessions recur across plots and plots sit in blocks. The
slope is fine; its *p*-value is not, and the rule leans on it. Only the trial
mean is removed, not block, so block × genotype confounding enters too.

**Fix.** Report both the conditional and the unconditional λ — both are
already in `crossval_summary.csv` — and label the first as conditional. Fit
the calibration at the accession level, or cluster the SE by accession, before
`lambda_p` is used to decide anything.

### H. The experiment validates a parameter; it does not test the breeding claim

The chain's objective is breeding oats for intercropping. What the trial tests
is: among high-producer accessions, do those predicted high-As raise partner
yield more than those predicted low-As? A pass licenses ranking on the As
BLUP. It does not establish that the index is worth selecting on, for two
reasons.

**H1 — total productivity is never tested, and it is free.** §8 item 5 records
monoculture checks as "omitted; would need extra plots" — correct, and so LER
is genuinely out of reach. But **total intercrop yield needs no extra plots at
all.** Both yields are already recorded on every plot. Pr and As are
negatively correlated (−0.22 oat, −0.37 pea in the current fit), so the As+
pool is a constrained selection, and whether it produces a better *system* is
exactly what a breeder will ask first. Oat + pea yield, and an
economically weighted index, should be pre-specified co-primary or at minimum
named secondary outcomes. Leaving them unstated risks the worst outcome: a
significant As contrast alongside a non-significant or negative total, decided
after the fact.

**H2 — there is no non-genomic benchmark.** The claim at stake is that
*genomic prediction* of associate effects works. What will be tested is that
*these BLUPs* rank accessions. If the As+ pool is simply the shorter, later,
lower-biomass oats, the trial will succeed and the method will have added
nothing over a tape measure. The fix is cheap: **record oat height, biomass,
flowering and maturity in the validation trial** and regress realised
associate effect on them. If the As BLUP retains signal after a phenotypic
proxy is accounted for, the result is a genomic-prediction result. If it does
not, that is worth knowing, and far better known from these plots than from
the next three years of selection. A third pool selected on a phenotypic proxy
would be stronger still, at the cost of ~25% more plots — and §5 shows plots
are the cheapest thing in this design.

### I. The more informative estimand is listed as secondary

§4 makes the pool contrast primary and the predicted-vs-realised correlation
secondary — while describing the latter as "a continuous validation, **more
informative** than the binary contrast." That is an internal contradiction,
and the continuous version is the better primary on three counts: it uses all
2n accessions rather than collapsing them to two means; its slope is λ
measured on new data, which is the quantity that feeds the next design and the
one every refresh is conditioned on; and it is directly comparable with the
cross-validation folds in §6. The pool contrast is a two-point summary of the
same regression.

**Fix.** Make the slope of realised on predicted the primary estimand, with
the pool contrast as its pre-registered two-point summary. This also dissolves
§8 item 4 (multiplicity across two primary tests): one estimand per species,
stated directionally.

---

## 4. NOTED

### J. Reliability is global, not per accession

`validation_inputs()` derives one reliability per species from
`Var(BLUP)/σ²` and gives every accession the same PEV. An accession with one
partner and one with twenty get the same prediction-error variance in
`σ²within`. A per-accession posterior SD would make the eligibility filter a
*reliability* filter rather than a partner count (which is the right form of
finding A1's fix, not just a recount), and would let `σ²within` be
pool-specific. The effects files save only four seed-level point estimates
per accession, so this needs BGLR's draws retained — a small change to
`BGLR_multi_trait_model.R`, worth making at the same time as A.

### K. Document drift

The `--refresh-doc` mechanism is sound and its limitation is honestly
documented. These are the sentences that have gone stale and now mislead:

| location | says | is |
|---|---|---|
| §5 prose | "reliabilities are low — oat `As` 0.31, pea `As` 0.24" | 0.35 and 0.27 (§5 table, same section) |
| §6 | "λ comes from four usable folds" | eight |
| §6 | "λ > 1 for oat means the realised contrast *exceeds* the predicted one" | λ_oat = 0.83; leftover from when only 2025 folds existed, and now the actively misleading reading (finding B) |
| §6 rehearsal table | five 2025 trials | nine trials now have folds |
| §7 | "Nothing in the scripts is hard-coded" | finding A |
| §7 | the 2026 trials moved the BLUPs | they did not (finding A3) |
| §5 vintage | 9 trials, 3,567 plots | the fit saw 5 and 1,985 |

### L. Two false comments in the code

Fixed in this commit rather than reported, per standing practice:
`validate_power.R`'s header and `simulate_validation()`'s docstring both claim
the simulation captures the attenuation from selecting on estimates. It does
not (finding D).

---

## 5. What I would do before ordering seed

In order. The first three are a day's work and change the numbers.

1. **Fix the trial set (A).** Derive `trials` from `output/trial_qc.csv`;
   write the fitted trial set out with the effects; make
   `validation_vintage()` read it from there and assert it matches the
   phenotype table. Re-run `validate_refresh.R`.
2. **Read `pool_diff.csv` first (A4).** This is the first honest test of the
   churn criterion. If retention is below 70% with three genuine new trials
   added under a fixed rule, §7's own standard says the effects are not stable
   enough to validate — and that is the finding, not a setback.
3. **Recount eligibility on the fitted data (A1).** Expect roughly 248 oat and
   244 pea candidates rather than 430 and 407, and expect the pools to change
   substantially. Better: filter on per-accession reliability (J).
4. **Re-estimate λ on a within-trial standardised scale (F)** and compare
   `interaction_frac`. Cheap, decisive, and it determines 65% of the oat
   variance.
5. **Quote λ with year in it (B)** and re-size on the forward-looking value,
   ~0.6. Settle §8 item 3 against that number.
6. **Pre-register one analysis (C).** Add `pool × location`; compute power for
   the model that will be fitted; extend the null simulation to generate
   As × location with zero mean.
7. **Add the free outcomes (H).** Total yield as co-primary; height, biomass,
   flowering and maturity recorded on every plot. Decide whether a
   phenotypic-proxy third pool is worth 25% more plots — §5 says plots are
   cheap, which makes this the best available use of them.
8. **Make the continuous slope primary (I).**

If 2 and 5 both come back badly — pools churning past 30% and λ near 0.4 —
the honest conclusion is that the 2026 trials have not settled the effects
enough to validate, and the better use of a field season is a design that
*estimates* what the current one cannot: replicated oat × pea combinations, a
proper factorial, and the specific-combination variance the whole project has
concluded is missing. That option is not currently on the table anywhere in
the documentation, and it should be, because the validation trial as designed
spends 40 of 400 plots gesturing at it through anchors.

---

## 6. Severity, honestly

Of the findings above, **A is a defect** — the numbers do not describe the
fit. **B, C1, F and H are judgement calls** where I think the current
treatment is wrong but a reasonable person could argue. **D, E, G and I are
overstatements or omissions** in how the evidence is described, not errors in
what was computed. **J, K, L are housekeeping.**

What A costs is not mainly power. It is that **12 accessions in the pools have
their producer and associate effects perfectly aliased**, and the pools are
built on separating exactly those two things. That is the kind of error that
survives a passing test suite, a successful refresh run, and a careful
reading — because every component works and the composition is what is wrong.

---

## 7. What I think is right and should not be changed

A critique that finds everything wrong is not useful. These are the decisions
I checked hardest and would defend:

- **Estimating λ instead of assuming it, and running it before power.** This
  is the single best structural decision in the chain.
- **The accession as the experimental unit, with PEV in `σ²within`.** This is
  the correct broad-inference denominator, and the conclusion it produces — that
  most of the SE cannot be bought with plots — is both right and the most useful
  thing in VALIDATION_DESIGN.md. Applying λ to the contrast while leaving
  `σ²within` at full size is the conservative choice and is correctly justified.

  **One qualification, added 2026-10-03.** I wrote that the term was "correctly
  derived (`Var(true) = Var(BLUP) + PEV`)". The *structure* is right and I would
  still defend it. The *value* was not: PEV was being backed out of that identity
  rather than measured, and the identity fails badly in this fit — for oat,
  `Var(BLUP) = 46.9` plus a measured `E[PEV] = 49.0` gives 95.9 against a fitted
  component of 160.9, with `mean(diag(G)) = 1.00`. The backed-out PEV was 2.3×
  the posterior's own value, which inflated `σ²within` and *understated* power.
  Measuring it per accession from the streamed draws is the fix.
- **One index for both pools (`single-v1`).** The fix is right and the retained
  account of the old defect — a power curve shaped by a flaw in the selection
  rule will argue for changing the design when the rule is what needs
  changing — is the most transferable lesson in the repository.
- **The scan rather than the bisection.** ΔPr(θ) is a step function with
  multiple sign changes; scanning is correct and the tie-break on ΔAs is a
  real subtlety handled properly.
- **Both species from one set of plots.** Reading the associate effect on the
  partner's yield means the 2×2 factorial validates oat and pea at once. Clean
  and genuinely elegant.
- **Partner-pool balance at zero imbalance for every accession.** The
  load-bearing property, asserted in code and in the tests.
- **The anchors' two-part justification.** Within-location duplication for
  plot error, fixed across locations for combination × location. Both halves
  are needed, both are asserted, and getting this wrong is easy.
- **§7b.** The argument against judging a trial by its effect on power is
  correct, non-obvious, and important. My criticism (G) is about the
  replacement rule's conditioning, not about the argument.
- **Renaming rather than dropping collapsed accessions**, and §7c's warning
  that a pooled analysis name is not a seed lot. That warning will save
  somebody a seed order.
- **Recording retractions in place.** The §5 note that the section said the
  opposite for one day, and the ND Victory retraction in BACKGROUND.md, are
  why this review could reconstruct what happened. Keep doing it.

---

## 8. What I could not check

- **Whether re-fitting on nine trials changes the conclusions.** It needs a
  ~20-minute BGLR run per seed and I did not launch one. Everything in finding
  A is about the mismatch, which is established; the *magnitude* of the
  correction is not.
- **Whether the three 2026 trials are themselves sound.** The QC screen passes
  ND and NY and flags 2026_IA for low repeatability. Their oat folds give the
  lowest λs in the set (0.26–0.53), which is either the honest forward-looking
  signal or a sign of a problem in the trials. I cannot separate those.
- **The BGLR fit's internals.** Convergence, mixing, prior sensitivity, and
  whether the `BRR`-on-`ZL` construction returns exactly `Var(Lb) = G ⊗ Σ` in
  the implementation. The four-seed agreement in the variance components is
  reassuring about mixing and says nothing about prior sensitivity.
- **The curation thresholds.** CURATION.md's bimodality argument is persuasive
  and I did not re-derive it from markers.
- **Whether the GRMs are right.** Marker QC, imputation and allele-frequency
  centring all sit upstream of everything here and were taken on trust.
- **Seed availability** (§8 item 2), which §8 correctly calls the real
  constraint on *n* and which no analysis can settle.
