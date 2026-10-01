# B4I follow-ups

Answers where the data settled it, open items where it did not.

Items 1–7 were written 2026-09-19; item 3 has since been folded into the curation
scripts and the rest are open. **Item 8, added 2026-09-30, is the brief on where
the model comparison stands** — written to be read cold, so it recaps the
reasoning rather than only the conclusions.

---

## Confirmed: the pea merge

`NDP170084G` and `ND VICTORY` correlate above 0.99 and were merged.

- 10 plots recorded as `NDP170084G` now carry `ND VICTORY`; that accession has 32 plots in the
  analysis, up from 22.
- `output/pea_analysis_names.csv` records the mapping and the reason.
- The pea GRM was collapsed 435 -> 434 lines, averaging the two rows and columns, so the
  relationship matrix and the phenotypes agree on the name.
- No plot in the fitted table carries `NDP170084G`.

The oat side collapsed 508 -> 468 the same way.

---

## 1. Pea pedigrees

- [ ] **Find a source of pea pedigrees.** Every one of the 435 B4I pea accessions reads
      `NA/NA` on T3/Oat — not a parsing problem, the field is genuinely empty.
- [ ] Decide whether to load them into T3 or keep them beside the project.
- [ ] Once loaded, re-run `code/curate_pea_accessions.R`: it tests for pedigrees each run and
      will start reporting full-sib families and parent comparisons without modification.

Worth noting this is currently a gap rather than a problem: the pea marker profiles are clean
(one duplicate pair in 415 genotyped accessions), so there is no clonal-family puzzle waiting
for pedigrees to explain. Pedigrees would add power to the GRM, not fix a defect.

## 2. Duplicated plot records in IA and ND

**Affected:** `B4I_2025_IA` (trialDbId 6822) and `B4I_2025_ND` (6823).
**Not affected:** `B4I_2025_AL` (6821), `B4I_2025_NY` (6824), `B4I_2025_IL` (6881),
`B4I_2026_IL` (6997).

They are stored differently, and the difference is measurable:

| trial | units | levels present | plot records | distinct plots |
|---|---|---|---|---|
| B4I_2025_AL | 1600 | plot, subplot | 400 | 400 |
| B4I_2025_NY | 1600 | plot, subplot | 400 | 400 |
| **B4I_2025_IA** | **2800** | plot, subplot | **1600** | **400** |
| **B4I_2025_ND** | **2800** | plot, subplot | **1600** | **400** |
| B4I_2025_IL | 400 | plot only | 400 | 400 |
| B4I_2026_IL | 400 | plot only | 400 | 400 |

AL and NY have the same subplot structure (1200 subplots, 3 per plot) and do *not* duplicate,
so subplots alone are not the cause.

**Mechanism found.** In IA and ND each plot-level `observationUnitDbId` is returned four times,
and the copies are identical except for `positionCoordinateX` / `positionCoordinateY`, which
take the values (0,0), (2,2), (0,2) and (2,0). Each copy carries zero observations. So the same
plot is being emitted once per cell of a 2x2 within-plot grid. It is not a pagination artifact:
the counts are identical at pageSize 500 and 5000.

- [ ] Ask whoever uploaded IA and ND why plot positions are on a 2x2 grid when AL and NY are not.
- [ ] Check whether the 2x2 grid is meaningful (a real sub-plot sampling layout) or an upload
      artifact. If meaningful, the subplot level should carry it, not the plot level.
- [ ] Confirm no *observations* were duplicated — the copies carry none, and the yields joined
      by `observationUnitDbId`, so the fitted data is unaffected.

`code/assemble_B4I_phenotypes.R` deduplicates on `observationUnitDbId`, which is safe given the
copies hold no distinct data. Left undetected this would have quadrupled IA and ND in the fit.

## 3. IL18-8562/IL18-805 versus IL18-8562/IL10-9872

**They are not distinguishable as separate crosses.**

| comparison | n pairs | min | median | mean | max |
|---|---|---|---|---|---|
| within IL18-8562/IL18-805 (8 members) | 28 | 0.9744 | 0.9853 | **0.9850** | 0.9956 |
| within IL18-8562/IL10-9872 (3 members) | 3 | 0.9960 | 0.9965 | **0.9971** | 0.9988 |
| **between the two** | 24 | 0.9753 | 0.9799 | **0.9827** | 0.9933 |

Between-cross correlation (0.983) is essentially the same as within the larger cross (0.985).
The two "crosses" are one clonal pool descending from IL18-8562, not two.

**IL23-2963, now `IL18-8562_self`** (recorded pedigree IL18-8562/IL14-2453), sits inside that
pool rather than beside it: mean 0.9854 against the ×IL18-805 members and 0.9839 against the
×IL10-9872 members — indistinguishable from an ordinary member of either.

- [x] **Merged.** `pool_clonal_families()` in `code/curation_functions.R` now pools clonal
      families that share a seed parent when their between-family correlation is within
      `pool_tolerance` (0.01) of the weaker within-family correlation. For IL18-8562 that is
      0.9827 against 0.9850, so the two families and IL23-2963 collapse into a single
      `IL18-8562_self` covering 12 accessions and three recorded pedigrees. The rule is
      data-driven, not a hard-coded exception: no other seed parent heads more than one clonal
      family, and none would merge on these numbers.
- [ ] Look at the boundary. IL23-2963's highest correlation with an accession *outside* these
      groups is 0.982, about the same as its correlation with members. The clonal pool probably
      extends past the 0.99 cut; a threshold sweep would show how many more accessions belong.

## 4. Missing blocks in AL, IA, ND and NY

You are right that this is a problem. What T3 holds:

| trial | rep values | block values | blocks |
|---|---|---|---|
| B4I_2025_AL | 1 | **1** | 1 |
| B4I_2025_IA | 1 | **2** | 1 |
| B4I_2025_ND | 1 | **4** | 1 |
| B4I_2025_NY | 1 | **5** | 1 |
| B4I_2025_IL | 1 | 1–8 | 8 (50 plots each) |
| B4I_2026_IL | 1 | 1,2 | 2 (200 plots each) |

The block field is populated in all six trials, but in four it holds a single constant — and a
different constant each time (1, 2, 4, 5). That pattern looks less like "blocks were not
recorded" and more like a single block *number* was attached to the whole trial, possibly the
trial's position in a larger planting. `rep` is 1 everywhere and carries nothing.

- [ ] Ask the AL, IA, ND and NY cooperators whether the trials were blocked in the field.
- [ ] If they were, recover the block assignments and re-upload; a 400-plot trial with no
      blocking is losing real precision.
- [ ] Find out what the constants 1, 2, 4, 5 mean — they may be the key to recovering the design.

Until then `code/BGLR_multi_trait_model.R` fits block effects only for the 10 blocks in the two
IL trials. A constant block within a trial is perfectly aliased with that trial's fixed effect,
so including those four would add a redundant random effect and hurt mixing without adding
information.

## 5. Credible intervals — answered

Computed from the saved MCMC samples (`output/BGLR_seed_*_Omega_*.dat`, `*_R.dat`), 6800 draws
pooled over four chains after discarding burn-in.

**Variance components, median [95% credible interval]:**

| component | median | 95% CrI |
|---|---|---|
| var Pr, oat (on oat yield) | 390.0 | [281.4, 536.9] |
| var As, oat→pea (on pea yield) | 104.6 | [70.5, 152.4] |
| var Pr, pea (on pea yield) | 181.7 | [124.8, 264.1] |
| var As, pea→oat (on oat yield) | 295.2 | [215.4, 398.2] |

All four are comfortably away from zero.

**Producer–associate covariances — neither is distinguishable from zero:**

| covariance | median | 95% CrI | P(< 0) | excludes 0 |
|---|---|---|---|---|
| oat Pr–As | −12.3 | **[−69.8, +36.7]** | 0.69 | no |
| pea Pr–As | −53.0 | **[−118.4, +1.8]** | 0.97 | no (only just) |
| residual oat–pea | −99.8 | **[−137.2, −64.5]** | >0.999 | **yes** |

On the correlation scale: oat −0.062 [−0.304, +0.189], pea −0.231 [−0.434, +0.009], residual
−0.121 [−0.165, −0.079].

**This tempers the headline.** The negative producer–associate correlations reported from the
point estimates are not supported as non-zero. The pea one has 97% posterior probability of
being negative, which is suggestive and worth pursuing, but its interval touches zero. The only
cross-trait covariance the data establishes is the residual, which is clearly negative:
competition within a plot, after trial, block and genetics are accounted for.

- [ ] Treat the pea Pr–As correlation as a hypothesis, not a result, until the design supports it.
- [ ] Re-check after the combination replication and blocking issues above are resolved — both
      cost precision exactly where this parameter needs it.

## 6. Single-partner accessions — your hypothesis was right

It is almost entirely 2026, which is one trial.

| species | year | accessions | median partners | with 1 partner |
|---|---|---|---|---|
| oat | 2025 | 240 | 6 | 2 (0.8%) |
| oat | **2026** | 234 | **1** | **126 (53.8%)** |
| pea | 2025 | 238 | 7 | 6 (2.5%) |
| pea | **2026** | 227 | **2** | **111 (48.9%)** |

Distribution of partner counts (capped at 10):

```
oat   1   2   3   4   5   6   7   8   9  10+
2025  2   7  17  27  43  34  35  26  21  28
2026 126  74  30   4   0   0   0   0   0   0

pea   1   2   3   4   5   6   7   8   9  10+
2025  6  10   6  20  27  35  30  34  21  49
2026 111  82  31   3   0   0   0   0   0   0
```

Within 2025, connectivity is good: 5 trials, median 6–7 partners, almost nothing with a single
partner. The pooled counts of 99 oat and 84 pea are accessions that appear only in 2026.

- [ ] Revisit once the other three 2026 trials (IA, ND, NY) are harvested — they are already in
      T3 with accessions attached and will connect the 2026 entries.
- [ ] Consider fitting 2025 alone as a sensitivity check; it is a well-connected design and would
      show how much the 2026 singletons are dragging on the Pr–As estimates.

## 7. The bivariate model is not cross-validated

`code/BGLR_multi_trait_model.R` reports variance components, per-accession
effects and credible intervals, but nothing out of sample. The four-chain rank
agreement it prints is a convergence diagnostic, not an accuracy estimate.

Sub-objective 1.4 asks for cross-validated accuracy of intercrop-related
breeding values per crop, which needs whole **accessions** held out, not cells.

- [ ] Implement accession-wise cross-validation for the bivariate model.
- [ ] Score producer and associate effects separately; they differ in accuracy
      and the associate effect is the harder one.
- [ ] Fold by family rather than at random, so clonal groups and full-sib
      families do not leak across folds and flatter the accuracy.
- [ ] Report the 99 oat and 84 pea single-partner accessions separately; their
      producer and associate effects are aliased and will predict badly whatever
      the model does.

The design is sketched in [CROSS_VALIDATION.md](../CROSS_VALIDATION.md#what-a-cross-validation-of-it-would-look-like).

---

## 8. Where the model comparison stands — a brief, 2026-09-30

Written to be read cold. If you have forgotten everything else, read this
section and then
[docs/BGLR_vs_MegaLMM.md](BGLR_vs_MegaLMM.md).

### The question

Two frameworks can estimate what an accession does for itself (**producer**) and
what it does to its partner (**associate**) in an oat–pea intercrop:

- **Bivariate DGE-IGE**, `BGLR::Multitrait`. Both yields analysed jointly, four
  genetic effects, kinship on all four through `Σ ⊗ G`. This is
  `code/BGLR_multi_trait_model.R`, the production analysis.
- **The MegaLMM factor model.** Treats the oat × pea matrix as a multi-trait
  problem — one species on the rows with kinship, the other as "environments" —
  and finds low-rank latent structure.

The hope for MegaLMM was never better general mixing ability. It was that a
factor model might **decompose the oat × pea interaction** into a few
interpretable dimensions, which a Kronecker kernel cannot do: a kernel gives you
a covariance, not a set of axes you could look at biologically.

### Why answering it needed a rewrite

The original simulation drew only **two** effects and one response — the oat's
producer effect and the *pea's* associate effect, both on oat yield. The oat's
associate effect and the pea's producer effect were not simulated at all, so half
the producer–associate model was untested and the DGE-IGE comparator had to be
univariate. It also ran MegaLMM in one orientation, which cannot see an oat's
associate effect, because that effect is read on **pea** yield.

So the generator went bivariate (four effects, two responses, a low-rank
interaction per trait), DGE-IGE became the real `Multitrait` production model,
and MegaLMM is now fitted in **both** orientations with each species' GMA
assembled across the two.

### The structural fact that explains most of the results

**MegaLMM has a per-column intercept and no per-row one.** So in either
orientation:

| | row species | column species |
|---|---|---|
| which effect is the margin | **producer** | **associate** |
| how it is estimated | latent factors + `U_R`, **kinship-shrunk** | the per-column intercept, **fixed and unshrunk** |

Almost everything below follows from this. It is why `megalmm_U` — scored on `U`,
which excludes that intercept — recovers the associate effect at r ≈ 0 while
`megalmm` does not; why the Eta-over-U gap **widens** with more data (0.06 → 0.26
on oat GMA from 1.6% to 48% observed) rather than closing, since the intercept is
the part more data estimates best; and why pinning the main effect helps when data
is scarce.

`Multitrait` has no such asymmetry: kinship is on all four effects. When there is
little information per accession, shrinkage is worth a great deal — which is the
whole GMA result in one sentence.

### What is settled

**1. Sparsity is not one factor among several. It is the experiment.** Partial
ω² = 0.97 on GMA, and r runs from 0.15 at 1.6% observed to 0.91 at 48%. Nothing
else is within an order of magnitude. ([docs/MegaLMM_anova.md](MegaLMM_anova.md))

**2. For GMA, use the BGLR `Multitrait` model.** It beats MegaLMM in 95–99% of
scenarios, and the gap is **largest where B4I sits**: −0.24 at 1.6% observed,
−0.006 at 48%. The two BGLR models (`additive` and `dge_ige`) are
indistinguishable from each other on GMA, which is the expected sanity check —
adding a specific-combination term should not disturb the main effects.

**3. Of MegaLMM's three levers, one matters, one is marginal, one is inert.**
`fixed_main_effect` is real (ω² 0.14), `K` is small but consistent with K = 5
ahead of 10 (ω² 0.06–0.10), and **`eigen_variance` does nothing at all** —
negative partial ω² on all four responses, meaning it explains less than its one
degree of freedom would by chance. 0.25 and 0.75 give 0.705 and 0.704.

**4. Pin the main effect when sparse, not when dense — and not at all if you want
the interaction.** On GMA, pinning is worth +0.11 at 1.6% observed and +0.20 at
4.8%, then stops helping by 16%. On interaction recovery it *costs* accuracy at
every density (0.525 → 0.478). There is no single rule that is right for both
responses.

**5. GxE hurts; the number of environments does not.** Ten *stable* environments
cost essentially nothing against one, and ten with GxE cost about 0.06 — and this
holds for both frameworks, so it is a property of the data rather than of either
model. Raw means of `r_addsurf_oat`:

| | one | ten_stable | ten_gxe |
|---|---|---|---|
| `dge_ige` | 0.903 | 0.907 | 0.844 |
| `megalmm` | 0.673 | 0.670 | 0.602 |

**6. MegaLMM does win on the interaction — conditionally.** Above a density
threshold and at low rank:

| `r_int_oat` gap (MegaLMM − `dge_ige`) | 1.6% | 4.8% | 16% | 48% |
|---|---|---|---|---|
| **rank 1** | −0.19 | −0.19 | **+0.37** | **+0.57** |
| **rank 5** | −0.21 | −0.32 | +0.04 | +0.36 |

Density is a **gate, not a modifier**: below it MegaLMM is worse at either rank.
Above it, rank sets the size of the win.

**7. Therefore the biological-*decomposition* hope is not available from the
current data — but the interaction itself is not out of reach.** B4I observes
about 3% of the matrix, with 91% of observed combinations in a single plot, which
is below the gate. The two models fail differently there, and the difference is
the actionable part:

| observed | `dge_ige` `r_int_oat` | `megalmm` `r_int_oat` |
|---|---|---|
| 1.6% | **0.205** | −0.002 |
| 4.8% | **0.285** | 0.029 |
| 16% | 0.416 | **0.549** |
| 48% | 0.479 | **0.852** |

The Kronecker kernel **degrades gracefully** — still r ≈ 0.2–0.3 at a few percent
observed, because it borrows across relatives rather than needing structure
repeated. MegaLMM is at **zero** there and then overtakes sharply. So below the
gate the move is not "abandon the interaction" but "use the kernel, and do not
expect interpretable axes". Reaching roughly 16% observed is what would buy the
decomposition, and it would be MegaLMM specifically that delivered it — a
statement about designing the next round of trials.

**8. Rank costs about twice what variance share buys.** Going rank 1 → 5 costs
~0.15 in `r_int`; doubling the interaction's variance share 10% → 20% buys ~0.07,
additively. A rank-5 interaction at 20% is harder to recover than a rank-1 one at
10%. This matters because a programme can influence how much specific combining
ability there is far more easily than its rank.

### What is unsettled, worst first

**1. The interaction comparison has a live confound: `dge_ige`'s Kronecker term
is severely rank-truncated and MegaLMM's is not.**

`SIM_KRON_RANK = 30` takes the 30 leading eigenvectors of each GRM, giving a
900-column basis. Measured on a 400-accession panel, 30 directions capture 62.7%
of the oat GRM's trace and 54.8% of the pea's — so the term can represent about
**34% of the Kronecker covariance, and 900 of 160,000 dimensions**. Reaching 95%
of the trace would need ~187 and ~204 directions, a basis of ~38,000 columns.

So "MegaLMM beats a Kronecker kernel at recovering the interaction" is, strictly,
"beats a **rank-30-truncated** Kronecker kernel". The basis is chosen by GRM
eigenvalue, not by where the interaction actually lives, so a rank-1 or rank-5
true interaction may fall largely outside it.

**`code/sim_kron_study.R` settles this, and it needs no extrapolation.** The
trick is that "full rank" is reachable after all. The *exact* kernel cannot be
used — it is defined only over observed combinations, so it returns nothing for
the never-observed cells where `r_int` is scored (measured: it gives `NA`). But at
`kron_rank = n_acc` the two bases span their whole spaces, so their row-wise
Kronecker product spans the **entire** Kronecker space: trace retained = 1,
nothing truncated, and still a full surface. At n = 60 that is 3,600 columns and
takes 18 seconds.

So the decision rule is direct. Define `delta(rank) = r_int(megalmm) −
r_int(dge_ige at rank)`; the published finding is `delta > 0`. Run at
`rank = n_acc` and you have `delta` at the **ceiling** — no truncated basis can
beat it. If `delta(full) > 0`, **no value of `SIM_KRON_RANK` can flip that
cell.**

A single-replicate check at n = 60, 48% observed, rank-1 interaction already
points one way:

| `kron_rank` | trace kept | `r_int_oat` |
|---|---|---|
| 10 | 24.5% | 0.259 |
| **60 (full)** | **100%** | **0.343** |
| `megalmm` | — | **0.815** |

Going from a quarter of the trace to **all** of it gains `dge_ige` 0.084, against
a gap of 0.47–0.56. If that holds up across cells and replicates, the truncation
costs `dge_ige` something real but nowhere near enough to explain the gap, and
finding 6 stands as stated. **One replicate of one cell is not the answer** — the
full design sweeps n = 60/100/200, 16% and 48% observed, and interaction ranks 1
and 5.

Note this affects the **simulation only**. The real analysis uses `kron_rank =
NA`, the exact kernel `G_oat[i,i'] · G_pea[j,j']` over observed combinations, with
no truncation beyond numerically null directions — tractable only because there
are 2,059 combinations.

**2. The gate's location, and the rank boundary, are both unlocated.** Four
sparsity levels cannot say where between 4.8% and 16% the switch happens, and
rank 1 against rank 5 says nothing about what lies between.
`code/sim_int_run.R` is the full factorial built to settle both — it adds 9%
sparsity and rank 3, and crosses the pinning switch instead of fixing it.
**It has not been run.**

**3. Everything in finding 6 and in MegaLMM_anova.md was measured under the old
scoring scheme.** Until 30 September the simulation held out 20% of the
*observed* cells, so models were fitted at **0.8 × the labelled sparsity** and
the per-cell metrics were scored on a small held-out subset. Read every sparsity
level in those documents as 0.8 × its label: the gate quoted as "between 4.8% and
16%" is between **3.84% and 12.8%** actually fitted. The scheme has since changed
— all observed plots are fitted and the metrics are scored on the never-observed
cells — but **the general simulation has not been rerun and the ANOVA has not
been redone**, by choice. `sim_int` uses the new scheme, so its 4.8% and 9% levels
both fall inside the old gate interval, which is convenient.

**4. Why MegaLMM collapses to zero on the interaction below the gate is not
established.** It is not that the task is impossible there — the Kronecker kernel
manages r ≈ 0.2–0.3 at the same densities (finding 7). Something specific about
the factor model fails, and the candidates have not been separated: too few
observations per column to estimate the loadings; the factors being spent on the
main effects instead; or the per-column intercept absorbing signal that should
have gone to the factors. `r_mainfactor_*` against `r_rowmean_*` is the diagnostic
already in the output and has not been looked at for this purpose. Knowing which
it is would say whether the gate can be moved by changing the model rather than
the design.

**5. `cor(I_oat, I_pea)` is set to zero in the simulation, and is unknown in
reality.** The two interaction surfaces are drawn independently, by choice,
because the real data cannot identify the parameter: 1,869 of 2,059 combinations
occur once, so `fit_mix_term` is off in production and specific-combination
covariance cannot be separated from plot error. The fitted residual correlation
between the two yields, −0.122, is a mixture of the two. **The validation trial's
anchor plots are the only replicated combinations in the programme and would be
the first clean measurement** — which was not why they were put there.

**6. There is no longer any number comparable to a real-data cross-validation.**
Dropping the holdout turned `r_obs_*` into `r_fit_*`, a goodness of fit measured
on the cells the model was fitted to. If a bridge to real-data accuracy is
wanted, it needs a holdout, and a holdout costs what item 3 describes.

**7. Three-way interactions are not estimable** under the 150-run fraction, so
whether the pinning crossover itself moves with panel size or with GxE is unknown.
That was the right trade at the time; it is the boundary of what that design
supports.

**8. The tuning advantage was measured on the data that chose it.** Giving
MegaLMM its best settings is worth 0.026–0.030 on GMA and 0.053–0.060 on the
interaction against an arbitrary choice. Too small to overturn the GMA result —
BGLR's margin is 0.07–0.10 — but it is a tuned-on-test figure and should be
quoted as such.

**9. [BOTH_ORIENTATIONS.md](../BOTH_ORIENTATIONS.md) idea 1 is unimplemented**:
carrying each species' kinship-shrunk main effect from its own orientation into
the other as a covariate. The orientation asymmetry above gives it a measurable
target.

### Traps in the plumbing, so they are not rediscovered

Three cost real time and all three are now guarded, but the guards are only
obvious once you know why they exist:

- **`simulation_results.csv` is rebuilt by globbing the cache**, so results from
  an earlier grid join it silently. 8,400 of 12,998 rows in the September run were
  from the pre-bivariate design. Always run
  `Rscript code/sim_filter_results.R --infer-design` before reading results.
- **Seeds are part of the cache key, and `rep` is the slow index** so that
  `--reps` is additive. It was not always: adding replicates to a finished grid
  used to renumber every seed and leave orphans that still globbed in, which put
  138 duplicated replicate-1 rows in the September results.
- **The cache filename carries a scheme marker** (`v3`). A change to what is
  scored, or to a column name, has to move it, or `bind_rows()` over a mixed cache
  produces both columns with half the rows `NA` in each.

Also: `r_addsurf_oat` and `r_oat_gma` are **different quantities** — the first is
the additive part of oat *yield* (oat producer + **pea** associate), the second an
oat accession's total contribution (oat producer + **oat** associate). They were
`r_gma_oat` and `r_oat_gma` until the names were judged too close to be safe.

### What to run next, in order

- [ ] **`Rscript code/sim_kron_study.R`** — the full-rank comparison above. The
      cheapest thing that could overturn a headline finding. `--check` is one
      cell in about a minute; the full design is n = 60/100/200 x two densities x
      two interaction ranks x 3 replicates.
- [ ] **`Rscript code/cross_validate_combinations.R`** — a proper 5-fold run.
      The `--quick` pass only verified the plumbing, and one oddity in it needs
      resolving (see [CROSS_VALIDATION.md](../CROSS_VALIDATION.md#status)).
- [ ] **Run `sim_int`** (120 scenarios × 3 replicates, ~65 single-core hours,
      ~3.3 h per task at `--array=1-20`). Clear `output/simulation_int/` on Ceres
      first: the existing cache predates both the seed renumbering and the scoring
      change.
- [ ] **Decide whether to rerun the general simulation** under the new scoring
      scheme. The conclusions do not appear to turn on it, but every sparsity
      figure in two documents currently needs a mental 0.8 ×.
- [ ] **Re-run `validate_power.R`.** Its null rejection rates (0.028 oat, 0.022
      pea) predate the anchor fix, which made replication slightly less even.
- [ ] **Accession-wise cross-validation of the bivariate model** — still item 7
      above, still the thing sub-objective 1.4 actually asks for.
