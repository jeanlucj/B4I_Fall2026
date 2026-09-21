# B4I follow-ups

Answers where the data settled it, open items where it did not.
Written 2026-09-19. Item 3 has since been folded into the curation scripts; the rest are open.

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
