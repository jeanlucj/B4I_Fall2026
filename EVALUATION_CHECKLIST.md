# Analysis evaluation checklist

Tick these off over time. Each line names the `arm_evaluation()` group; the
step-by-step walkthrough is in **[EVALUATION.md](EVALUATION.md)** at the
matching level (A1…A12). The order is offline/fast → online/slow on purpose:
work top to bottom so a free bug surfaces before an expensive one.
`disarm_evaluation()` at the end of each level.

Bootstrap once (EVALUATION.md §3):

```r
library(tidyverse); here::i_am("code/evaluation.R")
source(here::here("code", "evaluation.R"))
eval_load("analysis"); eval_conflicts(); eval_groups("analysis")
```

## Offline — no T3 login, reads `data/` and `output/`

- [ ] **A0 `eval_conflicts()`** — **no FUNCTION defined by more than one
      loaded file** (the `connect_t3` / `collapse_grm` / `read_grm` / 
      `grm_factor` duplicates were unified 2026-09-23). The settings list
      below it — `acc_files`, `cache_dir`, `db_name`, `fills`,
      `monoculture_labels`, `out_dir`, `refresh` — is expected. Any *function*
      reappearing here is a finding.
- [ ] **A1 `grm`** — oat 508 → **466**, pea 435 → **434**; `check_grm()`
      symmetric, min eigenvalue > 0; `max(abs(tcrossprod(L) - G)) < 1e-8`;
      **8 oat / 20 pea prior-only rows** noted and accepted.
- [ ] **A2 `curation_groups`** — planted duplicate in a synthetic similarity
      matrix is found and chained correctly; representative changes when the
      observation file is supplied; real oat groups = 19 covering 90.
- [ ] **A3 `curation_family`** — planted clonal family flagged, segregating one
      not; oat 9 of 74 above 0.985 sit in a **gap** (0.992 vs 0.747 baseline);
      pea 0 of 36; every `analysis_name` resolves to a real or constructed
      name.
- [ ] **A4 plot table** — 2400 cached units → **2377** rows → **2371** with both
      yields, and the 23 lost are accounted for one by one;
      `n_distinct(observationUnitDbId) == nrow(units)` (the IA/ND
      quadruplication really is deduped).
- [ ] **A5 `matrix`** — `Y` is **237 × 222, 1682 cells (3.20%)**, rows oats and
      columns peas *verified*; centering is within trial; `X_Env` 222 × **61**
      in groups 1/2/58 with no `NaN`; exactly **6** `units_*` files in the
      cache.
- [ ] **A6 `megalmm`** — missing-data map **1 of N** chosen; `fit_X = FALSE` at
      construction; `rownames(Eta_mean)` matches `rownames(Y)` (nothing in the
      code checks this); `rowSums(Lambda^2)` tail near zero, so `K = 10` was
      generous enough.
- [ ] **A7 `cv`** — 336 of 1682 cells held out; row/column floors hold; masked
      cells really `NA` in `Y_train`; **baseline r identical between the
      covariates-on and covariates-off arms** (the leak/consistency test);
      explain why `U_noU_R` (0.072) beats `Eta_mean` (0.031) when the code
      comment predicts the opposite.
- [ ] **A8 `bglr`** — `colnames(fit$ETAHat)` is `peaYield, oatYield` **in that
      order**, confirmed against the data (`cor(ETAHat[,1], pea_yield)` = 0.91); oat `var_Pr` (394) ≫ `var_As` (105) and pea the reverse
      (186 / 299); four chains agree; 99 oat and 84 pea single-partner
      accessions noted; 10 of 14 block levels fitted.

## Live BrAPI — needs `.Renviron`; smallest first

- [ ] **A9 `t3`** — on cached trial **6821**: `intercropGermplasmName` not all
      `NA`, `rep`/`block` populated, one row per plot.
- [ ] **A10 `trials`** — 941 accession names in; both selection stages
      independent; 34 trials out.
- [ ] **A11 `geno`** — GenoPea 13K (67) covers 421 of 435, no other protocol
      above 5%; 415 carry dosages (covered ≠ genotyped);
      **`verify_parent_order(conn, peds, 25)` run for pea** — never has been;
      only 1 of 112 pea parents is genotyped.
- [ ] **A12 `CURATION.md`** — `file.info("CURATION.md")$mtime` is **newer** than
      both `*_curation_settings.rds` `run_at`. *Currently it is not*: re-run
      `code/curation_report.R`.

## Whole system

- [ ] **Clean re-run** — `assemble_B4I_phenotypes.R` → `BGLR_multi_trait_model.R`
      reproduces `BGLR_variance_components.csv` to three digits.
- [ ] **Shapes** — every artifact in EVALUATION.md §5 exists and has the stated
      shape.
- [ ] **Decide about the prior-only rows** — 16 of the 20 ungenotyped peas now
      have pedigrees and genotyped full sibs, so they need not stay unrelated to
      everything.

------------------------------------------------------------------------

*Notes / anomalies noticed (the "a little funny" observations — the point of
the exercise):*

-

-

-
