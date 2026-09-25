# B4I Fall 2026 — oat–pea intercrop analysis

Analysis of the B4I oat–pea intercrop experiment: pull the data from
**T3/Oat** over BrAPI, work out which accessions are genetically distinct, and
fit models that separate what a genotype does for itself from what it does to
its partner.

- **How to run it** → this file.
- **What it is and how it fits together** → [DESIGN.md](DESIGN.md).
- **Why it works this way (theory + decisions)** → [BACKGROUND.md](BACKGROUND.md).
- **How accuracy is measured** → [CROSS_VALIDATION.md](CROSS_VALIDATION.md).
- **When would either framework work?** → [SIMULATION.md](SIMULATION.md).
- **Why sparsity hurts MegaLMM so much** → [docs/MegaLMM_sparsity_challenge.md](docs/MegaLMM_sparsity_challenge.md).
- **How the specific-combination term is built** → [docs/specific-combination_kronecker.md](docs/specific-combination_kronecker.md).
- **Using both analysis orientations** → [BOTH_ORIENTATIONS.md](BOTH_ORIENTATIONS.md).
- **Running the simulation on SciNet** → [code/scinet/README.md](code/scinet/README.md).
- **What was collapsed, and why** → [CURATION.md](CURATION.md).
- **Validating the producer/associate effects** → [VALIDATION_DESIGN.md](VALIDATION_DESIGN.md).
- **Checking that it does what it says** → [EVALUATION.md](EVALUATION.md)
  and [EVALUATION_SIMULATION.md](EVALUATION_SIMULATION.md), with
  [their](EVALUATION_CHECKLIST.md)
  [checklists](EVALUATION_SIMULATION_CHECKLIST.md).
- **Open questions** → [docs/B4I_followups.md](docs/B4I_followups.md).

This is a [workflowr](https://github.com/workflowr/workflowr) project: runnable
scripts and shared functions in `code/`, inputs in `data/`, generated results in
`output/` (gitignored — everything there regenerates).

## 1. Install dependencies

```r
remotes::install_github("TriticeaeToolbox/BrAPI.R")   # BrAPI wrapper
remotes::install_github("jeanlucj/T3BrapiHelpers")
remotes::install_github("jeanlucj/T3GenoTools")       # GRMs + shared geno cache
remotes::install_github("deruncie/MegaLMM")           # factor model
install.packages(c("tidyverse", "here", "BGLR", "patchwork", "workflowr"))
```

R ≥ 4.5. `httr` must be attached in any session that talks to the server —
BrAPI.R calls `httr::timeout()` unqualified — which every script here does.

## 2. Credentials — `.Renviron`

T3/Oat does not serve this data anonymously. Put your account in `.Renviron` at
the project root:

```
T3_USERNAME=yourname
T3_PASSWORD=yourpassword
```

`.Renviron` is gitignored; never commit it. Each script calls `readRenviron()`
on it directly, so it is picked up whatever directory R started in, including
inside `wflow_build()`. Missing credentials stop the run with an explicit
message rather than returning empty results.

## 3. Run the pipeline

The scripts are independent and each reads what the previous one wrote, so they
can be run one at a time and re-run in isolation. In order:

```bash
Rscript code/create_GRMs_T3.R                   # ~5 min first time, seconds after
Rscript code/find_trials_with_B4I_accessions.R  # ~10 min first time
Rscript code/curate_oat_accessions.R            # ~2 min
Rscript code/curate_pea_accessions.R            # ~2 min
Rscript code/curation_report.R                  # seconds; writes CURATION.md
Rscript code/assemble_B4I_phenotypes.R          # ~1 min
Rscript code/BGLR_multi_trait_model.R           # ~4 min
Rscript code/megalmm_build_inputs.R             # ~1 min
Rscript code/megalmm_oat_pea.R                  # ~20 min
```

The simulation framework is separate from the pipeline and answers a different
question — when either framework *would* work. See [SIMULATION.md](SIMULATION.md):

```bash
Rscript code/sim_run.R --check                  # sanity check; run this first
Rscript code/sim_run.R                          # the 60-scenario grid
```

First runs download VCFs and phenotypes; later runs hit caches (see
[DESIGN.md](DESIGN.md#caching)) and are much faster. Nothing needs deleting to
pick up a settings change except the specific cache a setting feeds.

### What each step produces

| step | key outputs |
|---|---|
| `create_GRMs_T3.R` | `GRM_Avena.rds` (508×508, Oat 3K), `GRM_Pisum.rds` (435×435, GenoPea 13K) |
| `find_trials_with_B4I_accessions.R` | 34 selected trials, `B4I_trait_availability.csv` (34×38), 139,604 observations |
| `curate_*_accessions.R` | `*_analysis_names.csv` — which accessions are really one genotype |
| `curation_report.R` | `CURATION.md` — thresholds used and every accession collapsed |
| `assemble_B4I_phenotypes.R` | `B4I_intercrop_pheno.rds` — 2,371 plots, both yields |
| `BGLR_multi_trait_model.R` | `BGLR_variance_components.csv`, per-accession producer/associate effects |
| `megalmm_build_inputs.R` | `megalmm_inputs_{raw,centered,standardized}.rds` |
| `megalmm_oat_pea.R` | `megalmm_cv_results.csv`, `megalmm_fit_centered.rds` |
| `sim_run.R` | `simulation_results.csv`, `simulation_summary.png` |

## 4. Settings worth knowing

Each script has a settings block at the top. The ones that change conclusions:

**`curate_oat_accessions.R` / `curate_pea_accessions.R`**
- `identity_threshold` (0.99) — above this, two accessions are one line.
- `family_mean_r_flag` (0.985) — a full-sib family averaging above this did not
  segregate.
- `pool_tolerance` (0.01) — how close between-family and within-family
  correlation must be before two clonal families of one female are pooled.

**`find_trials_with_B4I_accessions.R`**
- `min_b4i_accessions` (20) and `selection_trait` (oat grain yield) — the two
  selection stages. Set `selection_trait` to `NA` for the count-only rule.

**`BGLR_multi_trait_model.R`**
- `trials`, `study_years` — which trials enter the fit.
- `fit_mix_term` (FALSE) — the specific-combination term. Off because 91% of
  oat–pea combinations occur in a single plot; the script reports the
  replication either way.

**`megalmm_build_inputs.R` / `megalmm_oat_pea.R`**
- `fills` — the three cell values, all built so they can be compared.
- `min_obs_per_env` (3), `min_obs_per_oat` (2) — trimming. Below this MegaLMM's
  sampler produces `NaN` rather than failing cleanly.
- `eigen_variance` (0.80) — how much pea genetic variance the GRM eigenvectors
  should cover. Costs 58 eigenvectors here.
- `final_fill` (`centered`) — chosen from the cross-validation, not in advance.

## 5. Results so far

**The bivariate DGE-IGE model fits.** Producer and associate variances are all
solidly non-zero, and the residual covariance between the two yields on a plot
is clearly negative — competition, once trial and block are removed. Neither
producer–associate covariance excludes zero, so the negative correlations in the
point estimates are a hypothesis rather than a result.

**The MegaLMM framing does not work on this data yet.** Predicting held-out
oat×pea cells, the best MegaLMM configuration reaches r = 0.07 against 0.28 for
simply predicting each oat's own mean. The pea covariates do help — they roughly
double the model's accuracy — but not nearly enough. With 91% of combinations
appearing in a single plot there is no specific-combination signal to find, which
is the same conclusion the bivariate model reached from the other direction.

Both point at the same design change: replicate specific oat×pea pairs rather
than spreading singletons thinly. See [BACKGROUND.md](BACKGROUND.md) for the
reasoning and [docs/B4I_followups.md](docs/B4I_followups.md) for what else is open.

## 6. Data notes

Two things about the source data that will bite anyone re-deriving this:

- **Trials IA and ND return each plot four times.** The copies differ only in
  their within-plot grid coordinates and carry no observations. Deduplicate on
  `observationUnitDbId` or those two trials quadruple.
- **A plot's pea partner is not a germplasm record.** It lives on the
  observation unit as `additionalInfo$intercropGermplasm`, so observation units
  must be fetched separately from observations and joined by plot.

Also: BrAPI's default `pageSize` is 10 records. Raising it takes one trial's
observations from ~165 s to ~4 s.
