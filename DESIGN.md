# DESIGN

What this project is and how it is put together. For usage see
[README.md](README.md); for the reasoning behind the methods see
[BACKGROUND.md](BACKGROUND.md); for how accuracy is measured see
[CROSS_VALIDATION.md](CROSS_VALIDATION.md); for the simulation framework see
[SIMULATION.md](SIMULATION.md); for the open questions see
[docs/B4I_followups.md](docs/B4I_followups.md).

## Purpose

Analyse the B4I oat–pea intercrop experiment: pull the data from T3/Oat, work
out which accessions are genetically distinct, and fit models that separate what
a genotype does for itself from what it does to its partner.

Two model families, asking different questions of the same plots:

- a **bivariate DGE-IGE model** (BGLR) — producer and associate effects for both
  species in one analysis, with the covariance between them;
- an **extended MegaLMM factor model** — the oat–pea matrix read as
  genotype-by-environment, with each pea accession an environment for the oats.

## Shape

Two independent chains that meet at the phenotype table. Every step writes to
`output/` and reads what the previous one wrote, so any step can be re-run alone.

```
                    .Renviron  (T3_USERNAME / T3_PASSWORD, gitignored)
                         |
  data/Acc_B4I_Avena.txt |  data/Acc_B4I_Pisum.txt
           |             |             |
           v             v             v
  create_GRMs_T3.R ................... GRM_Avena.rds, GRM_Pisum.rds
           |
           |   find_trials_with_B4I_accessions.R
           |     trials evaluating >= 20 B4I accessions AND carrying oat
           |     grain yield on >= 20 of them  ->  34 trials
           |          |
           |          +-> B4I_trial_search.csv, B4I_trials_selected.csv
           |          +-> B4I_trait_availability.csv/.png   (34 x 38)
           |          +-> data/B4I_observations.rds/.csv.gz (139604 rows)
           |          +-> trial_cache/obs_*.rds, units_*.rds
           v
  curation_functions.R  (shared)
     |                    |
  curate_oat_accessions.R    curate_pea_accessions.R
     pedigrees, full-sib families, clonal pools,
     near-identical groups, analysis names
     |                    |
     +-> oat_analysis_names.csv    +-> pea_analysis_names.csv
     +-> *_curation_settings.rds   (the thresholds each run used)
                         |
            curation_report.R -> CURATION.md
                         |
                         v
            assemble_B4I_phenotypes.R
              plot table: oat acc x pea acc x both yields
              -> B4I_intercrop_pheno.rds/.csv
                         |
           +-------------+--------------+
           v                            v
  BGLR_multi_trait_model.R     megalmm_build_inputs.R
    bivariate DGE-IGE            oat x pea matrix, 3 fillings,
    -> variance components,      pea covariates
       Pr/As effects             -> megalmm_inputs_<fill>.rds
                                            |
                                 megalmm_setup.R (shared)
                                 megalmm_oat_pea.R
                                   CV sweep + final fit
                                   -> megalmm_cv_results.csv,
                                      megalmm_fit_<fill>.rds
```

### The other two chains

The simulation and the validation trial are separate from the pipeline above and
share only the GRMs and the fitted variance components.

```
  data/GRM_Avena.rds, GRM_Pisum.rds
           |
           v
  sim_config.R       SIM_LEVELS -> 120 data scenarios
                     the MegaLMM levers -> 960 candidates
                     AlgDesign::optFederov -> a D-optimal 150-run fraction
           |
           v
  sim_generate.R     one experiment under a known truth:
                     four genetic effects, two yields, a low-rank
                     interaction per trait, optional GxE
           |
           v
  sim_fit.R          split -> {additive, dge_ige} via fit_producer_associate,
                     MegaLMM in BOTH orientations -> two oat x pea surfaces
                     per model -> score against the simulated truth
           |
           v
  sim_run.R          the driver: caching, job-array slicing, --check,
                     the design model fitted to the outcomes
                     -> simulation_results.csv, simulation_design.csv,
                        simulation_design_effects.csv, simulation_summary.png
                     on a cluster: code/scinet/sim_array.sbatch


  output/BGLR_*_effects_all_seeds.csv, variance components
           |
           v
  validation_functions.R  (shared, sourced not run)
     |              |                |
  validate_      validate_        validate_      validate_
  pool_selection  power            design        crossval
     As+/As-      analytic and     field book    attenuation
     pools        simulated        + balance     (lambda),
     + churn      power            checks        leave-one-trial-out
```

Both chains are entirely offline: no T3 login, nothing downloaded. The
simulation additionally knows its own truth, which is what makes it checkable
against something other than another model.

## The scripts

| script | does |
|---|---|
| `create_GRMs_T3.R` | T3GenoTools GRMs for both species; skips protocols with negligible coverage |
| `find_trials_with_B4I_accessions.R` | trial discovery, phenotype download, trait availability matrix |
| `curation_functions.R` | shared curation: protocol check, dosages, correlations, families, analysis names |
| `curate_oat_accessions.R` | oat: pedigrees, parents, clonal families, renaming |
| `curate_pea_accessions.R` | pea: the same, minus the pedigree half (there are no pea pedigrees) |
| `curation_report.R` | writes `CURATION.md` from what the curation runs recorded |
| `assemble_B4I_phenotypes.R` | plot-level table for the models |
| `BGLR_multi_trait_model.R` | bivariate DGE-IGE fit |
| `megalmm_build_inputs.R` | oat x pea matrix and pea environmental covariates |
| `megalmm_setup.R` | MegaLMM model construction, sampling, posterior extraction |
| `megalmm_oat_pea.R` | cross-validation sweep and final MegaLMM fit |
| `sim_config.R` | the simulation factorial and its parameters, grounded in the B4I data |
| `sim_generate.R` | simulates one experiment under a known truth |
| `sim_fit.R` | fits both frameworks to a simulated experiment and scores them |
| `sim_run.R` | drives the simulation grid, with per-scenario caching |
| `sim_filter_results.R` | drops rows in `simulation_results.csv` that are not from the current grid. Needed because the CSV is rebuilt by globbing the cache |
| `sim_anova.R` | split-plot ANOVA of the MegaLMM results: which factors move accuracy, with effect sizes and expected mean squares |
| `sim_compare.R` | paired head-to-head of the two BGLR models against MegaLMM at its best settings |
| `sim_kron_study.R` | is the interaction result an artefact of `SIM_KRON_RANK`? Runs dge_ige at full rank, where no truncated basis can beat it |
| `cross_validate_combinations.R` | cross-validation over held-out oat x pea COMBINATIONS, both components kept in training |
| `sim_int_config.R` | the interaction-focused design: a full factorial, where the main grid is a fraction |
| `sim_int_run.R` | drives it. Shares the generator, fitters and scorer with `sim_run.R`; only the design and the output names differ |
| `sim_int_anova.R` | split-plot ANOVA of the interaction design, with `model` as the subplot factor: which factors decide whether MegaLMM or dge_ige wins |
| `interaction_decomp.R` | rewrites a fitted interaction surface as orthonormal scores x loadings, inside each MCMC draw. Sourced, not run |
| `sim_decomp_run.R` | drives the decomposition sweep: does the recovered structure match the simulated factors, or only the leading kinship directions? |
| `validation_functions.R` | shared machinery for the validation trial: pool construction, the field design, power. Sourced, not run |
| `validate_refresh.R` | the single entry point: drives the whole chain from T3 discovery to the field book when new trials land |
| `validate_pool_selection.R` | builds the As+ / As- pools and diffs them against the previous vintage |
| `validate_power.R` | analytic and simulation power, plus the false-positive check |
| `validate_design.R` | field book for the validation trial, with balance checks |
| `evaluation.R` | console tooling for [EVALUATION.md](EVALUATION.md): `arm_evaluation()`, `peek()`, `eval_load()`, the independent checks. Sourced by hand, never by a pipeline script |
| `evaluation_snippets.R` | the paste-along companion to [EVALUATION_SIMULATION.md](EVALUATION_SIMULATION.md). Not a script — blocks to copy, level by level |

## The tests

`tests/` holds the unit tests. `Rscript tests/run_all.R` runs the fast tier in
about 45 seconds; `--all` adds `test_fits.R`, which fits short real BGLR and
MegaLMM chains. Each file runs in its own process, so no file's seeds or loaded
functions can decide another's result.

Nothing in `tests/` reads `output/` or needs credentials — the oracles are
algebraic identities, planted answers and values the caller requested — so the
suite runs in a fresh clone. See [tests/README.md](tests/README.md) for what
each file pins.

## Conventions

- workflowr project: `code/` runnable scripts and shared functions, `analysis/`
  notebooks, `data/` inputs, `output/` generated results, `tests/` unit tests.
- Every script starts with `library(tidyverse)` and `here::i_am(...)`; other
  packages are called as `package::function()`.
- `output/` is gitignored apart from its README: everything in it regenerates.
- Two derived things live in `data/` instead, and are versioned: the GRMs and
  the downloaded observations. They are inputs to everything downstream, they
  take hours to rebuild from T3, and having them in the repository means a
  clone can run the models and the simulation without a T3 login.
- `output/trial_cache/`, `output/megalmm_runs/`, `output/simulation/` and
  `output/simulation_runs/` are caches and run state, regenerable and never
  committed.
- Credentials live only in `.Renviron`, which is gitignored.

## Caching

Three layers, because the expensive things differ in kind:

- **T3GenoTools' shared cache** (outside the repo, `T3GenoTools::geno_cache_root()`)
  holds VCFs, dosages and per-protocol GRMs. Rebuilt in hours, so it lives in the
  user data directory rather than a cache directory that the OS may purge.
- **`output/trial_cache/`** holds one file per trial of BrAPI observations and
  observation units. Re-running a trial search re-downloads nothing.
- **`output/megalmm_runs/`** is MegaLMM's own run state, which it needs on disk.

## Naming

Accessions that markers show to be one genotype are renamed rather than dropped,
so their phenotypes stay in the analysis under one name:

- `<seed>_<pollen>_no_cross` — a full-sib family that did not segregate;
- `<seed_parent>_self` — a line matching such a family and sharing its female,
  or a set of such families of one female that cannot be told apart;
- `<seed_parent>` — a line confirmed identical to its genotyped parent;
- for pea, the representative's own name, there being no pedigrees to reason from.

`oat_analysis_names.csv` and `pea_analysis_names.csv` carry the original name,
the analysis name and the reason. Both the phenotype table and the GRMs are put
through the same mapping, the GRM by averaging the rows and columns of a group.
