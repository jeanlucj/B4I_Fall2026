# DESIGN

What this project is and how it is put together. For usage see
[README.md](README.md); for the reasoning behind the methods see
[BACKGROUND.md](BACKGROUND.md); for how accuracy is measured see
[CROSS_VALIDATION.md](CROSS_VALIDATION.md); for the open questions see
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
           |          +-> B4I_observations.rds/.csv.gz      (139604 rows)
           |          +-> trial_cache/obs_*.rds, units_*.rds
           v
  curation_functions.R  (shared)
     |                    |
  curate_oat_accessions.R    curate_pea_accessions.R
     pedigrees, full-sib families, clonal pools,
     near-identical groups, analysis names
     |                    |
     +-> oat_analysis_names.csv    +-> pea_analysis_names.csv
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

## The scripts

| script | does |
|---|---|
| `create_GRMs_T3.R` | T3GenoTools GRMs for both species; skips protocols with negligible coverage |
| `find_trials_with_B4I_accessions.R` | trial discovery, phenotype download, trait availability matrix |
| `curation_functions.R` | shared curation: protocol check, dosages, correlations, families, analysis names |
| `curate_oat_accessions.R` | oat: pedigrees, parents, clonal families, renaming |
| `curate_pea_accessions.R` | pea: the same, minus the pedigree half (there are no pea pedigrees) |
| `assemble_B4I_phenotypes.R` | plot-level table for the models |
| `BGLR_multi_trait_model.R` | bivariate DGE-IGE fit |
| `megalmm_build_inputs.R` | oat x pea matrix and pea environmental covariates |
| `megalmm_setup.R` | MegaLMM model construction, sampling, posterior extraction |
| `megalmm_oat_pea.R` | cross-validation sweep and final MegaLMM fit |

## Conventions

- workflowr project: `code/` runnable scripts and shared functions, `analysis/`
  notebooks, `data/` inputs, `output/` generated results.
- Every script starts with `library(tidyverse)` and `here::i_am(...)`; other
  packages are called as `package::function()`.
- `output/` is gitignored apart from its README: everything in it regenerates.
- `output/trial_cache/` and `output/megalmm_runs/` are caches and run state,
  regenerable and never committed.
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
