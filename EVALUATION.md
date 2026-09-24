# EVALUATION — stepping through the analysis

This is the document for **checking that the analysis pipeline does what it
says**, from the RStudio console, one function at a time. [README.md](README.md)
is for running it; [DESIGN.md](DESIGN.md) is the static architecture;
[BACKGROUND.md](BACKGROUND.md) is why the methods are what they are;
[CURATION.md](CURATION.md) is what the curation actually did. For the
simulation see [EVALUATION_SIMULATION.md](EVALUATION_SIMULATION.md). The tick
list that goes with this file is
[EVALUATION_CHECKLIST.md](EVALUATION_CHECKLIST.md).

The code was written fast. That is not an argument that it is wrong, but it
does mean nothing here has been audited by anybody, and the failures worth
worrying about are the silent ones: a name that did not join, a matrix that
got transposed, a mask that leaked. Those produce a plausible number, not an
error. **The goal of this document is that working through it is faster and
more reliable than re-reading the code**, and that it puts a hand on every
place where a wrong answer would still look right.

The two results to keep in view while doing it:

- the bivariate DGE-IGE model says producer and associate variances are
  non-zero for both species and the residual oat–pea covariance on a plot is
  negative;
- the MegaLMM framing does **not** work on this data: r = 0.07 against 0.28
  for predicting each oat's own mean.

A negative result needs the same scrutiny as a positive one, and more of a
particular kind: the failure mode that makes a working method look broken is a
wiring mistake, and it looks exactly like "the method does not suit this data".

---

## 1. The modules, in one screen

Eleven files in `code/`, four of which are shared function libraries and the
rest runnable drivers. The **group** column is the name you pass to
`arm_evaluation()`.

| File (`code/…`) | Owns | Off/online | Group |
|---|---|---|---|
| `dge_ige_functions.R` | GRM reading, collapsing, factoring; incidence matrices; pulling covariances and effects out of a BGLR fit | offline | `grm`, `bglr` |
| `curation_functions.R` | connect, dosages, marker correlations, identity groups, full-sib families, analysis names | mixed | `curation_groups`, `curation_family`, `geno` |
| `t3_functions.R` | BrAPI download of observations and observation units; GRM construction wrapper | **online** | `t3` |
| `megalmm_setup.R` | MegaLMM state construction, burn-in/sampling, posterior extraction | offline | `megalmm` |
| `create_GRMs_T3.R` | driver: both species' GRMs | **online** | `geno` |
| `find_trials_with_B4I_accessions.R` | driver: trial discovery, phenotype download, trait availability | **online** | `trials` |
| `curate_oat_accessions.R` | driver: oat pedigrees, families, clonal pools, renaming | **online** | `curation_*` |
| `curate_pea_accessions.R` | driver: the same for pea | **online** | `curation_*` |
| `curation_report.R` | driver: writes `CURATION.md` from the recorded settings | offline | — |
| `assemble_B4I_phenotypes.R` | driver: the plot-level table | mostly cached | `t3` |
| `megalmm_build_inputs.R` | driver + functions: the oat × pea matrix and pea covariates | offline | `matrix` |
| `megalmm_oat_pea.R` | driver + functions: CV sweep and final fit | offline | `cv`, `megalmm` |
| `BGLR_multi_trait_model.R` | driver: the bivariate DGE-IGE fit | offline | `bglr` |

**Most of the pipeline can be evaluated offline.** `data/` carries the two
GRMs and all 139,604 downloaded observations, and `output/` carries the
curation results and the fitted models, so levels A1–A8 need no T3 login at
all. Only A9–A11 do.

---

## 2. How to evaluate: the plan

Three principles, in order:

1. **Offline before online.** A1–A8 read `data/` and `output/`. They are where
   the subtle bugs are — name joins, collapses, masks, transposes — and they
   cost nothing.
2. **Fast before slow**, which is the order the levels are in.
3. **One group armed at a time.** Arming everything makes the debugger
   unusable.

### The five tools (all in `code/evaluation.R`)

| Tool | Use it to | When |
|---|---|---|
| `eval_load("analysis")` | source the four shared function files, then pull each driver's own functions **without running the driver** | once per session |
| `eval_defs("code/x.R")` | the same for one script: evaluates only top-level function definitions and literal settings, reports how many driver expressions it skipped | when you want one script's internals |
| `arm_evaluation(group)` / `disarm_evaluation()` | `debug()` one module and watch the Global Environment fill line by line | stepping through logic |
| `peek(x)` | one-line health summary of the object between two steps: shape, NA count, fill %, degeneracy, rowname overlap. Returns its input, so it sits in a pipe | inspecting an intermediate |
| `check_grm()` / `check_matrix()` / `check_alignment()` / `eval_conflicts()` | re-derive a property independently instead of asking the pipeline whether it is happy | confirming a real result is real |

Debugger keys, once: **`n`** next line, **`s`** step into a call, **`c`** run
to this function's return, **`Q`** quit. `eval_groups("analysis")` prints the
menu and marks with `*` anything not yet loaded.

### What `eval_defs()` actually does

Every driver in `code/` defines its functions and then immediately runs.
Sourcing `megalmm_oat_pea.R` to get at `mask_cells()` would launch an
eighteen-fit cross-validation. `eval_defs()` reads the file, walks its
**top-level expressions**, and evaluates only two kinds:

1. **a function definition** — `f <- function(...) { ... }`;
2. **a static value** — a literal, or something built only from pure
   constructors (`c`, `list`, `here::here`, `file.path`, `paste0`, …), which is
   what a settings block is made of.

Everything else is skipped: `conn <- connect_t3(db_name)`,
`pheno <- readRDS(pheno_file)`, `fits <- purrr::map(seeds, fit_one)`, the
pipeline itself. A whitelisted expression that still refers to something the
driver would have built by then (`Y0 <- inputs[[1]]$Y`) is attempted, fails
harmlessly, and is counted as skipped.

**The rule looks only at the shape of the top-level statement, never inside a
function body.** This is the part worth being clear about: *defining* a
function never executes anything in it. `cv_once()` fits a MegaLMM model every
time it is called, and `fit_one()` runs a 20,000-iteration BGLR chain — both
are recreated by `eval_defs()` exactly like any other function, in full, with
their fits intact. They simply are not called. So there is no such thing as a
function too expensive to define; the cost lives in calling it, which is yours
to decide.

What you get is every function in the file plus its settings. What you do
**not** get is anything the driver computed, which matters when a function
reads one. `build_one()` in `megalmm_build_inputs.R` uses the driver's `pheno`,
`G_oat_full` and `blues`; after `eval_defs()` those do not exist, so calling
`build_one("centered")` errors until you build them yourself — which the A5
walkthrough does, in three lines. The error is immediate and obvious, not
silent.

### What `eval_conflicts()` is for

`debug()` attaches to a **function object**, not to a name. If two files define
the same name, `arm_evaluation()` debugs whichever copy was loaded last — and
if the pipeline calls the other one, you will step carefully through code that
is not running. `eval_conflicts()` lists every name that more than one loaded
file defines, splitting them into:

- **functions** — should be empty, and is. A duplicated function is a refactor
  waiting to happen, not a setting;
- **settings** — expected. Each driver has its own settings block, so
  `db_name`, `out_dir`, `cache_dir`, `fills`, `refresh`, `acc_files` and
  `monoculture_labels` appear in two or three of them. They are only in one
  session together because `eval_load()` put them there; no pipeline run ever
  sees two at once.

Until 2026-09-23 this repo had three duplicated *functions* — `connect_t3`,
`collapse_grm` and `read_grm` — the last two with different bodies, so the
BGLR model and the MegaLMM model collapsed the GRM through different code.
They have since been unified: see §6.1.

---

## 3. Bootstrap (paste once per session)

```r
library(tidyverse)
here::i_am("code/evaluation.R")            # or any file in the project
source(here::here("code", "evaluation.R"))

eval_load("analysis")                      # 4 sourced files + 4 drivers' functions
eval_conflicts()                           # names two files both define -- read this
eval_groups("analysis")                    # the menu
```

`eval_load()` prints what each file contributed. Expect:

```
megalmm_build_inputs.R: 9 function(s), 12 setting(s); skipped 8 driver expression(s)
megalmm_oat_pea.R:      3 function(s), 15 setting(s); skipped 39 driver expression(s)
find_trials_…:          9 function(s), 11 setting(s); skipped 41 driver expression(s)
create_GRMs_T3.R:       2 function(s),  6 setting(s); skipped 15 driver expression(s)
```

`eval_conflicts()` should answer **"no FUNCTION is defined by more than one
loaded file"**, then list seven settings names that two or three drivers each
define in their own settings block — `acc_files`, `cache_dir`, `db_name`,
`fills`, `monoculture_labels`, `out_dir`, `refresh`. That is the convention and
is fine. A *function* appearing in that first list is a finding.

For the online levels (A9–A11) open one connection and reuse it:

```r
conn <- connect_t3("T3/Oat")               # reads .Renviron itself
cfg  <- T3GenoTools::geno_config("T3/Oat", impute = "mean", progress = FALSE)
```

---

## 4. Level-by-level walkthrough

Each level: **arm** the group, run the console lines, and check the four
annotations — **→ returns** (the shape you should see), **🔍 eyeball** (what
healthy looks like, with this data's actual numbers), **🚩 red flag** (the
silent failure to hunt for), and the *Assumption:* the code is making, which
is the thing to test rather than the output. `disarm_evaluation()` at the end
of each level.

Numbers quoted are what the committed `data/` and the current `output/`
produce as of 2026-09-23. If yours differ, that is the finding.

---

### A1 — `grm` (offline, seconds)

Everything downstream joins on accession names, and this is where names get
rewritten.

```r
arm_evaluation("grm")
G_oat_raw <- read_grm(here::here("data", "GRM_Avena.rds"))
G_pea_raw <- read_grm(here::here("data", "GRM_Pisum.rds"))
G_oat <- collapse_grm(G_oat_raw, here::here("output", "oat_analysis_names.csv"))
G_pea <- collapse_grm(G_pea_raw, here::here("output", "pea_analysis_names.csv"))
L_oat <- grm_factor(G_oat)
disarm_evaluation()

check_grm(G_oat); check_grm(G_pea)
max(abs(tcrossprod(L_oat) - G_oat))          # must be < 1e-8
```

- **→ returns** `G_oat_raw` 508 × 508, `G_pea_raw` 435 × 435; after collapsing,
  **466 × 466** and **434 × 434**. `L_oat` is 466 × (466 − dropped).
- **🔍 eyeball** the oat collapse is 508 → 466: 50 accessions folded into 8
  names, exactly what [CURATION.md](CURATION.md) reports. Pea is 435 → 434, the
  single `NDP170084G` → `ND VICTORY` merge. `check_grm()` should say symmetric,
  min eigenvalue 0.0015 (oat) / 0.0100 (pea), both positive.
- **🚩 red flag** `check_grm()` reports **8 prior-only lines in oat and 20 in
  pea** — rows with diagonal 1 and every off-diagonal exactly 0, the
  ungenotyped accessions that `build_grm()` injected. They are in the GRM and
  therefore in the models, related to nothing. Their producer and associate
  effects rest on their own plots alone. That is defensible but it is not what
  "genomic prediction borrowed strength" means, and nothing downstream
  distinguishes them.
- **🔍 eyeball** there is now one `collapse_grm()` for the whole project, in
  `dge_ige_functions.R`, sourced by both `BGLR_multi_trait_model.R` and
  `megalmm_build_inputs.R`. So this level's numbers are the numbers both models
  use; until 2026-09-23 they came from two different copies of the function.
- **🚩 red flag** step into `collapse_grm()` and watch `A <- sweep(A, 2, colSums(A), "/")`.
  If a group's members are not contiguous in `rownames(G)` the `factor(…,
  levels = unique(new_name))` still handles it — but confirm the output row
  order is the order you expect, because everything after this indexes by name
  *and* the models index by position.
- *Assumption:* averaging the GRM rows and columns of a group equals the GRM of
  their averaged marker profiles. True when the lines really are identical,
  which is what the 0.99 threshold is asserting.

---

### A2 — `curation_groups` (offline, instant — with a planted duplicate)

The grouping logic can be tested without T3 at all, by handing it a similarity
matrix you built yourself with a known answer in it. This is the single
highest-value check in the document: it is the one place where you can compare
the code against ground truth rather than against its own output.

```r
arm_evaluation("curation_groups")

set.seed(1)
n <- 12
S <- diag(n); rownames(S) <- colnames(S) <- paste0("acc", 1:n)
S[lower.tri(S)] <- runif(sum(lower.tri(S)), 0.5, 0.8)
S[upper.tri(S)] <- t(S)[upper.tri(S)]
S["acc3", "acc7"] <- S["acc7", "acc3"] <- 0.995        # a planted duplicate
S["acc7", "acc9"] <- S["acc9", "acc7"] <- 0.993        # chained to a third

grp <- identity_groups(S, 0.99)                        # step through cutree
ped <- tibble::tibble(germplasmName = rownames(S), pedigree = NA_character_,
                      seed_parent = NA_character_, pollen_parent = NA_character_)
cls <- classify_groups(grp, ped, 0.5)
rep <- choose_representatives(cls, NULL)
disarm_evaluation()
grp; dplyr::select(rep, germplasmName, group, action, representative)
```

- **→ returns** `grp` has **three** rows — `acc3`, `acc7`, `acc9` — all in
  group 1. `rep` marks one `keep` and two `drop`.
- **🔍 eyeball** `acc3`–`acc9` were never made similar to each other
  (`S["acc3","acc9"]` is 0.504), yet they are one group. That is **single
  linkage, and it is deliberate** — a chain of near-identical accessions is one
  line however it is labelled. Confirm you agree with that on the real data
  too: the largest real oat group has 12 members.
- **🚩 red flag** single linkage is also how two genuinely different lines get
  merged through a bridging accession. On the real oat groups, check the
  *minimum* pairwise r inside each group, not the maximum:
  ```r
  sim <- NULL   # requires the real dosage matrix -- see A11
  ```
  Without T3, use the recorded groups as a proxy: 19 groups covering 90
  accessions, largest 12. A group whose size is far larger than its pedigree
  structure would explain is the signature.
- **🚩 red flag** `choose_representatives()` breaks ties by phenotype count,
  then alphabetically. When `observations_file` is missing every count is 0 and
  the representative is simply the **alphabetically first** name — silently, no
  message. Pass the real file and confirm the representative changes (or
  confirm it does not, and know why).
- *Assumption:* `cutree(hclust(as.dist(1 - S), "single"), h = 1 - threshold)`
  cuts at exactly the stated correlation. Worth stepping once: the height is a
  *distance*, and the sign convention is easy to invert.

---

### A3 — `curation_family` (offline, instant — with a planted family)

The pedigree half. Same approach: build a family whose answer you know.

```r
arm_evaluation("curation_family")

nm <- paste0("k", 1:6)
S2 <- matrix(0.75, 6, 6); diag(S2) <- 1
dimnames(S2) <- list(nm, nm)
S2[1:3, 1:3] <- 0.995; diag(S2) <- 1                   # k1-k3 did not segregate
ped2 <- tibble::tibble(
  germplasmName = nm,
  seed_parent   = c("A", "A", "A", "B", "B", "B"),
  pollen_parent = c("X", "X", "X", "Y", "Y", "Y"),
  pedigree      = paste0(c("A","A","A","B","B","B"), "/", c("X","X","X","Y","Y","Y")),
  both_parents  = TRUE)

fams <- pedigree_families(ped2, nm, 2)
fp   <- family_pair_correlations(fams, S2)
fs   <- summarise_families(fp, 0.985)
an   <- resolve_analysis_names(ped2, fs, S2, 0.99, 0.985, 0.01)
disarm_evaluation()
fs; an
```

- **→ returns** `fs` has two rows; `A/X` has `mean_r` ≈ 0.995 and
  `clonal_family = TRUE`, `B/Y` ≈ 0.75 and `FALSE`. `an` renames k1–k3 to
  `A_X_no_cross` and leaves k4–k6 alone.
- **🔍 eyeball** on the **real** data: oat has 9 of 74 families above 0.985
  averaging 0.992 against a segregating baseline of 0.747 — a clean separation,
  not a threshold sitting in the middle of a continuum. Pea now has 36 families
  and **none** above 0.985; the highest is 0.898. Plot both distributions and
  satisfy yourself the oat cut is in a gap:
  ```r
  readr::read_csv(here::here("output", "oat_family_correlations.csv")) |>
    ggplot2::ggplot(ggplot2::aes(mean_r)) + ggplot2::geom_histogram(bins = 40) +
    ggplot2::geom_vline(xintercept = 0.985, colour = "red")
  ```
- **🚩 red flag** `resolve_analysis_names()` applies three rules in sequence and
  **rule 3 overrides the first two**: a line confirmed identical to its
  genotyped seed parent is renamed to the parent, discarding the `_no_cross` /
  `_self` label. If the parent is not itself in the accession list, this
  invents a GRM row name that exists nowhere else. Check that every
  `analysis_name` is either an existing accession or a constructed
  `*_no_cross` / `*_self` label:
  ```r
  an_oat <- readr::read_csv(here::here("output", "oat_analysis_names.csv"))
  setdiff(an_oat$analysis_name[!grepl("_no_cross$|_self$", an_oat$analysis_name)],
          rownames(G_oat_raw))                      # expect character(0)
  ```
- **🚩 red flag** `pool_clonal_families()` merges two clonal families of one
  female when `between_r >= within_r - tolerance`. With `tolerance = 0.01` and
  within-family r around 0.985, that fires at between-family r ≈ 0.975 —
  which is *below* the 0.99 identity threshold used everywhere else. Two
  different thresholds are deciding "same line" in the same script. Decide
  whether you are comfortable with that; the oat `IL18-8562_self` group of 12
  rests on it (between 0.9827 vs within 0.985).
- *Assumption:* the pedigree string is `seed/pollen` in that order.
  `verify_parent_order()` checks 25 of them against the structured BrAPI
  endpoint and is only run in the oat script (`n_verify <- 25`). The pea script
  never calls it, and pea pedigrees are new. **Run it for pea** (A11).

---

### A4 — the plot table (offline if `output/trial_cache/` is populated)

`assemble_B4I_phenotypes.R` is a straight-line driver with no functions worth
arming. Check the artifact instead, by re-deriving its counts from the cache it
was built from.

```r
pheno <- readRDS(here::here("output", "B4I_intercrop_pheno.rds"))
peek(pheno)

units <- list.files(here::here("output", "trial_cache"), "^units_",
                    full.names = TRUE) |> purrr::map(readRDS) |> purrr::list_rbind()
dplyr::count(units, studyName)                        # independent re-derivation
dplyr::count(pheno, studyName)
```

- **→ returns** `pheno` is **2377 × 13**; `units` is 2400 rows, 400 per trial
  for six trials.
- **🔍 eyeball** 2400 plots in, 2377 out. The 23 lost are plots with neither
  yield recorded, plus monoculture plots. Account for them explicitly rather
  than accepting the number:
  ```r
  dplyr::anti_join(units, pheno, by = "observationUnitDbId") |>
    dplyr::count(studyName, germplasmName, intercropGermplasmName)
  ```
  2371 rows have both yields; 6 have pea yield but no oat yield; 0 the other
  way round.
- **🔍 eyeball** the renaming actually happened: `sum(pheno$oat_renamed)` and
  `sum(pheno$pea_renamed)` should be non-zero and should match the number of
  plots whose accession appears in the analysis-names files.
- **🚩 red flag** README documents that trials IA and ND return each plot four
  times. The dedup is `dplyr::distinct(observationUnitDbId, .keep_all = TRUE)`
  inside `fetch_observation_units()`, i.e. in the **cache-writing** path. A
  cache file written before that line existed would still hold quadruplicated
  rows and nothing downstream would notice. Confirm directly:
  `nrow(units) == dplyr::n_distinct(units$observationUnitDbId)`.
- **🚩 red flag** `yields |> distinct(observationUnitDbId, trait, .keep_all = TRUE)`
  silently keeps the *first* of any duplicated observation. Count how many it
  dropped before trusting it:
  ```r
  obs <- readRDS(here::here("data", "B4I_observations.rds"))
  obs |> dplyr::filter(observationVariableName == "Grain yield - g/m2|CO_350:0000260",
                       !is.na(value), value != "") |>
    dplyr::count(observationUnitDbId) |> dplyr::count(n)
  ```
- *Assumption:* `studyYear` is parsed out of the trial name with
  `str_extract(studyName, "(?<=_)\\d{4}(?=_)")`. It works for `B4I_2025_AL`. A
  trial named otherwise gets `NA` and then silently fails the
  `studyYear %in% study_years` filter in the model script.

---

### A5 — `matrix` (offline)

The oat × pea matrix and the pea covariates. This is where a transpose or a
misaligned covariate would do the most damage, because MegaLMM would still fit.

```r
arm_evaluation("matrix")
cells <- scale_within_trial(pheno, "centered")         # step: watch the group_by
bm    <- build_matrix(pheno, "centered")
Y     <- trim_matrix(bm$Y, 3, 2)
disarm_evaluation()

peek(bm$Y); peek(Y); check_matrix(Y, 3, 2)
inp <- readRDS(here::here("output", "megalmm_inputs_centered.rds"))
identical(dim(Y), dim(inp$Y))
```

- **→ returns** `bm$Y` is 442 × 423 before trimming; `Y` is **237 × 222** with
  **1682 filled cells (3.20%)**, min 2 per row and 3 per column.
- **🔍 eyeball** rows are oats and columns are peas. Verify rather than assume,
  because the whole framing depends on it:
  ```r
  all(rownames(Y) %in% pheno$germplasmName)             # TRUE
  all(colnames(Y) %in% pheno$intercropGermplasmName)    # TRUE
  ```
- **🔍 eyeball** centering is **within trial, before averaging**. Check one
  trial sums to zero: `cells |> group_by(studyName) |> summarise(m = mean(cell_value))`
  → all ~0 for `centered`, all ~0 with sd 1 for `standardized`, untouched for `raw`.
- **🚩 red flag** `pea_trait_blues()` globs **every** `units_*.rds` in
  `output/trial_cache/`, not the six intercrop trials. Today that is the same
  set (6 files), so the BLUEs are correct. Cache observation units for any
  other trial and the pea BLUEs silently widen to include it. Check:
  `length(list.files(here::here("output","trial_cache"), "^units_"))` is 6.
- **🚩 red flag** `X_Env` is 222 × **61** = 1 intercept + 2 pea phenotypes + 58
  GRM eigenvectors. The comment in `megalmm_build_inputs.R` says 80% of the pea
  variance "costs about 99 eigenvectors"; the real figure is 58 (README has it
  right). A stale comment is harmless; a stale *number in your head* about how
  much of the pea genome is represented is not.
- **🚩 red flag** `scale()` on the two phenotypic covariates returns `NaN` for a
  constant column. `stopifnot(!anyNA(inp$X_Env))` — currently 0 NA.
- **🔍 eyeball** the group vector is what stops 58 eigenvectors swamping 2
  phenotypes under the ARD prior:
  `table(inp$X_Env_groups)` → `1 / 2 / 58`.
- *Assumption:* `stopifnot(identical(rownames(X_Env), colnames(Y)))` in the
  script is the only thing keeping the covariates attached to the right pea.
  It is there; confirm it fires by breaking it once on a copy.

---

### A6 — `megalmm` (offline; needs the MegaLMM package)

```r
arm_evaluation("megalmm")
inp <- readRDS(here::here("output", "megalmm_inputs_centered.rds"))
accNames <- tibble::tibble(germplasmName = rownames(inp$Y))
state <- setup_megalmm_state(
  accNames = accNames, wideData = inp$Y, kinMat = inp$G_oat,
  envCov = list(X_Env = inp$X_Env, X_Env_groups = inp$X_Env_groups,
                newEnv = character(0)),
  runID = file.path(tempdir(), "eval_mm"), K = 10)
disarm_evaluation()
```

- **→ returns** a `MegaLMM_state`. The message `missing-data map 1 of N` should
  appear — **map 1, the most grouped**, which is the opposite of the choice
  that suits a dense matrix and is correct here.
- **🔍 eyeball** step to the `Lambda_prior` block and confirm `fit_X = FALSE` at
  construction. The covariates are switched on partway through burn-in
  (`fit_X_from`), not at the start; if they were on from the start the factors
  would never settle first.
- **🚩 red flag** `.undecorate_rows()` checks that MegaLMM's `U` rows match the
  accessions passed in — but **`Eta_mean` is returned unchecked**, and
  `Eta_mean` is what the cross-validation scores. Verify it yourself:
  ```r
  post <- readRDS(here::here("output", "megalmm_fit_centered.rds"))$posterior
  identical(rownames(post$Eta_mean), rownames(inp$Y))   # TRUE today
  identical(colnames(post$Eta_mean), colnames(inp$Y))   # TRUE today
  ```
  Both hold. They are not guaranteed by anything in the code.
- **🚩 red flag** `fixed_main_effect = TRUE` (used by the simulation, not by
  this script) changes `scale_Y` to FALSE and loosens `tot_F_var` for factor 1
  only. If you turn it on here, Y must be scaled globally by the caller first,
  or "loadings fixed at 1" means something different in every column. The test
  is `sd(post$Lambda[1, ]) < 1e-8`, not `all(post$Lambda[1, ] == 1)`.
- *Assumption:* `K = 10` is generous and the ARD prior shrinks the surplus.
  Test it on the saved fit: `round(rowSums(post$Lambda^2), 3)` — the tail
  should be near zero. If the smallest factor is still substantial, K was too
  small and the comparison understates MegaLMM.

---

### A7 — `cv` (offline; the crux of the negative result)

If the masking is wrong, the headline conclusion is wrong.

```r
arm_evaluation("cv")
Y <- readRDS(here::here("output", "megalmm_inputs_centered.rds"))$Y
sp <- mask_cells(Y, 0.20, 20260921, floor_obs = 2)
bl <- margin_baselines(sp$Y_train, sp$held)
disarm_evaluation()

length(sp$held); sp$achieved
all(rowSums(!is.na(sp$Y_train)) > 0); all(colSums(!is.na(sp$Y_train)) > 0)
all(is.na(sp$Y_train[sp$held]))                        # masked cells really gone
identical(sp$truth, Y[sp$held])                        # truth taken before masking
```

- **→ returns** 336 held-out cells of 1682 (the CSV records `n_held = 336` for
  every rep), `achieved` ≈ 1.0 of the 20% asked for.
- **🔍 eyeball** the floor works: `min(rowSums(!is.na(sp$Y_train)))` ≥ 2 and
  `min(colSums(!is.na(sp$Y_train)))` ≥ 2. Without it the ARD sampler returns
  `NaN` from several frames deep, not an error.
- **🚩 red flag — the leak test.** The baselines must be *identical* between the
  covariates-on and covariates-off arms, because the mask is seeded and the
  baselines ignore covariates. Check it in the recorded results:
  ```r
  cv <- readr::read_csv(here::here("output", "megalmm_cv_results.csv"))
  cv |> dplyr::filter(!stringr::str_starts(predictor, "MegaLMM")) |>
    dplyr::group_by(fill, seed, predictor) |>
    dplyr::summarise(n_distinct_r = dplyr::n_distinct(r), .groups = "drop") |>
    dplyr::count(n_distinct_r)                         # must be all 1
  ```
  It is 1 throughout. If it were not, the two arms saw different data and the
  covariate comparison is meaningless.
- **🚩 red flag — the direction of the result.** The baselines beat MegaLMM by
  4×, so the interesting failure is one that *handicaps the model*, not one
  that flatters it. Two things to satisfy yourself about:
  1. `MegaLMM_Eta` is scored against `Y[held]` on the same scale as the
     baselines — yes, both are cell values of the same matrix.
  2. `U` and `U_noU_R` are genetic values and exclude the per-column intercept,
     so they are structurally at a disadvantage *across columns*. And yet the
     **best MegaLMM number (0.072) is `U_noU_R`, not `Eta_mean` (0.031)**. That
     is backwards from what the code comment predicts and is worth explaining
     before the negative result is written up.
- **🔍 eyeball** the covariates help `U`/`U_noU_R` (0.018 → 0.072) but do
  nothing at all for `Eta_mean` (0.0308 → 0.0310). README's "the pea covariates
  roughly double the model's accuracy" is true of one target and not the other.
- *Assumption:* masking filled cells at random within a fitted matrix is CV2
  (a new oat × pea cell, both seen before), **not** CV0 (a new pea). `newEnv`
  is `character(0)` throughout, so no new-environment prediction is being
  tested here at all. See [CROSS_VALIDATION.md](CROSS_VALIDATION.md).

---

### A8 — `bglr` (offline; the result that *is* positive)

```r
fit <- readRDS(here::here("output", "BGLR_multitrait_pea_oat_yield_Gpea_Goat_Gmix.rds"))

# L must come from the GRM SUBSET to the accessions in the model, in the same
# sorted order the script used -- not the full L_oat from A1. Getting this
# wrong is instructive: the dimensions simply will not conform.
mod <- pheno |>
  dplyr::filter(!is.na(oat_yield), !is.na(pea_yield),
                studyYear %in% c(2025L, 2026L),
                studyName %in% c("B4I_2025_AL", "B4I_2025_IA", "B4I_2025_IL",
                                 "B4I_2025_ND", "B4I_2025_NY", "B4I_2026_IL"))
oatAccs   <- sort(unique(mod$germplasmName))
L_oat_sub <- grm_factor(G_oat[oatAccs, oatAccs, drop = FALSE])

arm_evaluation("bglr")
vc <- covariance_components(fit, "G_oat", c("peaYield", "oatYield"),
                            c(Pr = "oatYield", As = "peaYield"))
ef <- accession_effects(fit, "G_oat", L_oat_sub, c("peaYield", "oatYield"),
                        c(Pr = "oatYield", As = "peaYield"))
disarm_evaluation()
vc; peek(setNames(ef$PrEff, ef$accession))
```

- **→ returns** `vc` is one row: `var_Pr` ≈ 394, `var_As` ≈ 105,
  `cov_PrAs` ≈ −13, `cor_PrAs` ≈ −0.061. For `G_pea`: 186 / 299 / −55 / −0.233.
- **🚩 red flag — the transpose that inverts the paper.** The response is built
  as `Y <- as.matrix(grainWgt[, c("peaYield", "oatYield")])` — **pea first**.
  `effect_roles` then says the oat kernel's producer effect is the one on
  `oatYield`. Swap the two columns and every producer becomes an associate,
  silently, with no dimension error. This is the single highest-consequence
  ordering in the project. Verify it against the data, not against the label:
  ```r
  colnames(fit$ETAHat)                           # "peaYield" "oatYield"
  round(colMeans(fit$ETAHat), 1)                 # 74.9  167.4
  round(c(mean(mod$pea_yield), mean(mod$oat_yield)), 1)   # 74.9  167.3
  round(cor(fit$ETAHat[, 1], mod$pea_yield), 3)  # 0.908 -- column 1 IS pea
  round(cor(fit$ETAHat[, 2], mod$oat_yield), 3)  # 0.967
  ```
  (`mod` is the filtered table built above, and its row order must match the
  fit's — same filters, no re-sorting.) Then confirm `var_Pr` for `G_oat`
  (394) is much larger than its `var_As` (105): an oat's effect on its own
  yield should exceed its effect on the pea's, and that asymmetry is the sanity
  check the labels have to reproduce. `fit$y` is **not** saved, so the column
  names on `ETAHat` are the only label the saved object carries.
- **🔍 eyeball** four chains agree to three digits on every variance component,
  and the between-chain Spearman correlation of GMA is high. Chain agreement is
  necessary, not sufficient: four chains with the same bug agree perfectly.
- **🚩 red flag** the design diagnostics warn about **99 oat and 84 pea
  accessions with a single partner**, out of 442 and 423. Their producer and
  associate effects are aliased and rest entirely on genomic covariance with
  better-connected relatives. Combined with the prior-only GRM rows from A1,
  ask how many accessions have *both* problems — a prior-only line with one
  partner has nothing at all behind its effect:
  ```r
  # prior-only names from check_grm(G_oat)$prior_only[[1]], partners from the pheno table
  ```
- **🚩 red flag** `G_mix` (2059 × 2059) and its eigendecomposition
  `L_mix <- grm_factor(G_mix)` plus a `stopifnot` are computed **even when
  `fit_mix_term` is FALSE**, which it is. Minutes of work and a large
  allocation for a term that is not fitted. Not wrong, but if the script is
  ever slow or memory-hungry, that is why.
- **🔍 eyeball** only the two IL trials have more than one block, so 10 of 14
  block levels are fitted and four single-block trials are aliased with their
  trial effect and correctly dropped. `blockNumber` is never `NA` here; if it
  ever were, `paste(studyYear, studyName, blockNumber)` would make `"… NA"` a
  real level rather than a missing one.
- *Assumption:* a BRR on `Z L` is the RKHS model with kernel `Z G Z'`, with
  coefficients on the accession scale. The `stopifnot(max(abs(tcrossprod(L) - G)) < 1e-8)`
  in the script is what backs it.

---

> **Everything below is online.** Open the connection from §3 first.

### A9 — `t3` (live BrAPI; start on a cached trial)

```r
arm_evaluation("t3")
u <- fetch_observation_units(conn, "6821", here::here("output", "trial_cache"))
disarm_evaluation()
peek(u)
```

Trial `6821` is `B4I_2025_AL` and is already in `output/trial_cache/`, so this
touches no network on the happy path — the cheapest possible first live check.
The six intercrop trials are 6821 (AL), 6822 (IA), 6881 (IL), 6823 (ND),
6824 (NY) and 6997 (IL 2026).

- **→ returns** a tibble with `intercropGermplasmName`, `repNumber`,
  `blockNumber` populated.
- **🔍 eyeball** `intercropGermplasmName` is **not** all `NA`. It lives in
  `additionalInfo$intercropGermplasm`, not on the germplasm record, and an API
  change would empty it silently — producing a phenotype table with no pea side
  and no error.
- **🚩 red flag** `pageSize` defaults to 10 in BrAPI. The functions pass 5000 /
  10000. A page-size regression does not lose data, it makes one trial take 165
  s instead of 4; if a download suddenly crawls, look here first.
- *Assumption:* `distinct(observationUnitDbId, .keep_all = TRUE)` is the right
  dedup key for the IA/ND quadruplication (the copies differ only in grid
  coordinates and carry no observations).

---

### A10 — `trials` (live BrAPI; the two-stage selection)

```r
arm_evaluation("trials")
b4i <- read_b4i_accessions(c(here::here("data", "Acc_B4I_Avena.txt"),
                             here::here("data", "Acc_B4I_Pisum.txt")))
cand <- candidate_trials(conn, b4i)
disarm_evaluation()
length(b4i); nrow(cand)
```

- **→ returns** 941 accession names (508 oat + 435 pea − 2 monoculture
  labels); `cand` is every trial evaluating any of them.
- **🔍 eyeball** the two selection stages are independent: at least
  `min_b4i_accessions` (20) accessions, **and** oat grain yield recorded on at
  least that many. 34 trials survive.
- **🚩 red flag** `select_trials()` defaults `min_n = min_b4i_accessions`, a
  **global** read from the script's settings block. Change the setting at the
  top and call the function directly and you get the new value; pass it
  explicitly and you get yours. Two ways to set one threshold.
- *Assumption:* `conn$wizard()` is asked in batches of 200 because the filter
  travels in the request. The results are `distinct()`ed afterwards, so a
  trial matched by two batches is counted once.

---

### A11 — `geno` (live BrAPI; the expensive one)

```r
arm_evaluation("geno")
peas <- read_accessions(here::here("data", "Acc_B4I_Pisum.txt"),
                        c("NO_OATS_PLANTED", "NO_PEAS_PLANTED"))
peds <- fetch_pedigrees(conn, peas)
protocols <- confirm_single_protocol(conn, peds$germplasmDbId, cfg, "67")
dos <- marker_dosage(conn, peds, cfg, "67")
disarm_evaluation()
peek(dos); sum(peds$both_parents)
```

- **→ returns** 435 pedigree rows, **233 with both parents** (the pea pedigrees
  uploaded recently), `dos` is 415 × 6164.
- **🔍 eyeball** `confirm_single_protocol()` must report GenoPea 13K (id 67)
  covering 421 of 435 with no other protocol above 5%. For oat it is Oat 3K
  (id 66); three multi-GB GBS archives each cover a single accession and are
  correctly skipped.
- **🚩 red flag** 421 accessions are *covered* by the protocol but only 415
  carry dosages. "Covered" and "genotyped" are not the same thing — T3's
  coverage table can list an accession against a protocol whose archived VCF
  holds no sample of that name. `parent_genotype_status()` encodes exactly this
  distinction; read its `status` column rather than the coverage count.
- **🚩 red flag — do this one.** `verify_parent_order()` is called only in the
  oat script. Pea pedigrees are new and have never been checked:
  ```r
  verify_parent_order(conn, peds, 25)
  ```
  It should report 25 of 25 agreeing that `pedigree = seed/pollen`. If it does
  not, every "selfed female" conclusion on the pea side is reversed.
- **🚩 red flag** only **1 of 112 named pea parents** is genotyped on the array,
  so the parent-identity rule has almost nothing to work with for pea, and the
  pea curation script does not fetch parents into the dosage matrix at all. The
  oat and pea curations are not doing the same work despite the shared
  functions.
- *Assumption:* one protocol supplies every marker. The script stops if a
  second covers more than 5% of accessions — a real guard, not a comment.

---

### A12 — `CURATION.md` is generated, and can go stale

```r
lapply(c("oat", "pea"), \(sp)
  readRDS(here::here("output", paste0(sp, "_curation_settings.rds")))$run_at)
file.info(here::here("CURATION.md"))$mtime
```

- **🚩 red flag — live example.** As of 2026-09-23 the pea curation last ran
  **2026-09-22 14:05** and `CURATION.md` was written **2026-09-21 14:21**. The
  document therefore describes a run that no longer exists: it says pea has no
  pedigrees and no families, when the current outputs have 233 pedigrees and 36
  families. `curation_report.R` reads the recorded settings and is the only
  thing that keeps the two in step, and nothing runs it automatically.
- **Rule:** run `Rscript code/curation_report.R` after every curation run, and
  check the two timestamps before quoting `CURATION.md`.

---

## 5. What each module should generate

The artifact side of the same story: after a clean run, this is what should
exist and roughly what it should say. Anything absent, empty, or of the wrong
shape is a finding.

| Script | Writes | Expected shape / value |
|---|---|---|
| `create_GRMs_T3.R` | `data/GRM_Avena.rds` | list with `$G` 508 × 508, protocols = Oat 3K only, 8 prior-only rows |
| | `data/GRM_Pisum.rds` | `$G` 435 × 435, GenoPea 13K only, 20 prior-only rows |
| `find_trials_…R` | `output/B4I_trial_search.csv` | every candidate trial with `n_b4i` |
| | `output/B4I_trials_selected.csv` | 34 trials, both selection stages passed |
| | `output/B4I_trait_availability.csv/.png` | 34 × 38 counts, 0 = not measured |
| | `data/B4I_observations.rds/.csv.gz` | 139,604 rows, long format |
| | `output/trial_cache/obs_*.rds`, `units_*.rds` | 42 obs files, **6** units files |
| `curate_oat_accessions.R` | `oat_pedigrees.csv` | 508 rows |
| | `oat_parent_genotypes.csv`, `oat_accession_vs_parents.csv` | parent status; r to each parent |
| | `oat_family_correlations.csv` | 74 families, 9 clonal (mean r > 0.985) |
| | `oat_identity_groups.csv` | 19 groups, 90 accessions, largest 12 |
| | `oat_analysis_names.csv` | **50 rows** → 8 distinct analysis names |
| | `oat_curation_settings.rds` | the thresholds this run used, plus `run_at` |
| `curate_pea_accessions.R` | `pea_pedigrees.csv` | 435 rows, 233 with both parents |
| | `pea_family_correlations.csv` | 36 families, **0** clonal |
| | `pea_identity_groups.csv` | 1 group, 2 accessions |
| | `pea_analysis_names.csv` | **1 row**: `NDP170084G` → `ND VICTORY` |
| | *(not written)* | no `pea_parent_genotypes.csv` / `pea_accession_vs_parents.csv` — the pea script never fetches parents |
| `curation_report.R` | `CURATION.md` | must be newer than both `*_curation_settings.rds` |
| `assemble_B4I_phenotypes.R` | `B4I_intercrop_pheno.rds/.csv` | 2377 × 13; 2371 with both yields; 442 oat × 423 pea; 2059 combinations, 90.8% of them in a single plot |
| `BGLR_multi_trait_model.R` | `BGLR_variance_components.csv` | 12 rows (4 seeds × {G_oat, G_pea, residual}); oat Pr 394 / As 105, pea Pr 186 / As 299, residual cor −0.12 |
| | `BGLR_{oat,pea}_effects_all_seeds.csv` | 4 seeds × accessions, `PrEff`/`AsEff`/`GMA` |
| | `BGLR_{oat,pea}_rank_stability.csv` | mean rank and sd across chains |
| | `BGLR_multitrait_*.rds` | the seed-1 fit |
| | `BGLR_seed_*_{R,Omega_*}.dat` | BGLR's own MCMC samples (gitignored) |
| `megalmm_build_inputs.R` | `megalmm_inputs_{raw,centered,standardized}.rds` | `$Y` 237 × 222 with 1682 cells (3.20%); `$X_Env` 222 × 61 in groups 1/2/58; `$G_oat` 237², `$G_pea` 222² |
| `megalmm_oat_pea.R` | `megalmm_cv_results.csv` | 18 fits × 6 predictors; `n_held` = 336 |
| | | best MegaLMM 0.072, best baseline 0.279 |
| | `megalmm_fit_centered.rds` | `$posterior$Lambda` 10 × 222, `$U`/`$Eta_mean` 237 × 222 |
| | `megalmm_cv_accuracy.png` | the comparison |

---

## 6. Subtle-bug catalogue

Ranked by how badly a wrong answer would mislead, not by likelihood.

1. **Duplicated helpers — found, and fixed 2026-09-23.** `collapse_grm` and
   `read_grm` were defined in both `dge_ige_functions.R` and
   `megalmm_build_inputs.R`, with *different bodies*: the MegaLMM copies had
   dropped the column-name guard, the "collapsed N → M" message and the "run
   `create_GRMs_T3.R` first" error. The BGLR model and the MegaLMM model were
   collapsing the GRM through different code. They happened to agree. Likewise
   `connect_t3` (identical copies in `curation_functions.R` and
   `find_trials_with_B4I_accessions.R`) and `grm_factor` (`dge_ige_functions.R`
   and `sim_generate.R`). Each now has exactly one definition —
   `connect_t3` in `t3_functions.R`, the three GRM helpers in
   `dge_ige_functions.R` — and every script sources it. The check that keeps it
   that way is `eval_conflicts()` reporting no duplicated function.
2. **Trait column order in `Y`.** `c("peaYield", "oatYield")` plus
   `effect_roles` is the only thing assigning producer vs associate. A swap
   inverts every conclusion with no error.
3. **`Eta_mean` row/column names are never checked** the way `U`'s are, and
   `Eta_mean` is what the CV scores. Currently correct.
4. **Two "same line" thresholds.** 0.99 for identity, but
   `pool_clonal_families()` merges at `within_r − 0.01`, i.e. about 0.975. The
   12-member `IL18-8562_self` group depends on which one you believe.
5. **`pea_trait_blues()` globs the cache directory** instead of filtering to
   the six intercrop trials. Correct only because the cache happens to hold
   exactly those six.
6. **Single linkage in `identity_groups()`.** Deliberate, but it means group
   membership can be chained through a bridge. Check minimum within-group r,
   not mean.
7. **Prior-only GRM rows** (8 oat, 20 pea) enter every model as lines related
   to nothing. 16 of the 20 pea ones now have pedigrees and genotyped full
   sibs, so this is fixable and currently unfixed.
8. **`CURATION.md` is generated but not regenerated automatically.** It is
   stale right now.
9. **`select_trials()` reads a global default.** Threshold settable two ways.
10. **`G_mix` is built and factored when `fit_mix_term` is FALSE.** Wasted
    work, not a wrong answer.

---

## 7. Where things live

- `data/` — versioned inputs, including the two GRMs and the downloaded
  observations. Rebuilding them from T3 takes hours, which is why they are not
  in `output/`.
- `output/` — everything derived, gitignored apart from its README. Deleting it
  loses nothing but time.
- `output/trial_cache/` — one file per trial of observations (`obs_*`) and
  observation units (`units_*`). 42 obs, 6 units.
- `output/megalmm_runs/` — MegaLMM's on-disk run state; needed while sampling,
  disposable afterwards.
- `T3GenoTools::geno_cache_root()` — outside the repo (about 9 GB), holding
  VCFs, dosages and per-protocol GRMs, shared by every project on the machine.
- `.Renviron` — credentials, gitignored, read by `readRenviron()` inside each
  script so it works whatever directory R started in.
