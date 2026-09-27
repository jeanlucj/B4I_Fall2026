# Running the simulation on SciNet (Ceres)

The grid is **120 data scenarios** and a **D-optimal fraction of 150 MegaLMM
runs**, each of those fitted in both orientations. Per data scenario that is two
bivariate `BGLR::Multitrait` fits (additive and DGE-IGE) plus whichever MegaLMM
settings the fraction assigned to it. It is many hours on a laptop and
embarrassingly parallel, so it belongs on the cluster. This directory holds a
SLURM job array that splits it across tasks.

**The cost is now in the BGLR half, not MegaLMM.** Multitrait is slower than the
univariate model it replaced, it runs once per scenario regardless of the
fraction, and its Kronecker interaction term is what needs the memory. MegaLMM is
flat in sparsity and was the half the fraction cut. Size the jobs for BGLR.

Conventions here follow `T3Predictathon2026/scripts/Analysis_Claude/optimizer`;
its `RUNBOOK_SLURM.md` is the fuller guide to Ceres itself and is worth reading
once. What follows is only what this repository needs.

## 1. Get the code and the R packages there

```bash
ssh <first.last>@ceres.scinet.usda.gov
sacctmgr -Pns show user format=account,defaultaccount   # your account, for -A
my_quotas                                                # check there is room

cd /project/<account>
git clone https://github.com/jeanlucj/B4I_Fall2026.git
cd B4I_Fall2026
```

### The `.Renviron` trap

**R reads exactly one `.Renviron`** — the working directory's if one exists,
otherwise `$HOME`'s. It does not merge them and it does not walk up.

This repository has its own `.Renviron` for `T3_USERNAME` and `T3_PASSWORD`, so
a `~/.Renviron` holding `R_LIBS_USER` is **never read** and packages land in
the default library inside your 30 GB home. Put all three in the repository's
`.Renviron`:

```
T3_USERNAME=...
T3_PASSWORD=...
R_LIBS_USER=/project/<account>/R_packages/%v
```

Then install, on a compute node rather than the login node:

```bash
salloc -N1 -n8 --mem=32G -t 2:00:00 -A <account>
module load r/4.5.3

# the library directory must exist before R will use it
mkdir -p /project/<account>/R_packages/4.5

Rscript -e 'install.packages(c("remotes","tidyverse","here","BGLR","withr","AlgDesign"), repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("deruncie/MegaLMM")'
```

**`AlgDesign` is not optional.** `sim_config.R` calls
`AlgDesign::optFederov()` to build the D-optimal fraction, and it is called
while the script is loading its design — so without it every task dies before
fitting anything, including under `--check`. Confirm all five before submitting:

```bash
Rscript -e 'for (p in c("tidyverse","here","BGLR","MegaLMM","AlgDesign","withr"))
             cat(sprintf("%-10s %s\n", p, requireNamespace(p, quietly=TRUE)))'
```

**Pass `repos=` explicitly.** Without it a batch `Rscript -e
install.packages(...)` has no mirror to look in and reports
`packages ... are not available for this version of R`, which reads like an
R-version incompatibility and is not one — see
[below](#if-installpackages-says-a-package-is-not-available). `remotes` is in
the list because the next line needs it.

The simulation needs `BGLR`, `MegaLMM`, `AlgDesign`, `tidyverse`, `here` and
`withr`. It
does **not** need `BrAPI.R`, `T3GenoTools` or a T3 login: the GRMs live in
`data/` and are versioned, so `git clone` brings them and there is nothing to
copy across. Only `code/create_GRMs_T3.R` and
`code/find_trials_with_B4I_accessions.R` need credentials, and only to refresh
what the repository already holds.

### If install.packages says a package is "not available"

```
Warning: packages 'tidyverse', 'BGLR' are not available for this version of R
```

Nearly always this means R could not see a repository holding them, not that
they are incompatible with R 4.5.3. Diagnose in this order:

```r
R.version.string            # 4.5.3 after module load
getOption("repos")          # "@CRAN@" unresolved, or empty => this is the cause
nrow(available.packages())  # 0, or an error => no usable mirror
```

`@CRAN@` is a placeholder that an interactive session resolves by asking which
mirror to use. A batch session has nobody to ask, so it stays unresolved and
every package looks missing. Fixes, in increasing permanence:

```r
install.packages("BGLR", repos = "https://cloud.r-project.org")   # per call
options(repos = c(CRAN = "https://cloud.r-project.org"))          # per session
```

For something permanent, put that `options()` line in an `.Rprofile`. But note
the same one-file-only rule that bites with `.Renviron` above: R reads the
working directory's `.Rprofile` if there is one and `$HOME`'s otherwise, never
both — and **this repository has its own**, from workflowr. So a line in
`~/.Rprofile` will be ignored whenever you run from the repository. Put it in
the repository's `.Rprofile`, or pass `repos=` on the command line, which
always works.

Other causes, once the repository is definitely set:

- **The name is wrong.** Case matters: `MegaLMM`, not `megalmm`.
- **It is not on CRAN.** `MegaLMM`, `BrAPI.R`, `T3GenoTools` and
  `T3BrapiHelpers` are GitHub-only and need `remotes::install_github()`.
- **It failed to compile and the warning names a dependency.** `MegaLMM` builds
  C++ via `Rcpp` and `RcppEigen`. Read the log for the first error, not the last
  warning.
- **`R_LIBS_USER` points somewhere that does not exist**, so R silently falls
  back to the default library inside your 30 GB home. Check with
  `Sys.getenv("R_LIBS_USER")` and `.libPaths()`.

## 2. Shake it out first

Three steps, cheapest first. All of them are worth doing on a fresh clone,
because each rules out a different kind of failure.

```bash
# 1. the unit tests: no cluster, no data, ~45 s on the login node
Rscript tests/run_all.R

# 2. the sanity check: one dense scenario, ~1 min
sbatch -A <account> --qos=debug --time=00:30:00 --array=1-1 \
       code/scinet/sim_array.sbatch --check

# 3. the pilot: four of the cheapest real runs, end to end including the cache
sbatch -A <account> --qos=debug --time=00:30:00 --array=1-1 \
       code/scinet/sim_array.sbatch --pilot
```

`tests/run_all.R` catches a bad install and a bad clone without queueing
anything. `--check` fits one dense scenario in which MegaLMM must recover the
simulated signal and fails loudly if it does not — the wiring mistake it guards
against looks exactly like "the method does not work here". `--pilot` runs four
real scenarios from the cheapest corner of the design, so it also exercises the
cache paths and the combine step. Read `logs/sim_<jobid>_<task>.log`.

## 3. Run the grid

```bash
mkdir -p logs                                    # SLURM will not create it
sbatch -A <account> code/scinet/sim_array.sbatch
```

**`mkdir -p logs` matters.** `--output=logs/sim_%A_%a.log` is resolved by SLURM
when the task starts, before the script's own `mkdir` runs, so on a fresh clone —
where `logs/` does not exist, because `output/` and `logs/` are gitignored — every
task fails with nothing written anywhere. It is the most confusing version of
section 5's "no output" problem, because there is no file to read at all.

Defaults in the script: `--array=1-20`, 4 cpus, **16 GB per cpu**, 12 hours,
partition `ceres`. Twenty tasks over 120 scenarios is **six scenarios each**.
The array slices by *data scenario*, not by run, so a scenario's BGLR and
MegaLMM halves always land in the same task and no experiment is simulated
twice.

**Twelve hours is far more than a task needs.** Measured on the worst cell
(400 × 400 at 48%, 61,440 training plots per trait): the additive fit is 3.9 min
and the DGE-IGE fit with the rank-30 Kronecker term is 9.0 min, both at the
production 6,000 iterations. Scaling that across the design — BGLR is close to
linear in the number of observations — the whole 120-scenario BGLR half is about
**6 single-core hours**, which at `--array=1-20` is **0.1–0.4 h per task**. See
[SIMULATION.md](../../SIMULATION.md#how-long-it-takes) for the table. MegaLMM adds
to it but is flat in sparsity.

The sparsity axis spans 1.6% to 48% — a factor of thirty in the number of
observations — so the slices are not equally expensive, and the slicing is by
position rather than by cost. That also turns out not to matter: the 15 worst
cells are regularly spaced in the grid ordering, so at `--array=1-20` **each task
gets exactly six scenarios and at most one of them**, and the same holds at
`1-40`. Measured, not assumed.

If a task does hit the wall clock it costs only its unfinished scenario —
everything finished is cached — so the remedy is to resubmit the same array.

Arguments after the script name are passed through to `sim_run.R`:

```bash
# throttle to ten concurrent tasks
sbatch -A <account> --array=1-20%10 code/scinet/sim_array.sbatch

# more replicates
sbatch -A <account> code/scinet/sim_array.sbatch --reps 5

# the cheap half first, to get an answer while the rest runs
sbatch -A <account> --array=1-10 code/scinet/sim_array.sbatch --filter "n_acc == 200"

# a bigger fraction, or none at all (960 candidate runs -- days, not hours)
sbatch -A <account> code/scinet/sim_array.sbatch --runs 200
sbatch -A <account> code/scinet/sim_array.sbatch --full
```

`--filter` takes an R expression over the **composite** design columns, which is
what `sim_config.R` recodes the nested axes into — `n_acc`, `sparsity`,
`interaction` (`none`, `f1_i10`, `f1_i20`, `f5_i10`, `f5_i20`), `environment`
(`one`, `ten_stable`, `ten_gxe`), `K`, `eigen_variance`, `fixed_main_effect`. So
`--filter "environment == 'one'"`, not `--filter "n_envs == 1"`.

The full list of `sim_run.R` flags is at the top of `code/sim_run.R` and in
[SIMULATION.md](../../SIMULATION.md#running-it).

## 4. Combine

No single task sees every scenario, so the combined table is rebuilt from the
cache afterwards:

```bash
Rscript code/sim_run.R --combine
```

This refits nothing. It reads `output/simulation/*.rds` and writes
`output/simulation_results.csv`, `output/simulation_design_effects.csv` and
`output/simulation_summary.png`. It is safe to run while tasks are still going,
to see partial results.

`simulation_design_effects.csv` is the one to read first. Because the fraction
gives most scenarios only some of the MegaLMM settings, a table of cell means
compares unlike with unlike; the design model — main effects and two-way
interactions in the composite factors — is what the fraction was built to
estimate. See [SIMULATION.md](../../SIMULATION.md#reading-the-results).

## 5. When a job produces no output

The commonest confusion, and it is a reading problem rather than a failure:

**R writes almost everything to stderr.** `message()`, warnings, errors and the
package startup banners all go there; only `cat()` and `print()` go to stdout.
So a job that dies during the fit leaves a stdout file containing nothing but
the shell's own `echo` lines, and the entire explanation somewhere else.

Older versions of `sim_array.sbatch` split the streams into `.out` and `.err`.
It now merges both into `logs/sim_<jobid>_<task>.log`. If you are looking at a
run from before that change:

```bash
cat logs/sim_*.err        # this is where R actually spoke
```

Then find out how the job ended, which distinguishes an R error from the two
things SLURM does on its own:

```bash
sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS,ReqMem
```

| State | means |
|---|---|
| `FAILED` with `ExitCode 1:0` | R errored — read the log |
| `TIMEOUT` | wall clock; `--qos=debug` gives only 30 minutes |
| `OUT_OF_MEMORY`, or `FAILED` with `ExitCode 0:125` | exceeded `--mem-per-cpu`; Ceres kills rather than throttles |
| `COMPLETED` but no results | check that the array slice had scenarios to run |

The script now reports R's exit status explicitly. Under `set -e` a bare
`Rscript` failure ends the job with nothing written about it, which is the
other half of why an empty log is easy to produce.

## 6. What to watch

- **Memory, not time, is the likely failure.** Ceres kills a job that exceeds
  its allocation rather than throttling it. The worst cell is 400 × 400 at 48%:
  76,800 plots per trait, 61,440 of them in training after the 20% held out, and
  a 61,440 × 900 Kronecker design matrix — about 440 MB, with two temporaries of
  that size built before they are multiplied, and `Multitrait` carrying two
  traits through it. The script asks for 4 × 16 GB, which is deliberately
  generous. If tasks are still killed there, drop `SIM_KRON_RANK` from 30 to 20
  in `code/sim_config.R`, which cuts that term from 900 columns to 400.
- **Restarting is free.** Every scenario is cached by name, so resubmitting
  the same array skips what finished and picks up what did not. A task that
  hits the wall clock costs only its unfinished scenario.
- **`sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS`** after the fact
  tells you what the array actually used, which is how to size the next one.

## 7. Bringing results back

```bash
# from the laptop
scp '<first.last>@ceres.scinet.usda.gov:/project/<account>/B4I_Fall2026/output/simulation_results.csv' output/
```

The per-scenario `.rds` files stay on the cluster; the combined CSV is all the
analysis needs.
