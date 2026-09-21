# Running the simulation on SciNet (Ceres)

The 360-fit grid — 60 data scenarios × 6 MegaLMM settings — is about 18 hours
on a laptop and embarrassingly parallel, so it belongs on the cluster. This
directory holds a SLURM job array that splits it across tasks.

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

Rscript -e 'install.packages(c("remotes","tidyverse","here","BGLR","withr","patchwork"), repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("deruncie/MegaLMM")'
```

**Pass `repos=` explicitly.** Without it a batch `Rscript -e
install.packages(...)` has no mirror to look in and reports
`packages ... are not available for this version of R`, which reads like an
R-version incompatibility and is not one — see
[below](#if-installpackages-says-a-package-is-not-available). `remotes` is in
the list because the next line needs it.

The simulation needs `BGLR`, `MegaLMM`, `tidyverse`, `here` and `withr`. It
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

```bash
sbatch -A <account> --qos=debug --time=00:30:00 --array=1-1 \
       code/scinet/sim_array.sbatch --check
```

`--check` fits one dense scenario where MegaLMM must recover the signal and
fails loudly if it does not. Thirty minutes is plenty. Read `logs/sim_*.out`.

## 3. Run the grid

```bash
sbatch -A <account> code/scinet/sim_array.sbatch
```

Defaults in the script: `--array=1-20`, 4 cpus, 8 GB per cpu, 12 hours,
partition `ceres`. Twenty tasks over 60 scenarios is three scenarios each,
which is comfortable inside twelve hours even for the 400 × 400 cells.

Useful variations:

```bash
# throttle to ten concurrent tasks
sbatch -A <account> --array=1-20%10 code/scinet/sim_array.sbatch

# more replicates -- arguments after the script go to sim_run.R
sbatch -A <account> code/scinet/sim_array.sbatch --reps 5

# only the cheap half, to get a first answer quickly
sbatch -A <account> --array=1-6 code/scinet/sim_array.sbatch --filter "n_acc == 200"
```

## 4. Combine

No single task sees every scenario, so the combined table is rebuilt from the
cache afterwards:

```bash
Rscript code/sim_run.R --combine
```

This refits nothing. It reads `output/simulation/*.rds` and writes
`output/simulation_results.csv` and `output/simulation_summary.png`. It is
also safe to run while tasks are still going, to see partial results.

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
  its allocation rather than throttling it. The 400 × 400 scenarios at 45%
  build a 72,000 × 900 Kronecker design matrix, about 518 MB, with two
  temporaries of that size before they are multiplied. If tasks die there,
  either raise `--mem-per-cpu` or drop `SIM_KRON_RANK` from 30 to 20 in
  `code/sim_config.R`, which cuts that term to 400 columns.
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
