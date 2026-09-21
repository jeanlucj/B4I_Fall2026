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
Rscript -e 'install.packages(c("tidyverse","here","BGLR","withr","patchwork"))'
Rscript -e 'remotes::install_github("deruncie/MegaLMM")'
```

The simulation needs `BGLR`, `MegaLMM`, `tidyverse`, `here` and `withr`. It
does **not** need `BrAPI.R`, `T3GenoTools` or a T3 login — it reads the GRMs
from `output/`, which is gitignored, so copy those two files across:

```bash
# from the laptop
scp output/GRM_Avena.rds output/GRM_Pisum.rds \
    <first.last>@ceres.scinet.usda.gov:/project/<account>/B4I_Fall2026/output/
```

Alternatively regenerate them there with `code/create_GRMs_T3.R`, which does
need the T3 credentials.

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

## 5. What to watch

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

## 6. Bringing results back

```bash
# from the laptop
scp '<first.last>@ceres.scinet.usda.gov:/project/<account>/B4I_Fall2026/output/simulation_results.csv' output/
```

The per-scenario `.rds` files stay on the cluster; the combined CSV is all the
analysis needs.
