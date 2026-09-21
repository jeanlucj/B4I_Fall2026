#!/bin/bash
# sim_array.sbatch -- run the simulation grid as a SLURM job array on Ceres.
#
#   sbatch -A <account> code/scinet/sim_array.sbatch
#   sbatch -A <account> --array=1-20%10 code/scinet/sim_array.sbatch   # 10 at a time
#   sbatch -A <account> --qos=debug --array=1-2 --time=00:30:00 code/scinet/sim_array.sbatch
#
# Each array task takes a slice of the scenarios: task k of N runs scenarios
# k, k+N, k+2N, ... Slicing by scenario keeps a scenario's BGLR and MegaLMM
# halves in the same task, so no simulation is generated twice.
#
# Tasks share output/simulation/ and never write the same file, because the
# cache path is built from the scenario name. When they are all done, combine
# with `Rscript code/sim_run.R --combine`, which rebuilds the CSV from the
# cache without refitting anything -- necessary because no single task sees
# every scenario.
#
# Find your account with:  sacctmgr -Pns show user format=account,defaultaccount

#SBATCH --job-name=b4i_sim
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=12:00:00
#SBATCH --partition=ceres
# One file, not two. R sends message(), warnings and errors to stderr and only
# cat()/print() to stdout, so splitting the streams puts most of what this
# pipeline says -- including why it died -- in the file nobody reads, and
# interleaves nothing. Omitting --error makes SLURM merge both into --output.
#SBATCH --output=logs/sim_%A_%a.log

set -euo pipefail

# Ceres kills a job that exceeds its memory rather than throttling it. The
# 400x400 scenarios at 45% build a 72,000 x 900 Kronecker design matrix --
# about 518 MB, with two temporaries of that size before they are multiplied.
# 4 cpus x 8G leaves room; drop SIM_KRON_RANK to 20 in code/sim_config.R if a
# task is still killed.

module load r/4.5.3

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

ntasks="${SLURM_ARRAY_TASK_COUNT:-1}"
task="${SLURM_ARRAY_TASK_ID:-1}"

echo "host $(hostname)  task ${task} of ${ntasks}  started $(date)"

# Report R's own exit status rather than letting set -e swallow it silently:
# a bare `Rscript` failure under `set -e` ends the script with no explanation
# in the log at all.
set +e
Rscript code/sim_run.R --task "${task}" --ntasks "${ntasks}" "$@"
status=$?
set -e

if [ "${status}" -ne 0 ]; then
  echo "task ${task} FAILED with exit status ${status} at $(date)"
  exit "${status}"
fi

echo "task ${task} finished $(date)"
