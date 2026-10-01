#!/bin/bash
# sim_decomp_array.sbatch -- the interaction DECOMPOSITION sweep as a SLURM array.
#
#   sbatch -A <account> code/scinet/sim_decomp_array.sbatch
#   sbatch -A <account> --qos=debug --time=00:30:00 --array=1-1 \
#          code/scinet/sim_decomp_array.sbatch --pilot
#
# Asks whether the fitted dge_ige interaction, rewritten as scores x loadings
# (code/interaction_decomp.R), recovers the simulated factor structure or only
# the leading eigen-directions of the two GRMs. See code/sim_decomp_run.R.
#
# WHAT IT COSTS. 120 data scenarios x 5 replicates = 600 cells, plus the
# interaction_pct = 0 null cells that calibrate the spectrum (10 more scenarios
# x 5 replicates). Each cell is one BGLR dge_ige fit with the coefficient draws
# streamed, and one MegaLMM pair for the head-to-head. From the measured
# dge_ige medians in output/simulation_int/simulation_int_results.csv --
# 20/52/98/169/498 s at n = 200 and 92/260/487/847/2748 s at n = 400 -- that is
# about 88 core-hours for the BGLR fits and 32 for MegaLMM, so roughly 120 in
# total and about 6 h per task over 20 tasks.
#
# The wall clock is 24 h rather than the 12 h sim_int uses. Four times the
# estimate: the n = 400 dense cells dominate and are not evenly spread over
# tasks, so a task that draws several of them runs well above the mean, and
# raising --rank later multiplies the cost by (rank/30)^2.
#
# AFTERWARDS, RUN --combine. Each task writes the two CSVs from whatever was in
# the cache when that task finished, so the files a 20-task job leaves behind are
# each a partial view. On the login node:
#
#   Rscript code/sim_decomp_run.R --combine
#
# fits nothing, globs the whole cache and rewrites both CSVs complete.
#
# Cells are cached per scenario-replicate and the job is resumable: a task that
# runs out of time can be resubmitted and will skip what it finished. Each cell
# streams its MCMC draws to output/simulation_decomp/draws/ under a prefix that
# carries the scenario, replicate and seed, so two tasks can never write the
# same file -- which they would if saveAt were a bare directory, because BGLR
# treats it as a filename PREFIX.
#
# logs/ must exist before submitting: SLURM resolves --output when the task
# starts, before the mkdir below.

#SBATCH --job-name=b4i_sim_decomp
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --partition=ceres
#SBATCH --output=logs/sim_decomp_%A_%a.log

set -euo pipefail

module load r/4.5.3

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

ntasks="${SLURM_ARRAY_TASK_COUNT:-1}"
task="${SLURM_ARRAY_TASK_ID:-1}"

echo "host $(hostname)  task ${task} of ${ntasks}  started $(date)"

set +e
Rscript code/sim_decomp_run.R --task "${task}" --ntasks "${ntasks}" "$@"
status=$?
set -e

if [ "${status}" -ne 0 ]; then
  echo "task ${task} FAILED with exit status ${status} at $(date)"
  exit "${status}"
fi

echo "task ${task} finished $(date)"
