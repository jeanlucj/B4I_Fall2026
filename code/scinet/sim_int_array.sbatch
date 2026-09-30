#!/bin/bash
# sim_int_array.sbatch -- the INTERACTION-FOCUSED simulation as a SLURM array.
#
#   sbatch -A <account> code/scinet/sim_int_array.sbatch
#   sbatch -A <account> --qos=debug --time=00:30:00 --array=1-1 \
#          code/scinet/sim_int_array.sbatch --check
#
# This is the second, smaller design (code/sim_int_config.R), not the main grid.
# It locates the two boundaries at which MegaLMM starts to beat the Kronecker
# kernel on the oat x pea interaction: a density threshold and a rank threshold.
#
# 120 data scenarios x 3 replicates = 360, each with one BGLR dge_ige fit and
# two MegaLMM runs (pinned and not). Over the default 20 tasks that is 18
# scenarios each. Estimated ~65 single-core hours in total, about 3.3 h per
# task -- comfortable inside the wall clock below. MegaLMM dominates the cost
# here, not BGLR, because it runs twice per scenario and in both orientations.
#
# Its cache and outputs are named separately from the main simulation
# (output/simulation_int/, simulation_int_results.csv) so the two can never be
# combined by accident -- which matters, because the two designs sweep
# different levels and a globbed CSV would silently mix them.
#
# logs/ must exist before submitting: SLURM resolves --output when the task
# starts, before the mkdir below.

#SBATCH --job-name=b4i_sim_int
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=12:00:00
#SBATCH --partition=ceres
#SBATCH --output=logs/sim_int_%A_%a.log

set -euo pipefail

module load r/4.5.3

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

ntasks="${SLURM_ARRAY_TASK_COUNT:-1}"
task="${SLURM_ARRAY_TASK_ID:-1}"

echo "host $(hostname)  task ${task} of ${ntasks}  started $(date)"

set +e
Rscript code/sim_int_run.R --task "${task}" --ntasks "${ntasks}" "$@"
status=$?
set -e

if [ "${status}" -ne 0 ]; then
  echo "task ${task} FAILED with exit status ${status} at $(date)"
  exit "${status}"
fi

echo "task ${task} finished $(date)"
