#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=user
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=32gb
#SBATCH --job-name="pbc"
#SBATCH --output=./array_%A_%a.o.out     # Standard output
#SBATCH --error=./array_%A_%a.e.out      # and error log
#SBATCH --array=1 # -3%3 # <- set array of tasks

if [ -z "$1" ]; then
    config_file=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID 'FNR==arrayID{print $0}' seedfile.txt)
else
    config_file=$1
echo $SLURM_ARRAY_TASK_ID $config_file

. ./compute_pbc.sh "$config_file"
