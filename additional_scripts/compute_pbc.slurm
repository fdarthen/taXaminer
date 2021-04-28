#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=freya
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=32gb
#SBATCH --job-name="pbc"
#SBATCH --output=/share/project/freya/analysis/array_%A_%a.o.out     # Standard output
#SBATCH --error=/share/project/freya/analysis/array_%A_%a.e.out      # and error log
#SBATCH --array=1-3%3

config_file=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID -F'\t' 'FNR==arrayID{print $1}' seedfile.txt)
echo $SLURM_ARRAY_TASK_ID $config_file

. ./compute_pbc.sh "$config_file"
