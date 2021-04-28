#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=freya
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=32gb
#SBATCH --job-name="milts"
#SBATCH --output=/share/project/freya/NematodeData/MILTS_reports/logs/array_%A_%a.o.out     # Standard output
#SBATCH --error=/share/project/freya/NematodeData/MILTS_reports/logs/array_%A_%a.e.out      # and error log
#SBATCH --array=1-18%4

config_file=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID -F'\t' 'FNR==arrayID{print $1}' seedfile.txt)
echo $SLURM_ARRAY_TASK_ID $config_file

./MILTS.sh "$config_file"
