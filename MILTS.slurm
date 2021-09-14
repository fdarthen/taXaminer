#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=freya
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=16gb
#SBATCH --job-name="milts"
#SBATCH --output=./%A_%a.o.out     # Standard output  /share/project/freya/T_draco/%A_%a.o.out
#SBATCH --error=./%A_%a.e.out      # and error log    /share/project/freya/invertebrates/logs
#SBATCH --array=1

PATH=$PWD/tools:$PATH

if [ -z "$1" ]; then
    config_path=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID 'FNR==arrayID{print $0}' seedfile.txt)
else
    config_path=$1
fi

echo $SLURM_ARRAY_TASK_ID $config_path

./MILTS.sh "$config_path"
