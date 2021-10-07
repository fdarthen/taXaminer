#!/bin/bash
#SBATCH --partition=pool
#SBATCH --account=user
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --job-name="milts"
#SBATCH --output=./%A_%a.o.out     # Standard output
#SBATCH --error=./%A_%a.e.out      # and error log
#SBATCH --array=1

PATH=$PWD/tools:$PATH

if [[ "$1" == '-s' || "$1" == '--seed' ]]; then
    config_path=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID 'FNR==arrayID{print $0}' $2)
elif [[ -f "${1}" ]]; then
    config_path=$1
else
    echo "Config specification unclear. Either state path to config file or path to slurm seedfile using the option -s/--seed"
    [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
fi

echo $SLURM_ARRAY_TASK_ID $config_path

./MILTS.sh "$config_path"
