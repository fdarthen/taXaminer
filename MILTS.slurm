#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=user
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8gb
#SBATCH --job-name="milts"
#SBATCH --output=./%A_%a.o.out     # Standard output
#SBATCH --error=./%A_%a.e.out      # and error log
#SBATCH --array=1

PATH=$PWD/tools:$PATH

CMD=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
SCRIPT_DIR=$(dirname $(realpath $(echo $CMD | cut -d' ' -f1)))

if [[ "$1" == '-s' || "$1" == '--seed' ]]; then
    config_path=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID 'FNR==arrayID{print $0}' $2)
elif [[ -f "${1}" ]]; then
    config_path=$1
else
    echo "Config specification unclear. Either state path to config file or path to slurm seedfile using the option -s/--seed"
    [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
fi

echo $SLURM_ARRAY_TASK_ID $config_path

$SCRIPT_DIR/MILTS.sh "$config_path"
