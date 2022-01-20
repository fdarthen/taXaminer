#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=user
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8gb
#SBATCH --job-name="milts"
#SBATCH --output=./%A_%a.o.out     # Standard output
#SBATCH --error=./%A_%a.e.out      # and error log
#SBATCH --array=1

# path to tools if conda installation wasn't used
PATH=$PWD/tools:$PATH

## get path where script is located
# to enable call from outside of MILTS directory
CMD=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
SCRIPT_DIR=$(dirname $(realpath $(echo $CMD | cut -d' ' -f1)))

## read command line arguments
# if flag -s or --seed is used, the next arguments is considered as path to
# a seedfile containing in each line one path to a config file
if [[ "$1" == '-s' || "$1" == '--seed' ]]; then
    config_path=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID 'FNR==arrayID{print $0}' $2)
elif [[ -f "${1}" ]]; then
    config_path=$1
else
    echo "Config specification unclear. Either state path to config file or path to slurm seedfile using the flags -s/--seed"
    [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
fi

# print Task ID and path to config file
echo $SLURM_ARRAY_TASK_ID $config_path

# call MILTS
python3 $SCRIPT_DIR/milts.py "$config_path"
