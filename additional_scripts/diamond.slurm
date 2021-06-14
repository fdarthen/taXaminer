#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=user
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=64gb
#SBATCH --job-name="diamond"
#SBATCH --output=./array_%A_%a.o.out     # Standard output
#SBATCH --error=./array_%A_%a.e.out      # and error log
#SBATCH --array=1                        # Array range. The %X specifies that max. X jobs run in parallel

PATH=$PWD/tools:$PATH

if [ -z "$1" ]; then
    config_file=$(awk -v arrayID=$SLURM_ARRAY_TASK_ID 'FNR==arrayID{print $0}' seedfile.txt)
else
    config_file=$1
echo TaskID $SLURM_ARRAY_TASK_ID


output_path=$(cat "$config_path" | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_path'])")
[[ ! -d "${output_path}" ]] && mkdir -p "${output_path}"
[[ ! -d "${output_path}/tmp" ]] && mkdir -p "${output_path}/tmp"

proteins_path=$(cat "$config_path" | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['proteins_path'])")
taxon_hits_lca_path=$(cat "$config_path" | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['taxon_hits_lca_path'])")
best_taxon_hit_path=$(cat "$config_path" | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['best_taxon_hit_path'])")

nr_db_path=$(cat "$config_path" | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['nr_db_path'])")
tax_id=$(cat "$config_path" | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['tax_id'])")
taxon_exclude=$(cat "$config_path" | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['taxon_exclude'])")

echo ${proteins_path}
echo ${taxon_hits_lca_path}
echo ${best_taxon_hit_path}
echo ${nr_db_path}
echo ${tax_id}

echo "assess LCA and best hit start:"
time4_1=`date +%s`
if [ "${taxon_exclude}" = "TRUE" ]; then
		diamond blastp -q "${output_path}tmp/tmp.longest_cds.protein.fasta" -o "${taxon_hits_lca_path}" -d "${nr_db_path}" -f 102 -b2.0 --tmpdir /dev/shm --sensitive --top 10 -c1 --taxon-exclude "$tax_id"
		diamond blastp -q "${output_path}tmp/tmp.longest_cds.protein.fasta" -o "${best_taxon_hit_path}" -d "${nr_db_path}" -f 6 qseqid sseqid evalue bitscore score pident staxids sscinames -b2.0 --tmpdir /dev/shm --sensitive -c1 -k 1 --taxon-exclude "$tax_id"
else
		diamond blastp -q "${output_path}tmp/tmp.longest_cds.protein.fasta" -o "${taxon_hits_lca_path}" -d "${nr_db_path}" -f 102 -b2.0 --tmpdir /dev/shm --sensitive --top 10 -c1
		diamond blastp -q "${output_path}tmp/tmp.longest_cds.protein.fasta" -o "${best_taxon_hit_path}" -d "${nr_db_path}" -f 6 qseqid sseqid evalue bitscore score pident staxids sscinames -b2.0 --tmpdir /dev/shm --sensitive -c1 -k 1
fi
time4_2=`date +%s`
echo "assess LCA and best hit end (time elapsed:" $(($time4_2-$time4_1)) "s)"
