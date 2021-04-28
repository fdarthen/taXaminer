#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=freya
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=64gb
#SBATCH --job-name="diamond"
#SBATCH --output=array_%A_%a.o.out     # Standard output
#SBATCH --error=array_%A_%a.e.out      # and error log
#SBATCH --array=1                                                    # Array range. The %X specifies that max. X jobs run in parallel


echo TaskID $SLURM_ARRAY_TASK_ID


proteins_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['proteins_path'])")
taxon_hits_lca_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['taxon_hits_lca_path'])")
best_taxon_hit_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['best_taxon_hit_path'])")
nr_db_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['nr_db_path'])")
tax_id=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['tax_id'])")

echo ${proteins_path}
echo ${taxon_hits_lca_path}
echo ${best_taxon_hit_path}
echo ${nr_db_path}
echo ${tax_id}

diamond blastp -q "${proteins_path}" -o "${taxon_hits_lca_path}" -d "${nr_db_path}" -f 102 -b2.0 --tmpdir /dev/shm --sensitive --taxon-exclude "$tax_id" --top 10 -c1

#diamond blastp -q "${proteins_path}" -o "${best_taxon_hit_path}" -d "${nr_db_path}" -f 6 qseqid sseqid evalue bitscore score pident staxids sscinames -b2.0 --tmpdir /dev/shm --sensitive --taxon-exclude "$tax_id" -c1 -k 1
