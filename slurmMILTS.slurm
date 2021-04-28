#!/bin/bash
#SBATCH --partition=all
#SBATCH --account=freya
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16gb
#SBATCH --job-name="milts"
#SBATCH --output=./array_%A_%a.o.out     # Standard output
#SBATCH --error=./array_%A_%a.e.out      # and error log


start=`date +%s`

# read necessary variables from config file
fasta_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['fasta_path'])")
gff_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['gff_path'])")
output_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_path'])")
proteins_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['proteins_path'])")
tax_id=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['tax_id'])")
taxon_hits_lca_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['taxon_hits_lca_path'])")
best_taxon_hit_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['best_taxon_hit_path'])")
nr_db_path=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['nr_db_path'])")
# BOOLS
only_plotting=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['only_plotting'])")
extract_proteins=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['extract_proteins'])")
compute_tax_assignment=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['compute_tax_assignment'])")
compute_pbc=$(cat ./config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['compute_pbc'])")


echo "FASTA: " $fasta_path
echo "GFF: " $gff_path
echo "Proteins: " $proteins_path
echo "Output directory: " $output_path
echo "Taxonomic assignment: "
echo $taxon_hits_lca_path
echo $best_taxon_hit_path
echo "NCBI Taxon ID: " $tax_id


if [ "${only_plotting}" = "FALSE" ]
then

  # check if proteins FASTA should be extracted but exists
  if [ "${extract_proteins}" = "TRUE" ]; then
    if [[ -f "${proteins_path}" ]]; then
      echo "Proteins FASTA file exists but it is set to be extracted. This process will overwrite it. Do you want to continue? (y/n)"
      read input
      if [ "$input" = "n" ]; then
        echo "Please change the option for 'extract_proteins' to 'FALSE'"
        [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
      else
        echo "Proteins FASTA file will be overwritten"
      fi
    fi
  fi

  # check if taxonomic assignment should be performed but files exist
  if [ "${compute_tax_assignment}" = "TRUE" ]; then
     if [[ -f "${taxon_hits_lca_path}" ]]; then
       echo "LCA hit file exists but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
       read input
       if [ "$input" = "n" ]; then
         echo "Please change the option for 'compute_tax_assignment' to 'FALSE'"
        [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
       else
         echo "LCA hit file will be overwritten"
       fi
     fi
     if [[ -f "${best_taxon_hit_path}" ]]; then
       echo "Best hits file exists but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
       read input
       if [ "$input" = "n" ]; then
         echo "Please change the option for 'compute_tax_assignment' to 'FALSE'"
         [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
       else
         echo "Best hits file will be overwritten"
       fi
    fi
  fi

  # check if PBC file should be generated but file exists
  if [ "${compute_pbc}" = "TRUE" ]; then
    if [[ -f "${pbc_path}" ]]; then
      echo "PBC file exists but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
      read input
      if [ "$input" = "n" ]; then
        echo "Please change the option for 'compute_pbc' to 'FALSE'"
        [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
      else
        echo "PBC file will be overwritten"
      fi
    fi
    echo "compute per pase coverage start:"
    time0_1=`date +%s`
    bash ./additional_scripts/compute_pbc.sh
    time0_2=`date +%s`
    echo "compute per pase coverage end (time elapsed:" $(($time0_2-$time0_1)) "s)"
  fi


  # 1.a) remove newlines from fasta
  awk '/^>/{if(NR==1){print}else{printf("\n%s\n",$0)}next} {printf("%s",$0)} END{printf("\n")}' $fasta_path >> tmp.MILTS.fasta

  # 2) start python script --> produces descriptive gene statistics
  echo -e "\nproduce gene info start:"
  time1_1=`date +%s`
  python produce_gene_info.py
  time1_2=`date +%s`
  echo "produce gene info end (time elapsed:" $(($time1_2-$time1_1)) "s)"


  # 3) start R script
  echo "perform PCA and clustering start:"
  time2_1=`date +%s`
  Rscript perform_PCA_and_clustering.R --verbose >> $output_path"R_log.out" # in working directory, should also contain config.R
  time2_2=`date +%s`
  echo "perform PCA and clustering end (time elapsed:" $(($time2_2-$time2_1)) "s)"


  # 4.a) create directory in MILTS results to store orphans analysis
  [[ ! -d "${output_path}taxonomic_assignment" ]] && mkdir -p "${output_path}taxonomic_assignment"


  if [ "${extract_proteins}" = "TRUE" ]
  then
    echo "retrieving peptide sequenes start:"
    time3_1=`date +%s`
    gffread -S --table "@geneid,@cdslen" -y ${proteins_path} -g ${fasta_path} ${gff_path}
    python ./additional_scripts/longest_cds.py
    time3_2=`date +%s`
    echo "retrieving peptide sequenes end (time elapsed:" $(($time3_2-$time3_1)) "s)"
  fi
  date



  # 4.b) run taxonomic assignment with Diamond

  if [ "${compute_tax_assignment}" = "TRUE" ]
  then
    echo "taxonomic assignment start:"
    time4_1=`date +%s`
    diamond blastp -q "${proteins_path}" -o "${taxon_hits_lca_path}" -d "${nr_db_path}" -f 102 -b2.0 --tmpdir /dev/shm --sensitive --taxon-exclude "$tax_id" --top 10 -c1
    diamond blastp -q "${proteins_path}" -o "${best_taxon_hit_path}" -d "${nr_db_path}" -f 6 qseqid sseqid evalue bitscore score pident staxids sscinames -b2.0 --tmpdir /dev/shm --sensitive --taxon-exclude "$tax_id" -c1 -k 1
    time4_2=`date +%s`
    echo "taxonomic assignment end (time elapsed:" $(($time4_2-$time4_1)) "s)"
  fi
fi


# 5) plot genes with PCA coordinates and taxonomic assignment
echo "plot taxonomic assignment start:"
time5_1=`date +%s`
# grepping the GFF to only relevant lines accelerates gffread on large files immensely
zgrep -E "CDS|gene" ${gff_path} | gffread - -o ${output_path}"tmp.prot_gene_matching.txt" --table "@id,@geneid"
Rscript plotting.R --verbose >> $output_path"R_log.out"
time5_2=`date +%s`
echo "plot taxonomic assignment end (time elapsed:" $(($time5_2-$time5_1)) "s)"

# 1.b) remove the temporarily created newline-free FASTA
[[ -e ${output_path}"tmp.prot_gene_matching.txt" ]] && rm ${output_path}"tmp.prot_gene_matching.txt"
[[ -e tmp.MILTS.fasta ]] && rm tmp.MILTS.fasta


end=`date +%s`
runtime=$(($end-$start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
