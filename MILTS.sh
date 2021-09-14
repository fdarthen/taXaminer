#!/bin/bash

PATH=$PWD/tools:$PATH

start=`date +%s`

config_path=$1


# read necessary variables from config file
fasta_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['fasta_path'])")
gff_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['gff_path'])")
output_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_path'])")
proteins_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['proteins_path'])")
tax_id=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['tax_id'])")
# taxon_hits_lca_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['taxon_hits_lca_path'])")
# best_taxon_hit_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['best_taxon_hit_path'])")
nr_db_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['nr_db_path'])")
pbc_paths_list=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['pbc_paths'])")
pbc_path=$(echo $pbc_paths_list | cut -d ',' -f1 | awk -F '[' '{print $2}' | awk -F ']' '{print $1}' | awk -F "'" '{print $2}')

# BOOLS
only_plotting=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['only_plotting'])")
extract_proteins=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['extract_proteins'])")
compute_tax_assignment=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['compute_tax_assignment'])")
taxon_exclude=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['taxon_exclude'])")
compute_pbc=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['compute_pbc'])")
output_pdf=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_pdf'])")
output_png=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_png'])")

# config=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin))")
# echo $config # -> this is string of dictionary; parsing?

source <(grep = test_config.ini | sed 's/ *= */=/g') # TODO: new config file?


echo "Config: " $config_path
echo "FASTA: " $fasta_path
echo "GFF: " $gff_path
echo "Proteins: " $proteins_path
echo "PBC(s): " $pbc_paths_list
echo "Output directory: " $output_path
# echo "Taxonomic assignment: "
# echo $taxon_hits_lca_path
# echo $best_taxon_hit_path
echo "NCBI Taxon ID: " $tax_id
echo -e "\n"

[[ ! -d "${output_path}" ]] && mkdir -p "${output_path}"
[[ ! -d "${output_path}tmp/" ]] && mkdir -p "${output_path}tmp/"

if [ "${only_plotting}" = "FALSE" ]; then

    # 1.a) remove newlines from fasta
    awk '/^>/{if(NR==1){print}else{printf("\n%s\n",$0)}next} {printf("%s",$0)} END{printf("\n")}' $fasta_path > "${output_path}tmp/tmp.MILTS.fasta"
    samtools faidx "${output_path}tmp/tmp.MILTS.fasta"
    # 1.b) creating a tabular file for protein to gene ID matching and finding the protein with longest CDS for each gene
    # grepping the GFF to only relevant lines accelerates gffread on large files immensely

    # join -t $'\t'  -j 1  <(grep -P "\tCDS\t|\tgene\t|\tmRNA\t" /share/gluster/GeneSets/NCBI-Genomes/InvertebratesRefSeq/raw_dir/active/GCF_000001215.4/genomic.gff | sed 's/\tCDS\t/\tmRNA\t/g' | gffread - --table "@geneid,@id" | sort) <(grep -P "\tCDS\t|\tgene\t|\tmRNA\t" /share/gluster/GeneSets/NCBI-Genomes/InvertebratesRefSeq/raw_dir/active/GCF_000001215.4/genomic.gff | gffread - --table "@geneid,@cdslen" | sort) > "${output_path}tmp/tmp.prot_gene_matching.txt"

    # check if protein FASTA should be extracted but exists
    if [ "${extract_proteins}" = "TRUE" ]; then
        if [[ -f "${proteins_path}" ]]; then
            echo "Proteins FASTA file exists but it is set to be extracted. This process will overwrite it. Do you want to continue? (y/n)"
            read input
            if [ "$input" = "y" ]; then
                echo "Proteins FASTA file will be overwritten"
            else
                echo "Please reconsider your option for 'extract_proteins' or delete/move your proteins FASTA file before restarting"
                [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
            fi
        fi
    fi

    # check if taxonomic assignment should be performed but files exist
    if [ "${compute_tax_assignment}" = "TRUE" ]; then
        if [[ -f "${taxon_hits_lca_path}" ]]; then
            echo "LCA hit file exists but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
            read input
            if [ "$input" = "y" ]; then
                echo "LCA hit file will be overwritten"
            else
                echo "Please reconsider your option for 'compute_tax_assignment' or delete/move your LCA hit file file before restarting"
                [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
            fi
        fi
        if [[ -f "${best_taxon_hit_path}" ]]; then
            echo "Best hits file exists but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
            read input
            if [ "$input" = "y" ]; then
                echo "Best hits file will be overwritten"
            else
                echo "Please reconsider your option for 'compute_tax_assignment' or delete/move your best hit file file before restarting"
                [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
            fi
        fi
    fi

    # check if PBC file should be generated but file exists
    if [ "${compute_pbc}" = "TRUE" ]; then
        if [[ -f "${pbc_path}" ]]; then
            echo "PBC file exists but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
            read input
            if [ "$input" = "y" ]; then
                echo "PBC file will be overwritten"
            else
                echo "Please reconsider your option for 'compute_pbc' or delete/move your PBC file file before restarting"
                [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
            fi
        fi
        echo "compute per pase coverage start:"
        time0_1=`date +%s`
        bash ./additional_scripts/compute_pbc.sh "$config_path"
        time0_2=`date +%s`
        echo "compute per pase coverage end (time elapsed:" $(($time0_2-$time0_1)) "s)"
    fi


    # 2) start python script --> produces descriptive gene statistics
    echo -e "produce gene info start:"
    time1_1=`date +%s`
    # python3 produce_gene_info.py "$config_path"
    time1_2=`date +%s`
    echo "produce gene info end (time elapsed:" $(($time1_2-$time1_1)) "s)"


    # 3) start R script for PCA and clustering
    echo "perform PCA and clustering start:"
    time2_1=`date +%s`
    Rscript perform_PCA_and_clustering.R "$config_path" --verbose >> "${output_path}R_log.out"
    time2_2=`date +%s`
    echo "perform PCA and clustering end (time elapsed:" $(($time2_2-$time2_1)) "s)"


    # 4.a) create directory in MILTS results to store orphans analysis
    [[ ! -d "${output_path}taxonomic_assignment" ]] && mkdir -p "${output_path}taxonomic_assignment"

    # 4.b) extract peptide sequences from GFF and FASTA
    if [ "${extract_proteins}" = "TRUE" ]; then
        echo "retrieving peptide sequenes start:"
        time3_1=`date +%s`
        # extraction of peptide sequence for each CDS of each gene with gffread
        # length of cds is written into the header in the fasta
        gffread -S --table "@id,@cdslen" -y ${proteins_path} -g ${fasta_path} ${gff_path}
        time3_2=`date +%s`
        echo "retrieving peptide sequenes end (time elapsed:" $(($time3_2-$time3_1)) "s)"

    fi
fi


samtools faidx "${proteins_path}" -o "${output_path}tmp/tmp.proteins.fa.fai"


# 5.a) deduce taxonomic assignment based on LCA and best hit
echo "compute taxonomic assignment start:"
time5_1=`date +%s`
python3 taxonomic_assignment.py "$config_path"
time5_2=`date +%s`
echo "compute taxonomic assignment end (time elapsed:" $(($time5_2-$time5_1)) "s)"

# 5.b) plot genes with PCA coordinates and taxonomic assignment
echo "plot taxonomic assignment start:"
time6_1=`date +%s`
Rscript plotting.R "$config_path" --verbose #>> $output_path"R_log.out"
time6_2=`date +%s`
echo "plot taxonomic assignment end (time elapsed:" $(($time6_2-$time6_1)) "s)"


# 5.c) create static plots from json files
# ${output_path}tmp/*.json not in "" so that filenames are preserved
[[ "${output_pdf}" = "TRUE" ]] && orca graph ${output_path}tmp/*.json -f "pdf" -d "${output_path}taxonomic_assignment/"
[[ "${output_png}" = "TRUE" ]] && orca graph ${output_path}tmp/*.json -f "png" -d "${output_path}taxonomic_assignment/"


# 6) remove the temporarily created files
[[ -e "${output_path}tmp/" ]] && rm -rf "${output_path}tmp/"


end=`date +%s`
runtime=$(($end-$start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Total runtime: $hours:$minutes:$seconds (hh:mm:ss)"
