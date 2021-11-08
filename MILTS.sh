#!/bin/bash

# add path to locally installed tools to PATH
PATH=$PWD/tools:$PATH

start=`date +%s`

user_cofig=$1
output_path=$(cat "$user_cofig" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_path'])")

[[ ! -d "${output_path}" ]] && mkdir -p "${output_path}"
[[ ! -d "${output_path}tmp/" ]] && mkdir -p "${output_path}tmp/"

echo "make data and user input validity check start:"
time0_1=`date +%s`
python3 prepare_and_check.py "$user_cofig"
time0_2=`date +%s`
echo "make data and user input validity check  end (time elapsed:" $(($time0_2-$time0_1)) "s)"

config_path="${output_path}tmp/tmp.cfg.yml"
# read variables from config file

source <(grep : $config_path | sed 's/ *: */="/' | sed 's/$/"/g')

#TODO: print run info from python

if [ "${update_plots}" = "FALSE" ]; then

    samtools faidx "${fasta_path}" -o "${output_path}tmp/tmp.MILTS.fasta.fai"

    if [[ "${compute_coverage}" = 'TRUE' ]]; then
        echo "compute per pase coverage start:"
        time0_1=`date +%s`
        python3 prepare_coverage.py "$config_path"
        time0_2=`date +%s`
        echo "compute per pase coverage end (time elapsed:" $(($time0_2-$time0_1)) "s)"
    fi

    # 2) start python script --> produces descriptive gene statistics
    echo -e "produce gene info start:"
    time1_1=`date +%s`
    python3 produce_gene_info.py "$config_path"
    time1_2=`date +%s`
    echo "produce gene info end (time elapsed:" $(($time1_2-$time1_1)) "s)"


    # 3) start R script for PCA and clustering
    echo "perform PCA and clustering start:"
    time2_1=`date +%s`
    Rscript perform_PCA_and_clustering.R "$config_path"
    time2_2=`date +%s`
    echo "perform PCA and clustering end (time elapsed:" $(($time2_2-$time2_1)) "s)"


    # 4.b) extract peptide sequences from GFF and FASTA
    if [[ "${extract_proteins}" = "TRUE" ]]; then
        echo "extract protein FASTA start"
        time3_1=`date +%s`
        #TODO: do for other GFF rules
        # extraction of peptide sequence for each CDS of each gene with gffread
        # length of cds is written into the header in the fasta
        gffread -S --table "@id" -y ${proteins_path} -g ${fasta_path} ${gff_path}
        time3_2=`date +%s`
        echo "extract protein FASTA end (time elapsed:" $(($time3_2-$time3_1)) "s)"
    fi
fi

samtools faidx "${proteins_path}" -o "${output_path}tmp/tmp.proteins.fa.fai"

# 5.a) run DIAMOND and compute taxonomic assignment from LCA and best hit
echo "compute taxonomic assignment start:"
time5_1=`date +%s`
python3 taxonomic_assignment.py "$config_path"
time5_2=`date +%s`
echo "compute taxonomic assignment end (time elapsed:" $(($time5_2-$time5_1)) "s)"

# 5.b) plot genes with PCA coordinates and taxonomic assignment
echo "plot taxonomic assignment start:"
time6_1=`date +%s`
Rscript plotting.R "$config_path"
time6_2=`date +%s`
echo "plot taxonomic assignment end (time elapsed:" $(($time6_2-$time6_1)) "s)"


# 5.c) make html self-contained, hoverwindow text selectable and change title
echo "modify HTML start:"
time6_1=`date +%s`
python3 modify_html.py "$config_path"
time6_2=`date +%s`
echo "modify HTML end (time elapsed:" $(($time6_2-$time6_1)) "s)"


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
