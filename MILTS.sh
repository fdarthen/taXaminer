#!/bin/bash

# add path to locally installed tools to PATH
PATH=$PWD/tools:$PATH

start=`date +%s`

config_path=$1

# read variables from config file
source <(grep : $config_path | sed 's/ *: */=/g')


# extract first pbc_path from list in $pbc_paths
pbc_path=$(echo $pbc_paths | cut -d ',' -f1 | awk -F '[' '{print $2}' | awk -F ']' '{print $1}')

echo "Config: " $config_path
echo "FASTA: " $fasta_path
echo "GFF: " $gff_path
echo "Proteins: " $proteins_path
echo "PBC(s): " $pbc_paths
echo "Output directory: " $output_path
echo "Taxonomic assignment: " $tax_assignment_path
if [ "${assignment_mode}" = "quick" ]; then
    echo "Quick taxonomic assignment mode (search rank: $quick_mode_search_rank; match rank: $quick_mode_match_rank)"
else
    echo "Exhaustive taxonomic assignment mode"
fi
echo "Database: " $database_path
echo "NCBI Taxon ID: " $taxon_id "(taxon excluded: $taxon_exclude)"
echo -e "\n"

[[ ! -d "${output_path}" ]] && mkdir -p "${output_path}"
[[ ! -d "${output_path}tmp/" ]] && mkdir -p "${output_path}tmp/"

if [ "${update_plots}" = "FALSE" ]; then

    samtools faidx "${fasta_path}" -o "${output_path}tmp/tmp.MILTS.fasta.fai"

    # check if protein FASTA should be extracted but exists
    if [ "${extract_proteins}" = "TRUE" ]; then
        if [[ -f "${proteins_path}" ]]; then
            echo "Proteins FASTA file exists but it is set to be created. This process will overwrite it. Do you want to continue? (y/n)"
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
        if [ "${assignment_mode}" = "quick" ]; then
            ta_path_q1=$(sed 's/txt/quick_1.txt/' <<< "$tax_assignment_path")
            ta_path_q2=$(sed 's/txt/quick_2.txt/' <<< "$tax_assignment_path")
            if [[ -f "${ta_path_q1}" && -f "${ta_path_q2}" ]]; then
                echo "Both files with quick taxonomic assignment hits exist but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
                read input
                if [ "$input" = "y" ]; then
                    echo "Quick taxonomic assignment hit files will be overwritten"
                else
                    echo "Please reconsider your option for 'compute_tax_assignment' or delete/move your quick taxonomic assignment hit files before restarting"
                    [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
                fi
            elif [[ -f "${ta_path_q1}" ]]; then
                echo "File with quick taxonomic assignment hits from first search exists but taxonomic assignment is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
                read input
                if [ "$input" = "y" ]; then
                    echo "Quick taxonomic assignment hit file 1 will be overwritten; hit file 2 is non-existant"
                else
                    echo "Please delete/move delete your quick taxonomic assignment hit file 1 before restarting to enable taxonomic assignment"
                    [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
                fi
            elif [[ -f "${ta_path_q2}" ]]; then
                echo "File with quick taxonomic assignment hits from second search exists but first one is missing. Do you wish to overwrite the second hit file? (y/n)"
                read input
                if [ "$input" = "y" ]; then
                    echo "Quick taxonomic assignment hit file 2 will be overwritten"
                else
                    echo "Please move your quick taxonomic assignment hit file 2 before restarting to enable taxonomic assignment"
                    [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
                fi
            fi
        elif [ "${assignment_mode}" = "exhaustive" ]; then
            if [[ -f "${tax_assignment_path}" ]]; then
                echo "File with taxonomic assignment hits exists but it is set to be computed. This process will overwrite it. Do you want to continue? (y/n)"
                read input
                if [ "$input" = "y" ]; then
                    echo "Taxonomic assignment hit file will be overwritten"
                else
                    echo "Please reconsider your option for 'compute_tax_assignment' or delete/move file stated in 'tax_assignment_path' before restarting"
                    [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
                fi
            fi
        else
            echo "Option for 'assignment_mode' seems to be invalid. Please reconsider you choice before restarting"
            [[ "$BASH_SOURCE" == "$0" ]] && exit 1 || return 1
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
                echo "Please reconsider your option for 'compute_pbc' or delete/move your PBC file before restarting"
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
    python3 produce_gene_info.py "$config_path"
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
        gffread -S --table "@id" -y ${proteins_path} -g ${fasta_path} ${gff_path}
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
