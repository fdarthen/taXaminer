#!/bin/bash

# needs: bowtie2 (or other mapping program), samtools, bedtools

# required input files:
# my_assembly.fna
# my_reads_1.fq and my_reads_2.fq (if the sequencing was single read, change bowtie2 command accordingly)
# all other files mentioned in this pipeline will be created by the respective programs
# change file names to match your file names

# call:
# ./compute_pbc.sh path/to/config.yml

PATH=$PWD/../tools:$PATH

config_path=$1
echo $config_path

# read variables from config file
source <(grep : $config_path | sed 's/ *: */=/g')

# take first pbc path in list of pbc paths
# extract first pbc_path from list in $pbc_paths
pbc_path=$(echo $pbc_paths | cut -d ',' -f1 | awk -F '[' '{print $2}' | awk -F ']' '{print $1}')
echo "PBC path: " $pbc_path

if [[ $reads == "["*","*"]" ]]; then
    reads_1=$(echo $reads | cut -d ',' -f1 | awk -F '[' '{print $2}')
    reads_2=$(echo $reads | cut -d ',' -f2 | awk -F ']' '{print $1}')
else
    reads_un=$(echo $reads | awk -F '[' '{print $2}' | awk -F '[' '{print $2}')
fi

[[ ! -d "${output_path}" ]] && mkdir -p "${output_path}"

if [[ -f "${bam}" ]]; then # if bam file exists
    echo ">>>PER BASE COVERAGE"
    echo "computing PBC from BAM file"
    # compute per base coverage
    touch ${pbc_path}
    bedtools genomecov -ibam ${bam} -d > $pbc_path
else
    cov_path="${output_path}/mapping_files/"
    [[ ! -d "${cov_path}" ]] && mkdir -p "${cov_path}"
    [[ ! -d "${output_path}/tmp" ]] && mkdir -p "${output_path}/tmp"

    echo ">>>BOWTIEBUILD"
    bowtie2-build $fasta_path ${cov_path}my_assembly

    echo ">>>MAPPING"
    if [[ -f "${reads_1}" &&  -f "${reads_2}" ]]; then
        echo "mapping paired end data"
        [ "$insert_size" == "" ] && insert_size=200
        echo "Insert size: " $((insert_size-200)) "-" $((insert_size+200))
        bowtie2 --sensitive -I $((insert_size-200)) -X $((insert_size+200)) -a -p 16 -x ${cov_path}my_assembly -1 ${reads_1} -2 ${reads_2} -S ${cov_path}my_mapping.sam
    elif [[ -f "${reads_un}" ]]; then
        echo "mapping unpaired data"
        bowtie2 --sensitive -a -p 16 -x ${cov_path}my_assembly -U ${reads_un} -S ${cov_path}my_mapping.sam
    else
        echo "Read/BAM files seem to be non-existent. Check file path(s) and restart."
    fi

    echo ">>>CONVERT"
    # convert to BAM file and sort BAM file
    samtools view -b ${cov_path}my_mapping.sam | samtools sort -o ${cov_path}my_mapping_sorted.bam

    echo ">>>INDEX"
    # index the sorted mapping
    samtools index ${cov_path}my_mapping_sorted.bam

    # index the assembly file that has been used as a reference
    samtools faidx $fasta_path -o "${output_path}tmp/tmp.MILTS.fasta.fai"

    echo ">>>PER BASE COVERAGE"
    # compute per base coverage
    touch ${pbc_path}
    bedtools genomecov -ibam ${cov_path}my_mapping_sorted.bam -d -g "${output_path}tmp/tmp.MILTS.fasta.fai" > $pbc_path

    # SAM file and temporary FASTA files will be deleted
    [[ -e "${output_path}tmp/tmp.MILTS.fasta.fai" ]] && rm "${output_path}tmp/tmp.MILTS.fasta.fai"
    [[ -e ${cov_path}"my_mapping.sam" ]] && rm ${cov_path}"my_mapping.sam"
fi
