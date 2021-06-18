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

reads_1=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['reads_1'])")
reads_2=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['reads_2'])")
reads_un=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['reads_un'])")
bam=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['bam'])")

# take first pbc path in list of pbc paths
pbc_paths_list=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['pbc_paths'])")
pbc_path=$(echo $pbc_paths_list | cut -d ',' -f1 | awk -F '[' '{print $2}' | awk -F ']' '{print $1}' | awk -F "'" '{print $2}')
echo "PBC path: " $pbc_path
output_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_path'])")
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

    fasta_in_path=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['fasta_path'])")
    awk '/^>/{if(NR==1){print}else{printf("\n%s\n",$0)}next} {printf("%s",$0)} END{printf("\n")}' $fasta_in_path >> "${output_path}tmp/tmp.MILTS.fasta"
    fasta_path="${output_path}tmp/tmp.MILTS.fasta"

    echo ">>>BOWTIEBUILD"
    bowtie2-build $fasta_path ${cov_path}my_assembly

    echo ">>>MAPPING"
    if [[ -f "${reads_1}" &&  -f "${reads_2}" ]]; then
        echo "mapping paired end data"
        ins_size=$(cat "$config_path" | python3 -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['ins_size'])")
        echo "Insert size: " $((ins_size-200)) "-" $((ins_size+200))
        bowtie2 --sensitive -I $((ins_size-200)) -X $((ins_size+200)) -a -p 16 -x ${cov_path}my_assembly -1 ${reads_1} -2 ${reads_2} -S ${cov_path}my_mapping.sam
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
    samtools faidx $fasta_path

    echo ">>>PER BASE COVERAGE"
    # compute per base coverage
    touch ${pbc_path}
    bedtools genomecov -ibam ${cov_path}my_mapping_sorted.bam -d -g $fasta_path".fai" > $pbc_path

    # SAM file and temporary FASTA files will be deleted
    [[ -e ${fasta_path} ]] && rm ${fasta_path}
    [[ -e ${fasta_path}".fai" ]] && rm ${fasta_path}".fai"
    [[ -e ${cov_path}"my_mapping.sam" ]] && rm ${cov_path}"my_mapping.sam"
fi
echo "Done!"
