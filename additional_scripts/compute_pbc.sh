#!/bin/bash

# needs: bowtie2 (or other mapping program), samtools, bedtools

# required input files:
# my_assembly.fna
# my_reads_1.fq and my_reads_2.fq (if the sequencing was single read, change bowtie2 command accordingly)
# all other files mentioned in this pipeline will be created by the respective programs
# change file names in the config file to match your file names

# call:
# ./compute_pbc.sh

fasta_path=$(cat ../config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['fasta_path'])")
reads_1=$(cat ../config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['reads_1'])")
reads_2=$(cat ../config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['reads_2'])")
reads_un=$(cat ../config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['reads_un'])")
pbc_paths_list=$(cat ../config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['pbc_paths'])")
pbc_paths=${pbc_paths_list//[[\]]}
pbc_path=${pbc_paths%%,*}
output_path=$(cat ../config.yml | python -c "import sys, yaml; print(yaml.safe_load(sys.stdin)['output_path'])")

echo ">>>BOWTIEBUILD"
bowtie2-build $fasta_path ${output_path}my_assembly


echo ">>>MAPPING"
bowtie2 --sensitive -a -p 8 -x ${output_path}my_assembly -1 ${reads_1} -2 ${reads_2} -U ${reads_un} -S ${output_path}my_mapping.sam


echo ">>>CONVERT"
# convert to BAM file and sort BAM file
samtools view -b ${output_path}my_mapping.sam | samtools sort -o ${output_path}my_mapping_sorted.bam

echo ">>>INDEX"
# index the sorted mapping
samtools index ${output_path}my_mapping_sorted.bam

# index the assembly file that has been used as a reference
samtools faidx $fasta_path

echo ">>>PER BASE COVERAGE"
# compute per base coverage
#TODO: remove paranthesis from pbc_paths
bedtools genomecov -ibam ${output_path}my_mapping_sorted.bam -d -g $fasta_path".fai" > $pbc_path
