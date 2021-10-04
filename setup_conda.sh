#!/bin/bash

conda create -n milts python=3.6 -y
conda activate milts

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
echo "channels added"

echo "python packages"
conda install biopython scipy pyyaml taxopy=0.8.0 -y

echo "R & R packages"
conda install r-base=4.0.5 r-factoextra r-htmlwidgets=1.5.2 r-mclust r-dbscan r-plotly=4.9.2.1 bioconductor-biostrings -y
conda install -c plotly plotly-orca -y
conda install -c hcc r-paran -y


# other tools
echo "additional software"
conda install samtools bedtools bowtie2 diamond=2.0.11 gffread=0.12.1 -y

# setting MILTS.sh as an executable
chmod +x MILTS.sh

# download taxdump files for taxopy package
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxf 'taxdump.tar.gz' nodes.dmp names.dmp
rm -rf taxdump.tar.gz

echo "Done installing dependencies for MILTS"
