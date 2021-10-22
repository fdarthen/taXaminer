#!/bin/bash

conda create -n milts python=3.8 r-base=4.0.5 -y
source ~/anaconda3/etc/profile.d/conda.sh
conda activate milts

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
echo "channels added"

echo "install mamba"
conda install -c conda-forge mamba -y

echo "python packages"
mamba install biopython scipy pyyaml taxopy -y

echo "R & R packages"
mamba install r-ggplot2 r-factoextra r-htmlwidgets r-mclust r-dbscan r-plotly bioconductor-biostrings r-viridis orca -y
mamba install -c plotly plotly-orca -y
mamba install -c hcc r-paran -y


# other tools
echo "additional software"
mamba install samtools bedtools bowtie2 diamond gffread -y

# setting MILTS.sh as an executable
chmod +x MILTS.sh

# download taxdump files for taxopy package
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxf 'taxdump.tar.gz' nodes.dmp names.dmp
rm -rf taxdump.tar.gz

echo "Done installing dependencies for MILTS"

