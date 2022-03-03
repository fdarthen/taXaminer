#!/bin/bash

conda create -n taxaminer python=3.8 r-base=4.0.5 -y
source ~/anaconda3/etc/profile.d/conda.sh
conda activate taxaminer

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
echo "channels added"

echo "python packages"
conda install biopython scipy pyyaml beautifulsoup4 jsmin taxopy=0.8.0 -y

echo "R & R packages"
conda install r-ggplot2 r-factoextra r-htmlwidgets r-mclust r-dbscan r-plotly bioconductor-biostrings r-viridis orca -y
conda install -c r r-mass
conda install -c plotly plotly-orca -y
conda install -c hcc r-paran -y

# other tools
echo "additional software"
conda install samtools -y
conda install bedtools -y
conda install bowtie2 -y
conda install diamond -y

# download taxdump files for taxopy package
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxf 'taxdump.tar.gz' nodes.dmp names.dmp
rm -rf taxdump.tar.gz

echo "Done installing dependencies for taXaminer"
