conda create -n milts python=3.7 -y
conda activate milts

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
echo "channels added"

echo "python packages"
conda install biopython scipy pyyaml -y

echo "R & R packages"
conda install r-base=4.0.3 r-factoextra r-htmlwidgets=1.5.2 bioconductor-phyloprofile r-mclust r-dbscan r-plotly=4.9.2.1 -y
conda install -c plotly plotly-orca -y
conda install -c hcc r-paran -y


# other tools
echo "additional software"
conda install samtools bedtools bowtie2 diamond gffread=0.12.1 -y

# setting MILTS.sh as an executable
chmod +x MILTS.sh

echo "Done installing dependencies for MILTS"
