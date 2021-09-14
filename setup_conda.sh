conda create -n diamond_test python=3.6 -y
conda activate diamond_test

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
echo "channels added"

echo "python packages"
conda install biopython scipy pyyaml taxopy -y



echo "R & R packages"
conda install r-base=4.0.3 r-factoextra r-htmlwidgets=1.5.2 r-mclust r-dbscan r-plotly=4.9.2.1 -y
conda install -c plotly plotly-orca -y
conda install -c hcc r-paran -y


# other tools
echo "additional software"
conda install samtools bedtools bowtie2 diamond=2.0.9 gffread=0.12.1 -y

# setting MILTS.sh as an executable
chmod +x MILTS.sh

echo "Done installing dependencies for MILTS"
