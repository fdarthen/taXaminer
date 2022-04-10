#!/bin/bash

# download taxdump files for taxopy package
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxf 'taxdump.tar.gz' nodes.dmp names.dmp
rm -rf taxdump.tar.gz

# make taXaminer script executable
[[ ! -d tools ]] && mkdir -p tools
cd tools

#samtools
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
tar xjf samtools-1.14.tar.bz2
cd samtools-1.14/
./configure --prefix=$PWD
make

cd ..
ln -s samtools-1.14/samtools samtools


# bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod a+x bedtools


#Python packages
#standard library includes: itertools + sys

python3 -m pip install biopython # numpy + Bio.Seq
python3 -m pip install scipy
python3 -m pip install pyyaml
python3 -m pip install beautifulsoup4
python3 -m pip install jsmin
python3 -m pip install taxopy

#R packages

Rscript -e 'install.packages("factoextra", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("plotly", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("paran", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("mclust", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("dbscan", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org",quiet=TRUE)'

# Bowtie2
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip
unzip bowtie2-2.4.2-linux-x86_64.zip
ln -s bowtie2-2.4.2-linux-x86_64/bowtie2 bowtie2

# Diamond
wget https://github.com/bbuchfink/diamond/releases/download/v2.0.13/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

# orca
wget https://github.com/plotly/orca/releases/download/v1.3.1/orca-1.3.1.AppImage
chmod +x orca-1.3.1.AppImage
mv orca-1.3.1.AppImage orca

rm -rf diamond-linux64.tar.gz samtools-1.14.tar.bz2 bowtie2-2.4.2-linux-x86_64.zip

cd ..
