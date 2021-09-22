#!/bin/bash

# download taxdump files for taxopy package
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxf 'taxdump.tar.gz' nodes.dmp names.dmp
rm -rf taxdump.tar.gz

# make MILTS script executable
chmod +x MILTS.sh
[[ ! -d tools ]] && mkdir -p tools
cd tools

#samtools
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar xjf samtools-1.11.tar.bz2
cd samtools-1.11/
./configure --prefix=$PWD
make

cd ..
ln -s samtools-1.11/samtools samtools


# bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod a+x bedtools


#Python packages
#standard library includes: itertools + sys

python3 -m pip install biopython # numpy + Bio.Seq
python3 -m pip install scipy
python3 -m pip install pyyaml


#R packages

Rscript -e 'install.packages("factoextra", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("plotly", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("paran", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("mclust", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("dbscan", repos="https://cloud.r-project.org",quiet=TRUE)'
Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org",quiet=TRUE)'


# Diamond
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.11/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz


# gffread
wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.1.Linux_x86_64.tar.gz
tar xzf gffread-0.12.1.Linux_x86_64.tar.gz
ln -s gffread-0.12.1.Linux_x86_64/gffread gffread

# orca
wget https://github.com/plotly/orca/releases/download/v1.3.1/orca-1.3.1.AppImage
chmod +x orca-1.3.1.AppImage
ln -s orca-1.3.1.AppImage orca

rm -rf gffread-0.12.1.Linux_x86_64.tar.gz diamond-linux64.tar.gz samtools-1.11.tar.bz2

cd ..
