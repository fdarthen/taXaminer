#!/bin/bash

PATH=$PWD/tools:$PATH

# Database
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

tar -zxf 'taxdump.tar.gz' nodes.dmp names.dmp
diamond makedb --in "nr.gz" -d "nr_taxonomy.dmnd" -p 8 --taxonmap "prot.accession2taxid.FULL.gz" --taxonnodes "nodes.dmp" --taxonnames "names.dmp"
rm -rf nr.gz prot.accession2taxid.gz nodes.dmp names.dmp taxdump.tar.gz
