# taXaminer

taXaminer enables the interactive exploration of taxonomic footprints in gene sets. The specific goal is to detect and differentiate contamination and horizontal gene transfer.

Besides the taxonomic assignment of genes, taXaminer uses a total of 16 further indicators to faciliate this. Among these indicators are read coverage, sequence composition, gene length and position of genes within their scaffold. To identify genes which deviate from the mean set of genes, a principal component analysis (PCA) is used as it condenses data to fewer dimensions. Genes with similar values for certain variables are thereby clustered together, so that deviations are made visible. The results can be interactively examined in a 3D scatterplot, where the dot position respresents a combination of coverage, sequence composition and spatial information provided by the PCA and the color the taxonomic assignment.

# Table of Contents
* [Installation](#installation)
* [Usage](#usage)
* [Bugs](#bugs)
* [Contributors](#contributors)
* [License](#license)
* [Contact](#contact)

# Installation

## Using Conda
Dependencies for *taXaminer* can be installed within a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) environment. For this, run:
```
./setup_conda.sh
```
The script will create the new environment 'taxaminer' and download and install all required dependencies - this step will take a while. *(Note: Depending on your system this sometimes fails, please check the console log for error messages concerning the dependency installation)*

## Preparation of reference database
To prepare the reference database for taxonomic assignment of the gene set, run:
```
./setup_db.sh
```
This downloads the most recent NCBI non-redundant protein database and all additional files into the working directory to create a DIAMOND formatted version of it. 

# Usage
1. Make a copy of the default config file 'config.yml' from the taXaminer directory
2. Adjust the following parameters to fit your data:
```
fasta_path: "path/to/assembly.fasta" # path to assembly FASTA
gff_path: "path/to/assembly.gff" # path to annotation in GFF3 format
output_path: "path/to/output_directory/" # directory to save results in
taxon_id: "<NCBI taxon ID>" # NCBI Taxon ID of query species
database_path: "path/to/database.dmnd" # path to database for taxonomic assignment
```
3. To include coverage information, adapt one of the following parameters (this is optional):
```
pbc_path_1: "path/to/pbc.txt" # path to PBC file; omit to use default location in output directory
bam_path_1: "path/to/mapping.bam" # path to BAM file; omit to use default location in output directory
read_paths_1: ["path/to/read_file_1.fa","path/to/read_file_2.fa"] # path to read file(s)
```
* Note: When using multiple coverage sets, duplicate the parameter you need and increase the number in the suffix

To run *taXaminer* on a SLURM cluster environment (recommended), run:
```
sbatch taxaminer.slurm <config.yml>
```
To run it locally, run:
```
python taxaminer.py <config.yml>
```

For details on additional options see [Configuration parameters](https://github.com/fdhubert/taXaminer/wiki/Configuration-parameters). 

# Bugs
Any bug reports, comments or suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/taXaminer/issues/new) or reach out via [email](mailto:freya.hubert@gmail.com).

# Contributors
* [Freya Arthen](https://github.com/fdhubert)
* Simonida Zehr

# License
*taXaminer* is released under [MIT license](https://github.com/BIONF/taXaminer/blob/master/LICENSE).

# Contact
Please contact us via [email](mailto:freya.hubert@gmail.com).
