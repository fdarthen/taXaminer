## Introduction
MILTS enables the interactive exploration of taxonomic footprints in gene sets. The specific goal is to detect and differentiate contamination and horizontal gene transfer.

Besides the taxonomic assignment of genes, MILTS uses a total of 16 further indicators to faciliate this. Among these indicators are read coverage, sequence composition, gene length and position of genes within their scaffold. To identify genes which deviate from the mean set of genes, a principal component analysis (PCA) is used as it condenses data to fewer dimensions. Genes with similar values for certain variables are thereby clustered together, so that deviations are made visible. The results can be interactively examined in a 3D scatterplot, where the dot position respresents a combination of coverage, sequence composition and spatial information provided by the PCA and the color the taxonomic assignment.

## Documentation
Please see the [GitHub Wiki](https://github.com/fdhubert/MILTS/wiki) for further information. 
