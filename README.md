# Identification of domestic dog ancestry in Mexican wolves (*Canis lupus baileyi*)
The repository contains code for the analysis of the Mexican wolf SNP data.  Please follow the links below to get the code and associated descriptions for various analyses and data processing.  The results of this study have been published in:  

Fitak RR , Rinkevich S, and Culver M. (2018) Genome-wide analyses of SNPs is consistent with no domestic dog ancestry in the endangered Mexican wolf (*Canis lupus baileyi*). *Journal of Heredity*. **109**(4): 373-383. doi: [10.1093/jhered/esy009](https://doi.org/10.1093/jhered/esy009)  

Additional data can be found in Dryad:  

Fitak RR, Rinkevich SE, Culver M (2018) Data from: Genome-wide analysis of SNPs is consistent with no domestic dog ancestry in the endangered Mexican wolf (*Canis lupus baileyi*). Dryad Digital Repository. [https://doi.org/10.5061/dryad.g68k008](https://doi.org/10.5061/dryad.g68k008)  

Some ancillary data files and custom scripts are available in the [Data](./Data) folder.

1. [Analyzing the initial Mexican wolf genotypes](./MW-data-processing.md)
    - processing the MW genotypes alone, including data cleaning, LD pruning, PCA, heterozygosity analysis, etc
2. [Merging together the Mexican wolf genotypes with 4 other canine genotyping studies](./data-prep.md)
    - an in-depth process of merging together the CanineHD genotypes from the Mexican wolves and 4 other studies of dogs and gray wolves.
3. [Global ancestry analysis in the merged dataset](./global.md)
    - LD pruning, PCA, Admixture, and treemix analysis from the merged genotyping across all datasets
4.  [Phasing](./phasing.md)
    - Phasing and imputing missing genotypes in the merged dataset using Beagle v3.3.2
5.  [Identifying local admixture in MW with domestic dogs](./Lamp-ld.md)
    - using LAMPLD to find local haplotype assignment to parental populations
6.  [Simulating admixture](./simulations.md)
    - A series of simulations to evaluate expected local haplotype structure under various migration schemes between Mexican wolves and domestic dogs.

## All code and content herein is licensed under:
## [GNU General Public License v2.0](./LICENSE)
