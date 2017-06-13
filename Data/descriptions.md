# Description of the files in this folder

## [chromosomes.txt](./chromosomes.txt)
This file contains one line for each dog chromosome, 1-38.  The columns are:
- length in base pairs ([dog genome browser at NCBI](https://www.ncbi.nlm.nih.gov/genome?term=canis%20lupus%20familiaris))
- number of SNPs from the observed dataset
- the recombination rate ([Wong et al. (2010, doi:10.1534/genetics.109.106831)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2828735/))

## [ascertainment.txt](./ascertainment.txt)
This file is used by MACS to filter SNPs occurring at various frequencies.  The two columns are:
- the frequency bin (starting from either 0 or the previous line if multiple lines)
- the proportion of SNPs to keep in this frequency bin (for example, in this we remove SNPs with MAF < 0.05)

## [process-macs.R](./process-macs.R)
This file takes the output files from macs, removes SNPs deviating from HWE (p<0.001), thins to the desired number of SNPs, and reformats the output into the tped format.  The R package "data.table" must be installed previously.
It requires the following arguments after "Rscript process-macs.R":
- the genotype matrix (reformatted from the macs output)
- list of the SNP positions (reformatted from the macs output)
- chromosome number
- chromosome length
- number of SNPs desired
- the iteration number
