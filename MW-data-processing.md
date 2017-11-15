# Processing and analyzing the Mexican Wolf genotype data.
1. Initially, all data was loaded into GENOMESTUDIO (see [here](https://www.illumina.com/techniques/microarrays/array-data-analysis-experimental-design/genomestudio.html))
2. Clustering was performed using the canine.egt file
3. We removed any SNP with a GenTrain Score of <0.35
    - More stringent than the recommended 0.15 by Illumina for Infinium
4. Exported the more or less raw data to PLINK format
    - MW.ped and MW.map
    - 172114 out of the 173662 SNPs remained
5. Saved "Reproducibility and Heritability" report
    - this analysed the errors among replicates and Mendelian inheritance errors among parent-child relationships.
    - [view the report here](./Data/MW-5-9-2014-Reproducibility_and_Heritability_Report.csv)
6. Converted to binary format using [PLINK v1.07](http://zzz.bwh.harvard.edu/plink/)
```bash
plink \
   --noweb \
   --nonfounders \
   --dog \
   --file MW \
   --make-bed \
   --out MW
```
7.  Modified the .fam file with the updated sample identifications
8.  Calculated Missingness statistics (genotyping rates by locus and individual
```bash
plink \
   --noweb \
   --nonfounders \
   --dog \
   --bfile MW \
   --missing \
   --out MW
	# Results:
	# MW.imiss - missingness statistics by individual
	# MW.lmiss - missingness statistics by loci
```
9.  Made a file `remove.list` containing the IDs of the 8 individuals to remove
    - 8 duplicates (the replicate with the lowest call rate for each pair was removed)
10. Made a cleaned file:
    - by removing the 8 individuals above
    - SNPs with a minimum allele frequency <0.05
    - SNPs with a genotyping rate <90%
    - Individuals with a genotyping rate <90%
    - Hardy-Weinberg Equilibrium p<0.001
    - Removed X and Y SNPs (XY\_exclude.list)
```bash
# Get list of X and Y SNPs
grep \
   -e "^39" \
   -e "^40" MW.bim | \
   cut -f2 > XY_exclude.list

# Process in PLINK
plink \
   --noweb \
   --nonfounders \
   --dog \
   --bfile MW \
   --remove remove.list \
   --exclude XY_exclude.list \
   --mind 0.1 \
   --maf 0.05 \
   --geno 0.1 \
   --hwe 0.001 \
   --make-bed \
   --out MW.clean
   # Results
	# 172114 markers to be included from [MW.bim ]
	# Reading pedigree information from [ MW.fam ] 
	# 96 individuals read from [ MW.fam ] 
	# 50 individuals with nonmissing phenotypes
	# Assuming a quantitative trait
	# Missing phenotype value is -9
	# 50 males, 46 females, and 0 of unspecified sex
	# Reading genotype bitfile from [ MW.bed ] 
	# Detected that binary PED file is v1.00 SNP-major mode
	# Reading list of SNPs to exclude [ XY_exclude.list ] ... 5532 read
	# Reading individuals to remove [ remove.list ] ... 8 read
	# 8 individuals removed with --remove option
	# Before frequency and genotyping pruning, there are 166582 SNPs
	# 29 founders and 59 non-founders found
	# Writing list of removed individuals to [ MW.clean.irem ]
	# 4 of 88 individuals removed for low genotyping ( MIND > 0.1 )
	# 1036 markers to be excluded based on HWE test ( p <= 0.001 )
	# Total genotyping rate in remaining individuals is 0.988855
	# 3432 SNPs failed missingness test ( GENO > 0.1 )
	# 101742 SNPs failed frequency test ( MAF < 0.05 )
	# After frequency and genotyping pruning, there are 62219 SNPs
	# After filtering, 48 individuals with non-missing status
	# After filtering, 43 males, 41 females, and 0 of unspecified sex
	# Writing pedigree information to [ MW.clean.fam ] 
	# Writing map (extended format) information to [ MW.clean.bim ] 
	# Writing genotype bitfile to [ MW.clean.bed ] 
	# Using (default) SNP-major mode
```


# These individuals had low call rates (<0.90)
	# WILD	1053
	# X	836
	# GR	418
	# X	1033

