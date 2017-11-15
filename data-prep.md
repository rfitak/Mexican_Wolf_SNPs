# Downloading and preparing the input data from other studies to merge with the Mexican wolf genotypes

## Download and prepare LUPA dataset
The LUPA datasets can be found here ([LUPA data](http://dogs.genouest.org/SWEEP.dir/Supplemental.html)) and are associated with the publication by [Vaysse et al. 2011 PLoS Genetics](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002316).  This publication also describes the development of the [CanineHD beadchip](https://www.illumina.com/products/by-type/microarray-kits/caninehd.html) from Illumina.  The complete LUPA dataset contains 547 individuals (15 gray wolves + 532 domestic dogs from 48 breeds) and 174810 loci.
```bash
# Download LUPA genotype data
curl -O http://dogs.genouest.org/SWEEP.dir/HDselection_updated_trees.ped
curl -O http://dogs.genouest.org/SWEEP.dir/HDselection_updated_trees.map
mv HDselection_updated_trees.ped LUPA.ped
mv HDselection_updated_trees.map LUPA.map
```

The SNP IDs are formatted using the GenBank or other database assession codes.  For example, "BICF2P18660", or "TIGRP2P2672".  For ease of use and comparison across studies, convert to the format "chromosome\_position" (e.g., chr1\_12345678).

```bash
# Convert SNP ID format
cat LUPA.map | \
   perl -ne 'chomp; @a=split(/\t/,$_); print "$a[0] chr$a[0]_$a[3] 0 $a[3]\n"' | \
   tr " " "\t" > tmp
mv tmp LUPA.map
```

At this point, there are 13 SNPs that belong to chromosome 40 (Y) at the end of the file.  I manually went in with VIM or any text editor and changed the SNP name and position to "40 chr40\_1 0 1", "40 chr40\_2 0 2"... "40 chr40\_13 0 13".
Now make a list of these X and Y SNPs to be excluded later in downstream analyses.

```bash
# Make a list of XY SNPs to exclude
grep "^39" LUPA.map | \
   cut -f2 | \
   cat - <(grep "^40" LUPA.map | cut -f2) > exclude.list
```
The file `exclude.list` should contain 5744 X and Y SNP IDs.

The LUPA dataset is now ready for downstream processing

## Download and prepare Stronen genotypes
The Stronen dataset is from [Stronen et al. 2015 Ecology & Evolution](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667828/).  The dataset contains "59 unrelated wolves from four previously identified population clusters (northcentral Europe n = 32, Carpathian Mountains n = 7, Dinaric‐Balkan n = 9, Ukrainian Steppe n = 11).".  The genotype data were made available in Dryad [here](https://dx.doi.org/10.5061/dryad.p6598). Unfortunately, only genotypes from SNPs after quality filtering were made available (137978 SNP loci)

```bash
# Download genotype and individual data
curl -O http://datadryad.org/bitstream/handle/10255/dryad.96783/59_EuropeanWolves.ped
curl -O http://datadryad.org/bitstream/handle/10255/dryad.96785/59_EuropeanWolves.map
curl -O http://datadryad.org/bitstream/handle/10255/dryad.96786/Sample%20ID%20and%20locations%20for%2059%20European%20wolves.xlsx

# Clean up files (remove Windows newline characters to format for Unix)
tr "\r" "\n" < 59_EuropeanWolves.map | \
   tr -s "\n" > Stronen.map
rm 59_EuropeanWolves.map
mv 59_EuropeanWolves.ped Stronen.ped
mv Sample%20ID%20and%20locations%20for%2059%20European%20wolves.xlsx Stronen.samples.xlsx

# Convert SNP ID format as described above for LUPA data
cat Stronen.map | \
   perl -ne 'chomp; @a=split(/\t/,$_); print "$a[0] chr$a[0]_$a[3] 0 $a[3]\n"' | \
   tr " " "\t" > tmp
mv tmp Stronen.map

# NOW DON'T FORGET TO CHANGE Y CHROMOSOME FORMATS IN VIM AS MENTIONED ABOVE!!!!
   # Only 1 Y chromsome SNP in the dataset
```
The Stronen dataset is now ready for downstream processing

## Download and prepare Cronin genotypes
The Cronin dataset is from [Cronin et al. 2015 Journal of Heredity](https://academic.oup.com/jhered/article/106/1/26/882754). This dataset contains 431 samples, including 35 coyotes, 91 dogs, and 305 gray wolves (including 8 Mexican wolves).  The genotype data are available in Dryad [here](https://doi.org/10.5061/dryad.284tf).  Unfortunately, the authors also did not provide raw data and only provided genotype data for 123801 SNP loci.  Furthermore, they oddly provide the data in a CSV format rather than traditional ped/map format and had to undergo some extra processing.

```bash
# Download data
curl -O http://datadryad.org/bitstream/handle/10255/dryad.73530/data_123801SNP.csv

# Build a space-delimited fam file
paste -d" " \
   <(cut -d"," -f2 data_123801SNP.csv) \
   <(cut -d"," -f1 data_123801SNP.csv) | \
   sed '1d' | \
   perl -ne 'chomp; print "$_ 0 0 0 -9\n"' \
   > cronin.fam

# Format genotypes and combine with fam file to build a ped file
sed '1d' data_123801SNP.csv | \
   cut -d"," -f3- | \
   sed -e 's/,/ /g' -e 's/_/ /g' -e 's/?/0/g' | \
   paste -d" " <(cat cronin.fam) - > cronin.ped
rm cronin.fam

# Make a list of the SNP loci
head -1 data_123801SNP.csv | \
   tr "," "\n" | \
   sed '1d' | \
   sed '1d' | \
   sed "s/_.*$//g" | \
   sed "s/f/F/g" > loci

# For each Cronin dataset SNP, find its match in the LUPA data
c=1
while read line
do 
a=$(grep -m 1 "$line\t" LUPA.map)
if [ "$a" == "" ]
   then
   echo "$line" >> no-matches
   echo "Finished $c ... no match"
   else
   echo "Finished $c"
   echo "$a" >> cronin.map
fi
c=$(( $c + 1 ))
done < loci
rm loci

# Convert SNP ID format as described above for LUPA data
cat cronin.map | \
   perl -ne 'chomp; @a=split(/\t/,$_); print "$a[0] chr$a[0]_$a[3] 0 $a[3]\n"' | \
   tr " " "\t" > tmp
mv tmp cronin.map
```

The Cronin dataset is now ready for downstream processing

## Download and prepare Alaskan husky genotypes
The Alaskan husky data is from an association study performed in [Vernau et al. 2013 PLoS One](https://doi.org/10.1371/journal.pone.0057195).  The genotype data were not made publicly available, but upon contacting the corresponding author was sent the genotype data for 28 Alaskan huskies for 172115 SNPs (karen28.tped & karen28.tfam).  The last column of the tfam file is the case-control status, and there were 8 cases (#2) and 20 controls (#1).  Controls were often matched to the cases, including close relatives, so we only retained control individual, then further limited these to the 10 with the lowest genome-wide mean identity by descent (IBD).
```bash
# Starting with karen28.tped and karen28.tfam

# First, cleanup windows newline characters to Unix, just in case
tr "\r" "\n" < karen28.tfam | tr -s "\n" > tmp
mv tmp karen28.tfam
tr "\r" "\n" < karen28.tped | tr -s "\n" > tmp
mv tmp karen28.tped

# Reduce to only "control" individuals
grep "1$" karen28.tfam | cut -d" " -f1-2 > karen.controls
plink \
   --noweb \
   --nonfounders \
   --dog \
   --tfile karen28 \
   --keep karen.controls \
   --make-bed \
   --out karen20
   # Results: 172115 SNPs, call rate 0.987572

# Convert SNP ID format as described above for LUPA data
   # Note, this is now a BIM file and not a MAP file
cat karen20.bim | \
   perl -ne 'chomp; @a=split(/\t/,$_); print "$a[0]\tchr$a[0]_$a[3]\t0\t$a[3]\t$a[4]\t$a[5]\n"' > tmp.bim
   mv tmp.bim karen20.bim

# NOW DON'T FORGET TO CHANGE Y CHROMOSOME FORMATS IN VIM AS MENTIONED ABOVE!!!!
   # There were 10 Y chromsome SNPs in this dataset
```

The Husky dataset is now ready for downstream processing

## Now we will prepare the raw Mexican wolf dataset
We begin with the MW SNP data after it has been processed in GenomeStudio.  There are 96 individuals and 172114 SNPs.

```bash
# Make new bim file
tr "\r" "\n" < MW.bim | \
   tr -s "\n" | \
   perl -ne 'chomp; @a=split(/\t/,$_); print "$a[0] chr$a[0]_$a[3] 0 $a[3] $a[4] $a[5]\n"' | \
   tr " " "\t" | \
   sed -e "s/X/39/g" -e "s/Y/40/g" > tmp
mv tmp MW.bim

# NOW DON'T FORGET TO CHANGE Y CHROMOSOME FORMATS IN VIM AS MENTIONED ABOVE!!!!
   # There were 10 Y chromsome SNPs in this dataset
```

The MW dataset is now ready for downstream processing


## Addtional cleanup: flipping SNPs and removing individuals
In this section, we will remove individuals from files that need to be excluded and omit SNPs from the X and Y chromosomes.  Also, Illumina sends out different versions of the CanineHD beadchip overtime, occassionally flipping the reference strand for certain SNPs.  This can cause errors when trying to merge files.  Luckily, Illumina outputs the exact strand for each SNP when you open the genotypes into GenomeStudio.  I saved a version of this for the MW data as [SNP\_Table2.txt](./Data/SNP_Table2.txt).  It is a tab delimited file with the columns: Name	SNP	ILMN Strand	Customer Strand.  The other genotype files may contain SNPs that need to have their strands flipped in order to match that of the MW dataset.  I found that the studies needed to be flipped if they were on the "bottom" strand.  This was further separated as shown using the code below.

```bash
# Make list of SNPs to flip that are "bottom" and "bottom"
grep "BOT.*BOT" SNP_Table2.txt | cut -f1 > stronen2flip.txt
   # 39323 SNPs
   
# Make list of SNPs to flip that are "top" and "bottom"   
grep "TOP.*BOT" SNP_Table2.txt | cut -f1 > cronin2flip.txt
   # 47298 SNPs
```

Now we can flip the strand of SNPs as needed in the genotype files to match the MW strandedness and exclude XY SNPs.  This will also convert each file to the PLINK binary format to save disk space

```bash
# Process the LUPA data
plink \
   --noweb \
   --nonfounders \
   --dog \
   --exclude exclude.list \
   --file LUPA \
   --flip cronin2flip.txt \
   --make-bed \
   --out LUPA
    # Results: 169066/174810 SNPs, 547 individuals

# Process the Stronen data
plink \
   --noweb \
   --nonfounders \
   --dog \
   --file Stronen \
   --exclude exclude.list \
   --flip stronen2flip.txt \
   --make-bed \
   --out Stronen
    # Results: 134548/137978 SNPs, 59 individuals
    
# Process the Cronin data
plink \
   --noweb \
   --nonfounders \
   --dog \
   --file cronin \
   --exclude exclude.list \
   --flip cronin2flip.txt \
   --make-bed \
   --out cronin
    # Results: 120671/123801 SNPs, 431 individuals

# Process the Husky data
plink \
   --noweb \
   --nonfounders \
   --dog \
   --bfile karen20 \
   --exclude exclude.list \
   --flip stronen2flip.txt \
   --make-bed \
   --out karen
    # Results:  166583/172115 SNPs, 20 individuals call rate: 0.987535
```

Next, we have to make list of individuals to exclude from the various studies.  For Mexican wolves, we need to exclude the 8 replicated samples with the lowest call rate (see XXXXXX) and the 4 additional samples with call rate below the threshold (90%).  To do this, we compare the list of 96 samples (MW.fam) with the list of 84 samples from the MW-only analysis from earlier (see [here](....)).

```bash
#  Compare the two sample sets to get a list of samples to exclude
comm \
   <(sort MW.fam) \
   <(sort ../MW-ANALYSES/MW.clean.fam) | \
   cut -f1 | \
   grep -v "^$" | \
   cut -d" " -f1-2 > remove.list

# Generate new MW file that flips some SNPs, excludes XY SNPs, and removes the individuals
plink \
   --noweb \
   --nonfounders \
   --dog \
   --bfile MW \
   --exclude exclude.list \
   --flip stronen2flip.txt \
   --remove remove.list \
   --make-bed \
   --out MW.clean

# Remove unnecessary files if desired
rm *.nosex *.nof *.hh
```

For the Stronen data, we will initially include all 59 wolves so no further cleanup is needed.  For the Cronin data, we will exclude all coyotes and two Mexican wolves which overlap with our MW dataset.  Here is a summary of the Mexican wolf sampels form the Cronin dataset.

| Cronin ID | MW Studbook ID | Overlap Notes |
| :---: | :---: | :---: |
| NK=108296 | 1033 | overlap # 1033 was already removed for low genotyping in MW |
| NK=108404 | 1139 | none |
| NK=108445 | 921 | none |
| NK=108446 | 1043 | none |
| NK=108448 | 1177 | overlap - checked genotypes, 0.9970188 match rate |
| NK=226615 | 1052 | none |
| NK=226616 | 1215 | none |
| NK=226618| 1133 | overlap - checked genotypes, 0.9955988 match rate (73% match between samples) |

```bash

echo "WO_NewMexico NK108448" >> cronin.matches
echo "WO_NewMexico NK226618" >> cronin.matches

# Also add coyotes to be removed
grep "^CO_" cronin.fam | \
   cut -d" " -f1-2 >> cronin.matches

```








For the Husky data, we have to calculate the IBD and keep the 10 most distantly related.  This was done using [PLINK v1.07](http://zzz.bwh.harvard.edu/plink/).

```bash
# Get IBD for the Husky data
plink \
   --noweb \
   --nonfounders \
   --dog \
   --bfile karen \
   --maf 0.05 \
   --geno 0.1 \
   --genome \
   --out karen.IBD

# Get a list of the Husky IDs
cut -d" " -f2 karen.fam > karen.ind

# Make the file tad-delimited
tr " " "\t" < karen.IBD.genome | tr -s "\t" > tmp
mv tmp karen.IBD.genome
```

now in R

```R
# Read in data table
a = read.table("karen.IBD.genome", sep = "\t", header = T)

# Read in individuals
b = scan("karen.ind", what = "character")

# Setup empty vector to store results
ibd = vector()

# Loop through each individual to get their mean IDB with each other individual
for (i in b){
   c = rbind(subset(a, IID1 == i), subset(a, IID2 == i))
   d = mean(c$PI_HAT)
   names(d) = i
   ibd = c(ibd, d)
}

#
names(ibd[order(ibd)[11:20]])
```









