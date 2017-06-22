# Lamp-LD
This code is for running and visualizing the various results using LAMP-LD v1.0.  Lamp-LD ([Baran et al. 2012; doi: 10.1093/bioinformatics/bts144](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts144)) is a fast and accurate method for inferring local ancestry in populations from SNP data. We ran Lamp-LD using three different combination of populations:
- Mexican wolves as the hybrid population, then dogs and all other gray wolves as the parental populations
- Mexican wolves as the hybrid population, then dogs and only European gray wolves as the parental populations
- North American gray wolves as the hybrid population, then dogs and only European gray wolves as the parental populations

## Code to run LAMP-LD with North American gray wolves as the hybrid population and dogs and European wolves as the parental populations
```bash
# Move into folder of interest
cd /wrk/rfitak/LAIT/NA_WOLF

# Set up a loop for each chromosome
for i in {1..38}; do

   # Make a folder for each chromosome
   mkdir CHR${i}
   cd CHR${i}
   
   # Copy raw phased genotype data
   mv ../chr${i}_phased.txt .

   # Get list of SNPs for the chromosome in two formats
   grep "chr${i}_" ../SNP-master.tsv > chr${i}.map
   cut -f2 chr${i}.map > chr${i}.snps

   # Run R code to format input files
   Rscript ../geno-format.R chr${i}_phased.txt ${i}

   # Make new ped file
   plink \
      --noweb \
      --dog \
      --nonfounders \
      --bfile ../MERGED.clean \
      --chr ${i} \
      --keep ../NA_gw.fam \
      --recode \
      --out chr${i}_NA_gw.unphased
   rm -rf chr${i}_NA_gw.unphased.map chr${i}_NA_gw.unphased.nosex chr${i}_NA_gw.unphased.log

   # Make a shortcut to the program directory
   bin=/wrk/rfitak/LAIT/bin

   # Run LAIT to generate input files for LAMP-LD
   mkdir LAMPLD
   $bin/lait.pl \
      lampld 2 \
      chr${i}.map \
      chr${i}_NA_gw.unphased.ped \
      chr${i}.snps \
      chr${i}_EU_gw.hap \
      chr${i}_dog.hap \
      LAMPLD
   
   # Run LAMP-LD
   cd LAMPLD
   $bin/unolanc2way \
      100 \
      50 \
      chr.pos \
      pop1.hap \
      pop2.hap \
      genofile.gen \
      chr${i}_lampld.out
      
   # Run various LAIT programs to summarize the output
   $bin/convertLAMPLDout.pl chr${i}_lampld.out chr${i}_lampld.converted.out
   $bin/standardizeOutput.pl lampld 2 chr${i}_lampld.converted.out chr${i}_ancestry.standardized.txt
   $bin/averageAncestry.pl phased 2 chr${i}_ancestry.standardized.txt chr${i}_avg.ancestry.txt  
   
   # Move back into main folder
   cd ../..

# Close loop
done
```
