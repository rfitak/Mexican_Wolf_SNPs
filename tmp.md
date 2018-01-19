# Lamp-LD
This code is for running and visualizing the various results using LAMP-LD v1.3.  Lamp-LD ([Baran et al. 2012; doi: 10.1093/bioinformatics/bts144](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts144)) is a fast and accurate method for inferring local ancestry in populations from SNP data. We ran Lamp-LD using a 3-way analysis (3 possible parental populations):
- Mexican wolves as the hybrid population, then dogs, European gray wolves (EUGW), and North American gray wolves (NAGW) as the parental populations

The requirements to run this script are:
- the phased output from BEAGLE for all samples (chr\#_phased.txt)
- an R script ([geno-format.R](./Data/geno-format.R)) to convert formatting, see the script in the DATA folder
- A PLINK-formatted bed/bim/fam file of the total SNP data (MERGED.clean.{bed|bim|fam})
- A PLINK-formatted fam file of just the Mexican wolves to keep ([MW.fam](./Data/MW.fam))
- The scripts available in the [LAIT](http://www.pitt.edu/~wec47/lait.html) package
- PLINK v1.07
- LAMP-LD v1.3
```bash
# Move into folder of interest
cd /wrk/rfitak/3WAY

# Set up a loop for each chromosome
for i in {1..38}; do

   # Make a folder for each chromosome
   mkdir CHR${i}
   cd CHR${i}
   
   # Copy raw phased genotype data
   mv ../chr${i}_phased.txt .

   # Get list of SNPs for the chromosome in two formats
   grep "chr${i}_" ../MERGED.clean.bim | \
      cut -f1-4 > chr${i}.map
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
      --keep ../MW.fam \
      --recode \
      --out chr${i}_mw.unphased
   rm -rf chr${i}_mw.unphased.map chr${i}_mw.unphased.nosex chr${i}_mw.unphased.log

   # Make input files for each program then run it.
   bin=/wrk/rfitak/LAIT/bin

   # LAMPLD
   $bin/lait.pl \
      lampld \
      3 \
      chr${i}.map \
      chr${i}_mw.unphased.ped \
      chr${i}.snps \
      chr${i}_NAgw.hap \
      chr${i}_EUgw.hap \
      chr${i}_dog.hap \
      .

   /wrk/rfitak/LAIT/3WAY/LAMPLD-v1.3/bin/unolanc \
      3 \
      100 \
      50 \
      chr.pos \
      pop1.hap \
      pop2.hap \
      pop3.hap \
      genofile.gen \
      chr${i}_lampld.out

   $bin/convertLAMPLDout.pl chr${i}_lampld.out chr${i}_lampld.converted.out
   $bin/standardizeOutput.pl lampld 3 chr${i}_lampld.converted.out chr${i}_ancestry.standardized.txt
   $bin/averageAncestry.pl phased 3 chr${i}_ancestry.standardized.txt chr${i}_avg.ancestry.txt 

   # Move back into main folder
   cd ..

# Close loop
done
```

## Convert to bed format
Here we convert the output to a bed-type format (Chromosome Start End   Individual  Genotype[0,1,2]) which we can use for downstream plotting
```bash
module load gcc/4.7.2
module load bedtools/2.26.0
for c in {1..38}
   do
   cd CHR${c}
   for i in {1..176}
      do
      sed -n "$i"p chr${c}_lampld.converted.out | fold -w1 | paste chr.pos - > ${i}.tmp
      start=0
      snp=1
      while read line
         do
         stop=$(echo "${line}" | cut -f1)
         anc=$(echo "${line}" | cut -f2)
         echo -e "${c}\t${start}\t${stop}\tInd_${i}\t${anc}" >> ${i}.bed.tmp
         start=$stop
         echo "Finished chromosome ${c} haplotype ${i} SNP ${snp}"
         snp=$(( $snp + 1 ))
      done < ${i}.tmp
      rm -rf ${i}.tmp
      grep "0$" ${i}.bed.tmp | bedtools merge -i - | sed "s/$/\tInd_${i}\t0/g" >> ../tracts.bed.${c}
      grep "1$" ${i}.bed.tmp | bedtools merge -i - | sed "s/$/\tInd_${i}\t1/g" >> ../tracts.bed.${c}
      grep "2$" ${i}.bed.tmp | bedtools merge -i - | sed "s/$/\tInd_${i}\t2/g" >> ../tracts.bed.${c}
      rm -rf ${i}.bed.tmp
   done
   cd ..
done
```











# Old BED-maker code
```bash
for c in {1..38}
   do
   for i in {1..88}
      do
      haps=$(sed -n "$i"p CHR${c}/chr${c}_lampld.out | sed "s/^ //" | tr " " "\n")
      nsegs=$(echo "$haps" | wc -l)
      if [ "$nsegs" == "1" ]
         then
         snpend=$(echo "$haps" | cut -d":" -f2)
         geno=$(echo "$haps" | cut -d":" -f1)
         stop=$(sed -n "$snpend"p CHR${c}/chr.pos)
         echo -ne "${c}\t0\t${stop}\tInd_${i}\t${geno}\n" >> out.bed
      else
         start=0
         for hap in `echo "$haps"`
            do
            snp=$(echo -n "$hap" | cut -d":" -f2)
            stop=$(sed -n "$snp"p CHR${c}/chr.pos)
            geno=$(echo -n "$hap" | cut -d":" -f1)
            echo -ne "${c}\t${start}\t${stop}\tInd_${i}\t${geno}\n" >> out.bed
            start=$(( $stop + 1 ))
         done
      fi
      echo "Finished chr $c ind $i"
   done
done
```
