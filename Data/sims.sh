#!/bin/bash -l
# author: rfitak
#SBATCH -J sims1.%a
#SBATCH -o sims1.%a.out
#SBATCH -e sims1.%a.err
#SBATCH -p longrun
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-type=END
#SBATCH --mail-user=rfitak9@gmail.com
#SBATCH --array=1-38

cd /wrk/rfitak/MW_SIM_WITH_MIG/MIG_1
module load r-env/3.3.2
module load plink

c=$SLURM_ARRAY_TASK_ID
mkdir CHR${c}
cd CHR${c}
cp ../ascertainment.txt .
cp ../process-macs.R .

len=$(sed -n "$c"p ../chromosomes.txt | cut -d" " -f1)
nsnps=$(sed -n "$c"p ../chromosomes.txt | cut -d" " -f2)
rec=$(sed -n "$c"p ../chromosomes.txt | cut -d" " -f3)
rec=$(awk "BEGIN {print "${rec}"*40000}")

for i in {1..10}
   do
   echo "Beginning chromosome $c iteration $i"
# Run MACS simulations
   START=$(date +%s)
   /wrk/rfitak/LAIT/macs/macs 2042 $len \
      -i 1 \
      -t 0.0004 \
      -r $rec \
      -I 4 176 460 138 1268 0 \
      -n 1 0.06 \
      -n 2 0.35 \
      -n 3 0.46 \
      -n 4 0.205 \
      -em 0.00005 1 4 2000 \
      -em 0.000075 1 4 0 \
      -ej 0.045 1 2 \
      -en 0.045025 2 1.73 \
      -ej 0.0975 3 4 \
      -en 0.097525 4 0.8 \
      -ej 0.104167 2 4 \
      -en 0.104175 4 4.51 \
      -F ascertainment.txt 1 \
      2> errors | \
      /wrk/rfitak/LAIT/macs/msformatter > chr${c}_${i}.macs
   END=$(date +%s)
   DIFF=$(( $END - $START ))
   echo "Simulation $i took $DIFF seconds"
   rm -rf errors haplo* tree*

# Get genotypes and position files
   sed -n '6p' chr${c}_${i}.macs | \
      tr " " "\n" | \
      sed '1d' > chr${c}_${i}.pos
   tail -n +7 chr${c}_${i}.macs | \
      sed 's/\(.\)/\1 /g' | \
      sed "s/ $//g" > chr${c}_${i}.geno

# Process data in R
   Rscript \
      process-macs.R \
      chr${c}_${i}.geno \
      chr${c}_${i}.pos \
      $c \
      $len \
      $nsnps \
      $i
   rm -rf chr${c}_${i}.geno chr${c}_${i}.pos chr${c}_${i}.macs
      echo "Finished chromosome $c iteration $i"

# Change all alleles: 1=A, 2=G
      cut -d" " -f5- chr${c}_${i}.tped | sed -e "s/1/A/g" -e "s/2/G/g" | paste -d" " <(cut -d" " -f1-4 chr${c}_${i}.tped) - > chr${c}_${i}.tped2
      mv -f chr${c}_${i}.tped2 chr${c}_${i}.tped

# Get the MW genotypes in tped format
      cut -f1-180 -d" " chr${c}_${i}.tped > MW_${i}.tped
      
# Make a general tfam file
      for z in {1..88}; do echo "MW $z 0 0 0 0" >> MW_${i}.tfam; done


# Convert MW tped into ped format and remove unneeded files
      plink \
         --noweb \
         --nonfounders \
         --dog \
         --tfile MW_${i} \
         --recode \
         --out MW_${i}
      rm -rf MW_${i}.nosex MW_${i}.log MW_${i}.tfam MW_${i}.tped MW_${i}.map
      
      # Transpose the genoype matrix
      transpose --fsep " " -t --limit 10000x10000 chr${c}_${i}.tped | tail -n +5 | sed "s/ //g" > tmp
      
      # Get the gray wolf (GW) and DOG haplotypes
      head -774 tmp | tail -n 598 > chr${c}_${i}_gw.hap
      tail -n 1268 tmp > chr${c}_${i}_dog.hap
      
      # Get the map file (first 4 columns of tped format) and list of SNPs
      cut -f1-4 -d" " chr${c}_${i}.tped > chr${c}_${i}.map
      cut -f2 -d" " chr${c}_${i}.tped > chr${c}_${i}.snps
      
      # Remove unneeded files
      rm -rf tmp
      
      # Run lait.pl to build files for LAMP-LD
      mkdir LAMPLD_${i}
      /wrk/rfitak/LAIT/bin/lait.pl \
         lampld 2 \
         chr${c}_${i}.map \
         MW_${i}.ped \
         chr${c}_${i}.snps \
         chr${c}_${i}_gw.hap \
         chr${c}_${i}_dog.hap \
         LAMPLD_${i}
# Run LAMP-LD
cd LAMPLD_${i}
/wrk/rfitak/LAIT/bin/unolanc2way \
         100 \
         50 \
         chr.pos \
         pop1.hap \
         pop2.hap \
         genofile.gen \
         chr${c}_lampld.out
      /wrk/rfitak/LAIT/bin/convertLAMPLDout.pl chr${c}_lampld.out chr${c}_lampld.converted.out
      /wrk/rfitak/LAIT/bin/standardizeOutput.pl lampld 2 chr${c}_lampld.converted.out chr${c}_ancestry.standardized.txt
      /wrk/rfitak/LAIT/bin/averageAncestry.pl phased 2 chr${c}_ancestry.standardized.txt chr${c}_avg.ancestry.txt  
cd ..
echo "Finished chromosome $c round $i"
done

rm -rf ascertainment.txt process-macs.R
