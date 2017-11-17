# Phasing the Merged genotypes
In order to run LAMP-LD, we first have to have phased genotypes.  We phased genotypes using [Beagle v3.3.2](http://faculty.washington.edu/browning/beagle/b3.html) separately for each chromosome.  The Beagle phasing algorithm is described in [Browning and Browning 2007, *Am. J. Hum. Genet.*](http://www.sciencedirect.com/science/article/pii/S0002929707638828?via%3Dihub).  The phasing also included the imputation of missing genotypes.

```bash
# Make BEAGLE header file
echo -n "I ID" > beagleheader-file.txt
cat MERGED.clean.fam | \
   perl -ne 'chomp; @a=split(/\s+/,$_); print " $a[1] $a[1]"' >> beagleheader-file.txt
echo "" >> beagleheader-file.txt
   # The header file is a single line, with a space-delimited list of individual IDs (2 columns per individual)

# Set file prefix
file=MERGED.clean

# Phase for separately for each chromosome
for i in {1..38}
do
   plink --bfile $file --dog --nonfounders --chr $i --recode --transpose --out chr$i --noweb;
   mkdir Chr${i}
   mv chr${i}.* Chr${i}
   cd Chr${i}
   mv chr${i}.log chr${i}_transposed.log
   cut -d " " -f1,2,5- chr${i}.tped > chr${i}_temp
   sed -e "s/$i/M/" chr${i}_temp > chr${i}_temp2
   cat ../beagleheader-file.txt chr${i}_temp2 > chr${i}_beagle.txt
   cp ../beagle.jar .
   java -Xmx1000m -jar beagle.jar unphased=chr${i}_beagle.txt missing=0 nsamples=20 out=chr${i}_phased_beagle.txt
   gunzip chr${i}_phased_beagle.txt.chr${i}_beagle.txt.phased.gz
   cut -d" " -f2- chr${i}_phased_beagle.txt.chr${i}_beagle.txt.phased > chr${i}_phased.txt
   rm beagle.jar chr${i}_temp chr${i}_temp2 chr${i}.nosex
   cd ..
done
```

The phasing included the parameters:
- missing=0  :: this set the character for missing alleles to "0", as output by plink
- nsamples=20 :: From the manual: "a positive integer giving the number of haplotype pairs to sample for each individual during each iteration of the phasing algorithm. The nsamples argument is optional. The default value is nsamples=4. If you are phasing an extremely large sample (say > 4000 individuals), you may want to use a smaller nsamples parameter (e.g. 1 or 2) to reduce computation time. If you are phasing a small sample (say < 200 individuals), you may want to use a larger nsamples parameter (say 10 or 20) to increase accuracy."

### Separate phased haplotypes by population
Here we separated out the phased haplotypes by the three populations coded in the [canine.cluster](./Data/canine.cluster) file.  The script below utilizes a short and simple perl script found here: [extract_id_beagle.pl](./Data/extract_id_beagle.pl).  Thank you Consuelo Quinto Cortes!!!!

```bash
# Make a file listing each individual for each population
grep "Dog$" canine.cluster | \
   cut -f1-2 -d" " | \
   tr " " "-" > Dog.txt
grep "GrayWolf$" canine.cluster | \
   cut -f1-2 -d" " | \
   tr " " "-" > GrayWolf.txt
grep "MexWolf$" canine.cluster | \
   cut -f1-2 -d" " | \
   tr " " "-" > MexWolf.txt
   
### Pull out Mexican wolf phased haplotypes ###
# Set file and remove suffix
INFILE=MexWolf.txt
POP=`echo $INFILE | cut -d"." -f1`

# Loop through chromosomes
for i in {1..38}
do 
   cd Chr${i}
   perl ../extract_id_beagle.pl chr${i}_phased.txt ../$INFILE $i
   cut -d' ' -f1 chr${i}_phased.txt > rs${i}.txt
   paste -d' ' rs${i}.txt ${POP}.${i} > ${POP}_beagle.${i}.2
   LIN=$(wc -l < ${POP}_beagle.${i}.2)
   LIN2=`expr $LIN - 1` 
   echo $LIN2
   echo 'I' > temp
   for j in $(eval echo {1..$LIN2})
   do
      echo 'M' >> temp
   done
   paste -d' ' temp ${POP}_beagle.${i}.2 > ${POP}_beagle.${i}
   rm rs${i}.txt ${POP}.${i} ${POP}_beagle.${i}.2 temp
   cd ..
done

### Pull out Gray wolf phased haplotypes ###
# Set file and remove suffix
INFILE=GrayWolf.txt
POP=`echo $INFILE | cut -d"." -f1`

# Loop through chromosomes
for i in {1..38}
do 
   cd Chr${i}
   perl ../extract_id_beagle.pl chr${i}_phased.txt ../$INFILE $i
   cut -d' ' -f1 chr${i}_phased.txt > rs${i}.txt
   paste -d' ' rs${i}.txt ${POP}.${i} > ${POP}_beagle.${i}.2
   LIN=$(wc -l < ${POP}_beagle.${i}.2)
   LIN2=`expr $LIN - 1` 
   echo $LIN2
   echo 'I' > temp
   for j in $(eval echo {1..$LIN2})
   do
      echo 'M' >> temp
   done
   paste -d' ' temp ${POP}_beagle.${i}.2 > ${POP}_beagle.${i}
   rm rs${i}.txt ${POP}.${i} ${POP}_beagle.${i}.2 temp
   cd ..
done

### Pull out Dog phased haplotypes ###
# Set file and remove suffix
INFILE=Dog.txt
POP=`echo $INFILE | cut -d"." -f1`

# Loop through chromosomes
for i in {1..38}
do 
   cd Chr${i}
   perl ../extract_id_beagle.pl chr${i}_phased.txt ../$INFILE $i
   cut -d' ' -f1 chr${i}_phased.txt > rs${i}.txt
   paste -d' ' rs${i}.txt ${POP}.${i} > ${POP}_beagle.${i}.2
   LIN=$(wc -l < ${POP}_beagle.${i}.2)
   LIN2=`expr $LIN - 1` 
   echo $LIN2
   echo 'I' > temp
   for j in $(eval echo {1..$LIN2})
   do
      echo 'M' >> temp
   done
   paste -d' ' temp ${POP}_beagle.${i}.2 > ${POP}_beagle.${i}
   rm rs${i}.txt ${POP}.${i} ${POP}_beagle.${i}.2 temp
   cd ..
done
```
The final, phased haplotypes by chromosome can be found inside each folder as:
- chr\*\_phased.txt  :: phased haplotypes for the entire dataset
- MexWolf\_beagle.1  :: phased haplotypes for the Mexican Wolves only
- GrayWolf\_beagle.1  :: phased haplotypes for the Gray Wolves only
- Dog\_beagle.1  :: phased haplotypes for the Dogs only
