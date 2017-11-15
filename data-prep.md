# Downloading and preparing the input data

### Download LUPA dataset
The LUPA datasets can be found here ([LUPA data](http://dogs.genouest.org/SWEEP.dir/Supplemental.html)) and are associated with the publication by [Vayssee et al. 2011 PLoS Genetics](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002316).  This publication also describes the development of the [CanineHD beadchip](https://www.illumina.com/products/by-type/microarray-kits/caninehd.html) from Illumina.  The complete LUPA dataset contains 547 individuals (15 gray wolves + 532 domestic dogs from 48 breeds) and 174810 loci.
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


