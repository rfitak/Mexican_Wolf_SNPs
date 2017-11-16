# Global ancestry of the merged dataset

## LD Prune dataset
Here we prune the 118287 SNP merged dataset for SNPs in high LD.  We use the R package (in bioconductor) [SNPRelate v0.9.18](http://bioconductor.org/packages/release/bioc/html/SNPRelate.html).  Pruning utilized a 1 MB window sliding window, and removing one from each pair of SNPs if the r > 0.5.

```R
# Load library
library(SNPRelate)

# Convet to gds file
bed.fn = "MERGED.clean.bed"
bim.fn = "MERGED.clean.bim"
fam.fn = "MERGED.clean.fam"
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "MERGED.clean.gds")

# Open GDS file for pruning
genofile = openfn.gds("MERGED.clean.gds")


# Get sample IDs
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# LD Prune 
snps2keep = snpgdsLDpruning(genofile, autosome.only = F, method = "r", slide.max.bp = 1000000, verbose = T, ld.threshold = 0.5, sample.id = samp.id)
   # 74876 SNPs are selected in total.

# Get a list of of the SNPs retained
snps.id = unlist(snps2keep)

# Save pruned dataset
snpgdsGDS2BED(genofile, bed.fn = "MERGED.clean.pruned", snp.id = snps.id)
```

## PCA Analysis
Here 
