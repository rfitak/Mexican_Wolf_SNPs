# Global ancestry of the merged dataset
The global ancestry analyses include a PCA, ADMIXTURE analysis, and three population test.  Prior to doing these analyses, I also built a "cluster" file that is needed, especially for plotting.  The format for this file is simply three columns: FID   IID   POP.  I assigned all samples to three possible 'pops', MexWolf, GrayWolf, and Dog.  
You can access this file here:
[Canine cluster file](./Data/canine.cluster)

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
After LD pruning, 74876 SNPs remained.  If repeated, this number does vary slightly due to stochasticity in the analysis. Also, SNPRelate alters the fam file output, so we overwrite it using the old fam file:

```bash
cp MERGED.clean.fam MERGED.clean.pruned.fam
```

## PCA Analysis
Here we construct a PCA for the Merged dataset again using SNPRelate.  The R code below continues from where we left off.

```R
# Make sure the library is loaded
library(SNPRelate)

# PCA (limited to the retained SNPs (not in LD)
pca = snpgdsPCA(genofile, autosome.only=  FALSE, snp.id = snps.id, sample.id = samp.id)
pc.percent <- 100 * pca$eigenval[1:32]/sum(pca$eigenval, na.rm = T)

# Plot Principle components (published in the article)
pdf("Merged-PCbars.pdf", width = 9, height = 7)
barplot(pc.percent[1:20], las = 1, ylab = "Percent Variation", xlab = "Principle Component", names.arg = c(1:20), ylim = c(0,15), axis.lty = 1)
dev.off()

# Read in the canine cluster file for plotting
pops = read.table("canine.cluster", header = F)

# Plot pairwise PCA for first 4 components, only save as a pdf if needed later
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits = 2), "%", sep = "")
palette(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"))
pairs(pca$eigenvect[,1:4], col = as.factor(pops$V3), labels = lbls, pch = as.numeric(pops$V3))
dev.copy(pdf, "Merged-pairs-PCA.pdf", height = 4.5, width = 3.5)
   # Dog = green = circles
   # GrayWolf = orange =triangles
   # MexWolf = blue = crosses

# Make PCA Plot (published in the article)
pdf("MERGED.pca.pdf")
plot.new()
plot.window(xlim = c(min(pca$eigenvect[,1]), max(pca$eigenvect[,1])), ylim = c(min(pca$eigenvect[,2]), max(pca$eigenvect[,2])))
axis(1)
axis(2, las = 2)
box()
points(pca$eigenvect[,1], pca$eigenvect[,2], pch = as.numeric(pops$V3), col = as.factor(pops$V3))
title(xlab = "PC 1 - 13.2%", ylab = "PC 2 - 3.8%")
legend("topright", legend = c("Dog", "Gray Wolf", "Mexican Wolf"), col = c(1:3), pch = c(1:3), bty = 'n')
dev.off()
```

The main PCA and barplot are published in the manuscript.  The paired PCA plot for the first 4 principle components can be found below:

    - Dog = green = circles
    - GrayWolf = orange = triangles
    - MexWolf = blue = crosses

![Paired-PCA](./Data/Merged-pairs-PCA.png)


## Admixture analysis
In this section we analysed the merged dataset using the program [ADMIXTURE v1.3](https://www.genetics.ucla.edu/software/admixture/).  The paper describing the software can be found in [Alexander et al. 2009 *Genome Research*](http://genome.cshlp.org/content/early/2009/07/31/gr.094052.109). Here are the general steps to run this analysis using the pruned dataset:

1.  Supervised Admixture Analysis
    - This analysis forces Mexican wolves to assign to either a Dog or Gray wolf population
    - Requires a '.pop' file (Assignments of known individuals.
2.  Unsupervised analysis
    - First, in R we generated 100 bootstrap replicates of the 88 Mexican wolves + 88 randomly sampled dogs + 88 randomly sampled gray wolves.
    - Then we ran Admixture for each dataset with *K* (number of populations) between 1 - 10.
    - Next, we determined the best values for *K* by finding the lowest cross validation (CV) error.
    - Ran admixture with *K* = 2... best estimate of *K*

```bash
# Make the .pop file for the supervised analysis in admixture
cut -f3 -d" " canine.cluster | \
   sed "s/MexWolf/-/g" > MERGED.clean.pruned.pop

# Run the supervised admxiture analysis with k=2
admixture \
   -j2 \
   --supervised \
   --cv=10 \
   -C 0.0001 \
   -c 0.0001 \
   MERGED.clean.pruned.bed \
   2
   ```
   
#### Next, in R we generated the 100 bootstrapped datasets.
   
```R
# Read in cluster file
cluster = read.table("canine.cluster", header = F, sep = " ")
fam = rep("0", 88)

# Run bootstrap loop process
for (i in 1:100){
   gw = subset(cluster, V3 == "GrayWolf")
   gw = gw[sample(1:nrow(gw), 88, replace = F), 1:2]
   gw = cbind(fam,paste(gw[,1],gw[,2],sep = "-"))
   mw = subset(cluster, V3 == "MexWolf")[,1:2]
   mw = cbind(fam,paste(mw[,1], mw[,2], sep = "-"))
   dog = subset(cluster, V3 == "Dog")
   dog = dog[sample(1:nrow(dog), 88, replace = F), 1:2]
   dog = cbind(fam, paste(dog[,1], dog[,2], sep = "-"))
   ind = rbind(gw, dog, mw)
   write.table(ind, file = paste("ind.", i, ".sample", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")
   print(paste0("Finished ", i))
}
```

#### Run each bootstrap replicate in ADMXITURE using a SLURM array job submission script

```bash
#!/bin/bash -l
# author: rfitak
#SBATCH -J admix.%a
#SBATCH -o admix.%a.out
#SBATCH -e admix.%a.err
#SBATCH -p serial
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=END
#SBATCH --mail-user=rfitak9@gmail.com
#SBATCH --array=1-100

# Change to working directory
cd ADMIXTURE

# Make a bed|bim|fam file for each bootstrap
plink \
   --noweb \
   --dog \
   --nonfounders \
   --bfile MERGED.clean.pruned \
   --keep ind.${SLURM_ARRAY_TASK_ID}.sample \
   --make-bed \
   --out admix-input.${SLURM_ARRAY_TASK_ID}

# Remove unnecessary files
rm -rf admix-input.${SLURM_ARRAY_TASK_ID}.{log,nosex}

# Run ADMIXTURE for K = 1-10
for i in {1..10}
do
admixture \
   -j1 \
   --cv=10 \
   -C 0.0001 \
   -c 0.0001 \
   admix-input.${SLURM_ARRAY_TASK_ID}.bed \
   $i | \
   tee log${SLURM_ARRAY_TASK_ID}.k${i}.out
done
```
   
   
   
   
   
   
   
   The Supervised admixture plot can then be plotted in R
   
   ```R
levels=c("PDL","DOG","GSL","HUS","MIXED","EURO","WO_LUPA",
   "WO_BC","WO_INTAK","WO_SEAK","WO_ID","WO_MN","WO_MAT",
   "WO_WO","MB","AR","GR","X","MEX_XX")
pops=factor(scan("pops.txt", what="character"),levels=levels)
at=cumsum(summary(pops))
for (i in 2:6){
file=paste0("MERGED.clean.pruned.",i,".Q")
tbl=t(as.matrix(read.table(file)))
tbl=tbl[,order(pops)]
colors=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")
dev.new(width=12, height=7)
barplot(tbl, ylab="Ancestry", border=NA, space=0, las=2, cex.lab=1.3, col=colors)
axis(1, at=at, labels=levels)
dev.copy(pdf,paste0("barplot.",i,".pdf"))
dev.off()
dev.off()
}

   ```

