# USAGE: Rscript geno-format.R <phased data file> <chromosome #>

# Grab command line arguments
args = commandArgs(trailingOnly=TRUE)

# Read input BEAGLE phased data file as a table
input = read.table(args[1], sep = " ", header=F, stringsAsFactors=F)

# Read in the PLINK-formatted fam file.
header = read.table("../MERGED.clean.fam", sep = " ", header=F, stringsAsFactors=F)

# Separate populations
genos = input[2:nrow(input),2:ncol(input)]

dog = c(113:294,767:1830,2007:2008,2023:2042)
gw = c(1:112,295:530,543:766,1831:1856)
NA_gw = c(295:530,543:766)
EU_gw = c(1:112,1831:1856)
mw = c(531:542, 1857:2006, 2009:2022)

dog2 = dog[seq(from = 2, to = length(dog), by = 2)] / 2
gw2 = gw[seq(from = 2, to = length(gw), by = 2)] / 2
NA_gw2 = NA_gw[seq(from = 2, to = length(NA_gw), by = 2)] / 2
EU_gw2 = EU_gw[seq(from = 2, to = length(EU_gw), by = 2)] / 2
mw2 = mw[seq(from = 2, to = length(mw), by = 2)] / 2

mw.geno = t(genos[, mw])
gw.geno = t(genos[, gw])
NA_gw.geno = t(genos[, NA_gw])
EU_gw.geno = t(genos[, EU_gw])
dog.geno = t(genos[, dog])
all.geno = t(genos[, c(mw, NA_gw, EU_gw, dog)])

# Make PLINK-formatted ped files
mw.ped = vector()
for (ind in seq(from = 1, to = nrow(mw.geno), by = 2)){
   a1 = mw.geno[ind,]
   a2 = mw.geno[ind + 1,]
   alleles = c(rbind(a1, a2))
   mw.ped = rbind(mw.ped, alleles)
   print(paste0("Finished Mexican wolf ", (ind + 1) / 2))
}
mw.ped = cbind(header[mw2,], mw.ped)

NA_gw.ped = vector()
for (ind in seq(from = 1, to = nrow(NA_gw.geno), by = 2)){
   a1 = NA_gw.geno[ind,]
   a2 = NA_gw.geno[ind + 1,]
   alleles = c(rbind(a1, a2))
   NA_gw.ped = rbind(NA_gw.ped, alleles)
   print(paste0("Finished NA gray wolf ", (ind + 1) / 2))
}
NA_gw.ped = cbind(header[NA_gw2,], NA_gw.ped)

all.ped = vector()
for (ind in seq(from = 1, to = nrow(all.geno), by = 2)){
   a1 = all.geno[ind,]
   a2 = all.geno[ind + 1,]
   alleles = c(rbind(a1, a2))
   all.ped = rbind(all.ped, alleles)
   print(paste0("Finished all canids individual ", (ind + 1) / 2))
}
all.ped = cbind(header[c(mw2, NA_gw2, EU_gw2, dog2),], all.ped)

# Make .hap files (1 row per haplotype)
dog.hap = apply(dog.geno, 1, paste, collapse = "") 
EU_gw.hap = apply(EU_gw.geno, 1, paste, collapse = "")
gw.hap = apply(gw.geno, 1, paste, collapse = "")

# Write output files
write.table(NA_gw.ped, file=paste0("chr", args[2], "_NA_gw.ped"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(mw.ped, file=paste0("chr", args[2], "_mw.ped"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(all.ped, file=paste0("chr", args[2], "_all.ped"), quote = F, sep = "\t", row.names = F, col.names = F)
write(gw.hap, file=paste0("chr", args[2], "_gw.hap"), ncolumns = 1)
write(dog.hap, file=paste0("chr", args[2], "_dog.hap"), ncolumns = 1)
write(EU_gw.hap, file=paste0("chr", args[2], "_EU_gw.hap"), ncolumns = 1)