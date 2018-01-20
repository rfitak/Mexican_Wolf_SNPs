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

# Merge results together from all chromosomes
cat tracts.bed.{1..38} > tracts.bed
rm -rf tracts.bed.{1..38}
```

## Visualizing the results using R

First, we create a summary table of the number and length of tracts from the parental populations in each individual
```R
# Read in data files
bed = read.table("tracts.bed", header = F, sep = "\t")
fam = read.table("MW.fam", sep = " ", header = F)
chr = scan("chr.lengths")

# Setup empty vectors to store results
counts.NAGW = vector()
counts.EUGW = vector()
counts.dog = vector()
mean.len.NAGW = vector()
mean.len.EUGW = vector()
mean.len.dog = vector()
q = vector()
c = 1

# Loop through each individual
for (n in seq(from = 1, to = 175, by = 2)){
   ind1 = paste0("Ind_", n)
   ind2 = paste0("Ind_", n + 1)
   frags = subset(bed, V4 == ind1 | V4 == ind2)
   frags.NAGW = subset(frags, V5 == "0")
   frags.EUGW = subset(frags, V5 == "1")
   frags.dog = subset(frags, V5 == "2")
   counts.NAGW = c(counts.NAGW, nrow(frags.NAGW))
   counts.EUGW = c(counts.EUGW, nrow(frags.EUGW))
   counts.dog = c(counts.dog, nrow(frags.dog))
   mean.len.NAGW = c(mean.len.NAGW, mean(frags.NAGW[,3] - frags.NAGW[,2]))
   mean.len.EUGW = c(mean.len.EUGW, mean(frags.EUGW[,3] - frags.EUGW[,2]))
   mean.len.dog = c(mean.len.dog, mean(frags.dog[,3] - frags.dog[,2]))
   q.NAGW = sum(as.numeric(frags.NAGW[,3] - frags.NAGW[,2])) / (2 * sum(chr))
   q.EUGW = sum(as.numeric(frags.EUGW[,3] - frags.EUGW[,2])) / (2 * sum(chr))
   q.dog = sum(as.numeric(frags.dog[,3] - frags.dog[,2])) / (2 * sum(chr))
   q = rbind(q, c(q.NAGW, q.EUGW, q.dog))
   c = c + 1
}

# Make summary data table
pops = c("CL", "CL", "MB", rep("CL",3), "MB","CL","CL","CL","CL","MB","CL","CL","CL","CL","MB","MB","CL","CL","MB","CL","CL","CL","CL","MB","CL","CL","CL","MB","MB","MB","CL","CL","CL","MB","MB","CL","CL","CL","MB","CL","CL","MB","CL","CL","CL","CL","MB","MB","MB","MB","CL","GR","MB","MB","CL","MB","MB","MB","GR","MB","MB","MB","GR","AG","GR","MB","GR","CL","CL","MB","MB","MB","CL","MB","CL","MB","CL","CL","CL","CL","MB","CL","AG","GR","CL","MB")
data = cbind(fam, pops, counts.NAGW, counts.EUGW, counts.dog, mean.len.NAGW, mean.len.EUGW, mean.len.dog, q)
colnames(data) = c("FID", "IID", "POP", "NUM_NAGW", "NUM_EUGW", "NUM_DOG", "MEAN_LENGTH_NAGW", "MEAN_LENGTH_EUGW", "MEAN_LENGTH_dog", "Q_NAGW", "Q_EUGW", "Q_dog")
write.table(data, file = "MW-tracts.summary.tsv", sep = "\t", quote = F, row.names = F)
```
### Next we can load the tab-delimited table S1 from the study and plot the distribution of dog fragments.
```R
library("Rgraphviz")
library(scales)

# Read in pedigree data from table S1
a=read.table("TableS1.tsv", header=T, sep = "\t")

# Rescale data so pie-charts are circular
new = rescale(a$MEAN_LENGTH_dog, to = c(min(a$NUM_DOG), max(a$NUM_DOG)))

# Begin Plot
pdf("Dog-segments.pdf")
plot.new()
plot.window(ylim=c(min(new),max(new)),xlim=c(min(a$NUM_DOG), max(a$NUM_DOG)))
box()
lab = rescale(c(new, 60,70,80,90,100), to = c(min(a$MEAN_LENGTH_dog), max(a$MEAN_LENGTH_dog)))[89:93]
axis(2,las=1, labels = c(3.4, 3.9, 4.3, 4.8, 5.3), at = c(60, 70, 80, 90, 100), cex.axis=1.2)
axis(1,cex.axis=1.2)
colors=c("#3771c8","#d40000","green")
for (g in 1:nrow(a)){
   pie=c(a[g,7],a[g,8],a[g,9])
   if (pie[1]==1){
      points(a[g,16],new[g], col="black", bg="#3771c8",pch=21, cex=2)
   } else if (pie[2]==1){
      points(a[g,16],new[g], col="black", bg="#d40000",pch=21,cex=2)
   } else if (pie[3]==1){
      points(a[g,16],new[g], col="black", bg="green",pch=21, cex=2)
   } else {
      pie.col=which(pie>0)
      pie=pie[pie>0]
      pieGlyph(pie, a[g,16], new[g], col=colors[pie.col], edges=200,radius = 0.75, labels=NA, border=T)
   }
}
title(xlab="Number of admixed segments", ylab="Mean length (Megabases)", cex.lab = 1.3)
dev.off()
```

### Next we will plot the overall proportion of dog ancestry per chromosome in each population
```R
# Load data files
bed = read.table("tracts.bed", header = F, sep = "\t")
fam = read.table("MW.fam", sep = " ", header = F)
chr = scan("chr.lengths")

# Set population names
pops = c("CL", "CL", "MB", rep("CL",3), "MB","CL","CL","CL","CL","MB","CL","CL","CL","CL","MB","MB","CL","CL","MB","CL","CL","CL","CL","MB","CL","CL","CL","MB","MB","MB","CL","CL","CL","MB","MB","CL","CL","CL","MB","CL","CL","MB","CL","CL","CL","CL","MB","MB","MB","MB","CL","GR","MB","MB","CL","MB","MB","MB","GR","MB","MB","MB","GR","AG","GR","MB","GR","CL","CL","MB","MB","MB","CL","MB","CL","MB","CL","CL","CL","CL","MB","CL","AG","GR","CL","MB")

# Setup empty result data frame
dog.prop = data.frame()
x = 1

# Loop through each individual
for (n in seq(from = 1, to = 175, by = 2)){
   ind1 = paste0("Ind_", n)
   ind2 = paste0("Ind_", n + 1)
   frags = subset(bed, V4 == ind1 | V4 == ind2)
   frags.dog = subset(frags, V5 == "2")
   frags.dog = cbind(frags.dog, len = frags.dog[,3] - frags.dog[,2])
   props = aggregate(frags.dog$len, list(frags.dog$V1), sum)
   for (c in 1:38){
      if (any(props$Group.1 == c)) dog.prop[x, c] = props[which(props$Group.1 == c),2] / (2 * chr[c])
      else dog.prop[x, c] = 0
   }
   x = x + 1
}

# Format final data file
data = cbind(fam, pops, dog.prop)
colnames(data) = c("FID", "IID", "POP", paste0("chr",1:38))

# Write output data
write.table(data, file = "Dog-ancestry-by-chr.tsv", sep = "\t", quote = F, row.names = F)

# Plot in ggplot
library(reshape)
library(ggplot2)

# Prep data table
data = data[,3:ncol(data)]
data = melt(data)
data2 = summarySE(data, measurevar="value", groupvars=c("POP","variable"))
data2$variable = as.factor(data2$variable)

# Subset data by population
AG.data = subset(data2, POP == "AG")
GR.data = subset(data2, POP == "GR")
MB.data = subset(data2, POP == "MB")
CL.data = subset(data2, POP == "CL")

# Create X-axis labels
lab=c(1,rep("",3),5,rep("",4), 10, rep("",4), 15,rep("",4), 20, rep("",4), 25, rep("",4), 30, rep("",4), 35, rep("",3), "Total")

# Build Plots
AG.plot = ggplot(AG.data, aes(x=variable, y=value)) + 
    geom_bar(position=position_dodge(), stat="identity", fill="#d40000") +
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) + 
    coord_cartesian(ylim=c(0,0.55)) + theme_bw() + 
    theme(panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")) +
    scale_y_continuous(breaks=c(0,0.15,0.25,0.35,0.45,0.55), labels= c(0,0.15,0.25,0.35,0.45,0.55)) +
    scale_x_discrete(labels=lab) +
    theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
    xlab("Chromosome")+ylab("Proportion dog assignment") + 
    theme(legend.position="none")
GR.plot = ggplot(GR.data, aes(x=variable, y=value)) + 
    geom_bar(position=position_dodge(), stat="identity", fill="#009E73") +
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) + 
    coord_cartesian(ylim=c(0,0.55)) + theme_bw() + 
    theme(panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")) +
    scale_y_continuous(breaks=c(0,0.15,0.25,0.35,0.45,0.55), labels= c(0,0.15,0.25,0.35,0.45,0.55)) +
    scale_x_discrete(labels=lab) +
    theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
    xlab("Chromosome")+ylab("Proportion dog assignment") + 
    theme(legend.position="none")
MB.plot = ggplot(MB.data, aes(x=variable, y=value)) + 
    geom_bar(position=position_dodge(), stat="identity", fill="#3771c8") +
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) + 
    coord_cartesian(ylim=c(0,0.55)) + theme_bw() + 
    theme(panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")) +
    scale_y_continuous(breaks=c(0,0.15,0.25,0.35,0.45,0.55), labels= c(0,0.15,0.25,0.35,0.45,0.55)) +
    scale_x_discrete(labels=lab) +
    theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
    xlab("Chromosome")+ylab("Proportion dog assignment") + 
    theme(legend.position="none")
CL.plot = ggplot(CL.data, aes(x=variable, y=value)) + 
    geom_bar(position=position_dodge(), stat="identity", fill="#E69F00") +
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) + 
    coord_cartesian(ylim=c(0,0.55)) + theme_bw() + 
    theme(panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")) +
    scale_y_continuous(breaks=c(0,0.15,0.25,0.35,0.45,0.55), labels= c(0,0.15,0.25,0.35,0.45,0.55)) +
    scale_x_discrete(labels=lab) +
    theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
    xlab("Chromosome")+ylab("Proportion dog assignment") + 
    theme(legend.position="none")

# Write all plots to file
pdf("barplots.pdf", width=7, height=10)
multiplot(MB.plot, GR.plot, AG.plot, CL.plot, cols=1)
dev.off()
```
## Here are the extra summarySE and multiplot functions, load if needed.
```R
#############################################################
################# Multiple plot function ####################
#############################################################
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#############################################################
################### SummarySE function ######################
#############################################################
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
```

## Finally, painting the chromosomes
We start by making a bed file in a slightly different way from the LAMPLD-converted output
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

Now we can paint the chromosomes for each individual.
```R
# Setup the colors
colors = c("white", "gray", "black")

# Load data
bed = read.table("out.bed", sep = "\t", header = F, colClasses = c(rep("integer", 3), rep("character", 2)))
fam = read.table("MW.fam", sep = " ", header = F)
chr = scan("chr.lengths")

# Define the populations and individuals
pops = c("CL", "CL", "MB", rep("CL",3), "MB","CL","CL","CL","CL","MB","CL","CL","CL","CL","MB","MB","CL","CL","MB","CL","CL","CL","CL","MB","CL","CL","CL","MB","MB","MB","CL","CL","CL","MB","MB","CL","CL","CL","MB","CL","CL","MB","CL","CL","CL","CL","MB","MB","MB","MB","CL","GR","MB","MB","CL","MB","MB","MB","GR","MB","MB","MB","GR","AG","GR","MB","GR","CL","CL","MB","MB","MB","CL","MB","CL","MB","CL","CL","CL","CL","MB","CL","AG","GR","CL","MB")

ind = cbind(as.character(fam$V2), pops)

# Plot for each indvidual
for (w in 1:nrow(ind)){
   wolf = ind[w,1]
   pdf(paste0("INDIVIDUAL_PLOTS/", ind[w,2], "_", wolf, ".pdf"))
   par(mar=c(2,2,1,1))
   plot.new()
   plot.window(xlim = c(0,max(chr)), ylim = c(0,76), xaxs="i", yaxs="i")
   axis(1, at=10000000*c(0,2,4,6,8,10,12), labels=c(0,2,4,6,8,10,12))
   lab=c(1,rep("",3),5,rep("",4), 10, rep("",4), 15,rep("",4), 20, rep("",4), 25, rep("",4), 30, rep("",4), 35, rep("",3))
   axis(2, las=1, labels=lab, at = seq(from = 1, to = 75, by = 2))
   x = c(0,0, rep(chr[38:1], each = 2))
   y = c(0,76,76, 2 * rep(37:1,each = 2),0)
   polygon(x,y,border = "black", lwd = 0.5, col = colors[1])
   #i1 = paste0("Ind_",2 * w - 1)
   #i2 = paste0("Ind_",2 * w)
   bed.wolf = subset(bed, V4 == paste0("Ind_", w))

bed11 = bed.wolf[bed.wolf$V5 == "11",]
for (r in 1:nrow(bed11)){
   xr = c(bed11[r,2], bed11[r,2], bed11[r,3], bed11[r,3])
   yr = c((2 * bed11[r,1] - 2), (2 * bed11[r,1]), (2 * bed11[r,1]), (2 * bed11[r,1] - 2))
   polygon(xr, yr, border=NA, col=colors[2])
}

bed01 = bed.wolf[bed.wolf$V5 == "01",]
for (r in 1:nrow(bed01)){
   xr = c(bed01[r,2], bed01[r,2], bed01[r,3], bed01[r,3])
   yr = c((2 * bed01[r,1] - 2), (2 * bed01[r,1] - 1), (2 * bed01[r,1] - 1), (2 * bed01[r,1] - 2))
   polygon(xr, yr, border=NA, col=colors[2])
}

bed22 = bed.wolf[bed.wolf$V5 == "22",]
for (r in 1:nrow(bed22)){
   xr = c(bed22[r,2], bed22[r,2], bed22[r,3], bed22[r,3])
   yr = c((2 * bed22[r,1] - 2), (2 * bed22[r,1]), (2 * bed22[r,1]), (2 * bed22[r,1] - 2))
   polygon(xr, yr, border=NA, col=colors[3])
}

bed02 = bed.wolf[bed.wolf$V5 == "02",]
for (r in 1:nrow(bed02)){
   xr = c(bed02[r,2], bed02[r,2], bed02[r,3], bed02[r,3])
   yr = c((2 * bed02[r,1] - 2), (2 * bed02[r,1] - 1), (2 * bed02[r,1] - 1), (2 * bed02[r,1] - 2))
   polygon(xr, yr, border=NA, col=colors[3])
}

bed12 = bed.wolf[bed.wolf$V5 == "12",]
for (r in 1:nrow(bed12)){
   xr = c(bed12[r,2], bed12[r,2], bed12[r,3], bed12[r,3])
   yr = c((2 * bed12[r,1] - 2), (2 * bed12[r,1] - 1), (2 * bed12[r,1] - 1), (2 * bed12[r,1] - 2))
   polygon(xr, yr, border=NA, col=colors[3])
}

bed12 = bed.wolf[bed.wolf$V5 == "12",]
for (r in 1:nrow(bed12)){
   xr = c(bed12[r,2], bed12[r,2], bed12[r,3], bed12[r,3])
   yr = c((2 * bed12[r,1] - 1), (2 * bed12[r,1]), (2 * bed12[r,1]), (2 * bed12[r,1] - 1))
   polygon(xr, yr, border=NA, col=colors[2])
}

segments(rep(0,75),1:75,rep(chr, each = 2)[1:75],1:75, lwd=0.5, col="black")
   wolf = sub("NK108296", "1033", wolf)
   wolf = sub("NK108404", "1139", wolf)
   wolf = sub("NK108445", "921", wolf)
   wolf = sub("NK108446", "1043", wolf)
   wolf = sub("NK226615", "1052", wolf)
   wolf = sub("NK226616", "1215", wolf)
   wolf = paste0(ind[w,2], "_", wolf)
text(90000000, 66, sub("_", " ", wolf), cex = 4)
dev.off()
print(paste0("Finished wolf ", wolf))
}

```
