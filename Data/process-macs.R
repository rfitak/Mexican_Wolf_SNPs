# Process in R
# args = genofile, posfile, chr_number, chr_length, nsnps, iteration
# e.g., args=c("genos.mat", "pos.txt", "1", "1.23e8", "6219", "1")
args = commandArgs(trailingOnly=TRUE)
library(data.table)

# Load data
data = fread(args[1])
pos = scan(args[2])


################################
################################
################################
# This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
# Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
# Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  

# NOTE: return code of -1.0 signals an error condition

SNPHWE <- function(obs_hets, obs_hom1, obs_hom2){
   if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
      return(-1.0)

   # total number of genotypes
   N <- obs_hom1 + obs_hom2 + obs_hets
   
   # rare homozygotes, common homozygotes
   obs_homr <- min(obs_hom1, obs_hom2)
   obs_homc <- max(obs_hom1, obs_hom2)

   # number of rare allele copies
   rare  <- obs_homr * 2 + obs_hets

   # Initialize probability array
   probs <- rep(0, 1 + rare)

   # Find midpoint of the distribution
   mid <- floor(rare * ( 2 * N - rare) / (2 * N))
   if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

   probs[mid + 1] <- 1.0
   mysum <- 1.0

   # Calculate probablities from midpoint down 
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr

   while ( curr_hets >=  2)
      {
      probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
      mysum <- mysum + probs[curr_hets - 1]

      # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
      curr_hets <- curr_hets - 2
      curr_homr <- curr_homr + 1
      curr_homc <- curr_homc + 1
      }    

   # Calculate probabilities from midpoint up
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr
   
   while ( curr_hets <= rare - 2)
      {
      probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
      mysum <- mysum + probs[curr_hets + 3]
         
      # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
      curr_hets <- curr_hets + 2
      curr_homr <- curr_homr - 1
      curr_homc <- curr_homc - 1
      }    
 
    # P-value calculation
    target <- probs[obs_hets + 1]

    #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

    #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

    # This assignment is the last statement in the fuction to ensure 
    # that it is used as the return value
    p <- min(1.0, sum(probs[probs <= target])/ mysum)
    }

################################
################################
################################

# Calc HWE for each column
HWE.Wigginton = function(x){
   A1 = x[seq(1, length(x), 2)]
   A2 = x[seq(2, length(x), 2)]
   geno = paste0(A1, A2)
   HOM1 = length(geno[geno == "00"])
   HOM2 = length(geno[geno == "11"])
   HET = length(geno[geno == "10" | geno == "01"])
   p = SNPHWE(HET, HOM1, HOM2)
   return(p)
}

hwe = apply(data, 2, HWE.Wigginton)
keep = which(hwe > 0.001)
print(paste0("A total of ",
   ncol(data) - length(keep),
   " SNPs were removed (HWE p-value < 0.001) and ",
   length(keep),
   " remain."))

# Make new datatable
data.hwe = data[, keep, with = FALSE]
pos.hwe = paste0("chr", as.numeric(args[3]), "_", round(as.numeric(args[4]) * pos[keep]))

# Thin to the given number of snps
thin = seq(1, length(pos.hwe), round(length(pos.hwe) / as.numeric(args[5])))
data.hwe.thin = data.hwe[, thin, with = FALSE]
data.hwe.thin[data.hwe.thin == 0] <- 2
pos.hwe.thin = pos.hwe[thin]

# Convert to tped
tped = cbind(rep(as.numeric(args[3]),
   length(pos.hwe.thin)),
   pos.hwe.thin,
   rep(0, length(pos.hwe.thin)),
   as.numeric(sub(paste0("chr", as.numeric(args[3]), "_"), "", pos.hwe.thin)),
   t(data.hwe.thin))
write.table(tped, file = paste0("chr", as.numeric(args[3]), "_", as.numeric(args[6]), ".tped"), quote = F, row.names = F, col.names = F, sep = " ")
