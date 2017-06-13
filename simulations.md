# Simulating Canine Evolutionary History
This set of code contains the scripts for simulating the demographic history of Mexican wolves, North American gray wolves, European gray wolves, and domestic dogs.  The demographic model used is that from [Fan et al. (2016, doi:10.1101/gr.197517.115)](http://genome.cshlp.org/content/26/2/163). Below is Figure 6 from Fan et al. (2016):  
![Figure 6](./Images/Fan-fig6.gif)  
Caption: Demographic model inferred using G-PhoCS. Estimates of divergence times and effective population sizes (Ne) inferred by applying a Bayesian demography inference method (G-PhoCS) to sequence data from 13,647 putative neutral loci in a subset of 22 canid genomes (because of limitations in computational power). Estimates were obtained in four separate analyses (Methods; Supplemental Table 6). Ranges of Ne are shown and correspond to 95% Bayesian credible intervals. Estimates are calibrated by assuming a per-generation mutation rate of μ = 10−8. Mean estimates (vertical lines) and ranges corresponding to 95% Bayesian credible intervals are provided at select nodes. Scales are given in units of years by assuming an average generation time of 3 yr and two different mutation rates: μ = 10−8 (dark blue) and μ = 4 × 10−9 (brown). The model also considered gene flow between different population groups (see Table 1).

After limiting this to just the four groups of interest, the tree looks like this:  
((MW),(NA_GW)),((EU_GW),(DOGS))  
- MW = Mexican wolves  
- NA_GW = North American gray wolves  
- EU_GW = European gray wolves  
- DOGS = domestic dogs

The divergence times are thus (years before present with confidence intervals):
- MW vs NA_GW = 5400 (4000 - 6600)
- EU_GW vs DOGS: 11700 (11100 - 12300)
- (MW + NA_GW) vs (EU_GW + DOGS): 12500

And ancestral population sizes are (mean and confidence interval):
- Ancestral: 45,100 (44,900 - 45,900)
- MW + NA_GW: 17,300 (13,000 - 21,700)
- MW: 600 (400-700)
- NA_GW: 3,500 (2,600 - 4,300)
- EU_GW + DOG: 8,000 (3,400 - 16,100)
- EU_GW: 3,900 - 5,300
- DOG: 1,400 - 2,700

## Step 1:  Simulate SNP data using MACS
In this step we will simulate SNP data for the 38 chromosomes in the canid genome and 50 replicates.  We will sample 176 MW, 460 NA_GW, 138 EU_GW, and 1268 DOGS chromosomes for the final output.  These sample sizes match exactly those reported in our study.  The lengths of each chromosome were retrieved from the [dog genome browser at NCBI](https://www.ncbi.nlm.nih.gov/genome?term=canis%20lupus%20familiaris).  The recombination rates for each chromosome were from the dog linkage map reported in [Wong et al. (2010, doi:10.1534/genetics.109.106831)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2828735/).  
All population sizes are scaled relative to Ne=10000.  
All times are scales to 4Ne generations, with a generation tme of 3 years (Fan et al. 2016).  
The mutation rate is set to 1x10^-8 changes/site/generation (Fan et al. 2016).  
```
   /wrk/rfitak/LAIT/macs/macs 2042 $len \
      -i 1 \
      -t 0.0004 \
      -r $rec \
      -I 4 176 460 138 1268 0 \
      -n 1 0.06 \
      -n 2 0.35 \
      -n 3 0.46 \
      -n 4 0.205 \
      -ej 0.045 1 2 \
      -en 0.045025 2 1.73 \
      -ej 0.0975 3 4 \
      -en 0.097525 4 0.8 \
      -ej 0.104167 2 4 \
      -en 0.104175 4 4.51 \
      -F ascertainment.txt 1 \
      2> errors | \
      /wrk/rfitak/LAIT/macs/msformatter > chr${chr}_${i}.macs
  
```
