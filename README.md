# veQTL
Variable expression quantitative trait loci, an R script which identifies genes that are variably expressed between two or more genotypes. The veQTL_engine.R requires multiple edits to make it usable. Firstly, the expr (line 8) needs to be mapped to the expression data (IDs as rownames), secondly the thresholds can be estimated to limit the p-value computating. 

## Input data
Samples must be matched in order between two matrices (genotype and expression).
Requires genotype to be a 012 matrix where genotypes are represented as the count of the minor allele (-1 = no call), genotypes must be rownames. Currently can only handle rownames with unique SNP IDs and no other annotations.
For expression data this needs to be a matrix of expression values transcript ID (or other unique identifiers) as rownames.

## veQTL_engine
Currently, this computes only the Brown-Forsythe test (a robust levenes tests) on all genotypes compared to all transcripts. Any statistic meet the threshold (lines 85-90) will be retained and have p-values calcualated.
### Estimate W threshold
Use `1-pf(1:100, 3-1, N-1)` to estimate p-value for W statistic from 1-100, set N at the smallest number of called genotypes. Repeat for two genotypes e.g( `1-pf(1:100,2-1,N-1)` ), find the smallest W the meets your desired p-value.

## Output
Data is output as a list of SNP with transcript statistics and pvalues, this is messy but can use veQTL_wrangler.R to tidy up and add annotation as desired.  
