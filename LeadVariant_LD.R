library(httr)
library(jsonlite)
library(xml2)
library(openxlsx)
library(tidyr)

## Using API to query LD from NCI. LDlink predicts LD in 1000 genomes. Using 6 population (all european)
## LDproxy calculates query variant (LeadSNP) with all variants +/- 500kb and returns them if r2 > 0.01 
server <- "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var="
r2 <- 0.5

## Remove non-biallelic varaint, variants not in 1000G, without rsNumber and disputes over allele
Rs <- c("rs4562056","rs3819405","rs2223621","rs71557345","rs12207986","rs7971")

## Population used CEU, TSI, GBR, IBS and FIN
ext <- paste(server, Rs, "&pop=CEU%2BTSI%2BFIN%2BGBR%2BIBS&r2_d=r2&token=4d405a2829c8", sep="")

results= data.frame(NULL)
for(i in 1:length(ext)){
  r <- GET(ext[i])
  #cat(content(r, type ='text/read_html', as='text', encoding ="utf-8"),"\n", file="variantsLD.txt", append=T)
  tmp = capture.output(cat(content(r, type ='text/read_html', as='text', encoding ="utf-8"),"\n"))
  tmp <- as.data.frame(tmp[-c(1,length(tmp))])
  colnames(tmp) <- "x"
  tmp = separate(tmp,x, into = c( 'RS_Number', 'Coord','Alleles','MAF','Distance','Dprime','R2','Correlated_Alleles','RegulomeDB','Function'), sep="\t")
  tmp$LeadSNP <- Rs[i]
  results = rbind(results, tmp)
}
## Remove LD calc with self
results <- results[results$Distance != 0,]
results <- results[results$R2 >= R2,]
