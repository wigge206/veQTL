library(data.table)
library(Biobase)
args = commandArgs(trailingOnly=TRUE)

## This code performs the Brown-Forsythe Levene-type test for x SNPs for all genes
geno = read.table(args[1], header =T, row.names =1)
outfile = paste( args[1], "veQTL.R", sep ="_") 
expr = read.table("Gene_exprssion.txt", header =T, row.names =1)


## SNP and geno must be matrix with samples (in same order) as cols, genes/SNPs as rows and IDs as rownames

results<-apply(snps,1, function(group){

  group <- as.character(group)

  ## exclude uncalled genotype (Denoted with -1) and samples that are in a genotype with a least than 10 samples in total.
  gr<-names(table(group)[table(group) > 10])
  gr <- group %in% gr[gr!=-1]

  group <- group[gr]
  y <- expr[,gr]
  #tmp = rownames(group)
  #tmp = c(tmp,rownames(y))

  N <- length(group) # total obs.
  reorder <- order(group)
  group <- group[reorder]
  y <- y[,reorder]
  group <- as.factor(group)

  k <- length(levels(group)) # number of groups
  n <- tapply(group,group, FUN = length) # number of obs. per group

  ## Calc Yi_bar (rowMedians) for each transcript.
  if (length(levels(group)) == 2){
  Yi_bar1 <- rowMedians(as.matrix(y[,group==levels(group)[1]]))
  Yi_bar2 <- rowMedians(as.matrix(y[,group==levels(group)[2]]))
  Yi_bar <- rbind(Yi_bar1, Yi_bar2)
  rownames(Yi_bar) <- levels(group)
  colnames(Yi_bar) <- rownames(y)
  } else {
  Yi_bar1 <- rowMedians(as.matrix(y[,group==levels(group)[1]]))
  Yi_bar2 <- rowMedians(as.matrix(y[,group==levels(group)[2]]))
  Yi_bar3 <- rowMedians(as.matrix(y[,group==levels(group)[3]]))
  Yi_bar <- rbind(Yi_bar1, Yi_bar2,Yi_bar3)
  rownames(Yi_bar) <- levels(group)
  colnames(Yi_bar) <- rownames(y)
  }

  Zij <-abs(y - t(apply(Yi_bar,2, function(x) rep(x, n)))) # maybe imporve?

  ## Calc Zi. (rowMeans) for each abs deviation from medians for each transcript.
  if (length(levels(group)) == 2){
    Zi.1 <- rowMeans(as.matrix(Zij[,group==levels(group)[1]]))
    Zi.2 <- rowMeans(as.matrix(Zij[,group==levels(group)[2]]))
    Zi. <- rbind(Zi.1, Zi.2)
    rownames(Zi.) <- levels(group)
    colnames(Zi.) <- rownames(Zij)
  } else {
    Zi.1 <- rowMeans(as.matrix(Zij[,group==levels(group)[1]]))
    Zi.2 <- rowMeans(as.matrix(Zij[,group==levels(group)[2]]))
    Zi.3 <- rowMeans(as.matrix(Zij[,group==levels(group)[3]]))
    Zi. <- rbind(Zi.1, Zi.2,Zi.3)
    rownames(Zi.) <- levels(group)
    colnames(Zi.) <- rownames(y)
  }

  Z.. <- rowMeans(Zij)

  
  # calc test statistics for each transcript -- > W = ((N-k)/(k-1)) * (sum(n*(Zi. - Z..)^2)/sum((Zij-rep(Zi.,n))^2))
  x <- ((N-k)/(k-1))
  x2 <- rep(n, ncol(Zi.))* (Zi. - rep(Z.., each=k))^2
  x3 <-(Zij - t(apply(Zi. , 2, rep, n)))^2 ## maybe improve

  z1 <- colSums(x2)
  z2 <- rowSums(x3)

  W <-  x * z1/z2 

  ### A filer to only return W the are greater than threshold. Use 1- pf(W, k-1, N-k) to test what W (for a give k and N) would return a significant pvalue, 
  

  W_k2 <- W[k==2]
  k2=2
  p_k2 <- ifelse(W_k2 > 52, (1 - pf(W_k2, k2-1, N-k2)), NA) ##current threshold for two group comparison = 52
  W_k3 <- W[k ==3]
  k3=3
  p_k3 <- ifelse(W_k3 > 28, (1 - pf(W_k3, k3-1, N-k3)), NA) ##current threshold for two group comparison = 28
 
  k2_g <- !is.na(p_k2)
  k3_g <- !is.na(p_k3)

  p <- c(p_k2[k2_g], p_k3[k3_g])
  W <- c(W_k2[k2_g], W_k3[k3_g])

  return(c(statistic = W, pvalue = p))
})

save(results, file = outfile)
