## this is used to deal with the split veQTL outputs.
files=list.files(pattern="veQTL")
for(i in 1:length(files)){
	
    load(files[i])
	results=results[lapply(results, length)> 0]
	assign(files[i],results) 

}
rm(files, i, results)

ls()
tmp <- mget(ls())
tmp <- do.call(c,tmp)
names(tmp) <- gsub("\\w+\\.","", names(tmp))

tmp<-lapply(tmp, function(x) data.frame(matrix(x, ncol=2), ID = gsub('statistic\\.',"",names(x)[seq(1, length(x)/2)])))

veQTL_df<-do.call(rbind, tmp)
veQTL_df$SNP <- gsub("\\.\\w+", "",rownames(veQTL_df))

colnames(veQTL_df) <- c("Statistics", "p.value", "GeneID", "SNP")
rownames(veQTL_df) <- 1:nrow(veQTL_df)
veQTL_df <- veQTL_df[,c("GeneID", "SNP","Statistics", "p.value")]

## cahnge file name to respective veQTL
write.table(veQTL_df, file="veQTL_output.txt", row.names=F, quote=F)
