library(tidyverse)
library(DESeq2)
countMat <- read.table("Peak_count.txt",header=T)
col<-data.frame(condition=c("HN","LN","HN","LN"),
cultivar=c("KN","KN","J4","J4"))
rownames(col)<-colnames(countMat)
dds<-DESeqDataSetFromMatrix(countData=countMat,
 colData=col,
 design=~cultivar+condition)
dds = DESeq(dds)
vsd<-vst(dds, blind=FALSE)
count_norm<-assay(vsd)
count_norm <- cbind(merge_peak,count_norm)
write.table(count_norm,"peak_vsd.txt",sep="\t",quote=F,row.names=F)
