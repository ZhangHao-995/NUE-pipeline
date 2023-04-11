####1. spurious associations
no_diff <- read.table("/home/zhangh/n_project/atac/diffpeak/combine/all_non_diff_peak.bed")
colnames(no_diff) <- c("chr","start","end")
no_diff_peak<-makeGRangesFromDataFrame(no_diff)
no_diff_ano<-annotatePeak(no_diff_peak,tssRegion=c(-3000,0),TxDb=txdb.cs,level="gene",assignGenomicAnnotation=T,
  genomicAnnotationPriority=c("Promoter","Exon","Intron","5UTR","3UTR","Downstream","Intergenic"),
  addFlankGeneInfo=T,flankDistance=2500,overlap="all")
no_diff_distal <- no_diff_ano %>% as.data.frame() %>% filter(annotation=="Distal Intergenic")
  %>% filter(distanceToTSS > 500000 | distanceToTSS < -500000) ###distance to TSS > 500k
ran_peak <- no_diff_distal %>% group_by(seqnames) %>% sample_n(1000) ###random 1000 for test
ran_peak <- ran_peak %>% mutate(peak_id=paste("peak_id",1:nrow(.),sep="_"))

exp2 <- exp2 %>% mutate(m=apply(.[,2:9],1,max),sd=apply(.[,2:9],1,sd))
exp2 <- exp2 %>% arrange(sd)
exp_bottom <- exp2[1:53945,]
gene_bottom <- merge(tss,exp_bottom[,1,drop=F],by="gene")
gene_ran <- gene_bottom %>% group_by(chr) %>% sample_n(1000)###random 1000 for test
spurious_pair <- crossing(gene_ran,ran_peak)
spurious_pair <- spurious_pair %>% filter(peak_chr!=gene_chr)

spurious_cpm <- spurious_count %>% transmute(peak_id,KCC=(KCC*1000000)/57773906,
  KLC=(KLC*1000000)/59412484,JCC=(JCC*1000000)/62939983,JLC=(JLC*1000000)/68390349,
  KCD=(KCD*1000000)/65684888,KLD=(KLD*1000000)/63990148,JCD=(JCD*1000000)/65745751,
  JLD=(JLD*1000000)/66374569)
spurious_cpm_log <- spurious_cpm %>% transmute(peak_id,across(KCC:JLD,~log2(.+1)))
spurious_pair_exp <- merge(ran_pair,exp2,by="gene")
  %>% pivot_longer(KCC:JLD,names_to="stage",values_to="tpm")
  %>% transmute(gene,peak_id,stage,tpm_norm=log2(tpm+1))
spurious_cpm_log <- spurious_cpm_log %>% pivot_longer(KCC:JLD,names_to="stage",values_to="cpm_norm")
spurious_pair_peak_exp <- merge(spurious_pair_exp,spurious_cpm_log,by=c("peak_id","stage"))
  %>% arrange(peak_id,gene)
spurious_pair_co <- spurious_pair_peak_exp %>% group_by(peak_id,gene) %>% summarise(cor=cor(tpm_norm,cpm_norm))%>% filter(cor!="NA")

####2.cor test
dis_exp_cor2 <- dis_exp_cor %>% transmute(peak_id,Geneid,co,co_test=pnorm(co,mean=0.01714292,sd=0.3802186,lower.tail=F))
dis_exp_cor_filter <- dis_exp_cor2 %>% filter(co_test<0.05)
