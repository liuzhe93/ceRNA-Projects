setwd("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA")
remove(list=ls())
miRDB<-read.csv("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA/miRDB/inter_miRNAmRNA_miRDB.csv")
miRTarBase<-read.csv("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA/miRTarBase/inter_miRNAmRNA_miRTarBase.csv")
targetScan<-read.csv("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA/TargetScan/inter_miRNAmRNA_TargetScan.csv")
dim(miRDB)
#[1] 18098     2
dim(miRTarBase)
#[1] 191   2
dim(targetScan)
#[1] 535   2

library(tidyr)
miRDB_merged<-tidyr::unite(miRDB, "miRDB", miRNA, mRNA)
unique(miRDB_merged)
miRDB_merged<-unique(miRDB_merged)
dim(miRDB_merged)
#[1] 9606    1

miRTarBase_merged<-tidyr::unite(miRTarBase, "miRTarBase", miRNA, mRNA)
unique(miRTarBase_merged)

targetScan_merged<-tidyr::unite(targetScan, "targetScan", miRNA, mRNA)
unique(targetScan_merged)

all_inter<-c(miRDB_merged$miRDB,miRTarBase_merged$miRTarBase,targetScan_merged$targetScan)
all_inter_freq<-as.data.frame(table(all_inter))
sorted_inter<-all_inter_freq[order(all_inter_freq[,2]),]
inter_selected<-sorted_inter[sorted_inter$Freq>=2,]
dim(inter_selected)
#[1] 92  2
library(dplyr)
mydata<-inter_selected %>% separate(all_inter, c("miRNA", "mRNA"), "_")

demRNA<-read.csv("/Users/liuzhe/Desktop/ceRNA/01_DEG/mRNA/case-vs-control-diff-pval-0.05-FC-2-edgeR.gene.csv")
demRNA$mRNA<-demRNA$gene_id

results_final<-merge(demRNA,mydata,by="mRNA")
results_final<-subset(results_final,select=c("miRNA","log2FoldChange","FDR","mRNA"))
write.csv(results_final,"inter_miRNAmRNA.csv",quote=F,row.names = F)


