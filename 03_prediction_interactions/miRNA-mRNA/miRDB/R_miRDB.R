setwd("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA/miRDB")
remove(list=ls())

miRDB<-read.table("miRDB_v6.0_prediction_result.txt",sep="\t")
miRDB_human<-miRDB[grep("^hsa",miRDB[,1]),]

dim(miRDB)
#[1] 6831595       3
dim(miRDB_human)
#[1] 3375741       3

miRNA<-read.table("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/lncRNA-miRNA/miRNAs.txt")
dim(miRNA)
for(i in 1:nrow(miRNA)){
  temp <- miRNA$V1[i]
  substring(temp,1:7)<-"hsa-miR"
  miRNA$miRNA_name[i]<-temp
}

miRDB_human$miRNA_name<-miRDB_human$V1
results_miRNAmRNA<-merge(miRNA,miRDB_human,by="miRNA_name")
dim(results_miRNAmRNA)


library("org.Hs.eg.db")
library("clusterProfiler")

ENTREZID<-bitr(results_miRNAmRNA[,4], fromType = "REFSEQ", toType = c("SYMBOL","ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)
dim(ENTREZID)
#[1] 45216     4
miRDB_human$REFSEQ<-miRDB_human$V2
merged<-merge(miRDB_human,ENTREZID,by="REFSEQ")
dim(merged)
#[1] 18098     8

inter_miRNAmRNA<-subset(merged,select=c("V1","SYMBOL"))
colnames(inter_miRNAmRNA)<-c("miRNA","mRNA")
write.csv(inter_miRNAmRNA,"inter_miRNAmRNA_miRDB.csv",quote=F,row.names = F)


