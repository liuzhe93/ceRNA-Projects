setwd("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA/miRTarBase")
remove(list=ls())

library("tidyverse")
library("readxl")
MySheet<-read_excel("hsa_MTI.xlsx")
dim(MySheet)
#[1] 502652      9
miRNA<-read.table("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/lncRNA-miRNA/miRNAs.txt")
dim(miRNA)
#[1] 44  1
for(i in 1:nrow(miRNA)){
  temp<-miRNA$V1[i]
  substring(temp,1,7)<-'hsa-miR'
  miRNA$miRNA_name[i]<-temp
}
MySheet$miRNA_name<-MySheet$miRNA
merged<-merge(miRNA,MySheet,by="miRNA_name")
inter_miRNAmRNA<-subset(merged,select=c("miRNA_name","Target Gene"))
colnames(inter_miRNAmRNA)<-c("miRNA","mRNA")
write.csv(inter_miRNAmRNA,"inter_miRNAmRNA_miRTarBase.csv",quote=F,row.names = F)


