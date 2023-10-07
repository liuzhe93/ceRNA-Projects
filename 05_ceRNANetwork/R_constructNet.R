setwd("/Users/liuzhe/Desktop/ceRNA/05_ceRNANetwork")
remove(list=ls())

library("tidyverse")
library("readxl")
MySheet<-read_excel("/Users/liuzhe/Desktop/ceRNA/04_lncRNALoalization/Supplementary Table 3. Intracellular localization of lncRNAs.xlsx")
dim(MySheet)
lncRNA_all<-subset(MySheet,select=c("LncRNA_GeneSymbol","Location"))
lncRNA_all<-as.data.frame(lncRNA_all)
library(dplyr)
lncRNA<-subset(lncRNA_all,lncRNA_all$Location=="Cytoplasm")
lncRNA_kept <- lncRNA$LncRNA_GeneSymbol
lncRNA_miRNA_selected <- data.frame()
lncRNA_miRNA_all<-read.csv("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/lncRNA-miRNA/SuppTab1.csv")
for(i in 1:nrow(lncRNA_miRNA_all)){
  if(lncRNA_miRNA_all[i,1]  %in% lncRNA_kept){
    temp <- lncRNA_miRNA_all[i,]
    lncRNA_miRNA_selected<-rbind(lncRNA_miRNA_selected,temp)
  }
}
dim(lncRNA_miRNA_all)
#[1] 267   5
dim(lncRNA_miRNA_selected)
#[1] 89  5
miRNA_kept<-unique(lncRNA_miRNA_selected$miRNA)
length(miRNA_kept)
#[1] 25
for(i in 1:length(miRNA_kept)){
  temp <- miRNA_kept[i]
  substring(temp,1:7) <- "hsa-miR"
  miRNA_kept[i] <- temp
}
miRNA_mRNA_selected<-data.frame()
miRNA_mRNA_all<-read.csv("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA/inter_miRNAmRNA.csv")
for(i in 1:nrow(miRNA_mRNA_all)){
  if(miRNA_mRNA_all[i,1]  %in% miRNA_kept){
    temp <- miRNA_mRNA_all[i,]
    miRNA_mRNA_selected<-rbind(miRNA_mRNA_selected,temp)
  }
}
dim(miRNA_mRNA_all)
#[1] 16   4
dim(miRNA_mRNA_selected)
#[1] 9  4
mRNA_kept<- miRNA_mRNA_selected$mRNA
length(mRNA_kept)
#[1] 9
length(sort(mRNA_kept))
#[1] 9

lncRNA_kept
#[1] "DIAPH3-AS1" "ERVMER61-1" "KCNA3"      "LINC00114"  "LINC00462"  "LINC00486"  "LINC00501"  "LMO7-AS1"  
#[9] "LSAMP-AS1"  "MIR205HG"   "MRVI1-AS1"  "PART1"      "RBMS3-AS3" 
miRNA_kept
#[1] "hsa-mir-216a" "hsa-mir-182"  "hsa-mir-96"   "hsa-mir-210"  "hsa-mir-205"  "hsa-mir-204"  "hsa-mir-21"  
#[8] "hsa-mir-217"  "hsa-mir-195"  "hsa-mir-141"  "hsa-mir-503"  "hsa-mir-122"  "hsa-mir-200a" "hsa-mir-143" 
#[15] "hsa-mir-183"  "hsa-mir-508"  "hsa-mir-363"  "hsa-mir-17"   "hsa-mir-372"  "hsa-mir-93"   "hsa-mir-145" 
#[22] "hsa-mir-187"  "hsa-mir-222"  "hsa-mir-301b" "hsa-mir-100" 
mRNA_kept
#[1] "AFF3"  "ANLN"  "E2F3"  "EZH2"  "FRMD5" "NOVA1" "PROX1" "RUNX1" "SOX11"
length(lncRNA_kept)
#[1] 13
length(miRNA_kept)
#[1] 25
length(mRNA_kept)
#[1] 9

head(lncRNA_miRNA_selected)
for(i in 1:nrow(lncRNA_miRNA_selected)){
  temp <- lncRNA_miRNA_selected[i,5]
  substring(temp,1:7) <- "hsa-miR"
  lncRNA_miRNA_selected[i,5] <- temp
}
head(lncRNA_miRNA_selected)
network_lncRNA_miRNA<-subset(lncRNA_miRNA_selected,select=c("LncRNA_GeneSymbol","miRNA"))
write.csv(network_lncRNA_miRNA,"network_lncRNA_miRNA.csv",row.names=F)
head(miRNA_mRNA_selected)
network_miRNA_mRNA<-subset(miRNA_mRNA_selected,select=c("miRNA","mRNA"))
write.csv(network_miRNA_mRNA,"network_miRNA_mRNA.csv",row.names=F)

network1<-network_lncRNA_miRNA
colnames(network1)<-c("source","targert")
network2<-network_miRNA_mRNA
colnames(network2)<-c("source","targert")
cyto_input<-rbind(network1,network2)
write.csv(cyto_input,"cyto_input.csv",row.names = F,quote=F)

lncRNA_kept<-as.data.frame(lncRNA_kept)
att_lncRNA<-as.data.frame(rep("lncRNA",dim(lncRNA_kept)[1]))
node1<-cbind(lncRNA_kept,att_lncRNA)
colnames(node1)<-c("node","attr")

miRNA_kept<-as.data.frame(miRNA_kept)
att_miRNA<-as.data.frame(rep("miRNA",dim(miRNA_kept)[1]))
node2<-cbind(miRNA_kept,att_miRNA)
colnames(node2)<-c("node","attr")

mRNA_kept<-as.data.frame(mRNA_kept)
att_mRNA<-as.data.frame(rep("mRNA",dim(mRNA_kept)[1]))
node3<-cbind(mRNA_kept,att_mRNA)
colnames(node3)<-c("node","attr")

node<-rbind(node1,node2,node3)
write.csv(node,"node_attr.csv",row.names = F,quote=F)

