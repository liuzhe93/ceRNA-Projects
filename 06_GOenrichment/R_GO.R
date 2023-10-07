setwd("/Users/liuzhe/Desktop/ceRNA/06_GOenrichment")
#参考资料：https://wencke.github.io/
remove(list=ls())
#install.packages("GOplot")
library("GOplot")
library("tidyverse")
library("readxl")
MySheet<-read_excel("DAVID_mRNA_annotation.xlsx")
dim(MySheet)
david<-head(MySheet,n=5)
ec_david<-subset(david,select=c("Category","Term","Count","Genes","PValue"))
library("tidyr")
library("dplyr")
ec_david <- ec_david %>% separate(Term, c("ID","Term"), "~")
for(i in 1:nrow(ec_david)){
  temp<-ec_david[i,1]
  temp<-substring(temp,8,9)
  ec_david[i,1]<-temp
}
colnames(ec_david)<-c("category","ID","term","count","genes","adj_pval")
write.csv(ec_david,"ec_david.csv",quote=F,row.names=F)

genes<-c("ANLN", "NOVA1", "SOX11", "E2F3", "PROX1", "AFF3", "EZH2", "RUNX1")
degenes<-read.table("/Users/liuzhe/Desktop/ceRNA/01_DEG/mRNA/case-vs-control-diff-pval-0.05-FC-2-edgeR.gene.xls",sep="\t",header=T,row.names=1)
genelist<-degenes[genes,]
genelist$ID<-rownames(genelist)
genelist$logFC<-genelist$log2FoldChange

circ<-circle_dat(ec_david,genelist)
# Create the plot
gene_logfc<-subset(genelist,select=c("ID","logFC"))
chord <- chord_dat(data = circ, genes = gene_logfc, process = ec_david$term)
pdf("GOenrichment.pdf",width=16,height = 18)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
dev.off()


