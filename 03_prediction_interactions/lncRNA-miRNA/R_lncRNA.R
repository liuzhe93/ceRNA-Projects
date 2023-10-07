setwd("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/lncRNA-miRNA")
remove(list=ls())

delncRNA<-read.table("/Users/liuzhe/Desktop/ceRNA/01_DEG/lncRNA/case-vs-control-diff-pval-0.05-FC-2-edgeR.gene.xls",header=T)
head(delncRNA)
mircode<-read.table("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/lncRNA-miRNA/miRcode/mircode_highconsfamilies.txt",header=T,sep="\t")
intersect_lncRNA<-intersect(delncRNA$gene_id,mircode$gene_symbol)
mircode_selected<-mircode[mircode$gene_symbol %in% intersect_lncRNA,]
demiRNA<-read.table("/Users/liuzhe/Desktop/ceRNA/01_DEG/miRNA/case-vs-control-diff-pval-0.05-FC-2-edgeR.gene.xls",header=T)

dim(demiRNA)
#[1] 143   6
dim(mircode_selected)
#[1] 847  12

library(tidyr)
df2 <- mircode_selected %>% as_tibble() %>% separate_rows(microrna, sep = "/")
df2$miRNA<-df2$microrna
for(i in 1:nrow(df2)){
  if(startsWith(df2$miRNA[i], "miR") == TRUE){
    temp <- df2$miRNA[i]
    substring(temp, 1, 4) <- "mir-"
    temp <- paste0("hsa-",temp)
    df2$miRNA[i]<-temp
  }else{
    temp <- paste0("hsa-mir-",df2$miRNA[i])
    df2$miRNA[i]<-temp
  }
}
demiRNA$miRNA<-demiRNA$gene_id
inter_LncRNAMiRNA<-merge(df2,demiRNA,by="miRNA")
dim(inter_LncRNAMiRNA)
# count the interactions between lncRNA and miRNA
# there are 267 interaction pairs
# the no. of lncRNA is 26
# the no. of miRNA is 28
length(inter_LncRNAMiRNA$miRNA)
#[1] 267
length(unique(inter_LncRNAMiRNA$miRNA))
#[1] 28
length(unique(inter_LncRNAMiRNA$gene_symbol))
#[1] 26
write.csv(inter_LncRNAMiRNA,"interactions_lncRNAmiRNA.csv",row.names=F)
delncRNA$gene_symbol<-delncRNA$gene_id
results<-merge(delncRNA,inter_LncRNAMiRNA,by="gene_symbol")
supple_table1<-subset(results,select=c("gene_symbol","log2FoldChange.x","log2CPM.x","PValue.x","FDR.x","up_down.x","miRNA","gene_id.x"))
library(stringr)
supple_table1$GeneID<-substring(supple_table1$'gene_id.x', 1, 15)
finalres<-cbind(supple_table1$gene_symbol,supple_table1$GeneID,supple_table1$'log2FoldChange.x',supple_table1$'FDR.x',supple_table1$miRNA)
colnames(finalres)<-c("LncRNA_GeneSymbol","LncRNA_GeneID","LncRNA_Log2FC","LncRNA_FDR","miRNA")
write.csv(finalres,"SuppTab1.csv",quote=F,row.names = F)
