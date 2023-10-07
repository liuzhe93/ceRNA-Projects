setwd("/Users/liuzhe/Desktop/ceRNA/07_KEGGenrichment")
remove(list=ls())

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
keytypes(org.Hs.eg.db) 

x <- c("AFF3", "ANLN", "E2F3", "EZH2", "FRMD5", "NOVA1", "PROX1", "RUNX1", "SOX11")
test = bitr(x, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
kk <- enrichKEGG(gene = test$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
write.csv(summary(kk),"KEGG-enrich.csv",row.names =FALSE)
pdf("KEGG_enrich.pdf")
barplot(kk,showCategory=14,title="Enrichment KEGG",color = "pvalue")
dev.off()
