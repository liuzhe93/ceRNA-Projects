setwd("/Users/liuzhe/Desktop/ceRNA/08_ceRNA_survival")

rm(list=ls())
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)

#读取表达矩阵
dataExp <- read.csv("/Users/liuzhe/Desktop/ceRNA/TCGA_data/mRNA_lncRNA/FPKM/TCGA_STAD_final.csv", header = T,check.names = FALSE,row.names = 1)

#读取表型信息
barcode_cancer<-read.csv("/Users/liuzhe/Desktop/ceRNA/TCGA_data/clinical/Cancer_barcode.csv")
barcode_normal<-read.csv("/Users/liuzhe/Desktop/ceRNA/TCGA_data/clinical/Normal_barcode.csv")
barcode<-rbind(barcode_cancer,barcode_normal)
dim(barcode)
#[1] 407   1
dim(barcode_cancer)
#[1] 375   1
dim(barcode_normal)
#[1] 32  1
barcode$y<-c(rep("Cancer",375),rep("Normal",32))
colnames(barcode)<-c("SampleID","SampleType")
head(barcode)
tail(barcode)
#处理两个文件
#logExp <- log2(dataExp+1)#把表达矩阵log
logExp <- dataExp
logExp<- logExp[,order(colnames(logExp))]#根据行名排序表达矩阵
AllSamples<-colnames(logExp)
length(AllSamples)
#[1] 407
allSamples_t<-as.data.frame(AllSamples)
dim(allSamples_t)
colnames(allSamples_t)<-"x"
library(dplyr)
intersect_cancer<-dplyr::intersect(allSamples_t,barcode_cancer)
intersect_normal<-dplyr::intersect(allSamples_t,barcode_normal)
dim(intersect_cancer)
#[1] 375   1
dim(intersect_normal)
#[1] 32  1
logExp_cancer <- subset(logExp, select=c(intersect_cancer$x))
logExp_normal <- subset(logExp, select=c(intersect_normal$x))
dim(logExp_cancer)
#[1] 12978   375
dim(logExp_normal)
#[1] 12978    32

logExp_process<-cbind(logExp_cancer,logExp_normal)

save(logExp_process,barcode,file = 'uni_matrix.Rdata')##把处理好的数据存好

rm(list=ls())
##载入数据
load('uni_matrix.Rdata')
lncRNA_kept<-read.csv("lncRNA_kept.csv",header=T)
miRNA_kept<-read.csv("miRNA_kept.csv",header=T)
mRNA_kept<-read.csv("mRNA_kept.csv",header=T)
selected_genes<-c(lncRNA_kept$x,miRNA_kept$x,mRNA_kept$x)
for
logExp_process["DIAPH3-AS1",]




gene_set<-read.csv("/Users/liuzhe/Desktop/Xudong/omics/01_data/MSigDB/gly_genes.csv")
gene_set$y<-rep("KEGG_GLYCEROLIPID_METABOLISM",49)
colnames(gene_set)<-c("Metagene","KeggItem")
list<- split(as.matrix(gene_set)[,1], gene_set[,2])

#uniSigExp<-read.table("/Users/liuzhe/Desktop/Xudong/omics/01_data/TCGA-LIHC/uniCox/uniSigExp.txt",header=T)
#gene_set<-as.data.frame(t(t(colnames(uniSigExp)[4:ncol(uniSigExp)])))
#gene_set$V2<-rep("KEGG_GLYCEROLIPID_METABOLISM",9)
#colnames(gene_set)<-c("Metagene","KeggItem")
#list<- split(as.matrix(gene_set)[,1], gene_set[,2])

#geneSet<-getGmt("/Users/liuzhe/Desktop/Xudong/omics/01_data/MSigDB/c2.all.v7.5.1.symbols.gmt",geneIdType = SymbolIdentifier())
gsva_matrix<- gsva(as.matrix(logExp_process), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.csv(gsva_matrix,"noNormalize.csv",quote=F)

normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}#设定normalization函数
nor_gsva_matrix1 <- normalization(gsva_matrix)
write.csv(nor_gsva_matrix1,"Normalize.csv",quote=F)


dim(nor_gsva_matrix1)
#[1]   1 390
GLMS_data<-t(nor_gsva_matrix1)
cutoff=median(nor_gsva_matrix1)
summary(GLMS_data)
#KEGG_GLYCEROLIPID_METABOLISM
#Min.   :0.0000              
#1st Qu.:0.6373              
#Median :0.7266              
#Mean   :0.7107              
#3rd Qu.:0.8150              
#Max.   :1.0000
cutoff<-0.6373

GLMS_data<-as.data.frame(GLMS_data)
for(i in 1:390){
  if(GLMS_data$KEGG_GLYCEROLIPID_METABOLISM[i]>cutoff){
    GLMS_data$GLMS[i]<-"high"
  }else{
    GLMS_data$GLMS[i]<-"low"
  }
}
table(GLMS_data$GLMS)
write.csv(GLMS_data,"GLMS_group.csv",quote=F)


