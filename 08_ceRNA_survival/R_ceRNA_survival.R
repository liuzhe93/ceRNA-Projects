setwd("/Users/liuzhe/Desktop/ceRNA/08_ceRNA_survival")
remove(list=ls())
#参考资料：https://blog.csdn.net/weixin_45822007/article/details/120916062

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

write.csv(lncRNA_kept,"lncRNA_kept.csv",quote=F,row.names = F)
write.csv(miRNA_kept,"miRNA_kept.csv",quote=F,row.names = F)
write.csv(mRNA_kept,"mRNA_kept.csv",quote=F,row.names = F)
selected_genes<-c(lncRNA_kept$x,miRNA_kept$x,mRNA_kept$x)
write.csv(selected_genes,"genes_ceRNA.csv",quote = F,row.names = F)


lncRNA_kept<-read.csv("lncRNA_kept.csv")
miRNA_kept<-read.csv("miRNA_kept.csv")
mRNA_kept<-read.csv("mRNA_kept.csv")

for(i in 1:nrow(miRNA_kept)){
  temp<-miRNA_kept$x[i]
  substring(temp,1:7)<-"hsa-mir"
  miRNA_kept$x[i]<-temp
}
miRNA_exp<-read.csv("/Users/liuzhe/Desktop/ceRNA/TCGA_data/miRNA/RPM_STAD_miRNA.csv")
miRNA_exp$x<-miRNA_exp$miRNA_ID
miRNA_export<-merge(miRNA_kept,miRNA_exp,by="x")
miRNA_export<-miRNA_export[,-2]
row.names(miRNA_export)<-miRNA_export$x
miRNA_export<-miRNA_export[,-1]
miRNA_export_t<-t(miRNA_export)
dim(miRNA_export_t)
#[1] 491  25
miRNA_export_t<-as.data.frame(miRNA_export_t)
miRNA_export_t$Sample_Name<-row.names(miRNA_export_t)

lncRNA_exp<-read.csv("/Users/liuzhe/Desktop/ceRNA/01_DEG/lncRNA/lncRNA_exp_uniq.csv")
lncRNA_exp$x<-lncRNA_exp$X
lncRNA_export<-merge(lncRNA_kept,lncRNA_exp,by="x")
lncRNA_export<-lncRNA_export[,-2]
row.names(lncRNA_export)<-lncRNA_export$x
lncRNA_export<-lncRNA_export[,-1]
lncRNA_export_t<-t(lncRNA_export)
dim(lncRNA_export_t)
#[1] 407  13
lncRNA_export_t<-as.data.frame(lncRNA_export_t)
lncRNA_export_t$Sample_Name<-row.names(lncRNA_export_t)

mRNA_exp<-read.csv("/Users/liuzhe/Desktop/ceRNA/01_DEG/mRNA/mRNA_exp_uniq.csv")
mRNA_exp$x<-mRNA_exp$X
mRNA_export<-merge(mRNA_kept,mRNA_exp,by="x")
mRNA_export<-mRNA_export[,-2]
row.names(mRNA_export)<-mRNA_export$x
mRNA_export<-mRNA_export[,-1]
mRNA_export_t<-t(mRNA_export)
dim(mRNA_export_t)
#[1] 407   9
mRNA_export_t<-as.data.frame(mRNA_export_t)
mRNA_export_t$Sample_Name<-row.names(mRNA_export_t)

merged1<-merge(lncRNA_export_t,mRNA_export_t,by="Sample_Name")
merged2<-merge(merged1,miRNA_export_t,by="Sample_Name")

for(i in 1:nrow(merged2)){
  temp<-merged2$Sample_Name[i]
  temp<-substr(temp,1,12)
  temp<-gsub("\\.","-",temp)
  merged2$Sample_Name[i]<-temp
}

clinical<-read.table("/Users/liuzhe/Desktop/ceRNA/TCGA_data/clinical/clinical_rmNA.txt",header=T)
clinical$Sample_Name<-clinical$Id

matrix<-merge(merged2,clinical,by="Sample_Name")
dim(matrix)
#[1] 160  61

names_change<-colnames(matrix)
names_change<-gsub("-","_",names_change)
colnames(matrix)<-names_change
uniSigExp<-subset(matrix,select=c("Id","futime","fustat","LMO7_AS1","MRVI1_AS1","hsa_mir_216a","hsa_mir_96",
                                  "hsa_mir_217","hsa_mir_122","hsa_mir_508","hsa_mir_363","hsa_mir_17",
                                  "hsa_mir_372","hsa_mir_93","hsa_mir_145","hsa_mir_187","hsa_mir_100",
                                  "AFF3","EZH2","NOVA1","PROX1"))
write.csv(uniSigExp,"uniSigExp.csv",quote=F,row.names = F)

library("survival")
library("survminer")

allgenes<-c(lncRNA_kept$x,mRNA_kept$x,miRNA_kept$x)
allgenes<-gsub("-","_",allgenes)

#计算最佳截点
res.cut<-surv_cutpoint(matrix,time="futime",event="fustat",variables=allgenes)
summary(res.cut)
plot(res.cut,"KCNA3",palette="npg")
res.cat<-surv_categorize(res.cut)
head(res.cat)

allgenes
#[1] "DIAPH3_AS1"   "ERVMER61_1"   "KCNA3"        "LINC00114"    "LINC00462"    "LINC00486"    "LINC00501"   
#[8] "LMO7_AS1"     "LSAMP_AS1"    "MIR205HG"     "MRVI1_AS1"    "PART1"        "RBMS3_AS3"    "AFF3"        
#[15] "ANLN"         "E2F3"         "EZH2"         "FRMD5"        "NOVA1"        "PROX1"        "RUNX1"       
#[22] "SOX11"        "hsa_mir_216a" "hsa_mir_182"  "hsa_mir_96"   "hsa_mir_210"  "hsa_mir_205"  "hsa_mir_204" 
#[29] "hsa_mir_21"   "hsa_mir_217"  "hsa_mir_195"  "hsa_mir_141"  "hsa_mir_503"  "hsa_mir_122"  "hsa_mir_200a"
#[36] "hsa_mir_143"  "hsa_mir_183"  "hsa_mir_508"  "hsa_mir_363"  "hsa_mir_17"   "hsa_mir_372"  "hsa_mir_93"  
#[43] "hsa_mir_145"  "hsa_mir_187"  "hsa_mir_222"  "hsa_mir_301b" "hsa_mir_100"

pdf("LMO7_AS1.pdf")
fit<-survfit(Surv(futime,fustat)~LMO7_AS1, data=res.cat)
ggsurvplot(fit,
           data=res.cat,
           risk.table=TRUE,
           pval=T)
dev.off()
pdf("MRVI1_AS1.pdf")
fit<-survfit(Surv(futime,fustat)~MRVI1_AS1, data=res.cat)
ggsurvplot(fit,
           data=res.cat,
           risk.table=TRUE,
           pval=T)
dev.off()
pdf("hsa_mir_187.pdf")
fit<-survfit(Surv(futime,fustat)~hsa_mir_187, data=res.cat)
ggsurvplot(fit,
           data=res.cat,
           risk.table=TRUE,
           pval=T)
dev.off()
pdf("hsa_mir_100.pdf")
fit<-survfit(Surv(futime,fustat)~hsa_mir_100, data=res.cat)
ggsurvplot(fit,
           data=res.cat,
           risk.table=TRUE,
           pval=T)
dev.off()
pdf("AFF3.pdf")
fit<-survfit(Surv(futime,fustat)~AFF3, data=res.cat)
ggsurvplot(fit,
           data=res.cat,
           risk.table=TRUE,
           pval=T)
dev.off()
pdf("PROX1.pdf")
fit<-survfit(Surv(futime,fustat)~PROX1, data=res.cat)
ggsurvplot(fit,
           data=res.cat,
           risk.table=TRUE,
           pval=T)
dev.off()

