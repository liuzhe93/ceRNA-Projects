setwd("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/12_regression_DElncRNA_DEmRNA")
remove(list=ls())

library("tidyverse")
library("readxl")
ceRNAs<-read_excel("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/08_ceRNA_survival/Supplementary Table 5. Survival analysis of ceRNAs network-related genes.xlsx")
lncRNAs<-ceRNAs[1:26,1]
lncRNAs<-unique(lncRNAs)
lncRNAs$X<-lncRNAs$Gene_Name
mRNAs<-ceRNAs[77:94,1]
mRNAs<-unique(mRNAs)
mRNAs$X<-mRNAs$Gene_Name
miRNAs<-ceRNAs[27:76,1]
miRNAs<-unique(miRNAs)
miRNAs$X<-miRNAs$Gene_Name


lncRNA_exp<-read.csv("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/01_DEG/lncRNA/lncRNA_exp_uniq.csv")
mRNA_exp<-read.csv("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/01_DEG/mRNA/mRNA_exp_uniq.csv")
miRNA_exp<-read.csv("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/01_DEG/miRNA/RPM_STAD_miRNA.csv")

lncRNA_selected<-merge(lncRNAs, lncRNA_exp, by="X")
row.names(lncRNA_selected)<-lncRNA_selected$X
lncRNA_selected<-lncRNA_selected[,-(1:2)]
lncRNA_selected<-log2(lncRNA_selected+1)
lncRNA_selected<-t(lncRNA_selected)
mRNA_selected<-merge(mRNAs, mRNA_exp, by="X")
row.names(mRNA_selected)<-mRNA_selected$X
mRNA_selected<-mRNA_selected[,-(1:2)]
mRNA_selected<-log2(mRNA_selected+1)
mRNA_selected<-t(mRNA_selected)

miRNA_exp$X<-miRNA_exp$miRNA_ID
for(i in 1:nrow(miRNA_exp)){
  ss = miRNA_exp$miRNA_ID[i]
  if(grepl("-mir-", ss)){
    miRNA_exp$X[i]= gsub("mir","miR",ss) 
  }
}
miRNA_selected<-merge(miRNAs, miRNA_exp, by ="X")
row.names(miRNA_selected)<-miRNA_selected$X
miRNA_selected<-miRNA_selected[,-(1:3)]
miRNA_selected<-log2(miRNA_selected+1)
miRNA_selected<-t(miRNA_selected)

intersect_samples<-intersect(row.names(lncRNA_selected), row.names(miRNA_selected))
length(intersect_samples)

idx<-intersect_samples

lncRNA_intersect<-lncRNA_selected[idx,]
miRNA_intersect<-miRNA_selected[idx,]
mRNA_intersect<-mRNA_selected[idx,]

pval=data.frame()
corval=data.frame()
for(i in 1:ncol(lncRNA_intersect)){
  for(j in 1:ncol(mRNA_intersect)){
    temp<-cor.test(lncRNA_intersect[,i],mRNA_intersect[,j],method = "pearson")
    pval[i,j]<-temp$p.value
    corval[i,j]<-temp$estimate
#    pval<-temp$p.value
#    cor<-temp$estimate
  }
}
write.csv(pval,"Pvalue_lncRNA_mRNA.csv")
write.csv(corval,"CorValue_lncRNA_mRNA.csv")
pos=data.frame()
for(i in 1:nrow(pval)){
  for(j in 1:ncol(pval)){
    if(pval[i,j]<=0.05 & corval[i,j]>0.3){
      pos[i,j] = 1
    }else{
      pos[i,j] = 0
    }
  }
}
sum(pos)
#[1] 13
write.csv(pos,"Pos_lncRNA_mRNA.csv")

pval=data.frame()
corval=data.frame()
for(i in 1:ncol(lncRNA_intersect)){
  for(j in 1:ncol(miRNA_intersect)){
    temp<-cor.test(lncRNA_intersect[,i],miRNA_intersect[,j],method = "pearson")
    pval[i,j]<-temp$p.value
    corval[i,j]<-temp$estimate
    #    pval<-temp$p.value
    #    cor<-temp$estimate
  }
}
colnames(pval)<-colnames(miRNA_intersect)
rownames(pval)<-colnames(lncRNA_intersect)
colnames(corval)<-colnames(miRNA_intersect)
rownames(corval)<-colnames(lncRNA_intersect)
write.csv(pval,"Pvalue_lncRNA_miRNA.csv")
write.csv(corval,"CorValue_lncRNA_miRNA.csv")
neg=data.frame()
for(i in 1:nrow(pval)){
  for(j in 1:ncol(pval)){
    if(pval[i,j]<=0.05 & corval[i,j]< -0.3){
      neg[i,j] = 1
    }else{
      neg[i,j] = 0
    }
  }
}
sum(neg)
#[1] 34
colnames(neg)<-colnames(miRNA_intersect)
rownames(neg)<-colnames(lncRNA_intersect)
write.csv(neg,"Neg_lncRNA_miRNA.csv")


pval=data.frame()
corval=data.frame()
for(i in 1:ncol(mRNA_intersect)){
  for(j in 1:ncol(miRNA_intersect)){
    temp<-cor.test(mRNA_intersect[,i],miRNA_intersect[,j],method = "pearson")
    pval[i,j]<-temp$p.value
    corval[i,j]<-temp$estimate
    #    pval<-temp$p.value
    #    cor<-temp$estimate
  }
}
colnames(pval)<-colnames(miRNA_intersect)
rownames(pval)<-colnames(mRNA_intersect)
colnames(corval)<-colnames(miRNA_intersect)
rownames(corval)<-colnames(mRNA_intersect)
write.csv(pval,"Pvalue_mRNA_miRNA.csv")
write.csv(corval,"CorValue_mRNA_miRNA.csv")
neg=data.frame()
for(i in 1:nrow(pval)){
  for(j in 1:ncol(pval)){
    if(pval[i,j]<=0.05 & corval[i,j]< -0.3){
      neg[i,j] = 1
    }else{
      neg[i,j] = 0
    }
  }
}
sum(neg)
#[1] 36
colnames(neg)<-colnames(miRNA_intersect)
rownames(neg)<-colnames(mRNA_intersect)
write.csv(neg,"Neg_mRNA_miRNA.csv")


all_data<-cbind(lncRNA_selected,mRNA_selected)
all_data<-as.data.frame(all_data)

library("ggpubr")
pdf("ERVMER61-1_E2F3.pdf")
ggscatter(all_data,
          x="ERVMER61-1",
          y="E2F3",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="ERVMER61-1",
          ylab="E2F3")
dev.off()

pdf("KCNA3_AFF3.pdf")
ggscatter(all_data,
          x="KCNA3",
          y="AFF3",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="KCNA3",
          ylab="AFF3")
dev.off()

pdf("KCNA3_NOVA1.pdf")
ggscatter(all_data,
          x="KCNA3",
          y="NOVA1",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="KCNA3",
          ylab="NOVA1")
dev.off()

pdf("LINC00114_PROX1.pdf")
ggscatter(all_data,
          x="LINC00114",
          y="PROX1",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LINC00114",
          ylab="PROX1")
dev.off()


pdf("LINC00501_ANLN.pdf")
ggscatter(all_data,
          x="LINC00501",
          y="ANLN",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LINC00501",
          ylab="ANLN")
dev.off()

pdf("LINC00501_E2F3.pdf")
ggscatter(all_data,
          x="LINC00501",
          y="E2F3",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LINC00501",
          ylab="E2F3")
dev.off()

pdf("LINC00501_EZH2.pdf")
ggscatter(all_data,
          x="LINC00501",
          y="EZH2",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LINC00501",
          ylab="EZH2")
dev.off()

pdf("LINC00501_FRMD5.pdf")
ggscatter(all_data,
          x="LINC00501",
          y="FRMD5",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LINC00501",
          ylab="FRMD5")
dev.off()

pdf("LMO7-AS1_ANLN.pdf")
ggscatter(all_data,
          x="LMO7-AS1",
          y="ANLN",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LMO7-AS1",
          ylab="ANLN")
dev.off()

pdf("LMO7-AS1_E2F3.pdf")
ggscatter(all_data,
          x="LMO7-AS1",
          y="E2F3",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LMO7-AS1",
          ylab="E2F3")
dev.off()

pdf("LMO7-AS1_EZH2.pdf")
ggscatter(all_data,
          x="LMO7-AS1",
          y="EZH2",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LMO7-AS1",
          ylab="EZH2")
dev.off()

pdf("LMO7-AS1_FRMD5.pdf")
ggscatter(all_data,
          x="LMO7-AS1",
          y="FRMD5",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LMO7-AS1",
          ylab="FRMD5")
dev.off()

pdf("LMO7-AS1_PROX1.pdf")
ggscatter(all_data,
          x="LMO7-AS1",
          y="PROX1",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LMO7-AS1",
          ylab="PROX1")
dev.off()

pdf("LMO7-AS1_RUNX1.pdf")
ggscatter(all_data,
          x="LMO7-AS1",
          y="RUNX1",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LMO7-AS1",
          ylab="RUNX1")
dev.off()

pdf("LMO7-AS1_SOX11.pdf")
ggscatter(all_data,
          x="LMO7-AS1",
          y="SOX11",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="LMO7-AS1",
          ylab="SOX11")
dev.off()

pdf("MRVI1-AS1_AFF3.pdf")
ggscatter(all_data,
          x="MRVI1-AS1",
          y="AFF3",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="MRVI1-AS1",
          ylab="AFF3")
dev.off()

pdf("MRVI1-AS1_NOVA1.pdf")
ggscatter(all_data,
          x="MRVI1-AS1",
          y="NOVA1",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="MRVI1-AS1",
          ylab="NOVA1")
dev.off()

pdf("PART1_AFF3.pdf")
ggscatter(all_data,
          x="PART1",
          y="AFF3",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="PART1",
          ylab="AFF3")
dev.off()

pdf("PART1_NOVA1.pdf")
ggscatter(all_data,
          x="PART1",
          y="NOVA1",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="PART1",
          ylab="NOVA1")
dev.off()

pdf("RBMS3-AS3_AFF3.pdf")
ggscatter(all_data,
          x="RBMS3-AS3",
          y="AFF3",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="RBMS3-AS3",
          ylab="AFF3")
dev.off()

pdf("RBMS3-AS3_NOVA1.pdf")
ggscatter(all_data,
          x="RBMS3-AS3",
          y="NOVA1",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab="RBMS3-AS3",
          ylab="NOVA1")
dev.off()

library("ggalluvial")
mydata<-data.frame(
  LncRNA=c(rep("KCNA3",2)),
  miRNA=c(rep("hsa-miR-217",2)),
  mRNA=c("AFF3","NOVA1"),
  corval=c(0.552467925932138,0.382821522564157)
)
mydata

pdf("ggalluvial.pdf")
ggplot(mydata,
       aes(y = corval,
           axis1 = LncRNA, axis2 = miRNA, axis3 = mRNA)) +
  geom_alluvium(aes(fill=miRNA), width = 0, reverse = FALSE, discern = TRUE) +
  geom_stratum(width = 1/3, reverse = FALSE, discern = TRUE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = FALSE, size = 4, angle = 0, discern = TRUE) +
  scale_x_continuous(breaks = 1:3, labels = c("LncRNA", "miRNA", "mRNA"))+
  theme(legend.position = "none") +
  ggtitle("ceRNA")
dev.off()










