setwd("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.1_lasso_undo")
remove(list=ls())

library("glmnet")
library("survival")

lncRNA_kept<-read.csv("/Users/liuzhe/Desktop/ceRNA/08_ceRNA_survival/lncRNA_kept.csv")
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
dim(lncRNA_export_t)
#[1] 407  14
lncRNA_export_t[1:4,1:4]
for(i in 1:nrow(lncRNA_export_t)){
  temp<-lncRNA_export_t$Sample_Name[i]
  temp<-substr(temp,1,12)
  temp<-gsub("\\.","-",temp)
  lncRNA_export_t$Sample_Name[i]<-temp
}

clinical<-read.table("/Users/liuzhe/Desktop/ceRNA/TCGA_data/clinical/clinical_rmNA.txt",header=T)
clinical$Sample_Name<-clinical$Id

matrix<-merge(lncRNA_export_t,clinical,by="Sample_Name")
dim(matrix)
#[1] 317  27

names_change<-colnames(matrix)
names_change<-gsub("-","_",names_change)
colnames(matrix)<-names_change
uniSigExp<-subset(matrix,select=c("Id","futime","fustat","DIAPH3_AS1","ERVMER61_1","KCNA3","LINC00114",
                                  "LINC00462","LINC00486","LINC00501","LMO7_AS1","LSAMP_AS1",
                                  "MIR205HG","MRVI1_AS1","PART1","RBMS3_AS3"))
write.csv(uniSigExp,"uniSigExpID.csv",quote=F,row.names = F)
uniSigExp<-subset(matrix,select=c("futime","fustat","DIAPH3_AS1","ERVMER61_1","KCNA3","LINC00114",
                                  "LINC00462","LINC00486","LINC00501","LMO7_AS1","LSAMP_AS1",
                                  "MIR205HG","MRVI1_AS1","PART1","RBMS3_AS3"))
write.csv(uniSigExp,"uniSigExp.csv",quote=F,row.names = F)

rt<-read.csv("uniSigExpID.csv")
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/365

x<-as.matrix(rt[,c(3:ncol(rt))])
y<-data.matrix(Surv(rt$futime,rt$fustat))
fit<-glmnet(x,y,family="cox",maxit=1000)

pdf("lambda.pdf")
plot(fit,xvar="lambda",label=TRUE)
dev.off()

cvfit<-cv.glmnet(x,y,family="cox",maxit=1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef<-coef(fit,s=cvfit$lambda.min)
index<-which(coef!=0)
actCoef<-coef[index]
lassoGene<-row.names(coef)[index]
lassoGene<-c("futime","fustat",lassoGene)
lassoSigExp<-rt[,lassoGene]
lassoSigExp<-cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names = F,quote = F)





