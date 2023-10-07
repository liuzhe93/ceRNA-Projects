setwd("/Users/liuzhe/Desktop/ceRNA/10_uniIndep")
remove(list=ls())

riskClass<-read.table("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.2_multiCox/risk_class.txt")
riskClass$Id<-row.names(riskClass)
for(i in 1:nrow(riskClass)){
  temp<-riskClass$Id[i]
  temp<-strsplit(temp,split="_")[[1]][1]
  riskClass$Id[i]<-temp
}
clinical<-read.table("/Users/liuzhe/Desktop/ceRNA/TCGA_data/clinical/clinical_rmNA.txt",header=T)
matrix<-merge(clinical,riskClass,by="Id")
matrix_selected<-subset(matrix,select=c("Id","futime.x","fustat.x","age", "gender", "stage", "tissue_of_origin", 
                                        "primary_diagnosis","ajcc_T", "ajcc_M", "ajcc_N", 
                                        "race","riskScore"))
colnames(matrix_selected)<-c("Id","futime","fustat","age", "gender", "stage", "tissue_of_origin", 
                             "primary_diagnosis","ajcc_T", "ajcc_M", "ajcc_N", 
                             "race","riskScore")
 
library("survival")
library("forestplot")

clrs<-fpColors(box="green",line="darkblue",summary="royalblue")
rt<-matrix_selected

outTab=data.frame()
for(i in colnames(rt[,4:ncol(rt)])){
  cox<-coxph(Surv(futime,fustat)~rt[,i],data=rt)
  coxSummary=summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
               )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names = F,quote=F)

#绘制森林图
rt=read.table("uniCox.xls",header=T,sep="\t",row.names = 1,check.names = F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001,"<0.001",sprintf("%.3f",pVal))
tabletext<-
  list(c(NA,rownames(HR)),
       append("pvalue",pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")))
pdf(file="forest.pdf",onefile = FALSE,
    width = 8,
    height = 4)
forestplot(tabletext,
           rbind(rep(NA,3),HR),
           col=clrs,
           graphwidth=unit(50,"mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Harard ratio")
dev.off()





















