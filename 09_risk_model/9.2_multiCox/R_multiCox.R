setwd("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.2_multiCox")
remove(list=ls())

library("survival")
library("survminer")

rt<-read.csv("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.1_lasso_undo/uniSigExpID.csv")
#rt<-rt[,1:5]
for(i in 1:nrow(rt)){
  temp=rt$Id[i]
  temp<-paste0(temp,"_",i)
  rt$Id[i]<-temp
}
row.names(rt)<-rt$Id
rt<-rt[,-1]

multiCox<-coxph(Surv(futime,fustat) ~ ., data = rt)
multiCox<-step(multiCox,direction="both")
multiCoxSum<-summary(multiCox)

outTab<-data.frame()
outTab<-cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub('"',"",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names = F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=rt)
rt$riskScore<-riskScore

rt_saved<-rt
rt_saved$riskCLass<-ifelse(rt$riskScore>1.144427,"high","low")

write.table(rt_saved,"risk_class.txt",sep="\t",row.names = T, quote = F)
write.table(rt,"risk.txt",sep="\t",row.names = T, quote = F)
#计算最佳截点
rt=rt[,c(1,2,16)]
res.cut<-surv_cutpoint(rt,
                       time = "futime",
                       event = "fustat",
                       variables = "riskScore")
summary(res.cut)
pdf("riskScore_distribution.pdf",height = 10, width =6)
plot(res.cut,"riskScore",palette="npg")
dev.off()

res.cat<-surv_categorize(res.cut)
head(res.cat)

fit<-survfit(Surv(futime,fustat)~riskScore,data=res.cat)
pdf("survival_curve.pdf",height = 6, width = 6)
ggsurvplot(fit,
           data=res.cat,
           risk.table = TRUE,
           pval=T)
dev.off()

