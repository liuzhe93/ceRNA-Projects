setwd("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.3_ROC")
remove(list=ls())

library("survivalROC")

rt<-read.table("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.2_multiCox/risk.txt",row.names=1,sep="\t",header=T)
pdf(file="ROC.pdf",width = 6, height = 6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime,status=rt$fustat,marker = rt$riskScore,
                predict.time = 1, method = "KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1),ylim=c(0,1),col="red",
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ", round(roc$AUC,3),")"),
     lwd=2,cex.main=1.3,cex.lab=1.2,cex.axis=1.2,font=1.2)
abline(0,1)
dev.off()

