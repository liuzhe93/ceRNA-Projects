setwd("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.4_riskPlot")
remove(list=ls())

library("pheatmap")

rt<-read.table("/Users/liuzhe/Desktop/ceRNA/09_risk_model/9.2_multiCox/risk_class.txt",row.names = 1,sep="\t",header=T)
rt=rt[order(rt$riskScore),]

#绘制风险曲线
riskClass=rt[,"riskCLass"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
pdf(file="riskScore.pdf",width=12,height=5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk score)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=1.144427,v=lowLength,lty=2)
dev.off()

#绘制生存状态图
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat.pdf",width=12,height=5)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk score)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

#绘制风险热图
rt1=rt[,c(8,11)]
#rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",height = 2, width = 12)
pheatmap(rt1,
         annotaion=annotation,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row = 11,
         fontsize_col = 3,
         scale="row",
         color=colorRampPalette(c("blue","white","red"))(50))
dev.off()

