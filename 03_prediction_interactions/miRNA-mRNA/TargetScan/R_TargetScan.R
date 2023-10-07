setwd("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/miRNA-mRNA/TargetScan")
remove(list=ls())

targetScan<-read.table("Predicted_Targets_Context_Scores.default_predictions.txt",sep="\t",header=T)
dim(targetScan)
#[1] 1397978     13

targetScan_human<-targetScan[grep("^hsa",targetScan[,5]),]
dim(targetScan_human)
#[1] 228049     13
targetScan_human$miRNA_name<-targetScan_human$miRNA

miRNA<-read.table("/Users/liuzhe/Desktop/ceRNA/03_prediction_interactions/lncRNA-miRNA/miRNAs.txt")
dim(miRNA)
#[1] 44  1
for(i in 1:nrow(miRNA)){
  temp<-miRNA$V1[i]
  substring(temp,1,7)<-'hsa-miR'
  miRNA$miRNA_name[i]<-temp
}

merged<-merge(miRNA,targetScan_human,by="miRNA_name")
inter_miRNAmRNA<-subset(merged,select=c("miRNA_name","Gene.Symbol"))
colnames(inter_miRNAmRNA)<-c("miRNA","mRNA")
write.csv(inter_miRNAmRNA,"inter_miRNAmRNA_TargetScan.csv",quote=F,row.names = F)


