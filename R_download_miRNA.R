# 参考资料
# https://cloud.tencent.com/developer/article/1771685
# https://cloud.tencent.com/developer/article/1771576
setwd("/Users/liuzhe/Desktop/ceRNA/TCGA_data/miRNA")
#安装TCGAbiolinks
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("TCGAbiolinks")
#######################################一、数据下载阶段######################################
# 第一步：GDCquery（）筛选我们需要的数据，TCGAbiolinks包下载TCGA数据进行表达差异分析-肝癌案例
library(TCGAbiolinks)
library(DT)
library(dplyr)
library(SummarizedExperiment)
#选定要下载的cancer类型
TCGAbiolinks::getGDCprojects()$project_id
cancer_type="TCGA-STAD"
#选择下载你想要的数据类型
clinical<-GDCquery_clinic(project=cancer_type,type="clinical")
# 查看下载的数据
dim(clinical)
# 443  69
data_type <- "miRNA Expression Quantification"
data_category <- "Transcriptome Profiling"
query <- GDCquery(project = cancer_type,
                  data.category = data_category,
                  data.type = data_type,
                  workflow.type = "BCGSC miRNA Profiling")
samplesDown <- getResults(query,cols=c("cases"))  
#getResults(query, rows, cols)根据指定行名或列名从query中获取结果,此处用来获得样本的barcode
# 491
# 从samplesDown中筛选出TP（实体肿瘤）样本的barcodes
# TCGAquery_SampleTypes(barcode, typesample)
# TP代表PRIMARY SOLID TUMOR；NT-代表Solid Tissue Normal（其他组织样本可参考学习文档）
##此处共检索出446个TP样本barcodes
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
# 从samplesDown中筛选出NT(正常组织)样本的barcode
#此处共检索出45个NT样本barcodes
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
###设置barcodes参数，筛选符合要求的446个肿瘤样本数据和45正常组织数据
write.csv(dataSmTP,"Cancer_barcode_miRNA.csv",quote=F,row.names=F)
write.csv(dataSmNT,"Normal_barcode_miRNA.csv",quote=F,row.names=F)

queryDown <- GDCquery(project = "TCGA-STAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification", 
                      workflow.type = "BCGSC miRNA Profiling", 
                      barcode = c(dataSmTP, dataSmNT))

#barcode参数：根据传入barcodes进行数据过滤
######################################第二步：GDCdownload()下载GDCquery（）得到的结果#######################
# 下载数据，默认存放位置为当前工作目录下的GDCdata文件夹中。
setwd("/Users/liuzhe/Desktop/ceRNA/TCGA_data/miRNA")
GDCdownload(queryDown,method = "client", directory = "GDCdata",
            files.per.chunk = 10)
#method ；"API"或者"client"。"API"速度更快，但是容易下载中断。
#directory：下载文件的保存地址。Default: GDCdata。
#files.per.chunk = NULL:使用API下载大文件的时候，可以把文件分成几个小文件来下载，可以解决下载容易中断的问题。
GDCdownload(query = queryDown)
#读取下载的数据并将其准备到R对象中，在工作目录生成（save=TRUE）STAD_case.rda文件
# GDCprepare():Prepare GDC data,准备GDC数据，使其可用于R语言中进行分析
########第三步：GDCprepare()将前面GDCquery（）的结果准备成R语言可处理的SE（SummarizedExperiment）文件。######
dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename =
                          "STAD_case.rda")
########第四步：TCGAanalyze_Preprocessing()对数据进行预处理：使用spearman相关系数去除数据中的异常值##########
# 去除dataPrep1中的异常值，dataPrep1数据中含有肿瘤组织和正常组织的数据
testData<-load("STAD_case.rda")
testData
data
num_samples<-491
miRNA_name<-subset(data,select=miRNA_ID)
merged_count<-miRNA_name
for(i in 1:num_samples){
  temp<-data[,2+3*(i-1)]
  merged_count<-cbind(merged_count,temp)
}
colnames(merged_count)<-c("miRNA_ID",dataSmTP,dataSmNT)
counts_STAD<-merged_count

miRNA_name<-subset(data,select=miRNA_ID)
merged_FPKM<-miRNA_name
for(i in 1:num_samples){
  temp<-data[,3*i]
  merged_FPKM<-cbind(merged_FPKM,temp)
}
colnames(merged_FPKM)<-c("miRNA_ID",dataSmTP,dataSmNT)
FPKM_STAD<-merged_FPKM

write.csv(counts_STAD,"counts_STAD_miRNA.csv",row.names=F,quote=F)
write.csv(FPKM_STAD,"RPM_STAD_miRNA.csv",row.names=F,quote=F)



