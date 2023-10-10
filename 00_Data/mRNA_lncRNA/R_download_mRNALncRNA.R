# 参考资料
# https://cloud.tencent.com/developer/article/1771685
# https://cloud.tencent.com/developer/article/1771576
setwd("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/mRNA_lncRNA")
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
View(clinical)
setwd("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/mRNA_lncRNA/clinical")
save(clinical,file="STAD_clinical.Rdata")
write.csv(clinical, file="TCGAbiolinks-STAD-clinical.csv")
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
query <- GDCquery(project = cancer_type,
                  data.category = data_category,
                  data.type = data_type,
                  workflow.type = "STAR - Counts")
samplesDown <- getResults(query,cols=c("cases"))  
#getResults(query, rows, cols)根据指定行名或列名从query中获取结果,此处用来获得样本的barcode
# 407
# 从samplesDown中筛选出TP（实体肿瘤）样本的barcodes
# TCGAquery_SampleTypes(barcode, typesample)
# TP代表PRIMARY SOLID TUMOR；NT-代表Solid Tissue Normal（其他组织样本可参考学习文档）
##此处共检索出375个TP样本barcodes
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
# 从samplesDown中筛选出NT(正常组织)样本的barcode
#此处共检索出32个NT样本barcodes
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
###设置barcodes参数，筛选符合要求的375个肿瘤样本数据和32正常组织数据
write.csv(dataSmTP,"Cancer_barcode.csv",quote=F,row.names=F)
write.csv(dataSmNT,"Normal_barcode.csv",quote=F,row.names=F)

queryDown <- GDCquery(project = "TCGA-STAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "STAR - Counts", 
                      barcode = c(dataSmTP, dataSmNT))
#barcode参数：根据传入barcodes进行数据过滤
######################################第二步：GDCdownload()下载GDCquery（）得到的结果#######################
# 下载数据，默认存放位置为当前工作目录下的GDCdata文件夹中。
setwd("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/mRNA_lncRNA")
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
assays(data)
names(assays(data))[5]
dataPrep2 <- TCGAanalyze_Preprocessing(object = data,
                                       cor.cut = 0.6,
                                       datatype = "fpkm_unstrand")
#Number of outliers: 0
setwd("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/mRNA_lncRNA/FPKM")
#将预处理后的数据dataPrep2，写入新文件“STAD_dataPrep.csv”
write.csv(dataPrep2,file = "STAD_dataPrep.csv",quote = FALSE)
#########################第五步：TCGAtumor_purity（）筛选肿瘤纯度大于60%的肿瘤barcodes########################
# TCGAtumor_purity(barcodes, estimate, absolute, lump, ihc, cpe)，使用来自5种方法的5个估计值作为阈值对TCGA样本进行过滤，这5个值是estimate, absolute, lump, ihc, cpe，这里设置cpe=0.6（cpe是派生的共识度量，是将所有方法的标准含量归一化后的均值纯度水平，以使它们具有相等的均值和标准差）
#筛选肿瘤纯度大于等于60%的样本数据
purityDATA <- TCGAtumor_purity(colnames(data), 0, 0, 0, 0, 0.6)
# filtered 为被过滤的数据， pure_barcodes是我们要的肿瘤数据
Purity.STAD<-purityDATA$pure_barcodes
normal.STAD<-purityDATA$filtered
#########################第六步：将肿瘤表达矩阵与正常组织表达矩阵合并，进行基因注释############################
#获取肿瘤纯度大于60%的340个肿瘤组织样本+50个正常组织样本,共计390个样本
puried_data <-dataPrep2[,c(Purity.STAD,normal.STAD)]
###################################第七步：进行表达矩阵基因注释################################################
#基因注释,需要加载“SummarizedExperiment”包，“SummarizedExperiment container”每个由数字或其他模式的类似矩阵的对象表示。行通常表示感兴趣的基因组范围和列代表样品。
library("SummarizedExperiment")
rowData(data)   #传入数据dataPrep1必须为SummarizedExperiment对象
rownames(puried_data)<-rowData(data)$gene_name
write.csv(puried_data,file = "puried.STAD.csv",quote = FALSE)
#########################第八步：进行表达矩阵标准化和过滤，得到用于差异分析的表达矩阵##########################
library("EDASeq")
dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(tabDF = puried_data,
                                                    geneInfo = geneInfo,
                                                    method = "gcContent")
#将标准化后的数据再过滤，去除掉表达量较低（count较低）的基因，得到最终的数据
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
dim(dataNorm)
#[1] 17317   448
dim(dataFilt)
#[1] 12987   448
str(dataFilt)
#num [1:12925, 1:390] 0 30 803 0 0 2 0 263 6 4 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:12925] "A1BG" "A1CF" "A2M" "A4GALT" ...
#..$ : chr [1:390] "TCGA-ZP-A9CY-01A-11R-A38B-07" "TCGA-DD-AAVZ-01A-11R-A41C-07" "TCGA-DD-AADN-01A-11R-A41C-07" "TCGA-DD-A1EB-01A-11R-A131-07" ...
write.csv(dataFilt,file = "TCGA_STAD_final.csv",quote = FALSE)  
#保留的是407个样本（前375肿瘤，后32正常组织）

#读入表达矩阵文件
dataFilt_STAD_final <- read.csv("TCGA_STAD_final.csv", header = T,check.names = FALSE)
# 定义行名
rownames(dataFilt_STAD_final) <- dataFilt_STAD_final[,1]
dataFilt_STAD_final <- dataFilt_STAD_final[,-1]
# 先看一下矩阵长啥样，心里有个数：每一行是一个基因，每一列是一个样本
View(dataFilt_STAD_final)


