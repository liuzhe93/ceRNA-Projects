setwd("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/01_DEG/lncRNA")
library(SummarizedExperiment)

testData<-load("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/mRNA_lncRNA/STAD_case.rda")
testData
data
assays(data)
names(assays(data))[5]
library(TCGAbiolinks)
dataPrep2 <- TCGAanalyze_Preprocessing(object = data,
                                       cor.cut = 0.6,
                                       datatype = "unstranded")

purityDATA <- TCGAtumor_purity(colnames(dataPrep2), 0, 0, 0, 0, 0.6)
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
#[1] 12987   448purityDATA <- TCGAtumor_purity(colnames(data), 0, 0, 0, 0, 0.6)
# filtered 为被过滤的数据， pure_barcodes是我们要的肿瘤数据
Purity.STAD<-purityDATA$pure_barcodes
normal.STAD<-purityDATA$filtered
#########################第六步：将肿瘤表达矩阵与正常组织表达矩阵合并，进行基因注释############################
#获取肿瘤纯度大于60%的340个肿瘤组织样本+50个正常组织样本,共计390个样本
puried_data <-dataFilt[,c(Purity.STAD,normal.STAD)]
###################################第七步：进行表达矩阵基因注释################################################
#基因注释,需要加载“SummarizedExperiment”包，“SummarizedExperiment container”每个由数字或其他模式的类似矩阵的对象表示。行通常表示感兴趣的基因组范围和列代表样品。
library("SummarizedExperiment")
rowData(data)   #传入数据dataPrep1必须为SummarizedExperiment对象
write.csv(puried_data,file = "puried.STAD.csv",quote = FALSE)

puried_data<-as.data.frame(puried_data)
puried_data$geneName<-rownames(puried_data)
anno_info<-read.table("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/hg38_ref/geneInfo.txt",sep="\t")
colnames(anno_info)<-c("geneID","geneName","geneType")
pcg<-subset(anno_info,geneType == c("sense_overlapping","lincRNA","3prime_overlapping_ncrna","processed_transcript",
                                    "antisense","sense_intronic"))
pcg_exp<-merge(puried_data,pcg,by="geneName")
xxx<-aggregate(pcg_exp,by=list(geneIndex=pcg_exp$geneName),mean)
pcg_exp_uniq<-xxx[,-which(colnames(xxx) %in% c("geneID","geneName","geneType"))]
rownames(pcg_exp_uniq)<-pcg_exp_uniq$geneIndex
pcg_exp_uniq<-pcg_exp_uniq[,-which(colnames(pcg_exp_uniq) %in% c("geneIndex"))]
write.csv(pcg_exp_uniq,"lncRNA_exp_uniq.csv",quote=F)

library( "edgeR" )
group <- c(rep("T",375),rep("N",32))
y <- DGEList(counts=pcg_exp_uniq, group=group)
head(y$counts)
head(y$samples)
# 数据过滤
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
# 计算标准化因子
y <- calcNormFactors(y)
# 查看标准化因子
y$samples
# 计算离散度
y <- estimateDisp(y)
# 显著性检验
et <- exactTest(y)
# 获取排名靠前的基因，这里设置n=100000是为了输出所以基因
et <- topTags(et, n=100000)
# 转换为数据框类型
et <- as.data.frame(et)
head(et)
# 将行名粘贴为数据框的第一列
et <- cbind(rownames(et),et)
# 指定列名
colnames(et) <- c("gene_id", "log2FoldChange", "log2CPM", "PValue", "FDR")
# 保存数据到本地
write.table(et, "case-vs-control-all.gene.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
head(et)
# 差异基因筛选
etSig <- et[which(et$FDR < 0.05 & abs(et$log2FoldChange) > 1),]
# 加入一列，up_down 体现上下调信息
etSig[which(etSig$log2FoldChange > 0), "up_down"] <- "Up"
etSig[which(etSig$log2FoldChange < 0), "up_down"] <- "Down"
# 保存文件
write.table(etSig, "case-vs-control-diff-pval-0.05-FC-2-edgeR.gene.xls", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
head(etSig)



library(ggplot2)
head(et)
et$color <- ifelse(et$FDR<0.05 & abs(et$log2FoldChange)>= 1,ifelse(et$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")
p <- ggplot(et, aes(log2FoldChange, -log10(FDR), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (FDR)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pdf("DElncRNA_VP.pdf")
p
dev.off()

