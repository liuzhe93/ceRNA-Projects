setwd("/Users/liuzhe/Desktop/ceRNA/01_DEG/miRNA")
library(SummarizedExperiment)

counts_STAD<-read.csv("/Users/liuzhe/Desktop/ceRNA/TCGA_data/miRNA/counts_STAD_miRNA.csv")
counts_value<-counts_STAD[,2:492]
counts_value<-as.matrix(counts_value)
rownames(counts_value)<-counts_STAD$miRNA_ID

library( "edgeR" )
group <- c(rep("T",446),rep("N",45))

y <- DGEList(counts=counts_value, group=group)
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
pdf("DEmiRNA_VP.pdf")
p
dev.off()

