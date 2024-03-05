############################################################################===
#####  GSE153873_counts+DEG
############################################################################===

# 1 GSE153873_targets ----
library(GEOquery)
library(tidyverse)

GSE153873_sm <- getGEO(filename = "AD/GSE153873/GSE153873_series_matrix.txt.gz",
                       getGPL = F) 

GSE153873_pd <- pData(GSE153873_sm)

sample_id <- GSE153873_pd[,2]
sample_name <- GSE153873_pd[,1]
subtype <- str_split_fixed(sample_name,'-',3)[,3]
group <- str_replace_all(subtype, c('Old'='Ctrl','Young'='Ctrl'))

GSE153873_targets <- data.frame(sample_id=sample_id,
                                sample_name=sample_name,
                                subtype=subtype,
                                group=group,
                                row.names = sample_id)


# 2 GSE153873_counts ----
# counts清洗
library(tidyverse)

GSE153873_counts <- read_tsv('AD/GSE153873/GSE153873_summary_count.star.txt.gz') %>% column_to_rownames('refGene')

GSE153873_targets <- arrange(GSE153873_targets,sample_name)
GSE153873_counts <- GSE153873_counts[,sort(colnames(GSE153873_counts))]

identical(GSE153873_targets[,2],colnames(GSE153873_counts))

colnames(GSE153873_counts) <- rownames(GSE153873_targets)

# 剔除Young亚组
GSE153873_targets_old <- filter(GSE153873_targets,subtype!='Young')
GSE153873_counts_old <- GSE153873_counts[,rownames(GSE153873_targets_old)]
# 样本过滤
keep <- rowSums(GSE153873_counts_old> 0) >= 15  ##保留样本中表达量大于0的基因

head(keep)

GSE153873_counts_old <- GSE153873_counts_old[keep,]


# 3 vst ----
#首先看一下有没有缺失值
sum(is.na(GSE153873_counts_old))
# 无
GSE153873_targets_old <- arrange(GSE153873_targets_old,group)
library(DESeq2)  
group <- GSE153873_targets_old$group %>% as.factor()
colData = data.frame(sample_id = colnames(GSE153873_counts_old), group = group)
dds <- DESeq2::DESeqDataSetFromMatrix(GSE153873_counts_old, #不是整数用round近似取整
                                      colData = colData,
                                      design = ~ group)
# 设置被比较组别
dds$group<- relevel(dds$group, ref = "Ctrl") # 指定哪一组作为对照组


vst <- DESeq2::vst(dds)
GSE153873_counts_vst <- assay(vst)

# 4 质量评估可视化
## ~boxplot & density plot ----
#count raw文件虽然过滤了一大半为0的基因，但是仍然有些counts文件会残余部分0的值，log2时会报错，如果还有0的话可以整体+1，这次的过滤之后没有0.
boxplot(log2(GSE153873_counts_old), las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(log2(GSE153873_counts_old), legend = F)

boxplot(GSE153873_counts_vst, las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(GSE153873_counts_vst , legend = F)

## ~PCA plot ----
source("custom_functions.R")
# 未标准化
PCA_new(log2(GSE153873_counts_old), group = GSE153873_targets_old$group)
# VST标准化
PCA_new(GSE153873_counts_vst, group = GSE153873_targets_old$group,show_name = T)

# 加载包
library("factoextra")
library("FactoMineR")
# 准备数据
pca1 <- PCA(t(log2(GSE153873_counts_old+1)), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca1,
             label = "none", 
             habillage = as.factor(GSE153873_targets_old$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

pca2 <- PCA(t(GSE153873_counts_vst), graph = FALSE)
fviz_pca_ind(pca2,
             label = "none", 
             habillage = as.factor(GSE153873_targets_old$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)
## ~聚类树状图 ----
# 查看默认配色
palette.pals() #查看配色方案
palette() #查看当前配色方案
palette("Set 1") #选择一个配色方案
palette("default") # 重置为默认配色方案”R3“
# 查看brewer配色系统
display.brewer.all()

#未标准化
expr1 <- as.matrix(GSE153873_counts_old)
dist1 <- dist(t(expr1))
hc = hclust(dist1)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc, k = 3, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())

#vst标准化
expr2 <- as.matrix(GSE153873_counts_vst)
dist2 <- dist(t(expr2))
hc = hclust(dist2)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc, k = 4, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()


## ~基于arrayQulitMetrics ----

###  构建ExpressionSet ### 
library(arrayQualityMetrics)

#因为arrayQualityMetrics不支持lumi、ElistRaw格式的数据，此时需要先构建ExpressionSet
#因为arrayQualityMetrics出现NA会提示下标越界，GSE153873_sm里有NA值，此时需要利用处理过NA值的表达矩阵重新构建ExpressionSet

GSE153873_counts_old_QC <- ExpressionSet(assayData = as.matrix(GSE153873_counts_old), 
                                     phenoData = AnnotatedDataFrame(data = GSE153873_targets_old))
arrayQualityMetrics(GSE153873_counts_old_QC, 
                    outdir = "AD/GSE153873/GSE153873_counts_old_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)

GSE153873_counts_vst_QC <- ExpressionSet(assayData = as.matrix(GSE153873_counts_vst), 
                                         phenoData = AnnotatedDataFrame(data = GSE153873_targets_old))

arrayQualityMetrics(GSE153873_counts_vst_QC, 
                    outdir = "AD/GSE153873/GSE153873_counts_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)

dev.off()
# 4 DESeq2 差异表达分析 ----
## ~构建dds对象 ----
library(DESeq2)  
## ~DESeq一步完成差异分析 ----
time1 <- Sys.time()
dds2 <- DESeq(dds)
runtime1 <- Sys.time() - time1
runtime1   # Time difference of 4.238453 mins

time2 <- Sys.time()
dds2 <- DESeq(dds_filt, parallel = T)
runtime2 <- Sys.time() - time2
runtime2   # Time difference of 3.506404 mins

# save(dds2, file = "Covid/GSE153873/data/dds2.Rda")
load(file = "Covid/GSE153873/data/dds2.Rda")

## ~提取差异分析结果 ----
resultsNames(dds2)
res <- results(dds2)
head(res)

# res <- results(dds2, contrast=c("group","AD","Ctrl"))

# GSE153873_DEG_P <- as.data.frame(res) %>% 
#   na.omit() %>% 
#   arrange(padj) %>% 
#   dplyr::filter(abs(log2FoldChange) > 1, pvalue < 0.05)  #P值，不是adjP



GSE153873_DEG <- as.data.frame(res) %>%
  rownames_to_column("gene_name") %>% 
  na.omit() %>% dplyr::select(-2) %>% 
  dplyr::rename(., logFC = log2FoldChange )

GSE153873_DEG_sig <- as.data.frame(res) %>%
  rownames_to_column("gene_name") %>% 
  na.omit() %>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05) %>% 
  dplyr::select(-2) %>% 
  dplyr::rename(., logFC = log2FoldChange )

write.table(GSE153873_DEG,
            file = "AD/GSE153873/results/GSE153873_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(GSE153873_DEG_sig,
            file = "AD/GSE153873/results/GSE153873_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)



# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. 
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. 


## ~保存分析结果数据 ----
save(GSE153873_DESeq2_DEG,file = "AD/GSE153873/results/GSE153873_DESeq2_DEG.Rda")



length(intersect(res_deseq2$gene_name, GSE153873_DEG_sig[,1]))
library(VennDiagram)
venn.diagram(x = list("Covid-19" = res_deseq2$gene_name, 
                      "AD" = GSE153873_DEG_sig[,1]),
             filename = "AD/GSE153873/results/Venn2.jpeg",
             fill = c("green", "red"),
             scaled = F,
             cex = 1.5,
             cat.cex = 1.5,
             main = "Venn Diagram",
             main.cex = 2)
