############################################################################===
#####  GSE183533_targets+GSE183533_counts 
############################################################################===

# 1 GSE183533_targets ----
library(GEOquery)
library(tidyverse)
library(limma)

GSE183533_sm <- getGEO(filename = "Covid/GSE183533/GSE183533_series_matrix.txt.gz",
                       getGPL = F) 

GSE183533_pd <- pData(GSE183533_sm)

GSE183533_targets <- pData(GSE183533_sm) %>%
  dplyr::select(2, 1) %>%#选出两列，title 这一列其实包含了3列内容
  separate(title, into = paste0("x", 1:3), sep = "_") %>%#把title 这一列拆成3列，根据下划线来分割
  mutate(sample_id = geo_accession,
         group = str_replace_all(x1,c('Normal'='con','covid-19'='COV')),
         sample_name = x3,
         counts_name = str_replace_all(x3,c('4N'='N_4','17N'='N_17','7N'='N_7','19N'='N_19','21N'='N_21'))) %>%
  dplyr::select(sample_id:counts_name)  #'7N'='N_7','17N'='N_17'顺序反了17N会变成1N_7,需要注意。



# 写出targets文件
write.table(x = GSE183533_targets,
            file = "Covid/GSE183533/GSE183533_targets.txt",
            quote = F,
            sep = "\t",
            row.names = F)

# 读取targets文件
GSE183533_targets <- readTargets(file = "Covid/GSE183533/GSE183533_targets.txt",
                                 row.names = "sample_id")

# END
save(GSE183533_targets, file = "Covid/GSE183533/GSE183533_targets.Rdata")


# 2 GSE183533_counts ----
# counts清洗
library(tidyverse)
PCA_new <- function(expr, ntop = 500, group, show_name = F){
  library(ggplot2)
  library(ggrepel)
  rv <- genefilter::rowVars(expr)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(expr[select, ]))#最核心的代码
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1], 
                  PC2 = pca$x[, 2], 
                  group = group, 
                  name = colnames(expr))
  attr(d, "percentVar") <- percentVar[1:2]
  if (show_name) {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      geom_text_repel(aes(label = name),
                      size = 3,
                      segment.color = "black",
                      show.legend = FALSE )
  } else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2",color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
  }
}


GSE183533_counts <- read.csv('Covid/GSE183533/GSE183533_counts_cov_41CovarMatrix.csv.gz',header=T, row.names = 1)

identical(GSE183533_targets[,4],colnames(GSE183533_counts))

#将counts文件列名替换为GSM开头的
colnames(GSE183533_counts) <- rownames(GSE183533_targets)

# 样本过滤
keep <- rowSums(GSE183533_counts> 0)>= 21  ##保留至少在50%样本中表达量大于0的基因

head(keep)

GSE183533_counts <- GSE183533_counts[keep,]

# save
save(GSE183533_counts,file = "Covid/GSE183533/GSE183533_counts.Rdata")

# 3 vst ----
#首先看一下有没有缺失值
sum(is.na(GSE183533_counts))
# 无

library(DESeq2)  
group <- GSE183533_targets$group %>% as.factor()
colData = data.frame(sample_id = colnames(GSE183533_counts), group = group)
dds <- DESeq2::DESeqDataSetFromMatrix(GSE183533_counts, #不是整数用round近似取整
                                      colData = colData,
                                      design = ~ group)
vst <- DESeq2::vst(dds)
GSE183533_counts_vst <- assay(vst)

save(GSE183533_counts_vst, file = "Covid/GSE183533/GSE183533_counts_vst.Rdata")
# 4 质量评估可视化
## ~boxplot & density plot ----
#count raw文件虽然过滤了一大半为0的基因，但是仍然有些counts文件会残余部分0的值，log2时会报错，如果还有0的话可以整体+1，这次的过滤之后没有0.
boxplot(log2(GSE183533_counts), las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(log2(GSE183533_counts), legend = F)

boxplot(GSE183533_counts_vst, las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(GSE183533_counts_vst , legend = F)

## ggplot2 data   

# # 这是芯片数据里的ggplot箱线图数据准备代码，看看区别
# box_data <-  log2(GSE183533_counts + 1) %>% #对表达矩阵log2转换，整体+1避免log2(0)报错
#   as.data.frame() %>% #tideverse包不接受矩阵，只接受数据框
#   rownames_to_column("gene_id") %>% #行名变列，这一列叫probe_id
#   slice_sample(n = 50000) %>%  #因为运行速度很慢，所以只挑了50000个数据来运行
#   pivot_longer(cols = starts_with("GSM"),#对数据进行长宽转换，从宽数据到长数据，转换以gsm开头的列
#                names_to = "sample_id",#从gsm变成sample_id
#                values_to = "value") %>% #新建一列叫value
#   left_join(GSE183533_targets %>% dplyr::select(sample_id, group), #left_join相当于取交集，以左边的为主
#             by = "sample_id") #以sample_id为锚点取交集，多了一列就是group

# 测序数据代码
# 未标准化的counts数据准备
data_counts <- log2(GSE183533_counts) %>%  #对表达矩阵log2转换，如有0，就整体+1避免log2(0)报错
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE183533_targets %>% select(sample_id, group), by = "sample_id")

# vst标准化后的counts数据准备
data_vst <- GSE183533_counts_vst %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE183533_targets %>% select(sample_id, group), by = "sample_id")

## ggplot2 box plot 

# 未标准化箱线图
ggplot(data_counts, aes(x = sample_id, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() + 
  xlab("samples") + ylab("log2 count value") + #只做了log2转换，没有vst
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())

# vst标准化箱线图
ggplot(data_vst, aes(x = sample_id, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() + 
  xlab("samples") + ylab("vst transformed value") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())

## ggplot2 density plot

# 未标准化密度图
ggplot(data_counts, aes(x = value, col = sample_id)) +
  geom_density() +
  theme_bw() + 
  xlab("samples") +
  xlim(0, 20) +
  theme(legend.position = "none") 

# vst标准化密度图
ggplot(data_vst, aes(x = value, col = sample_id)) +
  geom_density() +
  theme_bw() + 
  xlim(0, 20) +
  theme(legend.position = "none") 

## ~PCA plot ----
source("custom_functions.R")
# 未标准化
PCA_new(log2(GSE183533_counts), group = GSE183533_targets$group)
# VST标准化
PCA_new(GSE183533_counts_vst, group = GSE183533_targets$group)

# 加载包
library("factoextra")
library("FactoMineR")
# 准备数据
pca1 <- PCA(t(log2(GSE183533_counts+1)), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca1,
             label = "none", 
             habillage = as.factor(GSE183533_targets$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

pca2 <- PCA(t(GSE183533_counts_vst), graph = FALSE)
fviz_pca_ind(pca2,
             label = "none", 
             habillage = as.factor(GSE183533_targets$group),
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
expr1 <- as.matrix(GSE183533_counts)
dist1 <- dist(t(expr1))
hc = hclust(dist1)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc, k = 2, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())

#vst标准化
expr2 <- as.matrix(GSE183533_counts_vst)
dist2 <- dist(t(expr2))
hc = hclust(dist2)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc, k = 2, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()


## ~基于arrayQulitMetrics ----

###  构建ExpressionSet ### 
library(arrayQualityMetrics)

#因为arrayQualityMetrics不支持lumi、ElistRaw格式的数据，此时需要先构建ExpressionSet
#因为arrayQualityMetrics出现NA会提示下标越界，GSE183533_sm里有NA值，此时需要利用处理过NA值的表达矩阵重新构建ExpressionSet

GSE183533_counts_QC <- ExpressionSet(assayData = as.matrix(GSE183533_counts), 
                                         phenoData = AnnotatedDataFrame(data = GSE183533_targets))
arrayQualityMetrics(GSE183533_counts_QC, 
                    outdir = "Covid/GSE183533/GSE183533_counts_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)

GSE183533_counts_vst_QC <- ExpressionSet(assayData = as.matrix(GSE183533_counts_vst), 
                                    phenoData = AnnotatedDataFrame(data = GSE183533_targets))

arrayQualityMetrics(GSE183533_counts_vst_QC, 
                    outdir = "Covid/GSE183533/GSE183533_counts_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)

dev.off()

# 4 第一次去除离群样本 ----
# 4.1 vst标准化前去除离群样本后再标准化，实际上和常规vst标准化之后直接去除离群样本区别不大----
GSE183533_counts_2 <- select(GSE183533_counts,-c(GSM5558935,GSM5558936,GSM5558943,GSM5558953,GSM5558954,GSM5558963,GSM5558972))	
GSE183533_targets_2 <- filter(GSE183533_targets,!rownames(GSE183533_targets)%in%c('GSM5558935','GSM5558936','GSM5558943','GSM5558953','GSM5558954','GSM5558963','GSM5558972'))
## ~vst ----
  #首先看一下有没有缺失值
  sum(is.na(GSE183533_counts_2))
# 无

library(DESeq2)  
group2 <- GSE183533_targets_2$group %>% as.factor()
colData2 = data.frame(sample_id = colnames(GSE183533_counts_2), group2 = group2)
dds2 <- DESeq2::DESeqDataSetFromMatrix(GSE183533_counts_2, #不是整数用round近似取整
                                      colData = colData2,
                                      design = ~ group2)
vst2 <- DESeq2::vst(dds2)
GSE183533_counts_2_vst <- assay(vst2)

save(GSE183533_counts_2_vst, file = "Covid/GSE183533/GSE183533_counts_2_vst.Rdata")
## ~质量评估可视化----

#箱线图、密度图
boxplot(GSE183533_counts_2_vst, las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(GSE183533_counts_2_vst , legend = F)

#PCA图
PCA_new(GSE183533_counts_2_vst, group = GSE183533_targets_2$group)

#fviz_PCA
library("factoextra")
library("FactoMineR")
# 准备数据
pca3 <- PCA(t(GSE183533_counts_2_vst), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca3,
             label = "none", 
             habillage = as.factor(GSE183533_targets_2$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

#聚类树状图
expr3 <- as.matrix(GSE183533_counts_2_vst)
dist3 <- dist(t(expr3))
hc3 = hclust(dist3)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc3, k = 2, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()


## ~基于arrayQulitMetrics----
GSE183533_counts_2_vst_QC <- ExpressionSet(assayData = as.matrix(GSE183533_counts_2_vst), 
                                         phenoData = AnnotatedDataFrame(data = GSE183533_targets_2))

arrayQualityMetrics(GSE183533_counts_2_vst_QC, 
                    outdir = "Covid/GSE183533/GSE183533_counts_2_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)



# 4.2 vst标准化之后去除离群样本 ----
GSE183533_counts_22_vst <- select(as.data.frame(GSE183533_counts_vst),-c(GSM5558935,GSM5558936,GSM5558943,GSM5558953,GSM5558954,GSM5558963,GSM5558972))	
GSE183533_targets_22 <- filter(GSE183533_targets,!rownames(GSE183533_targets)%in%c('GSM5558935','GSM5558936','GSM5558943','GSM5558953','GSM5558954','GSM5558963','GSM5558972'))

## ~质量评估可视化----

#箱线图、密度图
boxplot(GSE183533_counts_22_vst, las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(GSE183533_counts_22_vst , legend = F)

#PCA图
PCA_new(GSE183533_counts_22_vst, group = GSE183533_targets_22$group)

#fviz_PCA
library("factoextra")
library("FactoMineR")
# 准备数据
pca4 <- PCA(t(GSE183533_counts_22_vst), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca4,
             label = "none", 
             habillage = as.factor(GSE183533_targets_22$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

#聚类树状图
expr4 <- as.matrix(GSE183533_counts_22_vst)
dist4 <- dist(t(expr4))
hc4 = hclust(dist4)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc4, k = 2, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()


## ~基于arrayQulitMetrics----
GSE183533_counts_22_vst_QC <- ExpressionSet(assayData = as.matrix(GSE183533_counts_22_vst), 
                                                   phenoData = AnnotatedDataFrame(data = GSE183533_targets_22))

arrayQualityMetrics(GSE183533_counts_22_vst_QC, 
                    outdir = "Covid/GSE183533/GSE183533_counts_22_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)

# 5 第二次去除离群样本 ----
##还有一个离群样本，再删除一次
# 5.1 vst标准化前去除离群样本后再标准化 ----
GSE183533_counts_3 <- select(as.data.frame(GSE183533_counts),-c(GSM5558935,GSM5558936,GSM5558943,GSM5558953,GSM5558954,GSM5558963,GSM5558972,GSM5558937))	
GSE183533_targets_3 <- filter(GSE183533_targets,!rownames(GSE183533_targets)%in%c('GSM5558935','GSM5558936','GSM5558943','GSM5558953','GSM5558954','GSM5558963','GSM5558972','GSM5558937'))
## ~vst ----
#首先看一下有没有缺失值
sum(is.na(GSE183533_counts_3))
# 无

library(DESeq2)  
group3 <- GSE183533_targets_3$group %>% as.factor()
colData3 = data.frame(sample_id = colnames(GSE183533_counts_3), group3 = group3)
dds3 <- DESeq2::DESeqDataSetFromMatrix(GSE183533_counts_3, #不是整数用round近似取整
                                       colData = colData3,
                                       design = ~ group3)
vst3 <- DESeq2::vst(dds3)
GSE183533_counts_3_vst <- assay(vst3)

save(GSE183533_counts_3_vst, file = "Covid/GSE183533/GSE183533_counts_3_vst.Rdata")
## ~质量评估可视化----

#箱线图、密度图
boxplot(GSE183533_counts_3_vst, las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(GSE183533_counts_3_vst , legend = F)

#PCA图
PCA_new(GSE183533_counts_3_vst, group = GSE183533_targets_3$group)

#fviz_PCA
library("factoextra")
library("FactoMineR")
# 准备数据
pca4 <- PCA(t(GSE183533_counts_3_vst), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca4,
             label = "none", 
             habillage = as.factor(GSE183533_targets_3$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

#聚类树状图
expr4 <- as.matrix(GSE183533_counts_3_vst)
dist4 <- dist(t(expr4))
hc4 = hclust(dist4)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc4, k = 2, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()


## ~基于arrayQulitMetrics----
GSE183533_counts_3_vst_QC <- ExpressionSet(assayData = as.matrix(GSE183533_counts_3_vst), 
                                           phenoData = AnnotatedDataFrame(data = GSE183533_targets_3))

arrayQualityMetrics(GSE183533_counts_3_vst_QC, 
                    outdir = "Covid/GSE183533/GSE183533_counts_3_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)

# 5.2 vst标准化之后去除离群样本 ----
GSE183533_counts_33_vst <- select(as.data.frame(GSE183533_counts_vst),-c(GSM5558935,GSM5558936,GSM5558943,GSM5558953,GSM5558954,GSM5558963,GSM5558972,GSM5558937))	
GSE183533_targets_33 <- filter(GSE183533_targets,!rownames(GSE183533_targets)%in%c('GSM5558935','GSM5558936','GSM5558943','GSM5558953','GSM5558954','GSM5558963','GSM5558972','GSM5558937'))

## ~质量评估可视化----

#箱线图、密度图
boxplot(GSE183533_counts_33_vst, las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(GSE183533_counts_33_vst , legend = F)

#PCA图
PCA_new(GSE183533_counts_33_vst, group = GSE183533_targets_33$group)

#fviz_PCA
library("factoextra")
library("FactoMineR")
# 准备数据
pca5 <- PCA(t(GSE183533_counts_33_vst), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca5,
             label = "none", 
             habillage = as.factor(GSE183533_targets_33$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

#聚类树状图
expr5 <- as.matrix(GSE183533_counts_33_vst)
dist5 <- dist(t(expr5))
hc5 = hclust(dist5)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc5, k = 2, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set1"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()



## ~基于arrayQulitMetrics----

GSE183533_counts_33_vst_QC <- ExpressionSet(assayData = as.matrix(GSE183533_counts_33_vst), 
                                                    phenoData = AnnotatedDataFrame(data = GSE183533_targets_33))

arrayQualityMetrics(GSE183533_counts_33_vst_QC, 
                    outdir = "Covid/GSE183533/GSE183533_counts_33_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)
dev.off()







# 6 DESeq2 差异表达分析 ----

## ~构建dds对象 ----
library(DESeq2)  
counts <- select(as.data.frame(GSE183533_counts),-c(GSM5558935,GSM5558936,GSM5558943,GSM5558953,GSM5558954,GSM5558963,GSM5558972,GSM5558937))	
targets <- filter(GSE183533_targets,!rownames(GSE183533_targets)%in%c('GSM5558935','GSM5558936','GSM5558943','GSM5558953','GSM5558954','GSM5558963','GSM5558972','GSM5558937'))

group <- targets$group %>% as.factor()
colData = data.frame(sample_id = colnames(counts), group = group)
dds <- DESeq2::DESeqDataSetFromMatrix(counts, #不是整数用round近似取整
                                       colData = colData,
                                       design = ~ group)
## ~PCA图  ----
vsd <- vst(dds, blind = TRUE)
vsd_df <- assay(vsd) 
head(vsd_df, 5)[, 1:2]

DESeq2::plotPCA(vsd, intgroup = "group")

source("custom_functions.R")
PCA_new(expr = vsd_df, group = targets$group)

## ~remove batch effect---- 
# vsd_batch_rm <- limma::removeBatchEffect(x = vsd_df, batch = targets$batch)
# PCA_new(expr = vsd_batch_rm, group = targets$group)

## ~基因过滤 ----
table(targets$group)
keep <- rowSums(counts(dds) > 0) >= 33
dds_filt <- dds[keep, ]

## ~DESeq一步完成差异分析 ----
time1 <- Sys.time()
dds2 <- DESeq(dds_filt)
runtime1 <- Sys.time() - time1
runtime1   # Time difference of 4.238453 mins

time2 <- Sys.time()
dds2 <- DESeq(dds_filt, parallel = T)
runtime2 <- Sys.time() - time2
runtime2   # Time difference of 3.506404 mins

# save(dds2, file = "Covid/GSE183533/data/dds2.Rda")
load(file = "Covid/GSE183533/data/dds2.Rda")

## ~提取差异分析结果 ----
resultsNames(dds2)
res <- results(dds2)

#genecode中下载最新的v41gtf文件，提取里面gene_id和gene_name用于基因ID转换
#读取gtf文件
gtf_v41 <- rtracklayer::import('Covid/GSE183533/data/gencode.v41.annotation.gtf')
head(gtf_v41)
gtf_v41_df <- as.data.frame(gtf_v41)
colnames(gtf_v41_df)
#提取gene_id和gene_name
geneid_df <- dplyr::select(gtf_v41_df,c(gene_id,gene_name))#,gene_biotype
geneid_df<-unique(geneid_df) #删除列中所有重复向量
#类似ENSG00000284332.1这种末尾带小数点的gene_id是指版本号，有些表达矩阵没有添加这个版本号，删掉即可，这里对应的gtf注释也需要去除小数点
geneid_df$new <- geneid_df$gene_id
new_gtf_v41<-geneid_df
new_gtf_v41$new<-sapply(stringr::str_split(new_gtf_v41$new, "\\."), function(v)  return(v[1]))
new_gtf_v41<-unique(new_gtf_v41)

# 基因ID对应的基因名，开始基因注释
gtf_v41 <- geneid_df

GSE183533_DEG_sig <- as.data.frame(res) %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gtf_v41, by = "gene_id") %>% 
  relocate(gene_name, .after = "gene_id") %>% 
  na.omit() %>% 
  .[!duplicated(.$gene_name),] %>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05) %>% 
  dplyr::select(-c(1,3)) %>% 
  dplyr::rename(., logFC = log2FoldChange )

sum(duplicated(GSE183533_DEG_sig$gene_name))

GSE183533_DEG <- as.data.frame(res) %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gtf_v41, by = "gene_id") %>% 
  relocate(gene_name, .after = "gene_id") %>% 
  na.omit() %>% 
  .[!duplicated(.$gene_name),] %>% 
  dplyr::select(-c(1,3)) %>% 
  dplyr::rename(., logFC = log2FoldChange )

write.table(GSE183533_DEG,
            file = "Covid/GSE183533/results/GSE183533_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(GSE183533_DEG_sig,
            file = "Covid/GSE183533/results/GSE183533_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)




# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. 
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. 


## ~保存分析结果数据 ----
save(res_deseq2,file = "Covid/GSE183533/results/GSE183533_DESeq2_DEG.Rda")

#### End ----

