############################################################################===
#####  GSE163256_monos_targets+GSE163256_monos_counts
############################################################################===

# 1 GSE163256_targets ----
library(GEOquery)
library(tidyverse)
library(limma)

GSE163256_sm <- getGEO(filename = "Hemorrhagic stroke/GSE163256/GSE163256_series_matrix.txt.gz",
                       getGPL = F) 

GSE163256_pd <- pData(GSE163256_sm)

group <-  GSE163256_sm@phenoData@data[["characteristics_ch1.2"]] %>%
         str_replace_all(c("patientid: \\d+"="ICH","patientid: HC\\d+"="con"))

GSE163256_targets <- GSE163256_pd %>%
  dplyr::select(2) %>% 
  mutate(sample_id = geo_accession,
         sample_name = GSE163256_sm@phenoData@data[["description"]],
         group= group) %>%
  dplyr::select(sample_id,sample_name,group)


# 写出targets文件
write.table(x = GSE163256_targets,
            file = "Hemorrhagic stroke/GSE163256/GSE163256_targets.txt",
            quote = F,
            sep = "\t",
            row.names = F)

# 读取targets文件
# GSE163256_targets <- readTargets(file = "Hemorrhagic stroke/GSE163256/GSE163256_targets.txt",
#                               row.names = "sample_id")

# END
save(GSE163256_targets, file = "Hemorrhagic stroke/GSE163256/GSE163256_targets.Rdata")


# 2 GSE163256_counts ----
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

# counts清洗方法1：
# 读取counts文件并把第一行当做列名，第一列当做行名。
GSE163256_monos_counts<-read.csv("Hemorrhagic stroke/GSE163256/GSE163256_monos_counts.csv.gz",header=T, row.names = 1)
#转置一下把#列名#变成新的的#列#，且列名为samplename
tmp <-as.data.frame(t(GSE163256_monos_counts)) 
tmp2 <- rownames_to_column(tmp,'sample_name') 
#target文件是包含了所有样本信息，我们只要单核细胞相关的样本信息，取交集筛选出来
GSE163256_monos_targets<- semi_join(x=GSE163256_targets,y=tmp2,by='sample_name')
#筛选出来后排序方便后面匹配
GSE163256_monos_targets <- arrange(GSE163256_monos_targets,sample_name)
#将读取进来的counts文件排序，方便匹配后直接替换列名为gsm id
GSE163256_monos_counts <- GSE163256_monos_counts[,sort(names(GSE163256_monos_counts))]
#将counts文件列名替换为GSM开头的
colnames(GSE163256_monos_counts) <- rownames(GSE163256_monos_targets)

# counts清洗方法2：
# 读取counts文件并把第一行当做列名，第一列当做行名。
GSE163256_monos_counts_2<-read.csv("Hemorrhagic stroke/GSE163256/GSE163256_monos_counts.csv.gz",header=T, row.names = 1)

GSE163256_monos_targets_2 <- GSE163256_targets %>%
  dplyr::filter(sample_name %in% colnames(GSE163256_monos_counts_2)) 

GSE163256_monos_targets_2 <- arrange(GSE163256_monos_targets_2,sample_name)

GSE163256_monos_counts_2 <- GSE163256_monos_counts_2[,sort(names(GSE163256_monos_counts_2))]

colnames(GSE163256_monos_counts_2) <- rownames(GSE163256_monos_targets_2)

# 样本过滤
keep <- rowSums(GSE163256_monos_counts> 0)>= 120  ##保留至少在60%样本中表达量大于0的基因

head(keep)

GSE163256_monos_counts <- GSE163256_monos_counts[keep,]

# save
save(GSE163256_monos_targets, GSE163256_monos_counts,file = "Hemorrhagic stroke/GSE163256/GSE163256_monos_targets+counts.Rdata")

# 3 vst ----
#首先看一下有没有缺失值
sum(is.na(GSE163256_monos_counts))
# 无

library(DESeq2)  
group <- GSE163256_monos_targets$group %>% as.factor()
colData = data.frame(sample_id = colnames(GSE163256_monos_counts), group = group)
DDS <- DESeq2::DESeqDataSetFromMatrix(round(GSE163256_monos_counts), #不是整数用round近似取整
                                      colData = colData,
                                      design = ~ group)
vst <- DESeq2::vst(DDS)
GSE163256_monos_counts_vst <- assay(vst)

save(GSE163256_monos_counts_vst, file = "Hemorrhagic stroke/GSE163256/GSE163256_monos_counts_vst.Rdata")
# 4 质量评估可视化
## ~boxplot & density plot ----
#count raw文件虽然过滤了一大半为0的基因，但是仍然有些counts文件会残余部分0的值，log2时会报错，如果还有0的话可以整体+1，这次的过滤之后没有0.
boxplot(log2(GSE163256_monos_counts), las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(log2(GSE163256_monos_counts), legend = F)

boxplot(GSE163256_monos_counts_vst, las = 2, outline = F, col = group)
limma::plotDensities(GSE163256_monos_counts_vst , legend = F)

## ggplot2 data   

# # 这是芯片数据里的ggplot箱线图数据准备代码，看看区别
# box_data <-  log2(GSE163256_counts + 1) %>% #对表达矩阵log2转换，整体+1避免log2(0)报错
#   as.data.frame() %>% #tideverse包不接受矩阵，只接受数据框
#   rownames_to_column("gene_id") %>% #行名变列，这一列叫probe_id
#   slice_sample(n = 50000) %>%  #因为运行速度很慢，所以只挑了50000个数据来运行
#   pivot_longer(cols = starts_with("GSM"),#对数据进行长宽转换，从宽数据到长数据，转换以gsm开头的列
#                names_to = "sample_id",#从gsm变成sample_id
#                values_to = "value") %>% #新建一列叫value
#   left_join(GSE163256_targets %>% dplyr::select(sample_id, group), #left_join相当于取交集，以左边的为主
#             by = "sample_id") #以sample_id为锚点取交集，多了一列就是group

# 测序数据代码
# 未标准化的counts数据准备
data_counts <- log2(GSE163256_monos_counts) %>%  #对表达矩阵log2转换，如有0，就整体+1避免log2(0)报错
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE163256_monos_targets %>% select(sample_id, group), by = "sample_id")

# vst标准化后的counts数据准备
data_vst <- GSE163256_monos_counts_vst %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE163256_monos_targets %>% select(sample_id, group), by = "sample_id")

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
source("codes/custom_functions.R")
# 未标准化
PCA_new(log2(GSE163256_monos_counts), group = GSE163256_monos_targets$group)
# VST标准化
PCA_new(GSE163256_monos_counts_vst, group = GSE163256_monos_targets$group)

# 加载包
library("factoextra")
library("FactoMineR")
# 准备数据
pca <- PCA(t(log2(GSE163256_monos_counts+1)), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca,
             label = "none", 
             habillage = as.factor(GSE163256_monos_targets$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

## ~聚类树状图 ----
expr <- as.matrix(GSE163256_monos_counts_vst)
dist <- dist(t(expr))
hc = hclust(dist)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc, k = 6, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set2"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()


## ~基于arrayQulitMetrics ----

###  构建ExpressionSet ### 
library(arrayQualityMetrics)

GSE163256_monos_targets <- arrange(GSE163256_monos_targets,sample_id)

#因为arrayQualityMetrics不支持lumi、ElistRaw格式的数据，此时需要先构建ExpressionSet
#因为arrayQualityMetrics出现NA会提示下标越界，GSE163256_sm里有NA值，此时需要利用处理过NA值的表达矩阵重新构建ExpressionSet
GSE163256_monos_QC <- ExpressionSet(assayData = as.matrix(GSE163256_monos_counts_vst), 
                              phenoData = AnnotatedDataFrame(data = GSE163256_monos_targets))


arrayQualityMetrics(GSE163256_monos_QC, 
                    outdir = "Hemorrhagic stroke/GSE163256/GSE163256_QC_monos_counts", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)
dev.off()


