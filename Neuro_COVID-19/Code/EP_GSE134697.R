############################################################################===
#####  GSE134697_targets+GSE134697_counts
############################################################################===

# 1 GSE134697_targets ----
library(GEOquery)
library(tidyverse)
library(limma)

GSE134697_sm <- getGEO(filename = "Epilepsy/GSE134697/GSE134697_series_matrix.txt.gz",
                    getGPL = F) 

GSE134697_pd <- pData(GSE134697_sm)

GSE134697_targets <- GSE134697_pd %>%
  dplyr::select(2, 1) %>% 
  mutate(sample_id = geo_accession,
         sample_name = GSE134697_sm@phenoData@data[["description"]],
         group= str_extract(sample_name,"\\D+")) %>%
  dplyr::select(sample_id,sample_name,group)


# 写出targets文件
write.table(x = GSE134697_targets,
            file = "Epilepsy/GSE134697/GSE134697_targets.txt",
            quote = F,
            sep = "\t",
            row.names = F)

# 读取targets文件
# GSE134697_targets <- readTargets(file = "Epilepsy/GSE134697/GSE134697_targets.txt",
#                               row.names = "sample_id")

# END
save(GSE134697_targets, file = "Epilepsy/GSE134697/GSE134697_targets.Rdata")


# 2 GSE134697_counts ----
# counts简单清洗
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

GSE134697_counts<- read.table("Epilepsy/GSE134697/GSE134697_CountMatrix.tsv.gz",header=T,row.names = 1)

GSE134697_counts <- GSE134697_counts[,-(1:5)]


colnames(GSE134697_counts)

colnames(GSE134697_counts) <- colnames(GSE134697_counts)%>% 
  str_replace_all(c("X.home.projects.ku_10024.data.Human_data.BamFiles.CHKJ." = "",
                    "X.home.projects.ku_10024.data.Human_data.BamFiles.CTRL.0000" = "Cntrl1", 
                    ".Brain.RNA_D00140_0313.bam" = "")) #这一步其实可以省略，只是列名太长了，还是想先精简下再处理

GSE134697_counts <- GSE134697_counts[,sort(names(GSE134697_counts))]

GSE134697_targets <- arrange(GSE134697_targets,sample_name)

colnames(GSE134697_counts) <- rownames(GSE134697_targets)

# 重新按sample_id排列顺序
GSE134697_counts <- GSE134697_counts[,sort(names(GSE134697_counts))]
GSE134697_targets <- arrange(GSE134697_targets,sample_id)

#筛选非Hippocampus 的样本
GSE134697_targets <- filter(GSE134697_targets,group!= 'H') 
GSE134697_counts <- select(GSE134697_counts,rownames(GSE134697_targets))
identical( colnames(GSE134697_counts),rownames(GSE134697_targets))

GSE134697_targets$group <- str_replace_all(GSE134697_targets$group,c('Cntrl'='con','C'='epilepsy'))

# 样本过滤
keep <- rowSums(GSE134697_counts> 0)>= 10  ##保留至少在50%样本中表达量大于0的基因

head(keep)

GSE134697_counts <- GSE134697_counts[keep,]


# 3 vst ----
#首先看一下有没有缺失值
sum(is.na(GSE134697_counts))
# 无

library(DESeq2)  

GSE134697_targets <- arrange(GSE134697_targets,group)
group <- GSE134697_targets$group %>% as.factor()
colData = data.frame(sample_id = colnames(GSE134697_counts), group = group)
dds <- DESeq2::DESeqDataSetFromMatrix(GSE134697_counts,
                                      colData = colData,
                                      design = ~ group)
vst <- DESeq2::vst(dds)
GSE134697_counts_vst <- assay(vst)

save(GSE134697_counts_vst, file = "Epilepsy/GSE134697/GSE134697_counts_vst.Rdata")
# 4 质量评估可视化
# 查看默认配色
palette.pals() #查看配色方案
palette() #查看当前配色方案
palette("Set 1") #选择一个配色方案
# palette("default") # 重置为默认配色方案”R3“

## ~boxplot & density plot ----
#count raw文件虽然过滤了一大半为0的基因，但是仍然残余部分0的值，log2时会报错，所以整体+1.
boxplot(log2(GSE134697_counts), las = 2, outline = F, col = group)
limma::plotDensities(log2(GSE134697_counts), legend = F)

boxplot(GSE134697_counts_vst, las = 2, outline = F, col = group)
limma::plotDensities(GSE134697_counts_vst , legend = F)

## ggplot2 data   

# 这是芯片数据里的ggplot箱线图数据准备代码，看看区别
box_data <-  log2(GSE134697_counts) %>% #对表达矩阵log2转换，整体+1避免log2(0)报错
  as.data.frame() %>% #tideverse包不接受矩阵，只接受数据框
  rownames_to_column("gene_id") %>% #行名变列，这一列叫probe_id
  slice_sample(n = 50000) %>%  #因为运行速度很慢，所以只挑了50000个数据来运行
  pivot_longer(cols = starts_with("GSM"),#对数据进行长宽转换，从宽数据到长数据，转换以gsm开头的列
               names_to = "sample_id",#从gsm变成sample_id
               values_to = "value") %>% #新建一列叫value
  left_join(GSE134697_targets %>% dplyr::select(sample_id, group), #left_join相当于取交集，以左边的为主
            by = "sample_id") #以sample_id为锚点取交集，多了一列就是group

# 测序数据代码
# 未标准化的counts数据准备
data_counts <- log2(GSE134697_counts) %>%  #对表达矩阵log2转换，整体+1避免log2(0)报错
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE134697_targets %>% select(sample_id, group), by = "sample_id")

# vst标准化后的counts数据准备
data_vst <- GSE134697_counts_vst %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE134697_targets %>% select(sample_id, group), by = "sample_id")

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
PCA_new(log2(GSE134697_counts), group = GSE134697_targets$group)
# VST标准化
PCA_new(GSE134697_counts_vst, group = GSE134697_targets$group)

# 加载包
library("factoextra")
library("FactoMineR")
# 准备数据
pca <- PCA(t(log2(GSE134697_counts+1)), graph = FALSE)
# 作图
fviz_pca_ind(pca,
             label = "none", 
             habillage = as.factor(GSE134697_targets$group),
             palette = c("#00AFBB", "#E7B800","#E6B400"),
             addEllipses = TRUE)

## ~聚类树状图 ----
expr <- as.matrix(GSE134697_counts_vst)
dist <- dist(t(expr))
hc = hclust(dist)

library(factoextra)
library(RColorBrewer)

pdf("data/clust.pdf", height = 10, width = 20)
fviz_dend(hc, k = 2, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set2"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()


## ~基于arrayQulitMetrics ----

###  构建ExpressionSet ### 
library(arrayQualityMetrics)

GSE134697_targets <- arrange(GSE134697_targets,sample_id)

# rownames(GSE134697_targets) <- colnames(GSE134697_counts) 顺序不一致不要直接用，只是换列名，这列数据是没有变换的

#因为arrayQualityMetrics不支持lumi、ElistRaw格式的数据，此时需要先构建ExpressionSet
#因为arrayQualityMetrics出现NA会提示下标越界，GSE134697_sm里有NA值，此时需要利用处理过NA值的表达矩阵重新构建ExpressionSet
GSE134697_counts_vst_QC <- ExpressionSet(assayData = as.matrix(GSE134697_counts_vst), 
                           phenoData = AnnotatedDataFrame(data = GSE134697_targets))


arrayQualityMetrics(GSE134697_counts_vst_QC, 
                    outdir = "Epilepsy/GSE134697/GSE134697_counts_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)
dev.off()



# 4 DESeq2 差异表达分析 ----
## ~基因过滤 ----
table(GSE134697_targets$group)
keep <- rowSums(counts(dds) > 0) >= 19
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

save(dds2, file = "Epilepsy/GSE134697/data/dds2.Rda")
# load(file = "Epilepsy/GSE134697/data/dds2.Rda")

## ~提取差异分析结果 ----
resultsNames(dds2)
res <- results(dds2)

#genecode中下载最新的v19gtf文件，提取里面gene_id和gene_name用于基因ID转换
#读取gtf文件
gtf_v19 <- rtracklayer::import('Epilepsy/GSE134697/data/gencode.v19.annotation.gtf')
head(gtf_v19)
gtf_v19_df <- as.data.frame(gtf_v19)
colnames(gtf_v19_df)
#提取gene_id和gene_name
geneid_df <- dplyr::select(gtf_v19_df,c(gene_id,gene_name))#,gene_biotype
geneid_df<-unique(geneid_df) #删除列中所有重复向量
#删除gene_id的版本号小数点
geneid_df$gene_id <-  str_replace(geneid_df$gene_id,'\\.\\d+','')

gtf_v19 <- geneid_df

# 基因ID对应的基因名，开始基因注释并排序
GSE134697_DEG_sig <- as.data.frame(res) %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gtf_v19, by = "gene_id") %>% 
  relocate(gene_name, .after = "gene_id") %>% 
  na.omit() %>% 
  .[!duplicated(.$gene_name),] %>%
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05) %>% 
  dplyr::select(-c(1,3)) %>% 
  dplyr::rename(., logFC = log2FoldChange )

sum(duplicated(GSE134697_DEG_sig$gene_name))

GSE134697_DEG <- as.data.frame(res) %>% 
  rownames_to_column("gene_id") %>% 
  left_join(gtf_v19, by = "gene_id") %>% 
  relocate(gene_name, .after = "gene_id") %>% 
  na.omit() %>% 
  .[!duplicated(.$gene_name),] %>%
  dplyr::select(-c(1,3)) %>% 
  dplyr::rename(., logFC = log2FoldChange )

sum(duplicated(GSE134697_DEG$gene_name))

write.table(GSE134697_DEG,
            file = "Epilepsy/GSE134697/results/GSE134697_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(GSE134697_DEG_sig,
            file = "Epilepsy/GSE134697/results/GSE134697_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)

# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. 
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. 


## ~保存分析结果数据 ----
save(GSE134697_DESeq2_DEG,file = "Epilepsy/GSE134697/results/GSE134697_DESeq2_DEG.Rda")

#### End ----
length(intersect(res_deseq2$gene_name, GSE134697_DEG_sig[,1]))
library(VennDiagram)
venn.diagram(x = list("Covid-19" = res_deseq2$gene_name, 
                      "Epilepsy" = GSE134697_DEG_sig[,1]),
             filename = "Epilepsy/GSE134697/results/Venn.jpeg",
             fill = c("green", "red"),
             scaled = F,
             cex = 1.5,
             cat.cex = 1.5,
             main = "Venn Diagram",
             main.cex = 2)

