############################################################################===
#####  GSE163256_neuts_targets+GSE163256_neuts_counts+DEG
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

batch <- GSE163256_sm@phenoData@data[["characteristics_ch1.1"]] %>% str_replace(.,"batch: ","")

race <- GSE163256_sm@phenoData@data[["characteristics_ch1.7"]] %>% str_replace_all(c('race: '='','White'='W','Asian'='Y','African-American'='B'))

GSE163256_targets <- GSE163256_pd %>%
  dplyr::select(2) %>% 
  mutate(sample_id = geo_accession,
         sample_name = GSE163256_sm@phenoData@data[["description"]],
         group = group,
         batch = batch,
         race = race) %>%
  dplyr::select(sample_id:race)


# 写出targets文件
write.table(x = GSE163256_targets,
            file = "Hemorrhagic stroke/GSE163256/GSE163256_targets.txt",
            quote = F,
            sep = "\t",
            row.names = F)

# 读取targets文件
GSE163256_targets <- readTargets(file = "Hemorrhagic stroke/GSE163256/GSE163256_targets.txt",
                              row.names = "sample_id")

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
GSE163256_neuts_counts<-read.csv("Hemorrhagic stroke/GSE163256/GSE163256_neuts_counts.csv.gz",header=T, row.names = 1)
#转置一下把#列名#变成新的的#列#，且列名为samplename
tmp <-as.data.frame(t(GSE163256_neuts_counts)) 
tmp2 <- rownames_to_column(tmp,'sample_name') 
#target文件是包含了所有样本信息，我们只要中性粒细胞相关的样本信息，取交集筛选出来
GSE163256_neuts_targets<- semi_join(x=GSE163256_targets,y=tmp2,by='sample_name')
#筛选出来后排序方便后面匹配
GSE163256_neuts_targets <- arrange(GSE163256_neuts_targets,sample_name)
#将读取进来的counts文件排序，方便匹配后直接替换列名为gsm id
GSE163256_neuts_counts <- GSE163256_neuts_counts[,sort(names(GSE163256_neuts_counts))]

identical(GSE163256_neuts_targets[,2],colnames(GSE163256_neuts_counts))
#将counts文件列名替换为GSM开头的
colnames(GSE163256_neuts_counts) <- rownames(GSE163256_neuts_targets)

# counts清洗方法2：
# 读取counts文件并把第一行当做列名，第一列当做行名。
GSE163256_neuts_counts_2<-read.csv("Hemorrhagic stroke/GSE163256/GSE163256_neuts_counts.csv.gz",header=T, row.names = 1)

GSE163256_neuts_targets_2 <- GSE163256_targets %>%
  dplyr::filter(sample_name %in% colnames(GSE163256_neuts_counts_2)) 

GSE163256_neuts_targets_2 <- arrange(GSE163256_neuts_targets_2,sample_name)

GSE163256_neuts_counts_2 <- GSE163256_neuts_counts_2[,sort(names(GSE163256_neuts_counts_2))]

colnames(GSE163256_neuts_counts_2) <- rownames(GSE163256_neuts_targets_2)

# 样本过滤
keep <- rowSums(GSE163256_neuts_counts> 0)>= 120  ##保留至少在60%样本中表达量大于0的基因

head(keep)

GSE163256_neuts_counts <- GSE163256_neuts_counts[keep,]

# save
save(GSE163256_neuts_targets, GSE163256_neuts_counts,file = "Hemorrhagic stroke/GSE163256/GSE163256_neuts_targets+counts.Rdata")

# 3 vst ----
#首先看一下有没有缺失值
sum(is.na(GSE163256_neuts_counts))
# 无

library(DESeq2)  
group <- GSE163256_neuts_targets$group %>% as.factor()
colData = data.frame(sample_id = colnames(GSE163256_neuts_counts), group = group)
dds <- DESeq2::DESeqDataSetFromMatrix(round(GSE163256_neuts_counts), #不是整数用round近似取整
                                      colData = colData,
                                      design = ~ group)
vst <- DESeq2::vst(dds)
GSE163256_neuts_counts_vst <- assay(vst)

save(GSE163256_neuts_counts_vst, file = "Hemorrhagic stroke/GSE163256/GSE163256_neuts_counts_vst.Rdata")
# 4 质量评估可视化
## ~boxplot & density plot ----
#count raw文件虽然过滤了一大半为0的基因，但是仍然有些counts文件会残余部分0的值，log2时会报错，如果还有0的话可以整体+1，这次的过滤之后没有0.
boxplot(log2(GSE163256_neuts_counts), las = 2, outline = F, col = c("red", "blue"))
limma::plotDensities(log2(GSE163256_neuts_counts), legend = F)

boxplot(GSE163256_neuts_counts_vst, las = 2, outline = F, col = group)
limma::plotDensities(GSE163256_neuts_counts_vst , legend = F)

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
data_counts <- log2(GSE163256_neuts_counts) %>%  #对表达矩阵log2转换，如有0，就整体+1避免log2(0)报错
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE163256_neuts_targets %>% select(sample_id, group), by = "sample_id")

# vst标准化后的counts数据准备
data_vst <- GSE163256_neuts_counts_vst %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "sample_id") %>% 
  left_join(GSE163256_neuts_targets %>% select(sample_id, group), by = "sample_id")

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
PCA_new(log2(GSE163256_neuts_counts), group = GSE163256_neuts_targets$group)
# VST标准化
PCA_new(GSE163256_neuts_counts_vst, group = GSE163256_neuts_targets$group)

# 加载包
library("factoextra")
library("FactoMineR")
# 准备数据
pca <- PCA(t(log2(GSE163256_neuts_counts+1)), graph = FALSE) #这里+1，不+1提示无穷值
# 作图
fviz_pca_ind(pca,
             label = "none", 
             habillage = as.factor(GSE163256_neuts_targets$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)

pca2 <- PCA(t(GSE163256_neuts_counts_vst), graph = FALSE)
fviz_pca_ind(pca2,
             label = "none", 
             habillage = as.factor(GSE163256_neuts_targets$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)
## ~聚类树状图 ----
expr <- as.matrix(GSE163256_neuts_counts_vst)
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

# GSE163256_neuts_targets <- arrange(GSE163256_neuts_targets,sample_name)

#因为arrayQualityMetrics不支持lumi、ElistRaw格式的数据，此时需要先构建ExpressionSet
#因为arrayQualityMetrics出现NA会提示下标越界，GSE163256_sm里有NA值，此时需要利用处理过NA值的表达矩阵重新构建ExpressionSet
GSE163256_neuts_counts_vst_QC <- ExpressionSet(assayData = as.matrix(GSE163256_neuts_counts_vst), 
                                    phenoData = AnnotatedDataFrame(data = GSE163256_neuts_targets))


arrayQualityMetrics(GSE163256_neuts_counts_vst_QC, 
                    outdir = "Hemorrhagic stroke/GSE163256/GSE163256_neuts_counts_vst_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)
dev.off()


# 4 纳入 batch effect ----
group <- GSE163256_neuts_targets$group %>% as.factor()
batch <- GSE163256_neuts_targets$batch %>% as.factor()
race <- GSE163256_neuts_targets$race %>% as.factor()
colData_batch = data.frame(sample_id = colnames(GSE163256_neuts_counts), group = group, batch = batch, race = race)
dds_batch <- DESeqDataSetFromMatrix(countData = round(GSE163256_neuts_counts), 
                                    colData = colData_batch,
                                    design = ~ batch + race + group)

# vst 
vst_batch <- DESeq2::vst(dds_batch)
GSE163256_neuts_counts_vst_batch <- assay(vst_batch)

## ~质量评估可视化 （去不去批次样本箱线图、pca没有区别，可能只适合在做差异分析时候用）----
boxplot(GSE163256_neuts_counts_vst_batch, las = 2, outline = F, col = group)

limma::plotDensities(GSE163256_neuts_counts_vst_batch , legend = F)

PCA_new(GSE163256_neuts_counts_vst_batch, group = GSE163256_neuts_targets$group)

pca3 <- PCA(t(GSE163256_neuts_counts_vst), graph = FALSE)
fviz_pca_ind(pca3,
             label = "none", 
             habillage = as.factor(GSE163256_neuts_targets$group),
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE)


expr2 <- as.matrix(GSE163256_neuts_counts_vst_batch)
dist2 <- dist(t(expr2))
hc2 = hclust(dist2)

library(factoextra)
library(RColorBrewer)

fviz_dend(hc2, k = 6, 
          cex = 0.5,
          k_colors = brewer.pal(4, "Set2"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()

## ~基于arrayQulitMetrics ----

###  构建ExpressionSet ### 
library(arrayQualityMetrics)

GSE163256_neuts_counts_vst_batch_QC <- ExpressionSet(assayData = as.matrix(GSE163256_neuts_counts_vst_batch), 
                                           phenoData = AnnotatedDataFrame(data = GSE163256_neuts_targets))

arrayQualityMetrics(GSE163256_neuts_counts_vst_batch_QC, 
                    outdir = "Hemorrhagic stroke/GSE163256/GSE163256_neuts_counts_vst_batch_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)
dev.off()


# 5 去除离群样本----
# 用vst标准化之后的先看下去除离群值效果 ----
GSE163256_neuts_counts_vst_23 <- select(as.data.frame(GSE163256_neuts_counts_vst),-c(GSM4976148,GSM4976145,GSM4976149, 
                                                                                     GSM4976146,GSM4976157,GSM4976152, 
                                                                                     GSM4976153,GSM4976002,GSM4975998,
                                                                                     GSM4976004,GSM4976105,GSM4976106,
                                                                                     GSM4976109,GSM4976112,GSM4976134, 
                                                                                     GSM4976135,GSM4976137,GSM4976142,
                                                                                     GSM4975964,GSM4976082,GSM4975992,  
                                                                                     GSM4975991,GSM4976046,GSM4976015, 
                                                                                     GSM4976016,GSM4976069,GSM4976098,
                                                                                     GSM4976150,GSM4976131,GSM4975963,
                                                                                     GSM4976085,GSM4976021,GSM4976054,
                                                                                     GSM4976061,GSM4976155))	

GSE163256_neuts_targets_23 <- filter(GSE163256_neuts_targets,!rownames(GSE163256_neuts_targets)%in%outliers)
identical(colnames(GSE163256_neuts_counts_vst_23),rownames(GSE163256_neuts_targets_23))

library(arrayQualityMetrics)

GSE163256_neuts_counts_vst_23_QC <- ExpressionSet(assayData = as.matrix(GSE163256_neuts_counts_vst_23), 
                                                  phenoData = AnnotatedDataFrame(data = GSE163256_neuts_targets_23))

arrayQualityMetrics(GSE163256_neuts_counts_vst_23_QC, 
                    outdir = "Hemorrhagic stroke/GSE163256/GSE163256_neuts_counts_vst_23_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)
dev.off()

boxplot(GSE163256_neuts_counts_2_VST, las = 2, outline = F, col = group)

limma::plotDensities(GSE163256_neuts_counts_vst_23 , legend = F)

PCA_new(GSE163256_neuts_counts_vst_23, group = GSE163256_neuts_targets_23$group)
# 6 DESeq2 差异表达分析 ----
## ~构建DDS对象 （纳入批次）----
outliers <- c('GSM4976148','GSM4976145','GSM4976149','GSM4976146','GSM4976157',
              'GSM4976152','GSM4976153', 'GSM4976002','GSM4975998','GSM4976004',
              'GSM4976105','GSM4976106','GSM4976109', 'GSM4976112', 'GSM4976134', 
              'GSM4976135', 'GSM4976137','GSM4976142','GSM4975964','GSM4976082',
              'GSM4975992', 'GSM4975991','GSM4976046','GSM4976015','GSM4976016',
              'GSM4976069', 'GSM4976098','GSM4976150','GSM4976131','GSM4975963',
              'GSM4976085','GSM4976021','GSM4976054','GSM4976061','GSM4976155')

GSE163256_neuts_counts_2 <- select(GSE163256_neuts_counts,-c(GSM4976148,GSM4976145,GSM4976149, 
                                                             GSM4976146,GSM4976157,GSM4976152, 
                                                             GSM4976153,GSM4976002,GSM4975998,
                                                             GSM4976004,GSM4976105,GSM4976106,
                                                             GSM4976109,GSM4976112,GSM4976134, 
                                                             GSM4976135,GSM4976137,GSM4976142,
                                                             GSM4975964,GSM4976082,GSM4975992,  
                                                             GSM4975991,GSM4976046,GSM4976015, 
                                                             GSM4976016,GSM4976069,GSM4976098,
                                                             GSM4976150,GSM4976131,GSM4975963,
                                                             GSM4976085,GSM4976021,GSM4976054,
                                                             GSM4976061,GSM4976155))	

GSE163256_neuts_targets_2 <- filter(GSE163256_neuts_targets,!rownames(GSE163256_neuts_targets)%in%outliers)

identical(colnames(GSE163256_neuts_counts_2),rownames(GSE163256_neuts_targets_2))

GSE163256_neuts_targets_2 <- arrange(GSE163256_neuts_targets_2,group)

group <- GSE163256_neuts_targets_2$group %>% as.factor()
batch <- GSE163256_neuts_targets_2$batch %>% as.factor()
race <- GSE163256_neuts_targets_2$race %>% as.factor()
colData = data.frame(sample_id = colnames(GSE163256_neuts_counts_2), group = group, batch = batch, race = race)

#不纳入批次
dds <- DESeqDataSetFromMatrix(countData = round(GSE163256_neuts_counts_2), 
                              colData = colData,
                              design = ~ group)
# vst 
# PCA图
vsd <- vst(dds, blind = TRUE)
vsd_df <- assay(vsd) 
head(vsd_df, 5)[, 1:2]

DESeq2::plotPCA(vsd, intgroup = "group")

source("custom_functions.R")
PCA_new(expr = vsd_df, group = GSE163256_neuts_targets_2$group)

# remove batch effect 看下去除批次的效果如何，似乎比之前更好点，下面可以在DDS中纳入批次进行分析
vsd_batch_rm <- limma::removeBatchEffect(x = vsd_df, batch = GSE163256_neuts_targets_2$batch,batch2 = GSE163256_neuts_targets_2$race)
PCA_new(expr = vsd_batch_rm, group = GSE163256_neuts_targets_2$group)

boxplot(vsd_batch_rm, las = 2, outline = F, col = group)

#纳入批次
DDS <- DESeqDataSetFromMatrix(countData = round(GSE163256_neuts_counts_2), 
                                    colData = colData,
                                    design = ~ batch + race + group)

# vst 
# PCA图
VSD <- vst(DDS, blind = TRUE)
VSD_DF <- assay(VSD) 
head(vsd_df, 5)[, 1:2]

DESeq2::plotPCA(VSD, intgroup = "group")

source("custom_functions.R")
PCA_new(expr = VSD_DF, group = GSE163256_neuts_targets_2$group)

# remove batch effect 看下去除批次的效果如何，似乎比之前更好点
VSD_batch_rm <- limma::removeBatchEffect(x = VSD_DF, batch = GSE163256_neuts_targets_2$batch,batch2 = GSE163256_neuts_targets_2$race)
PCA_new(expr = VSD_batch_rm, group = GSE163256_neuts_targets_2$group)

boxplot(VSD_batch_rm, las = 2, outline = F, col = group)
## ~基因过滤 ----
table(GSE163256_neuts_targets_2$group)
keep <- rowSums(counts(DDS) > 0) >= 166
DDS_filt <- DDS[keep, ]

## ~DESeq一步完成差异分析 ----
time1 <- Sys.time()
DDS2 <- DESeq(DDS_filt)
runtime1 <- Sys.time() - time1
runtime1   # Time difference of 4.238453 mins

time2 <- Sys.time()
DDS2 <- DESeq(DDS_filt, parallel = T)
runtime2 <- Sys.time() - time2
runtime2   # Time difference of 3.506404 mins

save(DDS2, file = "Hemorrhagic stroke/GSE163256/data/dds2.Rda")
load(file = "Hemorrhagic stroke/GSE163256/data/DDS2.Rda")

## ~提取差异分析结果 ----
resultsNames(DDS2)
res <- results(DDS2)

GSE163256_DEG <- as.data.frame(res) %>%
  rownames_to_column("gene_name") %>% 
  na.omit() %>% dplyr::select(-2) %>% 
  dplyr::rename(., logFC = log2FoldChange )

GSE163256_DEG_sig <- as.data.frame(res) %>%
  rownames_to_column("gene_name") %>% 
  na.omit() %>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05) %>% 
  dplyr::select(-2) %>% 
  dplyr::rename(., logFC = log2FoldChange )

write.table(GSE163256_DEG,
            file = "Hemorrhagic stroke/GSE163256/results/GSE163256_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(GSE163256_DEG_sig,
            file = "Hemorrhagic stroke/GSE163256/results/GSE163256_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)


# #从临床信息中得知这里的人类基因版本是hg19，对应的是v19gtf文件，提取里面gene_id和gene_name用于基因ID转换
# #读取gtf文件
# gtf_v19 <- rtracklayer::import('Hemorrhagic stroke/GSE163256/data/gencode.v19.annotation.gtf')
# head(gtf_v19)
# gtf_v19_df <- as.data.frame(gtf_v19)
# colnames(gtf_v19_df)
# #提取gene_id和gene_name
# geneid_df <- dplyr::select(gtf_v19_df,c(gene_id,gene_name))#,gene_biotype
# geneid_df<-unique(geneid_df) #删除列中所有重复向量
# #类似ENSG00000284332.1这种末尾带小数点的gene_id是指版本号，有些表达矩阵没有添加这个版本号，删掉即可，这里对应的gtf注释也需要去除小数点
# geneid_df$new <- geneid_df$gene_id
# new_gtf_v19<-geneid_df
# new_gtf_v19$new<-sapply(stringr::str_split(new_gtf_v19$new, "\\."), function(v)  return(v[1]))
# new_gtf_v19<-unique(new_gtf_v19)
# # 基因ID对应的基因名，开始基因注释
# gtf_v19 <- geneid_df
# 
# res_deseq2 <- as.data.frame(res) %>% 
#   rownames_to_column("gene_id") %>% 
#   left_join(gtf_v19, by = "gene_id") %>% 
#   relocate(gene_name, .after = "gene_id") %>% 
#   na.omit() %>% 
#   arrange(padj) %>% 
#   dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05) 

# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. 
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. 


## ~保存分析结果数据 ----
save(GSE163256_DESeq2_DEG,file = "Hemorrhagic stroke/GSE163256/results/GSE163256_DESeq2_DEG.Rda")

#### End ----

# 新冠和脑出血的交集基因

length(intersect(GSE163256_DEG_sig$gene_name, res_deseq2$gene_name))

intersect(GSE163256_DEG_sig$gene_name, res_deseq2$gene_name)
library(VennDiagram)
venn.diagram(x = list("ICH" = GSE163256_DEG_sig$gene_name, 
                      "Covid" = res_deseq2$gene_name),
             filename = "Hemorrhagic stroke/GSE163256/results/Venn_ICH&Covid.jpeg",
             fill = c("blue", "red"),
             scaled = F,
             cex = 1.5,
             cat.cex = 1.5,
             main = "Venn Diagram",
             main.cex = 2)
 
