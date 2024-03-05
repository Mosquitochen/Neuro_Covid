
# 数据读取 ----
library(tidyverse)
library(readr)

expr <- read_tsv('Covid/E-MTAB-8871/nCOV2_Norm_log2_counts.txt') %>% column_to_rownames('Subject ID')%>% as.data.frame()
idf <- read_tsv('Covid/E-MTAB-8871/E-MTAB-8871/E-MTAB-8871.idf.txt') %>% as.data.frame()
sdrf <- read_tsv('Covid/E-MTAB-8871/E-MTAB-8871/E-MTAB-8871.sdrf.txt') %>% as.data.frame()

colnames(expr)
expr <- expr[-(1:2),] 
str(expr)
# 转换一下结构不然做不了后面的分析
expr[,1:32] <- as.numeric(unlist(expr[,1:32]))

sum(duplicated(rownames(expr)))

# targets ----
sample_name <- sdrf[,1]
group <- sdrf[,11] #作者上传的标准化的数据应该删掉一个qc不合格的所以少了一个样本
Filename <- sdrf[,25]

# 应该是过滤了一个质量不合格的样本
targets <- data.frame(sample_name = sample_name,
                     group = group,
                     Filename = Filename,
                     row.names = Filename ) 


# counts raw ----
counts_raw <- read_csv('Covid/E-MTAB-8871/counts_raw.csv') [-1,-2] %>% column_to_rownames('Probe Name')
colnames(counts_raw)
counts_raw <- head(counts_raw,-14)
# nanostring质量控制的时候剔除了一个
setdiff(rownames(targets),colnames(counts_raw))
targets <- filter(targets,rownames(targets)!='20200130_30102306760221-01_55-20200126_04.RCC') %>% 
     arrange(.,desc(group))

counts_raw <- counts_raw[,rownames(targets)]
identical(colnames(counts_raw),rownames(targets))
# 质量评估 ----
# 箱线图
palette("Set 1") 
boxplot(log2(counts_raw), las = 2, outline = F, col = as.factor(targets$group))


# PCA
PCA_new <- function(expr = expr, ntop = 500, group = targets$group, show_name = F){
  library(ggplot2)
  library(ggrepel)
  object <- expr
  rv <- genefilter::rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1], 
                  PC2 = pca$x[, 2], 
                  group = group, 
                  name = colnames(object))
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

PCA_new(counts_raw, 
        ntop = nrow(counts_raw),
        group = targets$group,
        show_name = F)

# 差异分析


# 6 DESeq2 差异表达分析 ----

## ~构建dds对象 ----
library(DESeq2)  

group <- targets$group %>% as.factor()
colData = data.frame(Filename = colnames(counts_raw), group = group)
dds <- DESeq2::DESeqDataSetFromMatrix(counts_raw, #不是整数用round近似取整
                                      colData = colData,
                                      design = ~ group)
## ~PCA图  ----
vsd <- vst(dds, blind = TRUE) # nsub	the number of genes to subset to (default 1000)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE) 
vsd_df <- assay(vsd) 
head(vsd_df, 5)[, 1:2]

boxplot(vsd_df, las = 2, outline = F, col = as.factor(targets$group))


DESeq2::plotPCA(vsd, intgroup = "group")

source("custom_functions.R")
PCA_new(expr = vsd_df, group = targets$group)

## ~DESeq一步完成差异分析 ----
dds$group<- relevel(dds$group, ref = "normal") # 指定哪一组作为对照组
dds2 <- DESeq(dds)

## ~提取差异分析结果 ----
resultsNames(dds2)
res <- results(dds2)



Imm_DEG_sig <- as.data.frame(res) %>% 
  rownames_to_column("gene_name") %>% 
  na.omit() %>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05) %>% 
  dplyr::select(-2) %>% 
  dplyr::rename(., logFC = log2FoldChange )


Imm_DEG <- as.data.frame(res) %>% 
  rownames_to_column("gene_name") %>% 
  na.omit() %>% 
  dplyr::select(-2) %>% 
  dplyr::rename(., logFC = log2FoldChange )

write.table(Imm_DEG ,
            file = "Covid/E-MTAB-8871/results/Imm_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(Imm_DEG_sig,
            file = "Covid/E-MTAB-8871/results/Imm_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)




# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. 
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. 


## ~保存分析结果数据 ----
save(Imm_DEG_sig,file = "Covid/GSE183533/results/Imm_DEG_sig.Rda")

#### End ----


