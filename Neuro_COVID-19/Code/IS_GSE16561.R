############################################################################===
#####  GSE16561 DGE
############################################################################===

# 1 GSE16561_targets ----
library(GEOquery)
library(tidyverse)

GSE16561_sm <- getGEO(filename = "Ischemic Stroke/GSE16561/GSE16561_series_matrix.txt.gz",
                      getGPL = F) 

GSE16561_pd <- pData(GSE16561_sm)

sample_name <- GSE16561_sm@phenoData@data[["title"]]
sample_id <- GSE16561_sm@phenoData@data[["geo_accession"]]
group <- GSE16561_sm@phenoData@data[["description"]]
age <- GSE16561_sm@phenoData@data[["age:ch1"]]
gender <- GSE16561_sm@phenoData@data[["gender:ch1"]]


GSE16561_targets <- data.frame(sample_name = sample_name,
                               sample_id = sample_id,
                               group = group,
                               age = age,
                               gender = gender,
                               row.names = sample_id ) %>% 
                               arrange(.,sample_name)

write.table(x = GSE16561_targets,
            file = "Ischemic Stroke/GSE16561/GSE16561_targets.txt",
            quote = F,
            sep = "\t",
            row.names = F)


# 2 GSE16561_expr_raw ----

GSE16561_expr_raw <- read_tsv('Ischemic Stroke/GSE16561/GSE16561_RAW.txt',) %>% 
                    column_to_rownames('ID_REF') %>% 
                   .[,sort(colnames(.))]

identical(colnames(GSE16561_expr_raw),GSE16561_targets[,'sample_name'])
colnames(GSE16561_expr_raw) <- rownames(GSE16561_targets)

#首先看一下有没有缺失值
sum(is.na(GSE16561_expr_raw))

# 3 芯片数据质量评估 ----
palette("Set 1")
#查看未标准化前的数据
GSE16561_pd$data_processing[1]
#箱线图
boxplot(GSE16561_expr_raw, 
        las = 2, 
        outline = FALSE,
        col = as.factor(GSE16561_targets$group))


# PCA_new
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

PCA_new(log2(GSE16561_expr_raw), 
        nrow(GSE16561_expr_raw), 
        group = GSE16561_targets$group,
        show_name = T)

# arrayQulitMetrics 
library(arrayQualityMetrics)

identical(colnames(GSE16561_expr_raw) ,rownames(GSE16561_targets))

GSE16561_QC_raw <- ExpressionSet(assayData = as.matrix(GSE16561_expr_raw), 
                             phenoData = AnnotatedDataFrame(data = GSE16561_targets))


arrayQualityMetrics(GSE16561_QC_raw, 
                    outdir = "Ischemic Stroke/GSE16561/GSE16561_QC_raw", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)
dev.off()

# 4 log2转换，标准化 ----
# 这里illumina芯片数据无需背景校正，先log2转换，再做标准化
# log2转换
GSE16561_expr_raw <- log2(GSE16561_expr_raw)

# 标准化
library(limma)
GSE16561_expr <- normalizeBetweenArrays(GSE16561_expr_raw)

boxplot(GSE16561_expr, 
        las = 2, 
        outline = FALSE,
        col = as.factor(GSE16561_targets$group))


PCA_new(GSE16561_expr, 
        nrow(GSE16561_expr), 
        group = GSE16561_targets$group,
        show_name = T)

# arrayQulitMetrics 
GSE16561_QC <- ExpressionSet(assayData = GSE16561_expr, 
                                 phenoData = AnnotatedDataFrame(data = GSE16561_targets))


arrayQualityMetrics(GSE16561_QC, 
                    outdir = "Ischemic Stroke/GSE16561/GSE16561_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)
dev.off()

# 5 注释+DGE ----
# 注释 ----
library(illuminaio)
GSE16561_anno  <- as.data.frame(readBGX('Ischemic Stroke/GSE16561/GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx')[1])

p2s <- GSE16561_anno[c('probes.Probe_Id','probes.Symbol')] %>% set_names("probe_id", "symbol")

## 芯片注释的自定义函数
annotate_expr <- function(expr, p2s, fun = mean){#自定义三个参数，p2s的列名是固定的
  library(dplyr)
  library(tibble)
  # expr
  expr <- expr %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("probe_id")#行名变成一列，叫probe_id
  # annotate
  expr_annotated <- p2s %>% 
    dplyr::inner_join(expr, by = "probe_id") %>% #expr和p2s取交集
    na.omit() %>% #去除NA
    dplyr::select(-probe_id) %>% #删掉probe_id这一列
    stats::aggregate(. ~ symbol, data =., FUN = fun) %>% #对于多个探针对应同一个symbol，对这些探针取均值
    tibble::column_to_rownames("symbol")#symbol这一列变成行名
  return(expr_annotated)#返回这个expr_annotated，再赋值给annotate_expr
}

## 芯片注释
GSE16561_expr_anno <- annotate_expr(GSE16561_expr, p2s)
head(GSE16561_expr_anno)[, 1:3]

# limma包差异分析 ---- 
library(limma)

design <- model.matrix(~ group, data = GSE16561_targets) #~ group等同于~ 1 + group，之前combat去除批次里有~ 1 + group，效果是一样的。
colnames(design)
colnames(design) <- c("Con", "ISvsCon")
head(design)

design2 <- model.matrix(~ 0 + group, data = GSE16561_targets) 
colnames(design2)
colnames(design2) <- c("Con", "IS")
head(design2)
contrasts <- makeContrasts("IS-Con", levels = design2)
contrasts

## arrayWeights
aw <- arrayWeights(GSE16561_expr_anno, design = design)
aw
barplot(aw)

## for design 
fit <- lmFit(GSE16561_expr_anno, design, weights = aw)#线性拟合
fit <- eBayes(fit)#贝叶斯
res <- topTable(fit, coef = "ISvsCon")#取最值，ISvsCon代表提取要对比的类型
res

## for design2  本质和for design 一模一样，只是用了管道符
res2 <- lmFit(GSE16561_expr_anno, design2, weights = aw) %>% 
  contrasts.fit(contrasts) %>% 
  eBayes() %>% 
  topTable(coef = "IS-Con")
res2

identical(res, res2)#小数点不完全一致，所以是false

library(dplyr)
sqrt(2) ^ 2 == 2 #2的平方根再平方，然而不等于2
near(sqrt(2) ^ 2, 2)
near(res, res2)#两个非常是接近了就认为相等
all(near(res, res2))

## 所有基因的差异分析结果
GSE16561_DEG <- topTable(fit, coef = "ISvsCon", number = Inf)
GSE16561_DEG_sig <- topTable(fit, coef = "ISvsCon", n = Inf, p.value = 0.05, lfc = 1)
#这里的p为矫正后的p

intersect(res_deseq2$gene_name, rownames(GSE16561_DEG_sig))
library(VennDiagram)
venn.diagram(x = list("Covid-19" = res_deseq2$gene_name, 
                      "IS" = rownames(GSE16561_DEG_sig)),
             filename = "Ischemic Stroke/GSE16561/results/Venn2.jpeg",
             fill = c("green", "red"),
             scaled = F,
             cex = 1.5,
             cat.cex = 1.5,
             main = "Venn Diagram",
             main.cex = 2)




## topTable & topTableF
topTable(fit)#两素比较用的多
topTableF(fit)#多组比较用这个
all(near(topTable(fit), topTableF(fit)))

#### 5 保存结果 ----
save(GSE16561_DEG,GSE16561_DEG_sig, file = "Ischemic Stroke/GSE16561/results/GSE16561_DEG_sig .Rda")

GSE16561_DEG <- GSE16561_DEG %>% rownames_to_column('gene_name')
GSE16561_DEG_sig <- GSE16561_DEG_sig %>% rownames_to_column('gene_name')

write.table(GSE16561_DEG,
            file = "Ischemic Stroke/GSE16561/results/GSE16561_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(GSE16561_DEG_sig,
            file = "Ischemic Stroke/GSE16561/results/GSE16561_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)

#### End ----


