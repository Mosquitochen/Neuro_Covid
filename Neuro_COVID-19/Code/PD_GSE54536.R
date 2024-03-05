############################################################################===
#####  GSE54536_targets+GSE54536_QC
############################################################################===

# 1 GSE54536_non_norm读取 ----
library(limma)
library(readr)
preview <- read_tsv("PD/GSE54536/GSE54536_Non-normalized_data.txt.gz")
colnames(preview)

GSE54536_non_norm <- read.ilmn(files = "PD/GSE54536/GSE54536_Non-normalized_data.txt.gz",
                               expr = "SAMPLE",#这里要根据数据实际的样子来改，有时候又变成了sample全
                               other.columns = "Detection",
                               probeid = "ID_REF")

head(GSE54536_non_norm$E)
dim(GSE54536_non_norm$E)

## save data 
save(GSE54536_non_norm, file = "PD/GSE54536/GSE54536_non-normalized.Rdata")

# 2 GSE54536质量评估及标准化 ----
library(limma)
library(tidyverse) 

## 2.1 GSE54536 data import ----
###  2.1.1 load data ----
load_input <- load("PD/GSE54536/GSE54536_non-normalized.Rdata")
load_input

###  2.1.2 GSE54536_targets ----
library(GEOquery)

GSE54536_sm <- getGEO(filename = "PD/GSE54536/GSE54536_series_matrix.txt.gz",
                      getGPL = F) 

GSE54536_pd <- pData(GSE54536_sm) 

group <- GSE54536_sm@phenoData@data[["diagnosis:ch1"]] %>% 
  gsub(pattern="Parkinson's disease",replacement="PD",.) %>% 
  str_replace_all(.,"healthy","con")

sample_name <- paste0(group, "_", 1:5) %>% str_replace_all(.,"5","1-5")

GSE54536_targets <- GSE54536_pd %>%
  dplyr::select(2) %>%
  mutate(sample_id = geo_accession,
         sample_name = sample_name,
         group = group) %>%
  dplyr::select(sample_id:group)

# 保存
write_tsv(x = GSE54536_targets, 
          file = "PD/GSE54536/GSE54536_targets.txt")

# END
save(GSE54536_targets, file = "PD/GSE54536/GSE54536_targets.Rdata")

rownames(GSE54536_targets) <- colnames(GSE54536_non_norm$E) <- GSE54536_targets$sample_id

## 2.2 Plots before QC ----


# 查看默认配色
palette.pals() #查看配色方案
palette() #查看当前配色方案
palette("Set 1") #选择一个配色方案
# palette("default") # 重置为默认配色方案”R3“
# 查看brewer配色系统
display.brewer.all()

boxplot(log2(GSE54536_non_norm$E),
        ylab = expression(log[2](intensity)),
        las = 2,
        col = rep(c("grey", "white"), each = 5), 
        outline = FALSE)

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

PCA_new(log2(GSE54536_non_norm$E), 
        ntop = nrow(GSE54536_non_norm$E),
        group = GSE54536_targets$group,
        show_name = T)


GSE54536_non_norm$E
## 2.3 data normalization ----
GSE54536_neqc <- neqc(GSE54536_non_norm, #neqc专门为illumina设计，又做背景矫正，又做标准化
                      detection.p="Detection")



GSE54536_expr_neqc <- GSE54536_neqc$E

sum(is.na(GSE54536_expr_neqc))

## 2.4 Plots after normalization ----
## boxplot after normalization
boxplot(GSE54536_expr_neqc,
        ylab = expression(log[2](intensity)),
        las = 2,
        col = rep(c("grey", "white"), each = 5),
        outline = FALSE)

## PCA after normalization
PCA_new(GSE54536_expr_neqc, 
        ntop = nrow(GSE54536_expr_neqc),
        group = GSE54536_targets$group,
        show_name = T)

## 2.5 QC using arrayQualityMrtics ----
library(arrayQualityMetrics)

rownames(GSE54536_targets) <- colnames(GSE54536_neqc$E) <- GSE54536_targets$sample_id

eset <- ExpressionSet(assayData = GSE54536_neqc$E, 
                      phenoData = AnnotatedDataFrame(data = GSE54536_targets))

arrayQualityMetrics(eset, 
                    outdir = "PD/GSE54536/GSE54536_QC", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()
##  GSE54536_sm
boxplot(exprs(GSE54536_sm),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
PCA_new(exprs(GSE54536_sm), 
        ntop = nrow(exprs(GSE54536_sm)),
        group = GSE54536_targets$group,
        show_name = T)


## 2.6 剔除离群样本 ----

GSE54536_targets <- filter(GSE54536_targets,sample_id!='GSM1318547')
GSE54536_non_norm$E <- GSE54536_non_norm$E[,rownames(GSE54536_targets)]
GSE54536_non_norm$other$Detection <- GSE54536_non_norm$other$Detection[,-1] #这里也要对应的剔除一列Pvalue，不然neqc会报错
# 剔除后再评估下
GSE54536_neqc2 <- neqc(GSE54536_non_norm, #neqc专门为illumina设计，又做背景矫正，又做标准化
                      detection.p="Detection")


GSE54536_expr_neqc2 <- GSE54536_neqc2$E


eset2 <- ExpressionSet(assayData = GSE54536_expr_neqc2, 
                      phenoData = AnnotatedDataFrame(data = GSE54536_targets))

arrayQualityMetrics(eset2, 
                    outdir = "PD/GSE54536/GSE54536_outliers_QC", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

boxplot(GSE54536_expr_neqc2,
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE,
        col=as.factor(group))
PCA_new(GSE54536_expr_neqc2, 
        ntop = nrow(GSE54536_expr_neqc2),
        group = GSE54536_targets$group,
        show_name = T)

####  save data
ls(pattern = "^GSE54536")
save(list = ls(pattern = "^GSE54536")[-3], 
     file = "PD/GSE54536/GSE54536_neqc_processed.Rda")



# 3 差异分析 ----

## 1 加载数据 ----
source("Custom_Functions.R")
load_input <- load("PD/GSE54536/GSE54536_non-normalized.Rdata")
load_input

## 2 整理数据 ----
## targets
targets <- GSE54536_targets

## 基因过滤
dim(GSE54536_expr_neqc2)

y <- GSE54536_neqc2
dim(y)                     #Detection就类似于P值，值越小，越可靠
expressed <- rowSums(y$other$Detection < 0.05) >= 5#筛选出y$other$Detection < 0.05，返回一个逻辑值，true就是1，且行之和大于等于5的
y <- y[expressed, ]#筛选出这些行
expr <- y$E#提取表达矩阵
dim(expr)

## 数据评估
boxplot(expr, las = 2, outline = F, col = as.factor(targets$group))

PCA_new(expr, show_name = T)
#可以看到左右两边分得很开，细胞实验批次效应相对来讲比较少，所以这是本身生物学效应的反应

## 3 芯片注释 ----
## 注释信息 GPL10558
library(illuminaHumanv4.db)
illuminaHumanv4()
p2s <- toTable(illuminaHumanv4SYMBOLREANNOTATED) %>% 
  set_names("probe_id", "symbol")#重新改个名字

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
expr_anno <- annotate_expr(expr, p2s)
head(expr_anno)[, 1:3]

head(targets)

# check表达矩阵的列名和分组信息中的样本名是否顺序一致
identical(colnames(expr_anno), targets$sample_id)#对，必须要保持一致

## 4 limma包差异分析 ----
 
library(limma)
#这里按照分组排序，让对照组在前面，实验组在后面
targets <- arrange(targets,group)
expr_anno <- expr_anno[,rownames(targets)]
identical(colnames(expr_anno), targets$sample_id)

## design & contrasts
design <- model.matrix(~ group, data = targets)
colnames(design)
colnames(design) <- c("Con", "PDvsCon")
head(design)

design2 <- model.matrix(~ 0 + group, data = targets) 
colnames(design2)
colnames(design2) <- c("Con", "PD")
head(design2)
contrasts <- makeContrasts("PD-Con", levels = design2)
contrasts

## arrayWeights
aw <- arrayWeights(expr_anno, design = design)
aw
barplot(aw)


## for design 
fit <- lmFit(expr_anno, design, weights = aw)#线性拟合
fit <- eBayes(fit)#贝叶斯
res <- topTable(fit, coef = "PDvsCon")#取最值，PDvsCon代表提取要对比的类型
res

results <- summary(decideTests(fit, lfc = 1))
results
## for design2  本质和for design 一模一样，只是用了管道符
res2 <- lmFit(expr_anno, design2) %>% 
  contrasts.fit(contrasts) %>% 
  eBayes() %>% 
  topTable(coef = "PD-Con")
res2


identical(res, res2)#小数点不完全一致，所以是false

library(dplyr)
sqrt(2) ^ 2 == 2 #2的平方根再平方，然而不等于2
near(sqrt(2) ^ 2, 2)
near(res, res2)#两个非常是接近了就认为相等
all(near(res, res2))

## 所有基因的差异分析结果
GSE54536_DEG <- topTable(fit, coef = "PDvsCon", number = Inf) %>% rownames_to_column('gene_name')
GSE54536_DEG_sig <- topTable(fit, coef = "PDvsCon", n = Inf, p.value = 0.05, lfc = 1)%>% rownames_to_column('gene_name')
#这里的p为矫正后的p

write.table(GSE54536_DEG,
            file = "PD/GSE54536/results/GSE54536_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(GSE54536_DEG_sig,
            file = "PD/GSE54536/results/GSE54536_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)


## topTable & topTableF
topTable(fit)#两素比较用的多
topTableF(fit)#多组比较用这个
all(near(topTable(fit), topTableF(fit)))

## 5 保存结果 ----
save(expr_anno, targets, GSE54536_res_DEG_sig, file = "PD/GSE54536/results/GSE54536_res_DEG.Rda")


length(intersect(res_deseq2$gene_name, GSE54536_DEG_sig$gene_name))
library(VennDiagram)
venn.diagram(x = list("Covid-19" = res_deseq2$gene_name, 
                      "PD" = GSE54536_DEG_sig$gene_name),
             filename = "PD/GSE54536/results/Venn.jpeg",
             fill = c("green", "red"),
             scaled = F,
             cex = 1.5,
             cat.cex = 1.5,
             main = "Venn Diagram",
             main.cex = 2)
