############################################################################===
#####  GSE37667原始数据清洗,Paired
############################################################################===

# 1 GSE37667_targets_unpaired ----
library(GEOquery)
library(tidyverse)

GSE37667_sm <- getGEO(filename = "Sleep Disorder/GSE37667/GSE37667_series_matrix.txt.gz",
                      getGPL = F) 
GSE37667_pd <- pData(GSE37667_sm)

sample_name <- GSE37667_pd[,1] %>% str_replace_all(c("Prolonged_wakefulness_" = "sleepless", "baseline_" = "ctrl",  "Sleep recovery_" = "excess","volunteer" = "_")) 

GSE37667_targets_paired <- GSE37667_pd %>% 
  dplyr::select(2, 1) %>% #取两列
  mutate(sample_id = geo_accession,#添加几列
         group = str_match(sample_name, "(.*)_\\d+")[,2],
         patient_id = as.numeric(str_extract(title, "\\d+")),
         sample_name = sample_name) %>% #把两个合并
  dplyr::select(sample_id:sample_name) %>%#选出这些列
  dplyr::filter(group != "excess") %>%
  add_count(patient_id, name = "IsPaired") %>% #对patient_id进行计数，就是每个id有几个，并且新列IsPaired
  dplyr::filter(IsPaired == 2) %>% #选出为2的
  arrange(patient_id)#按照patient_id重新排列


rownames(GSE37667_targets_paired) <- GSE37667_targets_paired[,1]

# 这里对GSE37667_targets_unpaired文件进行排序，按照sample_id的顺序升序排列，用于原始数据读取匹配文件名
library(dplyr)
GSE37667_targets_paired <- arrange(GSE37667_targets_paired,sample_id)


# 2 使用affy包读取CEL文件 ----

### 读入CEL
library(affy)

FileName <- list.files(path = "Sleep Disorder/GSE37667/GSE37667_RAW_Paired/")

FileName

## affy包中有2个读取CEL文件的函数，推荐ReadAffy(),更灵活易用
# cel <- read.affybatch(filenames = cel_files)
GSE37667_cel_paired <- ReadAffy(filenames = FileName,
                         celfile.path = "Sleep Disorder/GSE37667/GSE37667_RAW_Paired/",
                         phenoData = GSE37667_targets_paired) #注意phenoData的行名必须要确保和读进来的aff原始文件一致。



cel <- GSE37667_cel_paired 
# 表达矩阵 exprs
head(exprs(cel))[, 1:5]
# 样本信息 pData
GSE37667_targets_paired<- pData(cel)
# 样本名称 sampeNames     行名赋给了表达矩阵列名，即样本名
sampleNames(cel)
# ScanDate
cel@protocolData@data$ScanDate
# expression matrix
GSE37667_expr_paired <- exprs(cel)##

# 3 Affymetrix芯片数据标准化 ----

palette.pals() #查看配色方案
palette() #查看当前配色方案
palette("Set 1") #选择一个配色方案
# palette("default") # 重置为默认配色方案”R3“

# 芯片扫描日期可作为一种重要的批次效应因素
GSE37667_targets_paired$batch <- GSE37667_cel_paired@protocolData@data$ScanDate %>% 
  str_sub(1, 8)#提取1到8个字符

table(GSE37667_targets_paired$batch)

### *~RNA degradation  ---- 
#评估RNA的讲解程度
RNAdeg <- AffyRNAdeg(GSE37667_cel_paired)
summaryAffyRNAdeg(RNAdeg)
#RNA降解从5‘端开始降解，所5'信号值要比3'信号值低一些，所以假如低得多了，斜率就大了，说明降解得很厉害

#  RNA degradation plot
cols <- rainbow(nrow(GSE37667_targets_paired))
plotAffyRNAdeg(RNAdeg, cols = cols)
legend("topleft",
       ncol = 1,
       legend = sampleNames(GSE37667_cel_paired), 
       lty = 2, 
       lwd = 1,
       cex = 0.5,
       box.lty=0,
       bg = "transparent",
       col = cols)
box()#这里画图出了问题，legend的放置位置错了

### *~simpleaffy包的质量评估函数qc ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("simpleaffy")

library(simpleaffy)
GSE37667_qc <- simpleaffy::qc(GSE37667_cel_paired)
plot(GSE37667_qc)
# P/M/A:present/ marginal present/ absent
# actin 3/5 - 3, gapdh 3/5 - 1


### *~boxplot and histogram ----
#箱线图
boxplot(GSE37667_cel_paired, 
        las = 2, 
        outline = FALSE,
        col = as.factor(GSE37667_targets_paired$group))
#密度图
hist(GSE37667_cel_paired, lty = 1:3, col = cols)
legend("topright",
       legend = sampleNames(GSE37667_cel_paired),
       ncol = 2,
       lty = 1:3,
       cex = 0.8,
       col = cols,
       box.col = "transparent",
       xpd = TRUE)
box()

### *~PCA_new ----
## new function for PCA 
#参考了Deseq2包来写的PCA函数
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

PCA_new(log2(GSE37667_expr_paired), 
        nrow(GSE37667_expr_paired), 
        group = GSE37667_targets_paired$group, #筛选了18列的数据，用的是对应的targets文件
        show_name = T)

## ~芯片数据质量评估 基于arrayQulitMetrics ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("arrayQualityMetrics")

library(arrayQualityMetrics)
arrayQualityMetrics(GSE37667_cel_paired, 
                    outdir = "Sleep Disorder/GSE37667/GSE37667_cel_paired_QC", #通过该文件夹下的index来看结果
                    force = TRUE,#如果本身有个文件夹的话就把文件夹刷新了
                    intgroup = "group",#GSE37667_cel_paired@phenoData@data[["group"]]
                    do.logtransform = TRUE)
dev.off()#dev.off报错了也无妨，就是为了让它报错

## ~背景校正 ----
GSE37667_bgc <- affy::bg.correct(GSE37667_cel_paired, method = "rma")#最常用rma方法
bgcorrect.methods()

#背景矫正后作图
boxplot(GSE37667_bgc, 
        las = 2, 
        outline = FALSE,
        col = as.factor(GSE37667_targets_paired$group))

PCA_new(log2(exprs(GSE37667_bgc)), 
        nrow(GSE37667_bgc), 
        group = GSE37667_targets_paired$group,
        show_name = T)
#可以看到PCA图还是没什么改善

## ~芯片数据标准化 ----
###对于affy芯片来说，rma或者gcrma已经把背景矫正整合进去了，所以没有必要单独做背景矫正

### *~RMA ----

GSE37667_rma <- affy::rma(GSE37667_cel_paired)#用的是原始cel数据，不是背景矫正后的数据

boxplot(GSE37667_rma, 
        las = 2, 
        outline = FALSE,
        col = as.factor(GSE37667_targets_paired$group))
#看上去中位线就在一条线上了

PCA_new(exprs(GSE37667_rma), 
        nrow(GSE37667_rma), 
        group = GSE37667_targets_paired$group,
        show_name = T)

#芯片数据标准化之后看数据质量是否改善
library(arrayQualityMetrics)
arrayQualityMetrics(GSE37667_rma, 
                    outdir = "Sleep Disorder/GSE37667/GSE37667_cel_paired_rma_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)

### *~GCRMA ----
library(gcrma)
GSE37667_gcrma <- gcrma(GSE37667_cel_paired)
GSE37667_expr_paired_gcrma <- exprs(GSE37667_gcrma)
boxplot(GSE37667_expr_paired_gcrma, 
        las = 2, 
        outline = FALSE,
        col = as.factor(GSE37667_targets_paired$group))

dev.off()
PCA_new(exprs(GSE37667_gcrma), 
        nrow(GSE37667_gcrma), 
        group = GSE37667_targets_paired$group,
        show_name = F)

library(arrayQualityMetrics)
arrayQualityMetrics(GSE37667_gcrma, 
                    outdir = "Sleep Disorder/GSE37667/GSE37667_cel_paired_gcrma_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)
dev.off()




### *~缺失值补充 ----
#首先看一下有没有缺失值
sum(is.na(exprs(GSE37667_gcrma)))



### *~save data ----

save(GSE37667_expr_paired_gcrma, 
     file = "Sleep Disorder/GSE37667/GSE37667_expr_paired_gcrma.Rdata")

# 4 去除批次效应 ----
### *~加载数据 ----
load_input <- load("affymetrix/GSE37667/GSE37667_expr_paired_gcrma.Rda")
GSE37667_expr_paired_gcrma <- exprs(GSE37667_gcrma)

### *~批次效应去除 ----
# BiocManager::install("sva")
# devtools::install_github('zhangyuqing/sva-devel')
library(sva)
library(tidyverse)

table(GSE37667_targets_paired$batch)#看批次，分了5批

sum(is.na(GSE37667_expr_paired_gcrma))#看是否有缺失值

mod <-  model.matrix(~1 + group, data= GSE37667_targets_paired)#把分组group因素告诉mod
GSE37667_combat <- ComBat(dat = GSE37667_expr_paired_gcrma, 
                          batch = GSE37667_targets_paired$batch, 
                          mod = mod, 
                          par.prior=TRUE, 
                          prior.plots=FALSE)

boxplot(GSE37667_combat, las = 2, 
        col = as.factor(GSE37667_targets_paired$group), outline = F)

ntop <- nrow(GSE37667_combat)
p1 <- PCA_new(GSE37667_expr_paired_gcrma, 
              ntop = ntop,
              group = GSE37667_targets_paired$group,
              show_name = T) +
  ggtitle("PCA before batch") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- PCA_new(GSE37667_combat, 
              ntop = ntop,
              group = GSE37667_targets_paired$group,
              show_name = T) +
  ggtitle("PCA after batch") +
  theme(plot.title = element_text(hjust = 0.5))

p12 <- cowplot::plot_grid(p1, p2, nrow = 1)#组图函数
# pdf("batch/PCA_GSE37667_batch.pdf", height = 7, width = 15)
p12

dev.off()
### *~基于arrayQulitMetrics----

GSE37667_gcrma_batch_QC <- ExpressionSet(assayData = as.matrix(GSE37667_combat), 
                                            phenoData = AnnotatedDataFrame(data = GSE37667_targets_paired))

arrayQualityMetrics(GSE37667_gcrma_batch_QC, 
                    outdir = "Sleep Disorder/GSE37667/GSE37667_gcrma_batch_QC", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = F)
dev.off()


# sessionInfo()


# 5 limma包差异分析----
library(limma)
### *~芯片注释 ----
annotate_expr <- function(expr, p2s, fun = mean){
  # expr
  expr <- expr %>% 
    as.data.frame() %>% 
    rownames_to_column("probe_id")
  # annotate
  expr_annotated <- p2s %>% 
    inner_join(expr, by = "probe_id") %>% 
    na.omit() %>% 
    dplyr::select(-probe_id) %>% 
    aggregate(. ~ symbol, data =., FUN = fun) %>% 
    column_to_rownames("symbol") 
  return(expr_annotated)
}

library(hgu133plus2.db)
hgu133plus2()
p2s <- toTable(hgu133plus2SYMBOL) 
  
GSE37667_expr_anno <- annotate_expr(expr = GSE37667_combat,
                                    p2s = p2s,
                                    fun = mean)
# check表达矩阵的列名和分组信息中的样本名是否顺序一致
identical(colnames(GSE37667_expr_anno), GSE37667_targets_paired$sample_id)#对，必须要保持一致

### *~DEG ----
## 输入数据 expression data
expr <- GSE37667_expr_anno

## 输入数据 design
targets <- GSE37667_targets_paired

design <- model.matrix(~ group + patient_id, data = targets)#最大的区别就是patient_id
colnames(design)
colnames(design) <- c("Ctrl", "Sleepless vs Ctrl", "patient_id")
head(design)

## arrayWeights计算权重
aw <- arrayWeights(expr, design = design)
aw
barplot(aw)

## 差异分析
fit <- lmFit(expr, design, weights = aw)#添加一个权重信息，如果是人体样本处理差别大一些，细胞样本技术上差异小一些，所以是人体样本的话就计算weight
fit <- eBayes(fit)
#因为只有两组，所以不需要Contrasts

summary(decideTests(fit, lfc = 1))#提取结果
#            N   TvsN   patient_id T和N比较，有558个基因是下调的
# Down       0   558          0
# NotSig     0 12785      13680
# Up     13680   337          0

GSE37667_DEG <- topTable(fit, coef = "Sleepless vs Ctrl", n = Inf) %>% rownames_to_column('gene_name')

GSE37667_DEG_sig <- topTable(fit, coef = "Sleepless vs Ctrl", n = Inf, p = 0.05, lfc = 1) %>% rownames_to_column('gene_name') #lfc为logfc=1

nrow(GSE37667_DEG_sig)


write.table(GSE37667_DEG,
            file = "Sleep Disorder/GSE37667/results/GSE37667_DEG.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(GSE37667_DEG_sig,
            file = "Sleep Disorder/GSE37667/results/GSE37667_DEG_sig.txt",
            quote = F,
            sep = "\t",
            row.names = F)

length(intersect(res_deseq2$gene_name, GSE37667_DEG_sig$gene_name))
library(VennDiagram)
venn.diagram(x = list("Covid-19" = res_deseq2$gene_name, 
                      "Sleep Disorder" = GSE37667_DEG_sig$gene_name),
             filename = "Sleep Disorder/GSE37667/results/Venn.jpeg",
             fill = c("blue", "red"),
             scaled = F,
             cex = 1.5,
             cat.cex = 1.5,
             main = "Venn Diagram",
             main.cex = 2)



#p值越小，B值越大

15/(15+1)


### *~保存数据 ----
save(expr, targets, GSE37667_paired_res, GSE37667_paired_res_sig, file = "Sleep Disorder/GSE37667/results/GSE37667_paired_DEG.Rda")

#### End ----
