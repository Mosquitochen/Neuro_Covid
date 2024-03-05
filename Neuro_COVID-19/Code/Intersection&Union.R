
# 汇总显著性差异基因在一张表里 ----
Filename <- list.files('Result/Heatmap/DEGs_sig/',full.names = T)
Filename 

DEGs_sig <- read.table(Filename[1],header=T)[1]

for (i in 2:8) {
  DEGs_i <- read.table(Filename[i],header=T)[1] 
  DEGs_sig <- qpcR:::cbind.na(DEGs_sig,DEGs_i)
}

colnames(DEGs_sig) <- c('AD','Covid-19','EP','HS','Imm','IS','PD','SD')


c('Covid-19','SD','HS','AD','EP','PD','IS','Imm')

DEGs_sig <- DEGs_sig[,c('Covid-19','SD','HS','AD','EP','PD','IS','Imm')]

write.table(DEGs_sig,"Result/Heatmap/DEGs_sig.txt",sep="\t",row.names = F,na="")


# 分别计算交集基因 ，以上代码可以不用再看，直接从这里开始执行----
DEGs_sig <- read.delim("Result/Heatmap/DEGs_sig.txt")
library(tidyverse)

Covid_SD <- intersect(DEGs_sig[,1],DEGs_sig[,2])
Covid_HS <- intersect(DEGs_sig[,1],DEGs_sig[,3])
Covid_AD <- intersect(DEGs_sig[,1],DEGs_sig[,4])
Covid_EP <- intersect(DEGs_sig[,1],DEGs_sig[,5])
Covid_PD <- intersect(DEGs_sig[,1],DEGs_sig[,6])
Covid_IS <- intersect(DEGs_sig[,1],DEGs_sig[,7])

intersect_Covid <- data.frame(Covid_SD) %>% 
              qpcR:::cbind.na(Covid_HS,Covid_AD,Covid_EP,Covid_PD,Covid_IS)

write.table(intersect_Covid,"Result/Heatmap/Intersection&Union/intersect_Covid.txt",sep="\t",row.names = F,na="")

                            
Imm_SD <- head(intersect(DEGs_sig[,8],DEGs_sig[,2]),-1)
Imm_HS <- head(intersect(DEGs_sig[,8],DEGs_sig[,3]),-1)
Imm_AD <- head(intersect(DEGs_sig[,8],DEGs_sig[,4]),-1)
Imm_EP <- head(intersect(DEGs_sig[,8],DEGs_sig[,5]),-1)
Imm_PD <- head(intersect(DEGs_sig[,8],DEGs_sig[,6]),-1)
Imm_IS <- head(intersect(DEGs_sig[,8],DEGs_sig[,7]),-1)

intersect_Imm <- data.frame(Imm_SD) %>% 
  qpcR:::cbind.na(Imm_HS,Imm_AD,Imm_EP,Imm_PD,Imm_IS)

write.table(intersect_Imm,"Result/Heatmap/Intersection&Union/intersect_Imm.txt",sep="\t",row.names = F,na="")

# 热图 ----
# 热图输入数据准备
DEGs_file <- list.files('Result/Heatmap/DEGs/',full.names = T)
DEGs_file 

# 计算所有交集基因的并集
union_Covid <- union(union(union(union(union(Covid_SD,Covid_HS),Covid_AD),Covid_EP),Covid_PD),Covid_IS) %>% 
  data.frame(gene_name=.,row.names = .)

union_Imm <- union(union(union(union(union(Imm_SD,Imm_HS),Imm_AD),Imm_EP),Imm_PD),Imm_IS)%>% 
  data.frame(gene_name=.,row.names = .)

# Covid-19 热图数据 ----
# logFC热图数据
Heatmap_lfc_cov <- union_Covid 

for (i in c(1:4,6:8)) {
  DEGs_name <- str_match(DEGs_file[i],'//(.*)_')[2]
  Heatmap_lfc_cov <- read_tsv(DEGs_file[i]) %>% select(gene_name,logFC) %>% set_names('gene_name',DEGs_name) %>% 
  left_join(Heatmap_lfc_cov,.,by='gene_name')
}

colnames(Heatmap_lfc_cov) <- c('gene_name','AD','Covid-19','EP','HS','IS','PD','SD')


Heatmap_lfc_cov <- Heatmap_lfc_cov[,c('gene_name','Covid-19','SD','AD','EP','IS','PD','HS')]


# 矫正后P值热图数据
Heatmap_padj_cov <- union_Covid 

for (i in 1:4) {
  DEGs_name <- str_match(DEGs_file[i],'//(.*)_')[2]
  Heatmap_padj_cov <- read_tsv(DEGs_file[i]) %>% select(gene_name,padj) %>% set_names('gene_name',DEGs_name) %>% 
    left_join(Heatmap_padj_cov,.,by='gene_name')
}

for (i in 6:8) {
  DEGs_name <- str_match(DEGs_file[i],'//(.*)_')[2]
  Heatmap_padj_cov <- read_tsv(DEGs_file[i]) %>% select(gene_name,adj.P.Val) %>% set_names('gene_name',DEGs_name) %>% 
    left_join(Heatmap_padj_cov,.,by='gene_name')
}
  
colnames(Heatmap_padj_cov) <- c('gene_name','AD','Covid-19','EP','HS','IS','PD','SD')


Heatmap_padj_cov <- Heatmap_padj_cov[,c('gene_name','Covid-19','SD','AD','EP','IS','PD','HS')]

write.table(Heatmap_lfc_cov,"Result/Heatmap/Intersection&Union/Heatmap_lfc_cov.txt",sep="\t",row.names = F,na="")
write.table(Heatmap_padj_cov,"Result/Heatmap/Intersection&Union/Heatmap_padj_cov.txt",sep="\t",row.names = F,na="")

# Cov-Imm 热图数据 ----
# logFC热图数据
Heatmap_lfc_Imm <- union_Imm 

for (i in c(1,3:8)) {
  DEGs_name <- str_match(DEGs_file[i],'//(.*)_')[2]
  Heatmap_lfc_Imm <- read_tsv(DEGs_file[i]) %>% select(gene_name,logFC) %>% set_names('gene_name',DEGs_name) %>% 
    left_join(Heatmap_lfc_Imm,.,by='gene_name')
}

colnames(Heatmap_lfc_Imm) <- c('gene_name','AD','EP','HS','Cov-Imm','IS','PD','SD')


Heatmap_lfc_Imm <- Heatmap_lfc_Imm[,c('gene_name','Cov-Imm','IS','SD','AD','PD','EP','HS')]


# 矫正后P值热图数据
Heatmap_padj_Imm <- union_Imm 

for (i in c(1,3:5)){
  DEGs_name <- str_match(DEGs_file[i],'//(.*)_')[2]
  Heatmap_padj_Imm <- read_tsv(DEGs_file[i]) %>% select(gene_name,padj) %>% set_names('gene_name',DEGs_name) %>% 
    left_join(Heatmap_padj_Imm,.,by='gene_name')
}

for (i in 6:8) {
  DEGs_name <- str_match(DEGs_file[i],'//(.*)_')[2]
  Heatmap_padj_Imm <- read_tsv(DEGs_file[i]) %>% select(gene_name,adj.P.Val) %>% set_names('gene_name',DEGs_name) %>% 
    left_join(Heatmap_padj_Imm,.,by='gene_name')
}

colnames(Heatmap_padj_Imm) <- c('gene_name','AD','EP','HS','Cov-Imm','IS','PD','SD')

Heatmap_padj_Imm <- Heatmap_padj_Imm[,c('gene_name','Cov-Imm','IS','SD','AD','PD','EP','HS')]


write.table(Heatmap_lfc_Imm,"Result/Heatmap/Intersection&Union/Heatmap_lfc_Imm.txt",sep="\t",row.names = F,na="")
write.table(Heatmap_padj_Imm,"Result/Heatmap/Intersection&Union/Heatmap_padj_Imm.txt",sep="\t",row.names = F,na="")


