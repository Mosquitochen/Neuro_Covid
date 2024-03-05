# 合并covid肺组织和免疫来源的数据集
# 再分别和神经常见病取交集进行后续富集分析 
library(tidyverse)
DEGs_sig <- read.delim("Result/GO&Kegg/DEGs_sig.txt") %>% 
  qpcR:::cbind.na(covid= union(.[,1],.[,8])) %>% select(-1,-8) %>% 
    select(7,1:6) %>% .[1:3573,]


Covid_SD <- intersect(DEGs_sig[,1],DEGs_sig[,2])
Covid_HS <- intersect(DEGs_sig[,1],DEGs_sig[,3])
Covid_AD <- intersect(DEGs_sig[,1],DEGs_sig[,4])
Covid_EP <- intersect(DEGs_sig[,1],DEGs_sig[,5])
Covid_PD <- intersect(DEGs_sig[,1],DEGs_sig[,6])
Covid_IS <- intersect(DEGs_sig[,1],DEGs_sig[,7])

Neuro_Covid <- data.frame(Covid_SD) %>% 
  qpcR:::cbind.na(Covid_HS,Covid_AD,Covid_EP,Covid_PD,Covid_IS) 

write.table(Neuro_Covid,"Result/GO&Kegg/Neuro_Covid.txt",sep="\t",row.names = F,na="")

# 获取包含logFC值得差异基因表用于联合logFC的富集分析
filename <- list.files('Result/GO&Kegg/Neu_DEG',full.names = T)
filename[5]

Neuro_Covid <- read_tsv("Result/GO&Kegg/Neuro_Covid.txt")

AD_enrich <-  read.table(filename[1],header=T)[,1:2] %>% 
              set_names('Covid_AD','logFC') %>% semi_join(.,Neuro_Covid['Covid_AD'],'Covid_AD') %>% 
              set_names('id','logFC') 
write.table(AD_enrich,"Result/GO&Kegg/logFC-Enrichment/AD_enrich.txt",sep="\t",row.names = F,na="")

EP_enrich <-  read.table(filename[2],header=T)[,1:2] %>% 
  set_names('Covid_EP','logFC') %>% semi_join(.,Neuro_Covid['Covid_EP'],'Covid_EP') %>% 
  set_names('id','logFC') 
write.table(EP_enrich,"Result/GO&Kegg/logFC-Enrichment/EP_enrich.txt",sep="\t",row.names = F,na="")

HS_enrich <-  read.table(filename[3],header=T)[,1:2] %>% 
  set_names('Covid_HS','logFC') %>% semi_join(.,Neuro_Covid['Covid_HS'],'Covid_HS') %>% 
  set_names('id','logFC') 
write.table(HS_enrich,"Result/GO&Kegg/logFC-Enrichment/HS_enrich.txt",sep="\t",row.names = F,na="")

IS_enrich <-  read.table(filename[4],header=T)[,1:2] %>% 
  set_names('Covid_IS','logFC') %>% semi_join(.,Neuro_Covid['Covid_IS'],'Covid_IS') %>% 
  set_names('id','logFC') 
write.table(IS_enrich,"Result/GO&Kegg/logFC-Enrichment/IS_enrich.txt",sep="\t",row.names = F,na="")

PD_enrich <-  read.table(filename[5],header=T)[,1:2] %>% 
  set_names('Covid_PD','logFC') %>% semi_join(.,Neuro_Covid['Covid_PD'],'Covid_PD') %>% 
  set_names('id','logFC') 
write.table(PD_enrich,"Result/GO&Kegg/logFC-Enrichment/PD_enrich.txt",sep="\t",row.names = F,na="")

SD_enrich <-  read.table(filename[6],header=T)[,1:2] %>% 
  set_names('Covid_SD','logFC') %>% semi_join(.,Neuro_Covid['Covid_SD'],'Covid_SD') %>% 
  set_names('id','logFC') 
write.table(SD_enrich,"Result/GO&Kegg/logFC-Enrichment/SD_enrich.txt",sep="\t",row.names = F,na="")

