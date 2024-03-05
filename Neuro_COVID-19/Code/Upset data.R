
# Covid交集 ----
Filename <- list.files('Result/DEGs',full.names = T)
Filename

DEGs <- read.table(Filename[1],header=T)[1]

for (i in 2:7) {
  DEGs_i <- read.table(Filename[i],header=T)[1] 
  DEGs <- qpcR:::cbind.na(DEGs,DEGs_i)
}

colnames(DEGs) <- c('AD','Covid-19','EP','HS','IS','PD','SD')


# 按差异基因数量重新排列
#两种顺序
c('Covid-19','SD','HS','AD','EP','PD','IS')
c('IS','PD','EP','AD','HS','SD','Covid-19')
DEGs <- DEGs[,c('IS','PD','EP','AD','HS','SD','Covid-19')]

write.table(DEGs,"Result/Upset/Upset_data.txt",sep="\t",row.names = F,na="")

# Covid 免疫基因交集 ----
Filename <- list.files('Result/DEGs_Imm',full.names = T)
Filename

DEGs <- read.table(Filename[1],header=T)[1]

for (i in 2:7) {
  DEGs_i <- read.table(Filename[i],header=T)[1] 
  DEGs <- qpcR:::cbind.na(DEGs,DEGs_i)
}

colnames(DEGs) <- c('AD','EP','HS','COV_Imm','IS','PD','SD')


# 按差异基因数量重新排列
#两种顺序
c('SD','HS','AD','EP','PD','COV_Imm','IS')
c('IS','COV_Imm','PD','EP','AD','HS','SD')
DEGs <- DEGs[,c('IS','COV_Imm','PD','EP','AD','HS','SD')]

write.table(DEGs,"Result/Upset/Upset_Imm_data.txt",sep="\t",row.names = F,na="")
