rm(list=ls())
options(stringsAsFactors = F)
library(tximport) #Import transcript-level abundances and counts for transcript- and gene-level analysis packages
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table) #多核读取文件


Args <- commandArgs()
dir_name <- Args[6]

setwd(dir_name)
getwd()


####  salmon原始文件处理  ####
##载入transcript_id和symbol的对应关系文件
t2s <- fread("/share/database/openData/GRCh38_GENCODE/gencode.v35.annotation.TxToGene.csv", data.table = F, header = T); head(t2s)
#/datapool/yanzeqin/database/GRCm39_GENCODE/gencode.mus.TxToGene.csv
#小鼠的
#t2s <- fread("/datapool/yanzeqin/database/GRCm39_GENCODE/gencode.mus.TxToGene.csv", data.table = F, header = T); head(t2s)


##找到所有quant.sf文件所在路径  导入salmon文件处理汇总
files <- list.files(pattern="_quant.sf",recursive=T, full.names = T); files  #显示目录下所有符合要求的文件
txi <- tximport(files, type = "salmon", tx2gene = t2s)

##提取文件夹中的样品名作为counts行名
cn <- sapply(strsplit(files,'_quant.sf'), function(x) x[1])
cn2 <- sapply(strsplit(cn,'\\/'), function(x) x[length(x)])
colnames(txi$counts) <- gsub('_quant','',cn2); colnames(txi$counts)

##提取出counts/tpm表达矩阵
counts <- as.data.frame(apply(txi$counts,2,as.integer)) #将counts数取整
rownames(counts) <- rownames(txi$counts) 
dim(counts)
write.table(counts,file="exp_counts.csv",sep=",",col.names=T,row.names=T,quote=FALSE)

tpm <- as.data.frame(txi$abundance)  ###abundance为基因的Tpm值
colnames(tpm) <- colnames(txi$counts)
dim(tpm)

write.table(tpm,file="exp_tmp.csv",sep=",",col.names=T,row.names=T,quote=FALSE)
