rm(list=ls())
options(stringsAsFactors = FALSE)

library(tximport)  # Import transcript-level abundances and counts
library(tidyverse) # ggplot2, dplyr, tidyr, readr, purrr, tibble, forcats
library(data.table) # 高效读写

# 使用第一个命令行参数作为包含 *_quant.sf 的目录；若未提供，则使用当前目录
args <- commandArgs(trailingOnly = TRUE)
dir_name <- if (length(args) >= 1) args[1] else getwd()
setwd(dir_name)
print(getwd())

####  salmon 原始文件处理  ####
## 载入 transcript_id 和 symbol 的映射关系文件
# 人类：
t2s <- fread("/share/database/openData/GRCh38_GENCODE/gencode.v35.annotation.TxToGene.csv",
             data.table = FALSE, header = TRUE)
# 小鼠（如需小鼠，请切换下面此行并注释掉上面的人类）：
# t2s <- fread("/datapool/yanzeqin/database/GRCm39_GENCODE/gencode.mus.TxToGene.csv",
#              data.table = FALSE, header = TRUE)

## 找到所有 quant.sf 文件路径并导入处理
files <- list.files(pattern = "_quant.sf", recursive = TRUE, full.names = TRUE)
print(paste("Found", length(files), "quant.sf files"))
if (length(files) == 0) {
  stop("No quant.sf files found under: ", dir_name)
}

# 汇总为基因层面的计数与 TPM
txi <- tximport(files, type = "salmon", tx2gene = t2s)

## 提取样品名作为列名
cn <- sapply(strsplit(files, "_quant.sf"), function(x) x[1])
cn2 <- sapply(strsplit(cn, "/"), function(x) x[length(x)])
colnames(txi$counts) <- gsub("_quant", "", cn2)

## 导出 counts 与 tpm 表达矩阵
counts <- as.data.frame(apply(txi$counts, 2, as.integer)) # 将 counts 数取整
rownames(counts) <- rownames(txi$counts)
write.table(counts, file = "exp_counts.csv", sep = ",",
            col.names = TRUE, row.names = TRUE, quote = FALSE)

tpm <- as.data.frame(txi$abundance)  # abundance 为基因的 TPM 值
colnames(tpm) <- colnames(txi$counts)
write.table(tpm, file = "exp_tmp.csv", sep = ",",
            col.names = TRUE, row.names = TRUE, quote = FALSE)

cat("Saved exp_counts.csv and exp_tmp.csv under:", getwd(), "\n")
