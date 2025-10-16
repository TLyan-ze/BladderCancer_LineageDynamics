library(SeuratDisk)
library(Seurat)

setwd("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression")
getwd()
#LoadLoom路径不能有中文，不然报错
##直接h5ad转换成h5seurat 
##overwrite参数：覆盖源文件
#Convert("scRNA4.h5ad","h5seurat",assay = "RNA", overwrite = T)

library(anndata)
library(Seurat)
a <- read_h5ad("scRNA4_2000_all.h5ad")
b <- a$X

# 读取H5AD格式的单细胞数据
library(Seurat)
library(hdf5r)
library(SeuratDisk)
library(patchwork)
ct <- as.matrix(b)  #将稀疏矩阵转回普通文本矩阵
ct <- round(ct*100, digits = 0)
ct <- t(ct)
#genes.tsv
write.table(data.frame(rownames(ct),rownames(ct)),file = 'plot_all/genes.tsv',
            quote = F,sep = '\t',
            col.names = F,row.names = F)
#barcodes.tsv
write.table(colnames(ct),file = 'plot_all/barcodes.tsv',quote = F,
            col.names = F,row.names = F)


#matrix.mtx 文件是3列，第一列是行号，第二列是列号，第三列是基因表达量
#首先写一个头信息
file="plot_all/matrix.mtx"
sink(file)
cat("%%MatrixMarket matrix coordinate integer general\n")
cat("%\n")
cat(paste(nrow(ct),ncol(ct),sum(ct>0),"\n")) 
sink()

#再写入表达量信息
tmp=do.call(rbind,lapply(1:ncol(ct),function(i){
  return(data.frame(row=1:nrow(ct),
                    col=i,
                    exp=ct[,i]))
}) )
tmp=tmp[tmp$exp>0,]
head(tmp)
write.table(tmp,file = 'plot_all/matrix.mtx',quote = F,
            col.names = F,row.names = F,append = T )
#运行时间根据数据大小定。

##################################################################################################################################

