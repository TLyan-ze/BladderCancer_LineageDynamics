#设置工作目录 

setwd("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression")
getwd()

library(Seurat)

#构建对象
scrna_data1 <- Read10X_h5(
  "/WorkDir3/DewangZhou/scRNA-M/RC-3/outs/filtered_feature_bc_matrix.h5")
seob1 <- CreateSeuratObject(counts = scrna_data5, 
                            min.cells = 3, 
                            min.features = 200)

type1 ="RT"
seob1@meta.data$type = type1
sample1 = "CTR_RT"
seob1@meta.data$sample = sample1
seob1 <- RenameCells(seob1, new.names = paste0("CTR_RT", colnames(x = seob1)))

scrna_data2 <- Read10X_h5(
  "/WorkDir3/DewangZhou//scRNA-M/RTG/outs/filtered_feature_bc_matrix.h5")
seob2 <- CreateSeuratObject(counts = scrna_data6, 
                            min.cells = 3, 
                            min.features = 200)


type1 ="RT"
seob2@meta.data$type = type1
sample2 = "GEM_RT"
seob2@meta.data$sample = sample2
seob2 <- RenameCells(seob2, new.names = paste0("GEM_RT", colnames(x = seob2)))



##########################################################################

seob_list = list(CTR_RT = seob1, GEM_RT = seob2) # MGH 合并


seob <- merge(x = seob_list[[1]], 
              y = seob_list[-1])


#转换
library(SeuratDisk)
library(Seurat)      
#seurat2h5seurat中间过渡 
SaveH5Seurat(seob,filename="RT.h5seurat", overwrite = TRUE)
#数据转为最终h5ad格式
Convert("RT.h5seurat", dest = "h5ad", overwrite = TRUE) 
