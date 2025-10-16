#设置工作目录 

setwd("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression")
getwd()

library(Seurat)

#构建对象
scrna_data1 <- Read10X_h5(
  "/WorkDir3/DewangZhou/scRNA-M/MC-3/outs/MC_filtered_feature_bc_matrix.h5")
seob1 <- CreateSeuratObject(counts = scrna_data1, 
                            min.cells = 3, 
                            min.features = 200)

type1 ="MGH"
seob1@meta.data$type = type1
sample1 = "CTR_MGH"
seob1@meta.data$sample = sample1
seob1 <- RenameCells(seob1, new.names = paste0("CTR_MGH.", colnames(x = seob1)))

scrna_data2 <- Read10X_h5(
  "/WorkDir3/DewangZhou/scRNA-M/N-M/outs/N-M_filtered_feature_bc_matrix.h5"
)
seob2 <- CreateSeuratObject(counts = scrna_data2, 
                            min.cells = 3, 
                            min.features = 200)


type1 ="MGH"
seob2@meta.data$type = type1
sample2 = "NMI_MGH"
seob2@meta.data$sample = sample2
seob2 <- RenameCells(seob2, new.names = paste0("NMI_MGH.", colnames(x = seob2)))


scrna_data3 <- Read10X_h5(
  "/WorkDir3/DewangZhou/scRNA-M/M-M/outs/M-M_filtered_feature_bc_matrix.h5"
)
seob3 <- CreateSeuratObject(counts = scrna_data3, 
                            min.cells = 3, 
                            min.features = 200)

type1 ="MGH"
seob3@meta.data$type = type1
sample3 = "MI_MGH"
seob3@meta.data$sample = sample3
seob3 <- RenameCells(seob3, new.names = paste0("MI_MGH.", colnames(x = seob3)))

scrna_data4 <- Read10X_h5(
  "/WorkDir3/DewangZhou/scRNA-M/Me-M/outs/Me-M_filtered_feature_bc_matrix.h5"
)
seob4 <- CreateSeuratObject(counts = scrna_data4, 
                            min.cells = 3, 
                            min.features = 200)

type1 ="MGH"
seob4@meta.data$type = type1
sample4 = "Met_MGH"
seob4@meta.data$sample = sample4
seob4 <- RenameCells(seob4, new.names = paste0("Met_MGH.", colnames(x = seob4)))




#seob_list = list(CTR_MGH = seob1, NMI_MGH = seob2, 
#                 MI_MGH = seob3, Met_MGH = seob4
#) 


seob_list = list(CTR_MGH = seob1, NMI_MGH = seob2, 
MI_MGH = seob3, Met_MGH = seob4) # MGH 合并


seob <- merge(x = seob_list[[1]], 
              y = seob_list[-1])


#转换
library(SeuratDisk)
library(Seurat)      
#seurat2h5seurat中间过渡 
SaveH5Seurat(seob,filename="MGH.h5seurat", overwrite = TRUE)
#数据转为最终h5ad格式
Convert("MGH.h5seurat", dest = "h5ad", overwrite = TRUE) 
