

setwd("/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/allele_new/")
getwd()
mydata <- read.table("AP.alleleTable.unfiltered.txt", sep = "\t",header = T)

AP_info <- read.table("AP_info.csv", sep = ",",header = T)

AP_info2 <- read.table("sample_id", sep = "\t",header = F)






#mydata3 <- mydata[!duplicated(mydata[c("cellBC")]), ]

##########################################################
mydata_ab <- mydata[, c("Tumor", "sampleID")]
mydata_ab <- mydata_ab[!duplicated(mydata_ab[c("Tumor")]), ]

merged_trcr <- merge(mydata_ab, AP_info, by.x="sampleID", by.y="Tumor_combined")

trce <- merged_trcr[, c("Data_Inject", "Date_Sac", "Aging_day", "Cell_Clone", "Genotype", "MouseID", "Aging_Month", "Batch_Library", "Batch_Harvest", "Tumor_Name", "Sample_Name","Tumor", "Multi_BC", "X10x_Lane")]
write.table(trce, "trce.txt", sep="\t", row.names=F,quote=FALSE)
#############################################################

mydata2 <- mydata[!duplicated(mydata[c("cellBC")]), ]


merged_meta <- merge(mydata2, trce, by.x="Tumor", by.y="Tumor")

meta <- merged_meta[, c("cellBC", "X10x_Lane", "Tumor", "Cell_Clone", "Genotype", "MouseID", "Aging_Month", "Batch_Library", "Batch_Harvest", "Aging_day")]
rownames(meta) <- meta$cellBC
colnames(meta)[1] <- ",cellBC"

write.table(meta, "mete.csv", sep=",", row.names=T,quote=FALSE)


####################################
t1 <- trce[, c("Tumor", "Cell_Clone")]

file_t1 <- merge(mydata, t1, by.x="Tumor", by.y="Tumor")
write.table(file_t1, "AP.alleleTable.unfiltered.txt2", sep="\t", row.names=FALSE,quote=FALSE)
#######################
library(dplyr)
library(tidyr)

df <- file_t1[!duplicated(file_t1[c("cellBC")]), ]
# 按照克隆ID进行分组，统计每个克隆的细胞数量
df_count <- df %>%
  group_by(Tumor) %>%
  summarise(Cell_Count = n_distinct(cellBC))
df_count

# 将结果输出为数据框

