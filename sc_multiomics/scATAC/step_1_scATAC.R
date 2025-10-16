
#设置工作目录 
library(this.path)
cur_dir2 = dirname(this.path())
setwd(cur_dir2)
getwd()


httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(biomartr)
library(Seurat)
library(Signac)
library(ggplot2)
library(GenomicRanges)
library(future)

library("Matrix")
library("irlba")
#plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

#Merging objects
# read in peak sets
peaks.CTR <- read.table(
  file = "/datapool/yanzeqin/project/ZDW/scATAC/MC/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.NMI <- read.table(
  file = "/datapool/yanzeqin/project/ZDW/scATAC/NMI/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.MI <- read.table(
  file = "/datapool/yanzeqin/project/ZDW/scATAC/MI/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)


# convert to genomic ranges
gr.CTR <- makeGRangesFromDataFrame(peaks.CTR)
gr.NMI <- makeGRangesFromDataFrame(peaks.NMI)
gr.MI <- makeGRangesFromDataFrame(peaks.MI)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.CTR, gr.NMI, gr.MI))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#Create Fragment objects
# load metadata
md.CTR <- read.table(
  file = "/datapool/yanzeqin/project/ZDW/scATAC/MC/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.NMI <- read.table(
  file = "/datapool/yanzeqin/project/ZDW/scATAC/NMI/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.MI <- read.table(
  file = "/datapool/yanzeqin/project/ZDW/scATAC/MI/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


# perform an initial filtering of low count cells
md.CTR <- md.CTR[md.CTR$passed_filters > 500, ]
md.NMI <- md.NMI[md.NMI$passed_filters > 500, ]
md.MI <- md.MI[md.MI$passed_filters > 500, ]


# create fragment objects
frags.CTR <- CreateFragmentObject(
  path = "/datapool/yanzeqin/project/ZDW/scATAC/MC/outs/fragments.tsv.gz", 
  cells = rownames(md.CTR)
)

frags.NMI <- CreateFragmentObject(
  path = "/datapool/yanzeqin/project/ZDW/scATAC/NMI/outs/fragments.tsv.gz",
  cells = rownames(md.NMI)
)

frags.MI <- CreateFragmentObject(
  path = "/datapool/yanzeqin/project/ZDW/scATAC/MI/outs/fragments.tsv.gz",
  cells = rownames(md.MI)
)

print("start")
#Quantify peaks in each dataset
CTR.counts <- FeatureMatrix(
  fragments = frags.CTR,
  features = combined.peaks,
  cells = rownames(md.CTR)
)

NMI.counts <- FeatureMatrix(
  fragments = frags.NMI,
  features = combined.peaks,
  cells = rownames(md.NMI)
)

MI.counts <- FeatureMatrix(
  fragments = frags.MI,
  features = combined.peaks,
  cells = rownames(md.MI)
)

print("Finish:FeatureMatrix")

CTR_assay <- CreateChromatinAssay(CTR.counts, fragments = frags.CTR)
CTR <- CreateSeuratObject(CTR_assay, assay = "ATAC", meta.data=md.CTR)

NMI_assay <- CreateChromatinAssay(NMI.counts, fragments = frags.NMI)
NMI <- CreateSeuratObject(NMI_assay, assay = "ATAC", meta.data=md.NMI)

MI_assay <- CreateChromatinAssay(MI.counts, fragments = frags.MI)
MI <- CreateSeuratObject(MI_assay, assay = "ATAC", meta.data=md.MI)


print("Finish:CreateSeuratObject")
#Create the objects
print("Finish:")
#Merge objects
# add information to identify dataset of origin
CTR$dataset <- 'CTR'
NMI$dataset <- 'NMI'
MI$dataset <- 'MI'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = CTR,
  y = list(NMI, MI),
  add.cell.ids = c("CTR", "NMI", "MI")
)
print("Finish:merge")

combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)

saveRDS(combined, file = "combined_FindTopFeatures.rds")

combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

a1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
saveRDS(combined, file = "combined.rds")

ggsave(
  filename = "dimplot.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7,             # 宽
  height = 7,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)

print("Finish:all")
a2 <- CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr14-99700000-99760000"
)

ggsave(
  filename = "Finish.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7,             # 宽
  height = 7,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)


