
#设置工作目录 
library(this.path)
cur_dir2 = dirname(this.path())
setwd(cur_dir2)
Getwd()


httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(biomartr)
library(Seurat)
library(Signac)
library(ggplot2)
library(GenomicRanges)
library(future)


#plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

#Merging objects
# read in peak sets
peaks.CTR <- read.table(
  file = "/WorkDir3/DewangZhou/scATAC/MC/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.NMI <- read.table(
  file = "/WorkDir3/DewangZhou/scATAC/NMI/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.MI <- read.table(
  file = "/WorkDir3/DewangZhou/scATAC/MI/outs/filtered_peak_bc_matrix/peaks.bed",
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
  file = "/WorkDir3/DewangZhou/scATAC/MC/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.NMI <- read.table(
  file = "/WorkDir3/DewangZhou/scATAC/NMI/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.MI <- read.table(
  file = "/WorkDir3/DewangZhou/scATAC/MI/outs/singlecell.csv",
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
  path = "/WorkDir3/DewangZhou/scATAC/MC/outs/fragments.tsv.gz", 
  cells = rownames(md.CTR)
)

frags.NMI <- CreateFragmentObject(
  path = "/WorkDir3/DewangZhou/scATAC/NMI/outs/fragments.tsv.gz"
)

frags.MI <- CreateFragmentObject(
  path = "/WorkDir3/DewangZhou/scATAC/MI/outs/fragments.tsv.gz"
)


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

CTR_assay <- CreateChromatinAssay(CTR.counts, fragments = frags.CTR)
CTR <- CreateSeuratObject(CTR_assay, assay = "ATAC", meta.data=md.CTR)

NMI_assay <- CreateChromatinAssay(NMI.counts, fragments = frags.NMI)
NMI <- CreateSeuratObject(NMI_assay, assay = "ATAC", meta.data=md.NMI)

MI_assay <- CreateChromatinAssay(MI.counts, fragments = frags.MI)
MI <- CreateSeuratObject(MI_assay, assay = "ATAC", meta.data=md.MI)

#Create the objects
CTR_assay <- CreateChromatinAssay(CTR.counts, fragments = frags.CTR)
CTR <- CreateSeuratObject(CTR_assay, assay = "ATAC", meta.data=md.CTR)

NMI_assay <- CreateChromatinAssay(NMI.counts, fragments = frags.NMI)
pbmc1k <- CreateSeuratObject(NMI_assay, assay = "ATAC", meta.data=md.NMI)

MI_assay <- CreateChromatinAssay(MI.counts, fragments = frags.MI)
MI <- CreateSeuratObject(MI_assay, assay = "ATAC", meta.data=md.MI)

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
combined[["ATAC"]]


combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, group.by = 'dataset', pt.size = 0.1)

CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr14-99700000-99760000"
)


