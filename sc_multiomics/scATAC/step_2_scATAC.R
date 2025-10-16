#https://www.jianshu.com/p/c0edcd14ba66
#https://blog.csdn.net/ZIGRA/article/details/132312220
#https://stuartlab.org/signac/articles/overview # official tutorial
#设置工作目录 
setwd("/datapool/yanzeqin/project/ZDW/scATAC")
getwd()
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
#小鼠
#library(BSgenome.Mmusculus.UCSC.mm10)
#library(EnsDb.Mmusculus.v79)
#hg19
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(EnsDb.Hsapiens.v75)
#hg38版本
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(ggplot2)
#############################################################
#read rds
MGH <- readRDS("combined.rds")

############################################################
#更改路径
#old_path <- MGH@assays$ATAC@fragments[[1]]@path
# 新的路径
#new_path <- "/WorkDir3/DewangZhou/scATAC/MC/outs/fragments.tsv.gz"
# 更新路径
#MGH@assays$ATAC@fragments[[1]]@path <- new_path
# 验证更新
#print(MGH@assays$ATAC@fragments[[1]]@path)

#################################################################
# extract gene annotations from EnsDb.Hsapiens.v86
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(MGH) <- annotations
print(1)



############################################################
#计算QC的几个指标
#NC
MGH <- NucleosomeSignal(object = MGH)
print("finish")

#计算TSS
MGH <- TSSEnrichment(object = MGH , fast = FALSE)

MGH$pct_reads_in_peaks <- MGH$peak_region_fragments / MGH$passed_filters * 100
MGH$blacklist_ratio <- MGH$blacklist_region_fragments / MGH$peak_region_fragments

#TSS
MGH$high.tss <- ifelse(MGH$TSS.enrichment > 2, 'High', 'Low')
p <- TSSPlot(MGH, group.by = 'high.tss') + NoLegend()

ggsave("tss_1.png",width=6,height=6)

#NC
MGH$nucleosome_group <- ifelse(MGH$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = MGH, group.by = 'nucleosome_group')
ggsave("nc_1.png",width=6,height=6)

#
VlnPlot(
  object = MGH,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave("pct_1.png",width=6,height=6)
#根据不同指标过滤

MGH <- subset(
  x = MGH,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
MGH
print(2)
##################################################################
#归一化和线性降维
MGH <- RunTFIDF(MGH)
MGH <- FindTopFeatures(MGH, min.cutoff = 'q0')
MGH <- RunSVD(object = MGH)
DepthCor(MGH)
ggsave("depth_1.png",width=6,height=6)
###################################################################
#非线性降维和聚类
MGH <- RunUMAP(
  object = MGH,
  reduction = 'lsi',
  dims = 2:30
)
MGH <- FindNeighbors(
  object = MGH,
  reduction = 'lsi',
  dims = 2:30
)
MGH <- FindClusters(
  object = MGH,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

DimPlot(object = MGH, label = TRUE) + NoLegend()
ggsave("cluster_1.png",width=6,height=6)
####################################################################
#Create a gene activity matrix
# compute gene activities
gene.activities <- GeneActivity(MGH)

# add the gene activity matrix to the Seurat object as a new assay
MGH[['RNA']] <- CreateAssayObject(counts = gene.activities)
MGH <- NormalizeData(
  object = MGH,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(MGH$nCount_RNA)
)

print("feature")

DefaultAssay(MGH) <- 'RNA'
FeaturePlot(
  object = MGH,
  features = c('HMGA2','SNAI2',"VIM","IGFBP5","ZEB1","PDGFRA"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

ggsave("feature_1.png")

##################################################################
#Integrating with scRNA-seq data ?????\

#load the pre-processed scRNA-seq data
#allen_rna <- readRDS("MGH3.rds")

#allen_rna
#allen_rna <- UpdateSeuratObject(allen_rna)

# 提取表达矩阵
#data_matrix <- GetAssayData(allen_rna, slot = "data")

# 检查是否有NA/NaN/Inf值
#if (anyNA(data_matrix) || any(is.nan(data_matrix)) || any(is.infinite(data_matrix))) {
  # 将NA/NaN/Inf值替换为0
#  data_matrix[is.na(data_matrix)] <- 0
#  data_matrix[is.nan(data_matrix)] <- 0
#  data_matrix[is.infinite(data_matrix)] <- 0
  
  # 更新Seurat对象中的数据
#  allen_rna <- SetAssayData(allen_rna, slot = "data", new.data = data_matrix)
#}

# 归一化数据
#allen_rna <- NormalizeData(allen_rna)


#print("rna1")
#allen_rna <- FindVariableFeatures(
#  object = allen_rna,
#  nfeatures = 5000
#)

print("#######################################################\nrna")
#transfer.anchors <- FindTransferAnchors(
#  reference = allen_rna,
#  query = MGH,
#  reduction = 'cca',
#  dims = 1:30
#)

#predicted.labels <- TransferData(
#  anchorset = transfer.anchors,
#  refdata = allen_rna$subclass,
#  weight.reduction = MGH[['lsi']],
#  dims = 2:30
#)

#MGH <- AddMetaData(object = MGH, metadata = predicted.labels)

#plot1 <- DimPlot(allen_rna, group.by = 'subclass', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
#plot2 <- DimPlot(MGH, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
#plot1 + plot2

#ggsave("rna_MGH_1.png",width=6,height=6)
#print(3)
######################################################################
#查找clster间差异peaks区域
#switch back to working with peaks instead of gene activities
#DefaultAssay(MGH) <- 'peaks'
#Idents(MGH) <- "predicted.id"

#da_peaks <- FindMarkers(
#  object = MGH,
#  ident.1 = c("L2/3 IT"), 
#  ident.2 = c("L4", "L5 IT", "L6 IT"),
#  test.use = 'LR',
#  latent.vars = 'nCount_peaks'
#)

#Head(da_peaks)

 
#plot1 <- VlnPlot(
#  object = MGH,
#  features = rownames(da_peaks)[1],
#  pt.size = 0.1,
#  idents = c("L4","L5 IT","L2/3 IT")
#)
#plot2 <- FeaturePlot(
#  object = MGH,
#  features = rownames(da_peaks)[1],
#  pt.size = 0.1,
#  max.cutoff = 'q95'
#)
#plot1 | plot2
#ggsave("peaks_1.png",width=6,height=6)
#print(4)
#######################################################################
#基因组绘图
#how cell types with at least 50 cells
#idents.plot <- names(which(table(Idents(brain)) > 50))

#CoveragePlot(
#  object = brain,
#  region = c("Neurod6", "Gad2"),
#  idents = idents.plot,
#  extend.upstream = 1000,
#  extend.downstream = 1000,
#  ncol = 1
#)

####################################################################
#motif 分析
#加载数据，查看cluster
p1 <- DimPlot(MGH, label = TRUE, pt.size = 0.1) + NoLegend()
p1
ggsave("cheak_cluster_1.png",width=6,height=6)
#构建 motif 类

#get  list of motif position frequency matrices from the JASPAR database
#使用 getrMatrixSet函数从JASPAR数据库中提取Motif的PFM矩阵信息 
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
MGH <- AddMotifs(
  object = MGH,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
print("11")
#finding overrepresented motifs
da_peaks <- FindMarkers(
  object = MGH,
  ident.1 = 'Pvalb',
  ident.2 = 'Sst',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

#test enrichment
#使用findmotifs函数进行motif富集分析
enriched.motifs <- FindMotifs(
  object = MGH,
  features = top.da.peak
)

MotifPlot(
	object = MGH,
	motifs = head(rawnames(enriched.motifs))
)
ggsave("find_motifs_1.png",width=6,height=6)
#computing motif acitvities
MGH <- RunChromVAR(
  object = MGH,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(MGH) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = MGH,
  features = "HMGA2",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2
ggsave("activaty_motifs_1.png",width=6,height=6)

#
differential.activity <- FindMarkers(
  object = MGH,
  ident.1 = 'Pvalb',
  ident.2 = 'Sst',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = MGH,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)
ggsave("f_motifs_2.png",width=6,height=6)

#################################################################################
#visualization 

