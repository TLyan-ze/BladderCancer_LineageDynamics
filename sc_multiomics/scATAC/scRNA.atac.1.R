library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(readxl)



# hg38版本
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

library(patchwork)
library(ggplot2)



setwd("/datapool/yanzeqin/project/ZDW/scATAC")
getwd()

# reference --- https://stuartlab.org/signac/articles/pbmmc_multiomic.html
# reference --- https://stuartlab.org/signac/articles/snareseq.html
# reference --- https://stuartlab.org/signac/articles/mouse_brain_vignette.html

# get annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "Ensembl"
genome(annotations) <- "hg38"

# read RNA and ATAC data
samplenames=c("MC", "MI", "NMI")
counts=list()
for(i in samplenames)
{	
	dir.name = paste("./",i,"/",i,"_filtered_feature_bc_matrix.h5", sep="")
	counts[[i]] = Read10X_h5(dir.name)
}

project.list=list()
for(i in samplenames)
{
project.list[[i]]<- CreateSeuratObject(
  counts = counts[[i]]$Gene,
  assay = "RNA")
}

# scATAC
for(i in samplenames)
{
	fragpath <- paste("./", i, "/out/fragments.tsv.gz", sep="")
	project.list[[i]][["ATAC"]]<- CreateChromatinAssay( counts = counts[[i]]$Peaks, sep = c(":", "-"), fragments = fragpath, genome="hg38" )
	Annotation(project.list[[i]][["ATAC"]]) <- annotations
}
###############################################################
# 质量控制，包括TSS富集分析，生成小提起图
# check quality 
for(i in samplenames)
{
	DefaultAssay(project.list[[i]]) <- "ATAC"
	project.list[[i]] <- TSSEnrichment(project.list[[i]])
	project.list[[i]] <- NucleosomeSignal(project.list[[i]])
	# project.list[[i]]$blacklist_fraction <- FractionCountsInRegion(
	#   object = project.list[[i]],
	#   assay = 'ATAC',
	#   regions = blacklist_mm10
	# )
}

for(i in samplenames)
{
	pdf(paste("./result/",i, ".check.quality", ".pdf", sep=""), width=10, height=5)
	temp=VlnPlot(
	  object = project.list[[i]],
	  features = c("nCount_RNA", "nCount_ATAC", "nucleosome_signal", "TSS.enrichment"),
	  ncol = 4,
	  pt.size = 0
	)
	print(temp)
	# , "blacklist_fraction"
	dev.off()
}

# remove cells with potential low quality
for(i in samplenames)
{
	project.list[[i]] <- subset(
	  x = project.list[[i]],
	  subset = nCount_ATAC < 100000 &
	    nCount_RNA < 30000 &
	    nCount_ATAC > 1000 &
	    nCount_RNA > 1000 &
	    nucleosome_signal < 1.5 &
	    TSS.enrichment > 1
	)
	print(project.list[[i]])
}
# project.list.bak=project.list

######################################################################
# 使用MACS2调用峰值，并将其添加到seurat对象
# call peaks
# peaks.bak=list()
for(i in samplenames[1:2])
{
	# Peak calling
	# call peaks using MACS2
	peaks <- CallPeaks(project.list[[i]])
	# peaks.bak[[i]]=peaks

	# remove peaks on nonstandard chromosomes and in genomic blacklist regions
	peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
	peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38, invert = TRUE)
	# Warning message:
	# In .merge_two_Seqinfo_objects(x, y) :
	#   The 2 combined objects have no sequence levels in common. (Use
	#   suppressWarnings() to suppress this warning.)

	# quantify counts in each peak
	macs2_counts <- FeatureMatrix(
	  fragments = Fragments(project.list[[i]]),
	  features = peaks,
	  cells = colnames(project.list[[i]])
	)

	# create a new assay using the MACS2 peak set and add it to the Seurat object
	# CreateChromatinAssay-"fragments": Path to a tabix-indexed fragments file for the data
	#           contained in the input matrix. If multiple fragment files are
	#           required, you can add additional 'Fragment' object to the
	#           assay after it is created using the 'CreateFragmentObject'
	#           and 'Fragments' functions. Alternatively, a list of
	#           'Fragment' objects can be provided.

	fragpath <- paste("./", i, "/out/fragments.tsv.gz", sep="")
	project.list[[i]][["peaks"]] <- CreateChromatinAssay(
	  counts = macs2_counts,
	  fragments = fragpath,
	  annotation = annotations
	)
}

### correct sampleID "Sham-3m", "ADT-1m", "Enz-2m"
{
	project.list[["MC"]]@meta.data$sample.id="MC"
	project.list[["MI"]]@meta.data$sample.id="MI"
	project.list[["NMI"]]@meta.data$sample.id="NMI"
}

################################################################
# 合并样本
### merge 3 samples
project=merge(x = project.list[[1]], y = project.list[2:length(project.list)])
# project.bak=project
table(project@meta.data$sample.id)
saveRDS(project, "project.merged.rds")



project=readRDS("project.merged.rds")
############################################################################
# RNA数据预处理，归一化，高变基因，数据缩放和PCA降维
### preprocessing RNA data
{
	DefaultAssay(project)="RNA"
	project <- NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000)

	str(project@assays$RNA@var.features)
	project <- FindVariableFeatures(project, assay="RNA", nfeatures = 3000)
	str(project@assays$RNA@var.features)
	# chr [1:3000] "Sbp" "Spink1" "Tcaf3" "Pbsn" "Spink5" "Defb50" "Klk1" ...
	top10 <- head(VariableFeatures(project, assay="RNA"), 10)
	# plot variable features with and without labels
	plot1 <- VariableFeaturePlot(project, assay = "RNA")
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	toplot=plot1 + plot2
	png(paste("./result/temp.RNA.preprocess", ".png", sep=""), width=10*100, height=5*100)
	print(toplot)
	dev.off()

	project <- ScaleData(project, features = project@assays$RNA@var.features, assay="RNA" )

	project <- RunPCA(project, assay="RNA", verbose = FALSE, npcs = 100)
	ElbowPlot(project,  ndims = 100)
	dev.off()

}

################################################################################
# 对peaks 数据进行预处理，可视化，降维和聚类
### preprocessing Peaks data
{
	DefaultAssay(project) <- "peaks"
	project <- FindTopFeatures(project, min.cutoff = 5)
	project <- RunTFIDF(project)
	project <- RunSVD(project)
}

### Visualization based on RNA + peaks
{
	# Identify multimodal neighbors. These will be stored in the neighbors slot, 
	# and can be accessed using project[['weighted.nn']]
	# The WNN graph can be accessed at project[["wknn"]], 
	# and the SNN graph used for clustering at project[["wsnn"]]
	# Cell-specific modality weights can be accessed at project$RNA.weight
	project <- FindMultiModalNeighbors(
	  object = project,
	  reduction.list = list("pca", "lsi"), 
	  dims.list = list(1:50, 2:40),
	  modality.weight.name = c("RNA.weight", "peaks.weight"),
	  verbose = TRUE
	)

	project <- RunUMAP(
	  object = project,
	  nn.name = "weighted.nn",
	  assay = "RNA",
	  verbose = TRUE
	)

	DefaultAssay(project)="RNA"
	project <- FindClusters(project, verbose = FALSE, resolution = 0.8, graph.name = "wsnn", algorithm = 3)
	project@meta.data$fig.cluster=project@meta.data$seurat_clusters

	# saveRDS(project, "project.visualized.rds")
	# project=readRDS("project.visualized.rds")

	pdf(paste("./result/temp.project.umap.cluster", ".pdf", sep=""), height=10, width=12.5)
	DimPlot(project, label = TRUE, repel = TRUE, reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="fig.cluster")# + NoLegend()
	dev.off()

	pdf(paste("./result/temp.project.umap.sample", ".pdf", sep=""), height=10, width=12.5)
	DimPlot(project, label = TRUE, repel = TRUE, reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="sample.id")# + NoLegend()
	dev.off()

	# check.genes=c("Trp63", "Ar", "Krt8")
	check.genes=c("Ar")
	"Ar"%in%(rownames(project@assays$RNA@scale.data))
	DefaultAssay(project)="RNA"
	# SPG1 SPG2 SPG3(Sycp-low) SPG4(SYCP3-high Prdm9-high) Scytes CellCycle
	pdf(paste("./result/project.featureplot", ".pdf", sep=""), height=10, width=14)
	FeaturePlot(project, features = check.genes, raster=TRUE, pt.size=0.5, raster.dpi = c(1024, 1024))
	dev.off()

}
######################################################################################
# 细胞类型注释，和去除非上皮细胞，这一步是否要去掉

### find cell types (to remove remaining stromal and immune cells)
{
	project <- SetIdent(project, value = "fig.cluster")
	project.markers <- FindAllMarkers(project, min.pct=0.2, max.cells.per.ident=500, assay="RNA")
	# saveRDS(project.markers, "project.markers.by.fig.cluster.rds")
	# project.markers=readRDS("project.markers.by.fig.cluster.rds")

	### write markers of clusters as a table
	markers=project.markers
	markers=markers[order(markers$cluster, markers$avg_log2FC, decreasing = c(FALSE,TRUE)), ]
	markers$foldchange=2^(markers$avg_log2FC)
	markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.25) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
	write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("./result/markers.project.by.fig.cluster","csv",sep="."), row.name=T)

	### check celltype markers' expression to make sure every cluster's celltype
	markers=project.markers
	markers=markers[markers$p_val_adj<0.05, ]
	markers=markers[markers$avg_log2FC>log2(1.5) & markers$pct.1>0.2, ]
	### 注释另一个方法
	
	#library(Seurat)
	#library(tidyverse)
	#library(SingleR)
	#library(celldex)
	#library(RColorBrewer)
	#library(SingleCellExperiment)
	

	#sce <- as.SingleCellExperiment(DietSeurat(srat))
	#sce
	
	### celltypeAnno() need to be run in directory "processing"
	source("./functions.seurat.R")   #分细胞类型
	temp.celltype=celltypeAnno(project, markers, mouse=FALSE)
	# [1] "Decided cell type:"
	#  [1] "Luminal cell: 0.59"
	#  [2] "Basal cell_profRen: 0.95"
	#  [3] "Myofibroblast: -0.06"
	#  [4] "Basal cell_profRen: 0.32"
	#  [5] "Luminal cell: 0.42"
	#  [6] "MDSCfunction_profRen: 0.47"
	#  [7] "B/Plasma cell_cell201801: 0.59"
	#  [8] "Luminal cell: 0.24"
	#  [9] "Luminal cell: 0.26"
	# [10] "B/Plasma cell_cell201801: 0.7"
	# [11] "Basal cell_profRen: 1"
	# [12] "Neutrophil: 0.38"
	# [13] "Basal cell: 1.29"
	# [14] "Basal cell_profRen: 0.55"
	# [15] "Mast cell_cell201808: 0.48"
	# [16] "Mast cell_cell201808: 0.35"
	# [17] "Mast cell_cell201808: 0.3 Luminal cell: 0.29"
	# [18] "Fibroblast_cell201801: 4.45"
	# [19] "B cell Lymph node: 1.91 Dendritic cell_cell201801: 1.91"
	# [20] "B cell: 0.9"
	# [21] "B cell: 0.67"
	# [22] "Endothelial cell: 2.94"
	# [23] "T cell: 3.81"
	# [24] "T cell: 0.58"
	# [25] "Mast cell_cell201808: 0.89"
	# [26] "Mast cell_cell201808: 4.55"
	### cluster 12: Afap1l2+ Dcn+ Pdpn+ Apoe+ Col17a1+ Krt5+ --- under considering
	### cluster 13: cell-cycling
	### cluster 17: Dcn+ Pdpn+
	### cluster 21: Fap+ Krt.all- 
	### cluster 22: NE cells Ncam1+ Chga+ Ascl1+
	### cluster 23: Krt.all-
	### cluster 24: Ms4a8a+ Krt8+1.63fc Krt18+1.75fc
	### cluster 25: Kit+ Krt.all-

	### remove non-epithelial clusters and cell-cycling clusters
	project=subset(project, cells=colnames(project)[!project@meta.data$fig.cluster%in%c("13", "17", "21", "23", "24", "25")])
	### remove cluster 2 (mt genes and ribosome genes high expresed) and cluster 22 (NE but EPCAM-)
	project=subset(project, cells=colnames(project)[!project@meta.data$fig.cluster%in%c("2", "22")])

	pdf("./result/temp.project.umap.cluster.pdf", height=10, width=12.5)
	DimPlot(project, label = TRUE, repel = TRUE, reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="fig.cluster")# + NoLegend()
	dev.off()

	pdf("./result/temp.project.umap.sample.pdf", height=10, width=12.5)
	DimPlot(project, label = TRUE, repel = TRUE, reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="sample.id")# + NoLegend()
	dev.off()

	markers.to.plot <- c("Krt5","Krt6a","Ar", "Krt8")
	pdf("./result/temp.dotplot.pdf")
	DotPlot(project, features = markers.to.plot, cols = c("gray", "Red"), dot.scale = 6) + RotatedAxis()
	dev.off()

}

#######################################
#去除低质量细胞后，再重新归一化数据，聚类
### re-normalize, visualizing and clustering after removing cells
### preprocessing RNA data
{
	DefaultAssay(project)="RNA"
	project <- NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000)

	str(project@assays$RNA@var.features)
	# chr [1:3000] "Sbp" "Spink1" "Tcaf3" "Pbsn" "Spink5" "Defb50" "Klk1" ...
	project <- FindVariableFeatures(project, assay="RNA", nfeatures = 3000)
	str(project@assays$RNA@var.features)
	# chr [1:3000] "Sbp" "Spink1" "Tcaf3" "Pbsn" "Defb50" "Spink5" "Klk1" ...
	top10 <- head(VariableFeatures(project, assay="RNA"), 10)
	# plot variable features with and without labels
	plot1 <- VariableFeaturePlot(project, assay = "RNA")
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	toplot=plot1 + plot2
	png("./result/prepro_RNA_temp.png", width=10*100, height=5*100)
	print(toplot)
	dev.off()

	project <- ScaleData(project, features = project@assays$RNA@var.features, assay="RNA" )

	project <- RunPCA(project, assay="RNA", verbose = FALSE, npcs = 100)
	ElbowPlot(project,  ndims = 100)
	dev.off()

}



############################################
#可视化 RNA和peak

### preprocessing Peaks data
{
	DefaultAssay(project) <- "peaks"
	project <- FindTopFeatures(project, min.cutoff = 5)
	project <- RunTFIDF(project)
	project <- RunSVD(project)

}

### Visualization based on RNA + peaks
{
	require(harmony)
	#  theta: Diversity clustering penalty parameter. Specify for each
  #         variable in group.by.vars. Default theta=2. theta=0 does not
  #         encourage any diversity. Larger values of theta result in
  #         more diverse clusters.

  # lambda: Ridge regression penalty parameter. Specify for each variable
  #         in group.by.vars. Default lambda=1. Lambda must be strictly
  #         positive. Smaller values result in more aggressive
  #         correction.

  #  sigma: Width of soft kmeans clusters. Default sigma=0.1. Sigma
  #         scales the distance from a cell to cluster centroids. Larger
  #         values of sigma result in cells assigned to more clusters.
  #         Smaller values of sigma make soft kmeans cluster approach
  #         hard clustering.

  # nclust: Number of clusters in model. nclust=1 equivalent to simple
  #         linear regression.
	DefaultAssay(project)
	# [1] "RNA"
	project <- RunHarmony(project, group.by.vars = "sample.id", assay.use="RNA", reduction = 'pca', plot_convergence = TRUE, dims.use=1:30, reduction.save="pca.harmony", theta=2, lambda=2, sigma=0.01, nclust=5) 
	dev.off()

	project <- RunHarmony(project, group.by.vars = "sample.id", assay.use = 'peaks', reduction = 'lsi', plot_convergence = TRUE, dims.use=1:30, reduction.save="lsi.harmony", theta=2, lambda=2, sigma=0.01, nclust=5, project.dim = FALSE )
	dev.off()

	# Identify multimodal neighbors. These will be stored in the neighbors slot, 
	# and can be accessed using project[['weighted.nn']]
	# The WNN graph can be accessed at project[["wknn"]], 
	# and the SNN graph used for clustering at project[["wsnn"]]
	# Cell-specific modality weights can be accessed at project$RNA.weight
	project <- FindMultiModalNeighbors(
	  object = project,
	  reduction.list = list("pca.harmony", "lsi.harmony"), 
	  dims.list = list(1:30, 2:30),
	  modality.weight.name = c("RNA.weight", "peaks.weight"),
	  verbose = TRUE,
	)

	project <- RunUMAP(
	  object = project,
	  nn.name = "weighted.nn",
	  assay = "RNA",
	  verbose = TRUE
	)

	DefaultAssay(project)="RNA"
	project <- FindClusters(project, verbose = FALSE, resolution = 0.3, graph.name = "wsnn", algorithm = 3)
	project@meta.data$fig.cluster=project@meta.data$seurat_clusters

	project@meta.data$sample.id=factor(project@meta.data$sample.id, levels=c("MC", "MI", "NMI"))

	# saveRDS(project, "project.visualized.rds") # no-harmony version
	# project=readRDS("project.visualized.rds") # no-harmony version

	saveRDS(project, "project.visualized.harmony.rds")
	# project=readRDS("project.visualized.harmony.rds")


	pdf("./result/temp.project.umap.cluster.harmony.pdf", height=5, width=6)
	temp=DimPlot(project, label = TRUE, repel = TRUE, , reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="fig.cluster", label.size = 4)# + NoLegend()
	print(temp)
	dev.off()

	pdf("./result/temp.project.umap.sample.harmony.pdf", height=5, width=6)
	temp=DimPlot(project, label = FALSE, repel = TRUE, reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="sample.id", label.size = 4)# + NoLegend()
	print(temp)
	dev.off()

}

###################################################
#分别基于RNA和ATAC进行姜维和聚类

### get UMAP and clustering based on only RNA or ATAC
{
	project=readRDS("project.visualized.harmony.rds")
	# project.bak=project
	# project=project.bak

	# redu = "pca.harmony"  OR  "lsi.harmony"
	redu="lsi.harmony"
	if(redu=="pca.harmony")
	{
		assay="RNA"
	}else
	{
		assay="peaks"
	}
	project <- FindNeighbors(project, reduction = redu, dims = 1:30, verbose = FALSE)

	project <- RunUMAP(
	  object = project,
	  reduction = redu,
	  assay = assay,
	  dims = 1:30, 
	  verbose = FALSE
	)
	# nn.name = "weighted.nn",

	DefaultAssay(project)=assay
	project <- FindClusters(project, verbose = FALSE, resolution = 0.3, graph.name = paste(assay, "nn", sep="_"), algorithm = 3)
	project@meta.data$fig.cluster=project@meta.data$seurat_clusters

	project@meta.data$sample.id=factor(project@meta.data$sample.id, levels=c("MC", "MI", "NMI"))
	
	pdf(paste("./result/temp.project.umap.cluster.harmony", assay, "pdf", sep="."), height=5, width=6)
	temp=DimPlot(project, label = TRUE, repel = TRUE, , reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="fig.cluster", label.size = 4, pt.size=2)# + NoLegend()
	print(temp)
	dev.off()

	pdf(paste("./result/temp.project.umap.sample.harmony", assay, "pdf", sep="."), height=5, width=6)
	temp=DimPlot(project, label = FALSE, repel = TRUE, reduction = "umap", raster = TRUE, raster.dpi = c(1024, 1024), group.by="sample.id", label.size = 4, pt.size=2)# + NoLegend()
	print(temp)
	dev.off()

}
##############################################################
#检查特定基因的表达情况，并生成图

### check gene expresion
{
	DefaultAssay(project)="RNA"
	project <- SetIdent(project, value = "fig.cluster")
	markers.to.plot <- c("HMGA2","SNAI2","IGFBP5","ZEB1")
	pdf("./result/temp.dotplot.pdf")
	DotPlot(project, features = markers.to.plot, cols = c("gray", "Red"), dot.scale = 6) + RotatedAxis()
	dev.off()

	# check.genes=c("Trp63", "Ar",, "Krt5", "Fkbp5", "Acly", "Krt6a","Lypd3", "Krt8")
	# check.genes=c("nCount_RNA", "nFeature_RNA",  "percent.mt")
	check.genes=c("ZEB1")
	"Ar"%in%(rownames(project@assays$RNA@scale.data))
	DefaultAssay(project)="RNA"
	# SPG1 SPG2 SPG3(Sycp-low) SPG4(SYCP3-high Prdm9-high) Scytes CellCycle
	pdf("./result/project.featureplot.pdf", height=10, width=14)
	FeaturePlot(project, features = check.genes, raster=TRUE, pt.size=0.5, raster.dpi = c(1024, 1024))
	dev.off()

}

###############################################################
# monocle 轨迹推断
### monocle 
{
	library(monocle)
	# ### fig.cluster DEG
	# {
	# 	project <- SetIdent(project, value = "fig.cluster")
	# 	project.markers <- FindAllMarkers(project, min.pct=0.2, max.cells.per.ident=500, assay="RNA")
	# 	# saveRDS(project.markers, "project.markers.by.fig.cluster.harmony.rds")
	# 	# project.markers=readRDS("project.markers.by.fig.cluster.harmony.rds")

	# 	### write markers of clusters as a table
	# 	markers=project.markers
	# 	markers=markers[order(markers$cluster, markers$avg_log2FC, decreasing = c(FALSE,TRUE)), ]
	# 	markers$foldchange=2^(markers$avg_log2FC)
	# 	markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5), ]
	# 	write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.project.by.fig.cluster.harmony","csv",sep="."), row.name=T)
	# }

	### sample.id DEG
	{
		project <- SetIdent(project, value = "sample.id")
		project.markers <- FindAllMarkers(project, min.pct=0.2, max.cells.per.ident=500, assay="RNA")
		saveRDS(project.markers, "project.markers.by.sample.id.harmony.rds")
		# project.markers=readRDS("project.markers.by.sample.id.harmony.rds")

		### write markers of clusters as a table
		markers=project.markers
		markers=markers[order(markers$cluster, markers$avg_log2FC, decreasing = c(FALSE,TRUE)), ]
		markers$foldchange=2^(markers$avg_log2FC)
		markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5), ]
		write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.project.by.sample.id.harmony","csv",sep="."), row.name=T)
	}

	# turn to monocle2 conda environment to process data
	# conda activate monocle2
	library(monocle)
	library(Seurat)
	# project=readRDS("project.visualized.harmony.rds")
	temp=project
	temp=subset(project, cells=colnames(project)[project@meta.data$fig.cluster%in%c(3,2,6,5,1)])
	### use conda environment r423 to subset the "temp" object and save&read
	# saveRDS(temp, "temp.seurat.for.monocle.rds")
	# temp=readRDS("temp.seurat.for.monocle.rds")
	gene_metadata=data.frame(gene_short_name=rownames(temp@assays$RNA@counts))
	rownames(gene_metadata)=rownames(temp@assays$RNA@counts)  
	### https://github.com/cole-trapnell-lab/monocle-release/issues/262
	project.m <- newCellDataSet( temp@assays$RNA@counts,
                            phenoData = new("AnnotatedDataFrame", temp@meta.data[, c("fig.cluster", "sample.id","nFeature_RNA")]),
                            featureData = new("AnnotatedDataFrame", gene_metadata),
                            expressionFamily=negbinomial.size() )

	### Plan B to create monocle object --- https://github.com/cole-trapnell-lab/monocle-release/issues/262
	# {
	# 	data <- temp@assays$RNA@counts
	# 	pd <- new('AnnotatedDataFrame', data = temp@meta.data)
	# 	fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
	# 	fd <- new('AnnotatedDataFrame', data = fData)
	# 	#Construct monocle cds
	# 	project.m <- newCellDataSet(data,
	# 	                         phenoData = pd,
	# 	                         featureData = fd,
	# 	                         lowerDetectionLimit = 0.5,
	# 	                         expressionFamily = negbinomial.size())
	# }

	table(project.m@phenoData@data$fig.cluster)
	table(project.m@phenoData@data$sample.id)
	project.m <- estimateSizeFactors(project.m)
	project.m <- estimateDispersions(project.m)

	project.m <- detectGenes(project.m, min_expr = 0.1)
  	## only keep expressed genes
	# expressed_genes <- row.names(project.m)[project.m@featureData@data$num_cells_expressed>= 10]
	# project.m <- project.m[expressed_genes,]

	### use all significant markers of clusters as ordering genes
	# project.markers=readRDS("project.markers.by.fig.cluster.harmony.rds")
	project.markers=readRDS("project.markers.by.sample.id.harmony.rds")
	markers=project.markers
	markers=markers[markers$p_val_adj<0.05, ]
	markers=markers[(markers$avg_log2FC)>log2(1.5), ]
	markers$foldChange=2^(markers$avg_log2FC)
	# markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
	markers=markers[order(markers$foldChange, decreasing=TRUE), ]
	ordering.genes=unique(markers$gene)
	# ordering.genes=unique(rownames(markers))
	str(ordering.genes)
	# fig.cluster markers: foldchange>2
	# chr [1:162] "Sbp" "Spink1" "Defb50" "Lgals7" "Spink5" "Ltf" ...
	# sample.id markers: foldchange>1.5
	# chr [1:488] "Sbp" "Spink1" "Defb50" "Lgals7" "Spink5" "Ltf" ...

	# ### identify pair-wise DEGs for ordering cells
	# {
	# 	table(project@meta.data$fig.cluster)
	# 	#    0    1    2    3    4    5    6    7    8    9
	# 	# 4209 3568 2309 1975  744  510  412  326  280  180
	# 	DefaultAssay(project)="RNA"
	# 	### pairwise customized cluster's DEGs
	# 	markerlist=list()
	# 	fig.cluster=unique(project@meta.data$fig.cluster)
	# 	fig.cluster=sort(fig.cluster)
	# 	project = SetIdent(project, value = "fig.cluster")
	# 	for(i in fig.cluster)
	# 	{
	# 	  for(j in fig.cluster)
	# 	  {
	# 	      if(i!=j)
	# 	      {
	# 	        temp.name=paste(i, j, sep="__")
	# 	        markerlist[[temp.name]]=FindMarkers(project, ident.1 = i, ident.2 = j, slot="data", assay="RNA")
	# 	      }
	# 	  }
	# 	}
	# 	# saveRDS(markerlist, "project.fig.cluster.pairwise.rds")
	# 	# markerlist=readRDS("project.fig.cluster.pairwise.rds")

	# }

	# ### use all significant pair-wise epitype's DEGs as ordering genes
  # {
  #   unicluster=unique(project@meta.data$fig.cluster)
  #   markerlist.project=readRDS("project.fig.cluster.pairwise.rds")
  #   # str(project@assays$SCT@var.features)
  #   geneset=c()
  #   print("Including clusters:")
  #   print(unicluster)
  #   for(i in names(markerlist.project))
  #   {
  #     j=markerlist.project[[i]]
  #     mn=strsplit(i, split="__")[[1]]
  #     m=mn[1]
  #     n=mn[2]
  #     if( (m%in%unicluster)&(n%in%unicluster) )
  #     {
  #       # j=j[ (j$p_val_adj<0.05) & ( (j$avg_log2FC)>log2(1.5) ), ]
  #       j=j[ (j$p_val_adj<0.05) & ( (j$avg_log2FC)>log2(1.5) ), ]
  #       j=j[ order(j$avg_log2FC, decreasing=TRUE), ]
  #       # if(nrow(j)>5)
  #       # {j=j[1:5, ]}
  #       newset=rownames(j)
  #       geneset=union( geneset, newset )
  #       #  )
  #     }
  #   }
  #   print(str(geneset))
  #   # 
  #   # (j$p_val_adj<0.05) & ( (j$avg_log2FC)>log2(2) ) :
  #   # ..
  #   # (j$p_val_adj<0.05) & ( abs(j$avg_log2FC)>log2(2) ) :
  #   # 
  #   # (j$p_val_adj<0.05) & ( (j$avg_log2FC)>log2(1.5) ) :
  #   # chr [1:2955] "Clu" "Grid2" "Slc12a2" "Pde4d" "Cxcl2" "Pigr" "Basp1" ...
  #   ordering.genes=unique(geneset)
  # }

	project.m <- setOrderingFilter(project.m,  ordering.genes) # 
	pdf("plot_ordering_genes.pdf")
	plot_ordering_genes(project.m)
	dev.off()

	project.m <- reduceDimension(project.m, max_components = 2,method = 'DDRTree')

	project.m <- orderCells(project.m)
	pdf("plot_cell_trajectory_fig.cluster_project.m.pdf", width=5, height=5)
	plot_cell_trajectory(project.m, color_by = "fig.cluster")
	dev.off()

	pdf("plot_cell_trajectory_details_fig.cluster_project.m.pdf", width=30, height=6)
	temp=plot_cell_trajectory(project.m, color_by = "fig.cluster") + facet_wrap(~fig.cluster, nrow = 1)
	print(temp)
	dev.off()

	pdf("plot_cell_trajectory_State_project.m.pdf", width=5, height=5)
	# png("plot_cell_trajectory_byPseudotime_project.m.png")
	plot_cell_trajectory(project.m, color_by = "State", cell_size=0.05, cell_name_size = 8)
	dev.off()

	# set a indicator for time, and sort again
	project.m<-orderCells(project.m, root_state = 1)
	table(rownames(project.m@phenoData@data)==rownames(project@meta.data))
	# TRUE
	# 5255
	# project.m@phenoData@data$fig.celltype=project@meta.data$celltype


	pdf("plot_cell_trajectory_byPseudotime_project.m.pdf", width=5, height=5)
	# png("plot_cell_trajectory_byPseudotime_project.m.png")
	plot_cell_trajectory(project.m, color_by = "Pseudotime", cell_size=0.4)
	dev.off()

	pdf("plot_cell_trajectory_details_fig.cluster_project.m.pdf", width=30, height=6)
	temp=plot_cell_trajectory(project.m, color_by = "fig.cluster") + facet_wrap(~fig.cluster, nrow = 1)
	print(temp)
	dev.off()

	pdf("plot_cell_trajectory_fig.cluster_project.m.pdf", width=5, height=5)
	# png("plot_cell_trajectory_byPseudotime_project.m.png")
	plot_cell_trajectory(project.m, color_by = "epitype", cell_size=0.4)
	dev.off()

	pdf("plot_cell_trajectory_celltype_project.m.pdf", width=5, height=5)
	# png("plot_cell_trajectory_byPseudotime_project.m.png")
	plot_cell_trajectory(project.m, color_by = "celltype", cell_size=0.05, cell_name_size = 8)
	dev.off()

	pdf("plot_cell_trajectory_fig.celltype_project.m.pdf", width=5, height=5)
	# png("plot_cell_trajectory_byPseudotime_project.m.png")
	plot_cell_trajectory(project.m, color_by = "fig.celltype", cell_size=0.05, cell_name_size = 8)
	dev.off()

	pdf("plot_cell_trajectory_fig.zone_project.m.pdf", width=5, height=5)
	# png("plot_cell_trajectory_byPseudotime_project.m.png")
	plot_cell_trajectory(project.m, color_by = "fig.zone", cell_size=0.05, cell_name_size = 8)
	dev.off()

	pdf("plot_cell_trajectory_pheno_project.m.pdf", width=5, height=5)
	# png("plot_cell_trajectory_byPseudotime_project.m.png")
	plot_cell_trajectory(project.m, color_by = "pheno", cell_size=0.1, cell_name_size = 8)
	dev.off()

	temp.m.time=project.m@phenoData@data$Pseudotime
	names(temp.m.time)=rownames(project.m@phenoData@data)
	temp@meta.data[names(temp.m.time), "pseudotime"]=temp.m.time


}

############################################################
#不同样本的peak 分析，并找到最近的基因
# MI 比CTR； MI比NMI；NMI比CTR

# find differential peaks
{
	DefaultAssay(project) <- 'peaks'
	project <- SetIdent(project, value = "sample.id")
	enz.adt_peaks <- FindMarkers(
	  object = project,
	  ident.1 = c("MI"), 
	  ident.2 = c("NMI"),
	  min.pct = 0.05,
	  test.use = 'LR',
	  max.cells.per.ident=1000
	)
	head(enz.adt_peaks)
	saveRDS(enz.adt_peaks, "enz.adt_peaks.rds")
	# enz.adt_peaks=readRDS("enz.adt_peaks.rds")

	enz.adt_peaks=readRDS("enz.adt_peaks.rds")
	open_enz <- rownames(enz.adt_peaks[enz.adt_peaks$avg_log2FC > 0.25, ])
	open_adt <- rownames(enz.adt_peaks[enz.adt_peaks$avg_log2FC < -0.25, ])
	closest_enz <- ClosestFeature(project, open_enz)
	closest_adt <- ClosestFeature(project, open_adt)
	head(closest_enz)
}

###################################################
#绘制特定基因的图，展示不同样本中基因的表达情况

# plot peaks across clusters
{
	check.genes = c('HMGA2','SNAI2',"VIM","IGFBP5","ZEB1","PDGFRA")
	DefaultAssay(project) <- "peaks"
	project <- SetIdent(project, value = "sample.id")
	pdf("temp.coverage.plot.specified.markers.between.samples.pdf", height=70, width=10)
	temp=CoveragePlot(
	  object = project,
	  region = check.genes,
	  extend.upstream = 1000,
	  extend.downstream = 1000,
	  ncol = 1
	)
	print(temp)
	dev.off()
}

###################################################################
#计算每个peaks的GC含量，并与基因关联起来
### link peaks to genes
DefaultAssay(project) <- "peaks"

# first compute the GC content for each peak
project <- RegionStats(project, genome = BSgenome.Mmusculus.UCSC.mm10)

#genes.use = c("Trp63", "Ar", "Krt8")
genes.use = c('HMGA2','SNAI2',"VIM","IGFBP5","ZEB1","PDGFRA")

# link peaks to genes
project <- LinkPeaks(
  object = project,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = genes.use
)

# first compute the GC content for each peak
project <- RegionStats(project, genome = BSgenome.Mmusculus.UCSC.mm10)


#idents.plot <- c("B naive", "B intermediate", "B memory","CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive")
indents.plot <- c("CD4 Naive", "CD4 Memory", "CD8 Naive", "CD8 effector", "Double negative T cell", "NK dim", "pre-B cell", "B cell progenitor", "pDC", "Dendritic cell", "CD14+ Monocytes", "CD16+ Monocytes")

p1 <- CoveragePlot(
  object = project,
  region = "HMGA2",
  features = "HMGA2",
  expression.assay = "RNA",
  idents = project@meta.data$fig.cluster,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = project,
  region = "SNAI2",
  features = "SNAI2",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)

###############################################
#检查特定peaks的表达

### check peaks expression
# chr [1:234776] "1-3360935-3361229" "1-3398732-3399453" "1-3416119-3416839"
markers.to.plot <- c("IGFBP5")
pdf("./result/temp.featureplot.peaks.pdf")
plot2 <- FeaturePlot(
  object = project,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)
dev.off()


#################################################################################################
#使用tradeSeq包识别沿拟时间表达变化的基因和峰值
#使用Homer进行motif分析

### find TF along pseudotime
### refer to: "FOXA2 drives lineage plasticity and KIT pathway activation in neuroendocrine prostate cancer". Cancer Cell 2022
### section "Identification of lineage-pioneering TFs" in "METHOD DETAILS"
### tradeSeq tutorial (identify RNA or peaks along pseudotime) ---  https://www.bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/Monocle.html
### tradeSeq main tutorial --- https://bioconductor.org/packages/devel/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html
### Homer (identify motif) --- http://homer.ucsd.edu/homer/  
### http://homer.ucsd.edu/homer/introduction/install.html   ("Anaconda", "Installing the basic HOMER software", "Downloading Homer Packages")
### http://homer.ucsd.edu/homer/introduction/basics.html
### https://fonseca.lab.mcgill.ca/resources/20220810_ATAC_Analysis_Demo/guide.html
### "Motif analysis" --- findMotifsGenome.pl diff_peaks.txt mm10 ~/students/<your_name>/motifs -size given -nomotif
{
	library("tradeSeq")
	project.m=readRDS("Q_PC13_monocle0523.rds")
	project=readRDS("project.visualized.harmony_Q_PC13.rds")

	# check monocle object's PCA map
	pdf("temp.project.trajectory.sample.id.pdf")
  plot_cell_trajectory(project.m, color_by = "sample.id", cell_size=0.5)
  dev.off()

  pdf("temp.project.trajectory.State.pdf")
  plot_cell_trajectory(project.m, color_by = "State", cell_size=0.5)
  dev.off()
  ### considering cells only in State 3 may make sense

  ### identify gene with expression changed along pseudotime, only consider variable genes (or fitGAM can't take the too-big expression matrix)
	sce <- fitGAM(counts = project@assays$RNA@counts[project@assays$RNA@var.features, ],
              pseudotime = project.m@phenoData@data$Pseudotime,
              cellWeights = rep(1, ncol(project)) )
	saveRDS(sce, "sce.rds")

	project <- FindVariableFeatures(project, assay="peaks", nfeatures = 5000)
	sce.peaks <- fitGAM(counts = project@assays$peaks@counts[project@assays$peaks@var.features, ],
              pseudotime = project.m@phenoData@data$Pseudotime,
              cellWeights = rep(1, ncol(project)) )
	saveRDS(sce.peaks, "sce.peaks.rds")
	
	### https://statomics.github.io/tradeSeq/articles/tradeSeq.html
	### section "Association of gene expression with pseudotime"
	assoRes <- associationTest(sce)
	head(assoRes)
	str(assoRes)
	print("Trp63"%in%rownames(assoRes)[assoRes$pvalue<0.05])

	assoRes.peaks <- associationTest(sce.peaks)
	head(assoRes.peaks)
	#                        waldStat df       pvalue meanLogFC
	# 17-39842996-39848785 1391.91109  6 0.000000e+00 0.4273285
	# 16-9751716-9751926           NA NA           NA 2.1240975
	str(assoRes.peaks)
	print("x"%in%rownames(assoRes.peaks)[assoRes.peaks$pvalue<0.05])
	table(assoRes.peaks$pvalue<0.05)


	### http://homer.ucsd.edu/homer/ngs/peakMotifs.html
	### "Acceptable Input files"
	### make dynamic peaks information into format as input for Homer
	# HOMER peak files should have at minimum 5 columns (separated by TABs, additional columns will be ignored):
	# Column1: Unique Peak ID
	# Column2: chromosome
	# Column3: starting position
	# Column4: ending position
	# Column5: Strand (+/- or 0/1, where 0="+", 1="-")
	######################################
	# 将动态peaks信息转换为homer的输入
	peaks.inf=rownames(assoRes.peaks)
	peaks.inf=data.frame(id=peaks.inf)
	temp=strsplit(peaks.inf$id, split="-")
	peaks.inf$chromosome=paste("chr",sapply(temp, `[`, 1), sep="")
	peaks.inf$start=sapply(temp, `[`, 2)
	peaks.inf$end=sapply(temp, `[`, 3)
	peaks.inf$strand=0
	write.table(peaks.inf, "peaks.inf.txt", quote = FALSE, row.names=FALSE, col.names=FALSE, sep = "\t")

	### Run Homer to find enriched motifs
	### in linux, type: 
	#/datapool/yanzeqin/software/homer
	PATH=$PATH:/home/dingqiuxia/gjc/project.scrna.atac.20230501/homer
	PATH=$PATH:/home/dingqiuxia/gjc/project.scrna.atac.20230501/homer/bin
	cd /home/dingqiuxia/gjc/project.scrna.atac.20230501/homer
	findMotifsGenome.pl ../peaks.inf.txt ./data/genomes/mm10 ./outs -size 200  -nomotif 
	### if meeting  an ERROR: cannot access mm10 preparsed cgbins
	### https://www.biostars.org/p/67216/
	### "Hi, I met the same problem and I solve it by this: findMotifsGenome.pl my.bed /home/Software/Homer/data/genomes/hg19 <output directory=""> -preparsedDir"
	findMotifsGenome.pl ../peaks.inf.txt ./data/genomes/mm10 ./outs -size 200 -preparsedDir ./data/genomes/mm10/preparsed  -nomotif 
	### if meeting ERROR: !!Could not open file for 1 (.fa or .fa.masked)
	### https://stackoverflow.com/questions/30162494/homer-de-novo-motif-discovery-cannot-open-hg19-fasta-files
	### "One thing to check: if the chromosome naming in your bed file and consistent with the chrom naming in the genome your are using: for example you shouldnt have '12' for chromosome 12 in your bed file whereas in the genome of your interest it is 'chr12'"

	### after running findMotifsGenome.pl, ./outs includes the results
	### output format: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
	### section "findMotifsGenome.pl Output"

	### read Homer's prediction result
	temp=read.csv("/home/dingqiuxia/gjc/project.scrna.atac.20230501/homer/outs/knownResults.txt", sep="\t")
	temp$"TF"=sapply(strsplit(temp$Motif.Name, split="\\("), `[`, 1)
	temp$"NegLogPvalue"=-1*log2(temp$P.value)

	### check genes overlapped between Homer result and dynamic genes
	temp=temp[temp$TF%in%rownames(assoRes)[assoRes$pvalue<0.05], ] 
	temp=temp[order(temp$NegLogPvalue, decreasing=TRUE), ]
	temp=temp[, c("TF","NegLogPvalue")]
	temp$TF=factor(temp$TF, levels=temp$TF)

	# https://ggplot2.tidyverse.org/reference/annotate.html
	pdf("lineage.pioneering.TFs.pdf")
	# toprint=plot(temp)
	toprint=ggplot(temp, aes(x = TF, y = NegLogPvalue)) + geom_point(size=3, col = "red") + theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
	# + annotate("text", x = 4, y = 25, label = "Some text")
	print(toprint)
	dev.off()

}

#########################################################################
#使用AddModuleScore函数计算特定基因集的路径激活得分，并在UMAP上进行可视化
### plot pathway activation score on umap
{
	str(project@meta.data)
	str(project@assays$RNA@scale.data)

	# check pathway activity (by addModuleScore function)
	DefaultAssay(project)="RNA"
	project <- ScaleData(project, features = rownames(project), assay="RNA" )
	filenames=c("project_high200", "squamous.signature", "FGFR_SIGNATURE", "AR_signature", "stem_signature")
	for(i in filenames)
	{
		temp=paste(i, "xlsx", sep=".")
		signature<- readxl::read_xlsx(temp)
		gene <- signature$GENE
		gene=gene[gene%in%rownames(project)]
		temp<- AddModuleScore(
		  object = project,
		  features = gene,
		  ctrl = 100, #默认值是100
		  name = i
		)
		project@meta.data[, i]=temp@meta.data[, paste(i,1,sep="")]
		pdf(paste(i, "pdf", sep="."), width=10, height=4)
		temp=VlnPlot(project,features =i, pt.size=0)
		print(temp)
		dev.off()

	}
	project <- ScaleData(project, features = project@assays$RNA@var.features, assay="RNA" )

	pdf("project.featureplot.pathways.pdf", height=10, width=8)
	FeaturePlot(project, features = filenames, raster=TRUE, pt.size=4, raster.dpi = c(1024, 1024))
	dev.off()

}

############################################################################
### 查找跨簇的差异表达峰值，并绘制热图展示这些峰值的表达情况。
### find differentially expressed peaks across clusters, and plot heatmap
{
	project=readRDS("F1_2sample_ARC_project.rds")
	### find differentially expressed peaks
	DefaultAssay(project) <- 'peaks'
	project <- SetIdent(project, value = "new_cluster")
	dpeaks <- FindAllMarkers(
	  object = project,
	  min.pct = 0.1,
	  test.use = 'LR',
	  max.cells.per.ident=500
	)
	head(dpeaks)
	# saveRDS(dpeaks, "project.markers.peaks.by.new_cluster.rds")
	# dpeaks=readRDS("project.markers.peaks.by.new_cluster.rds")

	# https://stuartlab.org/signac/articles/mouse_brain_vignette.html
  dpeaks=dpeaks[dpeaks$p_val_adj<0.05, ]
  dpeaks=dpeaks[order(dpeaks$avg_log2FC, decreasing = TRUE), ]
  dpeaks=dpeaks[dpeaks$avg_log2FC>log2(4), ] # 
	dpeaks <- dpeaks$gene
	str(dpeaks)

	### summary each gene's peaks, add to project@assays$genepeaks
	str(project@assays$peaks@data[dpeaks, ])
	peaks.name <- ClosestFeature(project, rownames(project@assays$peaks@data[dpeaks, ]))
	peaks.by.gene=aggregate(project@assays$peaks@data[dpeaks, ], by = list(peaks.name$gene_name), FUN = mean, na.rm = TRUE)
	rownames(peaks.by.gene)=peaks.by.gene$Group.1
	peaks.by.gene=peaks.by.gene[, colnames(peaks.by.gene)!="Group.1"]
	peaks.by.gene=as.matrix(peaks.by.gene)

	project[["genepeaks"]] <- CreateAssayObject(counts = peaks.by.gene )
	project <- ScaleData(project, features = rownames(project@assays$genepeaks), assay="genepeaks" )

	### find differentially expressed gene-level peaks
	DefaultAssay(project) <- 'genepeaks'
	project <- SetIdent(project, value = "new_cluster")
	dpeaks <- FindAllMarkers(
	  object = project,
	  min.pct = 0.1,
	  test.use = 'LR',
	  max.cells.per.ident=500
	)
	head(dpeaks)
	# saveRDS(dpeaks, "project.markers.genepeaks.by.new_cluster.rds")
	# dpeaks=readRDS("project.markers.genepeaks.by.new_cluster.rds")

  dpeaks=dpeaks[dpeaks$p_val_adj<0.05, ]
  # dpeaks=dpeaks[order(dpeaks$avg_log2FC, decreasing = TRUE), ]
  # dpeaks=dpeaks[dpeaks$avg_log2FC>log2(1.5), ] # 
  str(dpeaks)
  table(dpeaks$cluster)

  # plot top 10 gene-level peaks for each cluster
  library(dplyr)
  top10 <- dpeaks %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
  temp <- SetIdent(temp, value = "new_cluster")
  temp=subset(project, downsample=300)
  tocheck=top10$gene
  pdf(paste("heatmap.project.genepeaks.pdf", sep=""),  width=(length(unique(tocheck))/4+5), height=length(tocheck)/8)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="genepeaks", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
	
}



###########################################################
### 将簇的顺序进行调整，以便在CoveragePlot中更好地展示不同簇之间的比较
### change cluster's order in CoveragePlot
{
	table(project@meta.data$new_cluster)
	project@meta.data$new_cluster=as.character(project@meta.data$new_cluster)
	project@meta.data$new_cluster=factor(project@meta.data$new_cluster, levels=c("3", "1", "2", "0"))
	table(project@meta.data$new_cluster)
}
### 绘制特定基因在不同簇中的覆盖图，以展示不同簇中基因的表达情况。
### plot peaks across clusters
{
	library(ggplot2)
	library(patchwork)
	require(RColorBrewer)

	check.genes = c("Trp63", "Ar", "Krt8")
	DefaultAssay(project) <- "peaks"
	project <- SetIdent(project, value = "new_cluster")
	pdf("temp.coverage.plot.specified.markers.between.samples.pdf", height=70, width=10)
	cell_type_cols <- c(brewer.pal(9, "Set1")[c(1,8,3,5)])  
	temp=CoveragePlot(
	  object = project,
	  region = check.genes,
	  extend.upstream = 1000,
	  extend.downstream = 1000,
	  ncol = 1
	)
	temp=temp & scale_fill_manual(values = cell_type_cols)
	print(temp)
	dev.off()

}
### 绘制不同簇在拟时间上的密度曲线，以展示不同簇在拟时间上的分布情况。
### plot density curve of pseudotime for clusters
{
	project.m=readRDS("F1_2Sample_4_monocle.rds")

	# give monocle's pseudotime to seurat object
	project.m.time=project.m@phenoData@data$Pseudotime
	names(project.m.time)=rownames(project.m@phenoData@data)
	project@meta.data[names(project.m.time), "pseudotime"]=project.m.time

	temp.df=project@meta.data[, c("new_cluster", "pseudotime")]
	# get colors
	require("scales")
	cell_type_cols <- c(brewer.pal(9, "Set1")[c(1,8,3,5)])  
	# color=color[length(oorder):1]# 
	# specify Clike's colour
	# color[length(oorder)]="#FF6600"
	pdf("CellTypeDensityAlongPseudotime.pdf", width=10, height=5)
	# transparancy desnity plot. alpha: transparent degree
	# {
	#   p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=fig.celltype, fill=fig.celltype)) +
	#     geom_density(adjust=1.5, alpha=.4)
	# }
	# stacked density plot
	# {
	  # p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=fig.celltype, fill=fig.celltype)) +
	  # geom_density(adjust=1.5, position="fill")
	# }
	{
	  p1<- ggplot(temp.df, aes(x=pseudotime, color=new_cluster)) +
	  geom_density() +
	  labs(title="Cell type density curve",x="Pseudotime", y = "Density") +
	  scale_color_manual(values=cell_type_cols)+
	  theme_classic()
	}
	print(p1)
	dev.off()
}

