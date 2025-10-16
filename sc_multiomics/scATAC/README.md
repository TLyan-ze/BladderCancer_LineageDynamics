# Single-cell ATAC-seq Analysis

This directory contains the single-cell ATAC-seq (scATAC-seq) analysis pipeline for chromatin accessibility profiling in bladder cancer samples.

## Overview

The scATAC-seq analysis pipeline processes chromatin accessibility data from three sample conditions:
- **CTR (Control)**: Control samples  
- **NMI (Non-Muscle Invasive)**: Non-muscle invasive bladder cancer
- **MI (Muscle Invasive)**: Muscle invasive bladder cancer

## Analysis Workflow

### Step 1: Data Preprocessing and Peak Calling (`step_1_scATAC.R`)

**Key Functions:**
- Read peak sets from Cell Ranger output
- Create unified peak set across all samples
- Filter peaks by length (20bp - 10kb)
- Create genomic ranges objects

### Step 2: Quality Control and Analysis (`step_2_scATAC.R`)

**Key Functions:**
- Add genome annotations (hg38)
- Compute quality control metrics
- Perform dimensionality reduction
- Cell type annotation
- Motif analysis

### Step 3: Integration and Visualization (`scATAC_seq.R`, `scATAC_seq_new.R`)

**Key Functions:**
- Merge multiple samples
- Create Fragment objects
- Quantify peaks across datasets
- Run TF-IDF normalization
- Perform SVD and UMAP
- Generate coverage plots

## Key Analysis Steps

### 1. Peak Processing
```r
# Read peak sets
peaks.CTR <- read.table("peaks.bed", col.names = c("chr", "start", "end"))

# Create unified peak set
combined.peaks <- reduce(x = c(gr.CTR, gr.NMI, gr.MI))

# Filter by peak width
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
```

### 2. Fragment Object Creation
```r
# Create fragment objects for each sample
frags.CTR <- CreateFragmentObject(path = "fragments.tsv.gz", cells = rownames(md.CTR))
```

### 3. Peak Quantification
```r
# Quantify peaks in each dataset
CTR.counts <- FeatureMatrix(fragments = frags.CTR, features = combined.peaks, cells = rownames(md.CTR))
```

### 4. Seurat Object Creation and Integration
```r
# Create Seurat objects
CTR_assay <- CreateChromatinAssay(CTR.counts, fragments = frags.CTR)
CTR <- CreateSeuratObject(CTR_assay, assay = "ATAC", meta.data = md.CTR)

# Merge datasets
combined <- merge(x = CTR, y = list(NMI, MI), add.cell.ids = c("CTR", "NMI", "MI"))
```

### 5. Dimensionality Reduction and Visualization
```r
# Standard scATAC-seq workflow
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

# Visualization
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
CoveragePlot(object = combined, group.by = 'dataset', region = "chr14-99700000-99760000")
```

## Key Scripts

### Main Analysis Scripts
- **`scATAC_seq.R`**: Main integration and analysis pipeline
- **`scATAC_seq_new.R`**: Updated version of main pipeline
- **`step_1_scATAC.R`**: Initial data processing and peak calling
- **`step_2_scATAC.R`**: Quality control and advanced analysis

### Utility Scripts
- **`functions.seurat.R`**: Cell type annotation functions
- **`scRNA.atac.1.R`**: scRNA-ATAC integration functions
- **`run_scATAC.sh`**: Batch processing script

## Dependencies

### R Packages
```r
library(Seurat)
library(Signac)
library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
```

## Sample Information

### Sample Types
- **MC (CTR)**: Control samples
- **NMI**: Non-muscle invasive bladder cancer
- **MI**: Muscle invasive bladder cancer

### Data Structure
Each sample contains:
- Peak coordinates (`peaks.bed`)
- Cell metadata (`singlecell.csv`)
- Fragment information (`fragments.tsv.gz`)
- Feature-barcode matrix (`filtered_feature_bc_matrix.h5`)

