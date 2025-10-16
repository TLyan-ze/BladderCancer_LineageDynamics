#!/usr/bin/env Rscript
# Single-cell ATAC-seq Analysis Pipeline Template
# Based on bladder cancer scATAC-seq analysis workflow

# Load required libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(GenomicRanges)
library(future)
library(Matrix)
library(irlba)

# Optional: for genome annotation and motif analysis
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(EnsDb.Hsapiens.v86)
# library(JASPAR2020)
# library(TFBSTools)

# Set up parallel processing
options(future.globals.maxSize = 50000 * 1024^2)

# Function to run complete scATAC-seq analysis
run_scATAC_analysis <- function(sample_paths, sample_names, output_dir = ".", 
                               min_fragments = 500, max_peak_width = 10000, 
                               min_peak_width = 20, min_cutoff = 20) {
  
  cat("Starting scATAC-seq analysis pipeline...\n")
  
  # Step 1: Read peak sets from all samples
  cat("Step 1: Reading peak sets...\n")
  peak_list <- list()
  
  for (i in seq_along(sample_paths)) {
    sample_name <- sample_names[i]
    peak_file <- file.path(sample_paths[i], "outs/filtered_peak_bc_matrix/peaks.bed")
    
    if (file.exists(peak_file)) {
      peaks <- read.table(peak_file, col.names = c("chr", "start", "end"))
      peak_list[[sample_name]] <- makeGRangesFromDataFrame(peaks)
      cat(sprintf("  Loaded %d peaks for %s\n", length(peak_list[[sample_name]]), sample_name))
    } else {
      stop(sprintf("Peak file not found: %s", peak_file))
    }
  }
  
  # Step 2: Create unified peak set
  cat("Step 2: Creating unified peak set...\n")
  combined.peaks <- reduce(x = do.call(c, peak_list))
  
  # Filter peaks by width
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths < max_peak_width & peakwidths > min_peak_width]
  cat(sprintf("  Unified peak set: %d peaks\n", length(combined.peaks)))
  
  # Step 3: Load metadata and create fragment objects
  cat("Step 3: Loading metadata and creating fragment objects...\n")
  seurat_objects <- list()
  
  for (i in seq_along(sample_paths)) {
    sample_name <- sample_names[i]
    sample_path <- sample_paths[i]
    
    # Load metadata
    metadata_file <- file.path(sample_path, "outs/singlecell.csv")
    if (file.exists(metadata_file)) {
      md <- read.table(metadata_file, stringsAsFactors = FALSE, sep = ",", 
                      header = TRUE, row.names = 1)[-1, ]  # remove first row
      
      # Filter cells by fragment count
      md <- md[md$passed_filters > min_fragments, ]
      cat(sprintf("  %s: %d cells after filtering\n", sample_name, nrow(md)))
      
      # Create fragment object
      fragment_file <- file.path(sample_path, "outs/fragments.tsv.gz")
      if (file.exists(fragment_file)) {
        frags <- CreateFragmentObject(path = fragment_file, cells = rownames(md))
        
        # Quantify peaks
        counts <- FeatureMatrix(fragments = frags, features = combined.peaks, 
                               cells = rownames(md))
        
        # Create Seurat object
        assay <- CreateChromatinAssay(counts, fragments = frags)
        seurat_obj <- CreateSeuratObject(assay, assay = "ATAC", meta.data = md)
        seurat_obj$dataset <- sample_name
        
        seurat_objects[[sample_name]] <- seurat_obj
        
      } else {
        stop(sprintf("Fragment file not found: %s", fragment_file))
      }
    } else {
      stop(sprintf("Metadata file not found: %s", metadata_file))
    }
  }
  
  # Step 4: Merge objects
  cat("Step 4: Merging Seurat objects...\n")
  if (length(seurat_objects) > 1) {
    combined <- merge(x = seurat_objects[[1]], 
                     y = seurat_objects[-1], 
                     add.cell.ids = sample_names)
  } else {
    combined <- seurat_objects[[1]]
  }
  
  cat(sprintf("  Combined object: %d cells, %d features\n", 
              ncol(combined), nrow(combined)))
  
  # Step 5: Standard scATAC-seq processing
  cat("Step 5: Running standard scATAC-seq workflow...\n")
  
  # TF-IDF normalization
  combined <- RunTFIDF(combined)
  
  # Find top features
  combined <- FindTopFeatures(combined, min.cutoff = min_cutoff)
  
  # SVD for dimensionality reduction
  combined <- RunSVD(combined)
  
  # UMAP embedding
  combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
  
  # Step 6: Generate basic plots
  cat("Step 6: Generating visualization plots...\n")
  
  # UMAP plot colored by dataset
  p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1) +
    ggtitle("UMAP by Dataset")
  
  ggsave(file.path(output_dir, "umap_by_dataset.png"), p1, 
         width = 8, height = 6, dpi = 300)
  
  # Coverage plot (example region)
  if (length(sample_names) > 1) {
    p2 <- CoveragePlot(object = combined, group.by = 'dataset', 
                      region = "chr1-1000000-2000000")
    
    ggsave(file.path(output_dir, "coverage_plot_example.png"), p2, 
           width = 10, height = 6, dpi = 300)
  }
  
  # Step 7: Save results
  cat("Step 7: Saving results...\n")
  saveRDS(combined, file.path(output_dir, "scATAC_combined_object.rds"))
  
  # Save peak information
  peak_df <- as.data.frame(combined.peaks)
  write.table(peak_df, file.path(output_dir, "unified_peaks.bed"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("Analysis completed successfully!\n")
  return(combined)
}

# Function to add genome annotations (optional)
add_genome_annotations <- function(seurat_obj, genome = "hg38") {
  
  if (genome == "hg38") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(EnsDb.Hsapiens.v86)
    
    # Extract gene annotations
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- "hg38"
    
    # Add annotations to object
    Annotation(seurat_obj) <- annotations
    
  } else if (genome == "mm10") {
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(EnsDb.Mmusculus.v79)
    
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- "mm10"
    
    Annotation(seurat_obj) <- annotations
  }
  
  return(seurat_obj)
}

# Example usage
if (FALSE) {
  # Define sample paths and names
  sample_paths <- c(
    "/path/to/CTR_sample",
    "/path/to/NMI_sample", 
    "/path/to/MI_sample"
  )
  
  sample_names <- c("CTR", "NMI", "MI")
  
  # Run analysis
  combined_obj <- run_scATAC_analysis(
    sample_paths = sample_paths,
    sample_names = sample_names,
    output_dir = "scATAC_results",
    min_fragments = 500,
    max_peak_width = 10000,
    min_peak_width = 20,
    min_cutoff = 20
  )
  
  # Optional: Add genome annotations
  combined_obj <- add_genome_annotations(combined_obj, genome = "hg38")
  
  # Save final object
  saveRDS(combined_obj, "scATAC_final_object.rds")
}