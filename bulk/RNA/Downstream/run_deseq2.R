#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

option_list <- list(
  make_option(c("-c", "--counts"), type = "character", help = "Path to counts matrix (TSV/CSV). Rows=genes, cols=samples."),
  make_option(c("-m", "--metadata"), type = "character", help = "Sample metadata TSV/CSV with columns: sample, condition."),
  make_option(c("--gene-col"), type = "character", default = NA, help = "Optional gene column name in counts file. If absent, rownames are used, or parsed from column 'X'."),
  make_option(c("--filter-min-mean"), type = "double", default = 1.0, help = "Filter genes with rowMeans <= threshold (default 1)."),
  make_option(c("--lfc"), type = "double", default = 1.0, help = "Absolute log2 fold change threshold (default 1)."),
  make_option(c("--p"), type = "double", default = 0.05, help = "P-value or adjusted P-value threshold (default 0.05)."),
  make_option(c("--use-padj"), action = "store_true", default = TRUE, help = "Use adjusted p-value (padj) for significance (default TRUE)."),
  make_option(c("--fitType"), type = "character", default = "mean", help = "DESeq fitType (parametric, local, mean). Default 'mean' to match legacy script."),
  make_option(c("-o", "--out-prefix"), type = "character", default = "deseq2", help = "Output prefix for results files and plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$counts) || is.na(opt$metadata)) {
  stop("--counts and --metadata are required.")
}

read_table_auto <- function(path) {
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("tsv", "txt")) readr::read_tsv(path) else readr::read_csv(path)
}

counts_raw <- read_table_auto(opt$counts)
md <- read_table_auto(opt$metadata)
colnames(md) <- tolower(colnames(md))
if (!all(c("sample", "condition") %in% colnames(md))) stop("Metadata must have columns: sample, condition")

# Prepare count matrix
if (!is.na(opt$`gene-col`) && opt$`gene-col` %in% colnames(counts_raw)) {
  gene_col <- opt$`gene-col`
  rownames(counts_raw) <- counts_raw[[gene_col]]
  counts_raw[[gene_col]] <- NULL
} else if ("X" %in% colnames(counts_raw)) {
  # Legacy format: parse gene name from 'X' split by '|'
  gn <- sapply(strsplit(counts_raw$X, "\\|"), function(x) if (length(x) >= 3) x[length(x)-2] else x[1])
  rownames(counts_raw) <- gn
  counts_raw$X <- NULL
} else if (!is.null(rownames(counts_raw)) && all(rownames(counts_raw) != "")) {
  # rownames already set
} else {
  stop("Cannot determine gene identifiers. Provide --gene-col or ensure rownames are genes.")
}

# Ensure samples match
samples <- colnames(counts_raw)
md <- md[match(samples, md$sample), , drop = FALSE]
if (any(is.na(md$sample))) stop("Metadata missing entries for some samples in counts matrix.")

# Filter and round counts
counts <- as.data.frame(counts_raw)
counts <- counts[rowMeans(counts) > opt$`filter-min-mean`, , drop = FALSE]
counts_int <- round(as.matrix(counts))

condition <- factor(md$condition)
colData <- data.frame(row.names = colnames(counts_int), condition = condition)

dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = colData, design = ~ condition)
dds <- DESeq(dds, fitType = opt$fitType, parallel = FALSE)
res <- results(dds)
res_df <- as.data.frame(res)

# Order and write all results
res_df <- res_df[order(res_df$pvalue, res_df$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
all_out <- paste0(opt$`out-prefix`, "_DEseq2_all.csv")
readr::write_csv(res_df, all_out)

# Significant genes
sig_metric <- if (opt$`use-padj`) "padj" else "pvalue"
padj_vec <- res_df[[sig_metric]]
sig_df <- res_df[!is.na(padj_vec) & (abs(res_df$log2FoldChange) >= opt$lfc) & (padj_vec < opt$p), ]
sig_out <- paste0(opt$`out-prefix`, "_DEseq2_sig.csv")
readr::write_csv(sig_df, sig_out)

# Volcano plot
res_df$color <- ifelse(!is.na(padj_vec) & (padj_vec < opt$p) & abs(res_df$log2FoldChange) >= opt$lfc,
                       ifelse(res_df$log2FoldChange > 0, 'red', 'blue'), 'gray')
color_map <- c(red = "red", gray = "gray", blue = "blue")
p <- ggplot(res_df, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  scale_color_manual(values = color_map) +
  labs(x = "log2 (fold change)", y = "-log10 (adj p-value)") +
  geom_hline(yintercept = -log10(opt$p), lty = 4, col = "grey", lwd = 0.6) +
  geom_vline(xintercept = c(-opt$lfc, opt$lfc), lty = 4, col = "grey", lwd = 0.6) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
plot_out <- paste0(opt$`out-prefix`, "_volcano.png")
ggsave(plot_out, p, width = 7, height = 5, dpi = 300)

message("[INFO] Wrote:")
message("  ", all_out)
message("  ", sig_out)
message("  ", plot_out)