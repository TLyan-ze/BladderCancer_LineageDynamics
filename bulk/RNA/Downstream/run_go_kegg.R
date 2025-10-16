#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(egg)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "CSV/TSV with differential results; must include gene symbols and p/adjusted p and log2FC."),
  make_option(c("--symbol-col"), type = "character", default = "SYMBOL", help = "Column name for gene symbols (default SYMBOL)."),
  make_option(c("--pcol"), type = "character", default = "pvalue", help = "Column name for p-value (default pvalue)."),
  make_option(c("--lfc-col"), type = "character", default = "log2FC", help = "Column name for log2 fold change (default log2FC)."),
  make_option(c("--p"), type = "double", default = 0.05, help = "Significance threshold for DE gene selection (default 0.05)."),
  make_option(c("--lfc"), type = "double", default = 1.0, help = "Absolute log2FC threshold (default 1)."),
  make_option(c("--organism"), type = "character", default = "hs", help = "Organism: hs (human) or mm (mouse)."),
  make_option(c("-o", "--out-prefix"), type = "character", default = "enrich", help = "Output prefix for results and plots."),
  make_option(c("--top"), type = "integer", default = 20, help = "Top N GO BP terms to plot (default 20).")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.na(opt$input)) stop("--input is required")

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt")) readr::read_tsv(path) else readr::read_csv(path)
}

x <- read_table_auto(opt$input)
if (!(opt$`symbol-col` %in% colnames(x))) stop("Symbol column not found: ", opt$`symbol-col`)
if (!(opt$`pcol` %in% colnames(x))) stop("p-value column not found: ", opt$`pcol`)
if (!(opt$`lfc-col` %in% colnames(x))) stop("log2FC column not found: ", opt$`lfc-col`)

# Select DE genes
x <- x %>% dplyr::filter(!is.na(.data[[opt$`symbol-col`]]))
diff_gene <- x %>% dplyr::filter(!is.na(.data[[opt$`pcol`]]), !is.na(.data[[opt$`lfc-col`]])) %>%
  dplyr::filter(.data[[opt$`pcol`]] < opt$p & abs(.data[[opt$`lfc-col`]]) > opt$lfc)

if (nrow(diff_gene) == 0) stop("No genes pass thresholds; adjust --p/--lfc")

orgDb <- switch(opt$organism,
  hs = org.Hs.eg.db,
  mm = org.Mm.eg.db,
  stop("Unsupported organism: ", opt$organism)
)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
})

# Map SYMBOL to ENTREZID
gene.diff_name <- bitr(diff_gene[[opt$`symbol-col`]], fromType = "SYMBOL",
                       toType = c("SYMBOL","ENTREZID"), OrgDb = orgDb)

if (!("ENTREZID" %in% colnames(gene.diff_name))) stop("ENTREZID mapping failed")

gene <- gene.diff_name$ENTREZID

# GO enrichment BP/CC/MF
ego_BP <- enrichGO(gene = gene, OrgDb = orgDb, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ego_CC <- enrichGO(gene = gene, OrgDb = orgDb, keyType = "ENTREZID",
                   ont = "CC", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ego_MF <- enrichGO(gene = gene, OrgDb = orgDb, keyType = "ENTREZID",
                   ont = "MF", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)

# Save tables
readr::write_csv(as.data.frame(ego_BP), paste0(opt$`out-prefix`, "_ego_BP.csv"))
readr::write_csv(as.data.frame(ego_CC), paste0(opt$`out-prefix`, "_ego_CC.csv"))
readr::write_csv(as.data.frame(ego_MF), paste0(opt$`out-prefix`, "_ego_MF.csv"))

# Plot top BP terms
ego_result_BP <- as.data.frame(ego_BP)
ego_result_BP <- ego_result_BP[order(ego_result_BP$Count, decreasing = TRUE), ]
if (nrow(ego_result_BP) > opt$top) ego_result_BP <- ego_result_BP[1:opt$top, ]
rownames(ego_result_BP) <- 1:nrow(ego_result_BP)
ego_result_BP$order <- factor(rev(as.integer(rownames(ego_result_BP))), labels = rev(ego_result_BP$Description))
ego_result_BP$p.adjust <- as.numeric(ego_result_BP$p.adjust)

p <- ggplot(ego_result_BP, aes(x = reorder(order, Count), y = Count, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradient(low = "#fee6ce", high = "#e6550d") +
  labs(title = "GO BP Enrichment", y = "Gene Count", x = "", fill = "-log10(p.adjust)") +
  theme_test() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  scale_y_continuous(limits = c(0, max(ego_result_BP$Count) + 1), expand = c(0, 0)) +
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        plot.title = element_text(size = 16, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text.x = element_text(angle = 0, size = 10, colour = 'black'),
        axis.text.y = element_text(angle = 0, size = 10, colour = 'black'))

plot_out <- paste0(opt$`out-prefix`, "_go_bp.pdf")
# Optional: p_modified <- egg::set_panel_size(p, width = unit(8, "in"), height = unit(7, "in"))
ggsave(filename = plot_out, plot = p, width = 10, height = 8, units = 'in', dpi = 300)

# KEGG enrichment
org <- if (opt$organism == "hs") "human" else "mouse"
kk <- enrichKEGG(gene = gene, keyType = "kegg", organism = org, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
readr::write_csv(as.data.frame(kk), paste0(opt$`out-prefix`, "_kegg.csv"))

hh <- as.data.frame(kk)
if (nrow(hh) > 0) {
  rownames(hh) <- 1:nrow(hh)
  hh$order <- factor(rev(as.integer(rownames(hh))), labels = rev(hh$Description))
  p2 <- ggplot(hh, aes(x = reorder(order, Count), y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_gradient(low = "#fee6ce", high = "#e6550d") +
    labs(title = "KEGG Pathways Enrichment", y = "Gene numbers", x = "") +
    theme_test() +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_y_continuous(limits = c(0, max(hh$Count) + 1), expand = c(0, 0)) +
    theme(legend.title = element_text(color = "black", size = 12),
          legend.text = element_text(color = "black", size = 12),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          plot.title = element_text(size = 16, colour = 'black'),
          axis.title = element_text(size = 14, colour = 'black'),
          axis.text.x = element_text(angle = 0, size = 10, colour = 'black'),
          axis.text.y = element_text(angle = 0, size = 10, colour = 'black'))
  kegg_out <- paste0(opt$`out-prefix`, "_kegg.pdf")
  ggsave(filename = kegg_out, plot = p2, width = 12, height = 10, units = 'in', dpi = 300)
}

message("[INFO] Wrote enrichment outputs with prefix ", opt$`out-prefix`)
