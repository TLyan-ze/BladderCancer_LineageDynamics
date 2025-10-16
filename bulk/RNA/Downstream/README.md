# Bulk RNA Downstream Analysis

This directory contains downstream analysis scripts for bulk RNA.

Contents:
- `run_deseq2.R`: Parameterized DESeq2 differential expression with volcano plot.
- `run_DEseq.R`: Legacy DESeq2 workflow tied to `GEMM mus_counts.tsv` format.
- `KEGG.R`: GO and KEGG enrichment visualization from a CSV of DE results.

## Dependencies
R packages: `DESeq2`, `optparse`, `dplyr`, `ggplot2`, `ggrepel`, `readr`, `clusterProfiler`, `AnnotationDbi`, `org.Hs.eg.db` or `org.Mm.eg.db` (choose organism), `enrichplot` (optional).

Install example:
```r
install.packages(c("optparse","ggplot2","ggrepel","readr","dplyr"))
BiocManager::install(c("DESeq2","clusterProfiler","AnnotationDbi","org.Hs.eg.db","org.Mm.eg.db","enrichplot"))
```

## Usage
### 1) DESeq2 (generalized)
Prepare counts matrix (`counts.tsv`) with rows=genes and columns=samples, and metadata (`samples.tsv`) with columns: `sample`, `condition` matching counts columns.

```bash
Rscript run_deseq2.R \
  --counts counts.tsv \
  --metadata samples.tsv \
  --filter-min-mean 1 \
  --lfc 1 \
  --p 0.05 \
  --out-prefix my_cohort
```
Outputs:
- `my_cohort_DEseq2_all.csv`: all results
- `my_cohort_DEseq2_sig.csv`: significant genes
- `my_cohort_volcano.png`: volcano plot

If your counts have a gene column (e.g. `Gene_Name`), add `--gene-col Gene_Name`. Legacy format with an `X` column will be auto-parsed.

### 2) GO/KEGG enrichment
Use `KEGG.R` by editing input lines to point to your DE CSV (e.g. `MI_IGFKO_2.csv`) and organism. It writes `ego_*.csv` and `kegg.pdf` plots.

Recommended: a future wrapper `run_go_kegg.R` will accept `--input`, `--organism`, and output prefix to avoid manual edits.

### 3) Legacy script
`run_DEseq.R` demonstrates a fixed-format pipeline including GO enrichment for mouse (`org.Mm.eg.db`). Use as reference; prefer `run_deseq2.R` for reproducibility.

## Tips
- Ensure sample order in metadata matches counts column order.
- Filter low-expression genes with `--filter-min-mean` to stabilize dispersion.
- Set `--use-padj` to use adjusted `padj` (default TRUE).
- For human vs mouse, load `org.Hs.eg.db` vs `org.Mm.eg.db` in enrichment.
