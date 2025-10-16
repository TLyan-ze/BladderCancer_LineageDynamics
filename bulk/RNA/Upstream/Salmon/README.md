# Salmon 上游流程（Bulk RNA）

本目录用于 bulk RNA 的上游定量（Salmon）。原 `sc_multiomics/scRNA/Salmon` 的内容属于 bulkRNA 上游；scRNA 的分析流程与结果位于谱系相关目录（例如 `Script/AP-data/expression`），请勿混淆。

- 目录位置：`bulk/RNA/Upstream/Salmon`
- 包含文件：`run_salmon.sh`、`get_salmon.R`、`samples_info`、`samples_list`、`sample1`、`sample2`、`sample_finish`、`readme`、`result/`

## 依赖
- Salmon 二进制（示例路径：`/share/Data01/yanzeqin/software/snakemake/conda/envs/salmon/bin/salmon`）
- R 包：`tximport`、`tidyverse`、`data.table`
- 转录本到基因映射 CSV（任选其一并在 `get_salmon.R` 中设置）：
  - 人：`/share/database/openData/GRCh38_GENCODE/gencode.v35.annotation.TxToGene.csv`
  - 鼠：`/datapool/yanzeqin/database/GRCm39_GENCODE/gencode.mus.TxToGene.csv`

## 输入准备
- `samples_info` 四列：`sample_id`、`replicate`、`fastq1`、`fastq2`
- 参考示例（存于 `readme`）：
  - `cat sample1 | while read l1; do fq1=\`ls /datapool/yanzeqin/project/ZDW/RNA/MI-hIGFBP5/00.CleanData/$l1/${l1}_1.clean.fq.gz\`; fq2=\`ls /datapool/yanzeqin/project/ZDW/RNA/MI-hIGFBP5/00.CleanData/$l1/${l1}_2.clean.fq.gz\`; echo -e "$l1\t$l1\t$fq1\t$fq2"; done`
  - `cat sample2 | while read l1; do fq1=\`ls /datapool/yanzeqin/project/ZDW/RNA/MI-2nd/00.CleanData/$l1/${l1}_1.clean.fq.gz\`; fq2=\`ls /datapool/yanzeqin/project/ZDW/RNA/MI-2nd/00.CleanData/$l1/${l1}_2.clean.fq.gz\`; echo -e "$l1\t$l1\t$fq1\t$fq2"; done`

## 运行 Salmon 定量
1. 编辑 `run_salmon.sh` 中的变量：
   - `input`：原始 FASTQ 所在目录
   - `output`：输出目录（示例：`${input}/03.salmon`）
   - `index`：Salmon 索引（人或鼠）
2. 运行：`bash run_salmon.sh`（脚本会为每个样本生成量化目录并通过 `qsub` 提交任务）
3. 完成后，在 `result/` 下可见各样本的 `quant.sf` 文件。

## 生成基因表达矩阵
- 在 Salmon 定量完成后运行：
  - `Rscript bulk/RNA/Upstream/Salmon/get_salmon.R bulk/RNA/Upstream/Salmon/result`
- 功能：递归查找 `result/` 下的 `*_quant.sf`，用 `tximport` 汇总为基因层面的表达矩阵。
- 输出：
  - `bulk/RNA/Upstream/Salmon/result/exp_counts.csv`（整数化的 counts）
  - `bulk/RNA/Upstream/Salmon/result/exp_tmp.csv`（TPM）

## 下游分析衔接
- 差异分析：`bulk/RNA/Downstream/run_deseq2.R` 或 `run_DEseq.R`
- 富集分析：`bulk/RNA/Downstream/run_go_kegg.R`（GO/KEGG）

## 约定与说明
- 统一将 Salmon 的输入/输出放在 `bulk/RNA/Upstream/Salmon/result/`
- `sc_multiomics/scRNA/Salmon` 目录仅作为旧位置标识，不再使用；scRNA 的流程与结果位于 `Script/AP-data/expression`。
- 如需在鼠数据上运行，请在 `get_salmon.R` 中切换映射 CSV 路径（鼠）。
