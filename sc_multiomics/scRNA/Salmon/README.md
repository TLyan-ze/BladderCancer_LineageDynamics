# Salmon RNA-seq 上游流程整理

本目录整理了常规 RNA 的上游分析以生成表达矩阵（counts 与 TPM）的脚本与使用说明，包含批量运行 Salmon 和结果聚合。

## 目录与文件
- `run_salmon.sh`：批量提交 Salmon `quant` 任务（SGE 队列环境），默认人类索引。
- `get_salmon.R`：聚合 `quant.sf`，生成 `exp_counts.csv` 与 `exp_tmp.csv`。
- `prepare_samples_info.sh`：根据样本列表自动生成 `samples_info`（四列：样本ID、样本ID、FQ1、FQ2）。
- `run_get_salmon.sh`：封装调用 `get_salmon.R` 的脚本（无需改 R 文件参数）。
- `sample1`、`sample2`：两批样本名列表（每行一个样本名）。
- `samples_info`、`samples_list`：原始清单（可由 `prepare_samples_info.sh` 重新生成）。
- `readme`：原始记录文件。

## 快速开始
1) 准备样本清单与原始数据路径
- 将样本名分别写入 `sample1`、`sample2`（每行一个，例如 `MGHU3_S1`）。
- 原始数据目录遵循以下布局（与历史脚本一致）：
  - 批次1：`/datapool/yanzeqin/project/ZDW/RNA/MI-hIGFBP5/00.CleanData/<样本名>/<样本名>_1.clean.fq.gz`
            与同目录下的 `<样本名>_2.clean.fq.gz`
  - 批次2：`/datapool/yanzeqin/project/ZDW/RNA/MI-2nd/00.CleanData/<样本名>/<样本名>_1.clean.fq.gz`
            与同目录下的 `<样本名>_2.clean.fq.gz`
  - 如路径不同，可编辑 `prepare_samples_info.sh` 中 `BASE1`、`BASE2` 值。

2) 生成 `samples_info`
- 运行：
  - `bash prepare_samples_info.sh`
- 输出：在当前目录生成/更新 `samples_info`（四列，以制表符分隔）：
  - `样本ID\t样本ID\tFQ1绝对路径\tFQ2绝对路径`
- 注：该脚本是对如下手动命令的自动化封装（保持同样输出）：
  - `cat sample1 | while read l1; do fq1=...; fq2=...; echo -e "$l1\t$l1\t$fq1\t$fq2"; done`
  - `cat sample2 | while read l1; do fq1=...; fq2=...; echo -e "$l1\t$l1\t$fq1\t$fq2"; done`

3) 批量运行 Salmon quant
- 编辑 `run_salmon.sh` 中以下变量（若需要）：
  - `input=/datapool/yanzeqin/project/ZDW/RNA/merge/MI`（项目根目录）
  - `output=${input}/03.salmon`（Salmon输出目录）
  - `index=/share/database/openData/GRCh38_GENCODE/transcripts_index_salmon`（人类索引）
    - 小鼠可改为：`/datapool/yanzeqin/database/GRCm39_GENCODE/salmon`
- 运行：`bash run_salmon.sh`
- 说明：脚本会读取 `samples_info`，为每个样本提交 SGE 任务（`qsub`）。
- 完成后，`03.salmon/<样本名>/quant.sf` 应存在。

4) 聚合 `quant.sf` 为表达矩阵
- 运行：`bash run_get_salmon.sh -d /datapool/yanzeqin/project/ZDW/RNA/merge/MI/03.salmon`
- 该脚本会调用 `get_salmon.R` 并在指定目录生成：
  - `exp_counts.csv`（整数化的 counts）
  - `exp_tmp.csv`（TPM 矩阵；R脚本历史命名为 tmp）
- `get_salmon.R` 依赖映射文件：
  - 人类：`/share/database/openData/GRCh38_GENCODE/gencode.v35.annotation.TxToGene.csv`
  - 小鼠示例（注释行）：`/datapool/yanzeqin/database/GRCm39_GENCODE/gencode.mus.TxToGene.csv`

## 依赖与环境
- Salmon 可执行：`/share/Data01/yanzeqin/software/snakemake/conda/envs/salmon/bin/salmon`
- SGE 队列环境（`qsub`、`qstat`、`qdel`）
- R 包：`tximport`、`tidyverse`、`data.table`

## 本地运行（无队列）参考
- 若无 SGE，可参考以下循环（需设置好 `index`）：
```
while read -r sid sname fq1 fq2; do 
  outdir="/path/to/03.salmon/${sid}"; mkdir -p "$outdir"
  salmon quant -i "$index" --gcBias -l A --thread 8 \
    -1 "$fq1" -2 "$fq2" -o "$outdir"
done < samples_info
```

## 常见问题
- `quant.sf` 不存在：检查 Salmon 是否成功运行、输出目录是否正确。
- R 脚本报错找不到 `TxToGene`：确认映射文件路径存在且可读。
- `Eqw` 队列状态：脚本已包含 `qdel` 处理；仍异常时请联系集群管理员。

## 结果文件
- `03.salmon/<样本>/quant.sf`：Salmon 原始结果。
- `exp_counts.csv`：基因层面的 counts 表达矩阵（样本为列）。
- `exp_tmp.csv`：TPM 表达矩阵（样本为列）。
