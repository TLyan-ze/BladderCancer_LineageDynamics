# BladderCancer_LineageDynamics

本仓库用于整理并复现与膀胱癌谱系追踪相关的分析代码与流程，覆盖三大模块：
- 谱系分析（上游）：从 FASTQ 到 UMI/allele 表的处理与对齐
- 谱系分析（下游）：树重构、扩增/适应度分析、模块评分与图形化
- 多组学与常规：scRNA+scATAC 的联动分析，以及常规 RNA-seq 与 WES 的基础流程

> 说明：当前仓库正从工作目录中的散落脚本整合迁移（如 `Script/`、`RNA/`、`ATAC/` 等）。本 README 给出统一的代码框架与运行逻辑，随迁移进度持续完善。

---

## 1. 目录结构（规划）

```
BladderCancer_LineageDynamics/
├── config/                      # 全局配置与样本信息
│   ├── samples.list             # 标准化样本列表（待校正样本名）
│   └── references/
│       └── PCT48.ref.fa        # 参考序列与资源（已复制）
├── lineage/                     # 谱系分析
│   ├── upstream/                # 上游：FASTQ → UMI/allele 表
│   └── downstream/              # 下游：树重构、扩增/适应度、图形化
│       ├── tree_reconstruction/
│       ├── expansions/
│       └── fitness/
├── sc_multiomics/               # scRNA 与 scATAC 及其整合
│   ├── scRNA/
│   ├── scATAC/
│   └── integrate/
├── bulk/                        # 常规 RNA-seq 与 WES
│   ├── RNA/
│   └── WES/
├── scripts/                     # 可执行脚本与工具（逐步迁移中）
│   └── utils/
├── resources/                   # 基因集、外部资源（如 .gmt）
│   ├── gene_sets/
│   └── references/
├── notebooks/                   # 复现论文图的 Notebook（可选）
└── results/                     # 统一输出目录
    ├── lineage/
    │   ├── upstream/
    │   └── downstream/
    ├── sc_multiomics/
    └── bulk/
```

---

## 2. 配置与样本命名

- `config/samples.list` 为统一的样本列表文件，用于 Cassiopeia 等上游处理。
- 工作目录中已有的 `samplelist`（如：`/WorkDir4/.../samplelist`）与 `Script/samples.list` 将被整理为标准格式。
- 样本名若有变动（例如 `MP-M-scATAC`、`M-AP_L1_S0004B0004` 等），请统一在 `config/samples.list` 中维护映射，保证下游可复现。

示例（制表符分隔，建议三列最简格式）：
```
<sample_id>  <R1_fastq_path>  <R2_fastq_path>
MP-M-scATAC  scRNA-seq/N-M_S3_L001_R1_001.fastq.gz  scRNA-seq/N-M_S3_L001_R2_001.fastq.gz
M-AP         AP/M-AP_L1_S0004B0004.R1.fastq.gz      AP/M-AP_L1_S0004B0004.R2.fastq.gz
```
> 若需要保留更多元数据（lane、批次、组别），可在首列拼接或扩展为 5 列，但上游脚本需同步支持。

---

## 3. 谱系分析（上游）

目标：从原始 FASTQ 生成 Cassiopeia 可用的 UMI 表与 allele 表，产出统一的中间结果供下游计算。

- 典型流程：
  - `fastp` 质控与转换 → `metapi`/`cassiopeia` 读入 → UMI collapse → 对齐与 allele call → molecule table 过滤
- 现有脚本来源：
  - `Script/step1_fastp.py`（单脚本串流程）
  - `Script/Cassiopeia.smk`（Snakemake 版本）
- 迁移目标：将上游脚本统一到 `lineage/upstream/`，对接 `config/samples.list` 与 `config/references/`。
- 输出位置：`results/lineage/upstream/<sample_id>/...`

运行示例（完成迁移后，示意命令）：
```
python lineage/upstream/run_cassiopeia_pipeline.py \
  --units config/samples.list \
  --ref   config/references/PCT48.ref.fa \
  --out   results/lineage/upstream \
  --threads 16
```

---

## 4. 谱系分析（下游）

目标：在上游结果基础上完成树重构、扩增/适应度、模块评分与图形化，复现论文关键图。

- 子模块与来源：
  - 树重构：`Script/ReconstructingTrees/run_M.py`、`tree.py` 等 → `lineage/downstream/tree_reconstruction/`
  - 扩增与适应度：`Script/call_expanions.py`、`fitness_scores.py` → `lineage/downstream/{expansions,fitness}/`
  - 资源与基因集：`Script/AP-data/.../*.gmt`、`AP-data_MGH/...` → `resources/gene_sets/`
- 输出位置：`results/lineage/downstream/<analysis_name>/...`

示意命令（完成迁移后）：
```
python lineage/downstream/tree_reconstruction/run_tree.py \
  --molecule_table results/lineage/upstream/<sample_id>/molecule_table.tsv \
  --out results/lineage/downstream/tree_reconstruction/<sample_id>
```

---

## 5. scRNA + scATAC 联动

目标：完成 scRNA 与 scATAC 的基本 QC、聚类、注释与差异分析，并在需要时做跨模态整合（如 WNN/Seurat v4）。

- scRNA：`scanpy` 读取 `.h5ad`、质控/聚类/marker/模块评分（现有：`Script/AP_SCRNA.py`） → `sc_multiomics/scRNA/`
- scATAC：常用框架包括 `ArchR`、`Signac`、`cisTopic`，后续迁移至 `sc_multiomics/scATAC/`
- 整合：`sc_multiomics/integrate/` 中实现跨模态整合与可视化（UMAP、热图、river plot 等）。

示意命令（完成迁移后）：
```
python sc_multiomics/scRNA/basic_qc.py --in data/AP_scRNA.h5ad --out results/sc_multiomics/scRNA
python sc_multiomics/integrate/wnn.py --rna results/sc_multiomics/scRNA/adata.h5ad --atac results/sc_multiomics/scATAC/atac.h5ad
```

---

## 6. 常规 RNA-seq 与 WES 基础分析

- RNA-seq：
  - 来源脚本：`RNA/shell_script/run_salmon.sh`、`get_salmon.R`（工作目录中）
  - 迁移目标：将 `salmon` 定量与后续整合整理至 `bulk/RNA/`
  - 输出位置：`results/bulk/RNA/<sample_id>/quant.sf`

- WES：
  - 建议流程：`bwa mem` → `samtools sort` → 标记重复 → `BQSR` → `HaplotypeCaller`
  - 迁移目标：统一到 `bulk/WES/`，按 GATK Best Practices 组织
  - 输出位置：`results/bulk/WES/<sample_id>/vcf/*.vcf.gz`

---

## 7. 环境与依赖

- Python：`>=3.9`，建议使用 Conda 管理环境
- 主要 Python 包：`cassiopeia`、`scanpy`、`anndata`、`pandas`、`numpy`、`scikit-learn`、`pysam`、`seaborn`、`leidenalg`
- 外部工具：`fastp`、`samtools`、`bwa`、`gatk`（如涉及 WES）
- R 生态（可选）：`ArchR` / `Seurat` / `Signac`

> 后续将提供 `environment.yml`/`requirements.txt` 以便一键部署。

---

## 8. 快速开始（占位示例）

1) 准备配置：
```
# 校正或生成统一样本列表
vim config/samples.list

# 放置参考序列与基因集资源
ls config/references/
ls resources/gene_sets/
```

2) 跑谱系上游（待迁移脚本）：
```
python lineage/upstream/run_cassiopeia_pipeline.py \
  --units config/samples.list \
  --ref   config/references/PCT48.ref.fa \
  --out   results/lineage/upstream
```

3) 跑下游树重构与扩增：
```
python lineage/downstream/tree_reconstruction/run_tree.py --molecule_table ...
python lineage/downstream/expansions/call_expansions.py --input ...
```

4) 跑 scRNA 基础分析：
```
python sc_multiomics/scRNA/basic_qc.py --in <adata.h5ad> --out results/sc_multiomics/scRNA
```

5) 跑 RNA-seq `salmon`：
```
bash bulk/RNA/run_salmon.sh --samples config/samples.list --out results/bulk/RNA
```

---

## 9. 迁移与映射计划

- 将以下现有代码逐步迁移至上述结构，并统一入口与输出：
  - `Script/step1_fastp.py` → `lineage/upstream/`
  - `Script/Cassiopeia.smk` → `lineage/upstream/`（可保留 Snakemake 版本）
  - `Script/ReconstructingTrees/*`、`call_expanions.py`、`fitness_scores.py` → `lineage/downstream/`
  - `Script/AP_SCRNA.py` → `sc_multiomics/scRNA/`
  - `RNA/shell_script/*` → `bulk/RNA/`
  - WES 流程脚本（待整理） → `bulk/WES/`
- 基因集与外部资源从 `Script/AP-data*/` 统一至 `resources/gene_sets/`

---

## 10. 提交与复现

- 本仓库将作为 GitHub 公开代码的主体；待迁移完成后执行：
```
git init && git add . && git commit -m "init: code structure and README"
git branch -M main
# 远端地址：https://github.com/TLyan-ze/BladderCancer_LineageDynamics.git
# 添加远端与推送（将于代码审核后执行）
# git remote add origin https://github.com/TLyan-ze/BladderCancer_LineageDynamics.git
# git push -u origin main
```

---

## 11. 备注

- 样本名变更会影响上下游一致性，请以 `config/samples.list` 为唯一真源（single source of truth）。
- 大文件（原始 FASTQ、对齐中间结果）不建议入仓；统一写到 `results/`，并在 README 记录复现步骤。
- 后续将补充脚本入口、参数说明与示例数据，以便审稿复现。
