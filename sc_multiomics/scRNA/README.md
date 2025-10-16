# scRNA 分析流程（谱系结果）

本目录收纳来自谱系结果目录的 scRNA 分析脚本，目录参考：
`/WorkDir4/yanzeqin/220823_A00403_0863_AHMNGLDSX3/Script/AP-data/expression`

- 仓库位置：`sc_multiomics/scRNA`
- 脚本位置：`sc_multiomics/scRNA/scripts`
- 说明：此处仅保存流程脚本；大体量结果文件（如 h5ad / h5seurat / mtx / plots）不纳入版本控制。

## 包含脚本
- `convert_MGH.R`、`convert_RT.R`、`convert_all.R`：数据转换与对象构建（如 h5seurat / h5ad）。
- `step1_h5ad.r`、`step_h5ad_MGH.r`、`step_h5ad_RT.r`：预处理与对象生成步骤。
- `step2_MGH.analysis.py`、`step2_RT.analysis.py`、`step2_all.analysis.py`：下游分析与绘图。

## 依赖建议
- R：`Seurat`、`tidyverse` 等（依具体脚本而定）
- Python：`scanpy`、`anndata`、`pandas` 等（依具体脚本而定）
- 输入数据：参考谱系结果目录中的 `matrix.mtx`、`genes.tsv`、`barcodes.tsv` 等。

## 快速使用
1. 在工作目录准备输入数据（参考谱系目录），并创建结果输出文件夹。
2. 运行 R 预处理脚本（示例）：
   - `Rscript sc_multiomics/scRNA/scripts/step1_h5ad.r`
   - `Rscript sc_multiomics/scRNA/scripts/step_h5ad_MGH.r`
3. 运行 Python 分析脚本（示例）：
   - `python sc_multiomics/scRNA/scripts/step2_MGH.analysis.py`
4. 结果输出建议存放于项目的结果目录，不纳入仓库。

## 与 bulk RNA 的关系
- bulk RNA 上游（Salmon）位于：`bulk/RNA/Upstream/Salmon`
- bulk RNA 下游差异与富集：`bulk/RNA/Downstream`
- scRNA 与 bulk RNA 的流程分开维护，避免混淆。

