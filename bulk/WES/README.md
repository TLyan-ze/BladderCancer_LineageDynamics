# WES 全外显子上游流程（Snakemake）

本目录整理了 WES 上游流程与注释，包括对照/肿瘤样本比对、重复标记、Mutect2肿瘤突变检测，以及 ANNOVAR 注释示例。

## 目录结构
- `new.smk`：Snakemake 工作流（fastp → BWA → MarkDuplicates → Mutect2 等）
- `config/`：配置目录
  - `config.yaml`：项目路径、输出子目录、软件路径、资源数据库、PON等
  - `cluster.yaml`：集群资源设置（队列、内存、线程）
  - `sample_now.tsv`：当前样本清单（样本ID、R1、R2；或按workflow所需格式）
  - `samplelist`：历史样本清单
- `mutect.sh`、`MI_MGH.mutect.sh`、`NMI_MGH.mutect.sh`、`MGH_U3.mutect.sh`：单样本/成对样本 Mutect2 调用示例
- `annovar/`：注释示例与说明
  - `readme`：`convert2annovar.pl` 与 `table_annovar.pl` 的示例调用
- `logs/`：示例日志与报错记录
- `run_snakemake_log`：一次完整运行的DAG输出示例

## 运行前准备
- 校验 `config/config.yaml`：
  - `units` 指向样本表：默认 `/datapool/yanzeqin/project/ZDW/Script/config/sample_now.tsv`
  - `output.relative` 作为项目输出根：默认 `/datapool/yanzeqin/project/ZDW`
  - `softwares.fastp`、`softwares.samtools` 等可执行路径存在
  - `databases.Human`、`databases.pon` 指向 hg38 参考与PON资源（示例为 GATK Best Practices）
- 样本清单 `sample_now.tsv`：保证R1/R2路径有效，样本ID与对照关系在workflow内或配置中定义（如正常样本ID：`CTR_MGH`）。

## 运行 Snakemake（集群示例）
建议使用封装脚本：`run_wes_snakemake.sh`
- 内容：设定工作流文件、配置、输出日志与集群参数，然后调用 snakemake。
- 基本命令（参考）：
```
snakemake -s new.smk \
  --configfile config/config.yaml \
  --cluster-config config/cluster.yaml \
  --cluster "qsub -clear -cwd -q all_el6.q -l vf={cluster.mem_mb}G,num_proc={threads}" \
  --jobs 20 --latency-wait 120 --rerun-incomplete --printshellcmds
```

## Mutect2 单独运行示例
以 `MI_MGH.mutect.sh` 为例（对照样本 `CTR_MGH`）：
```
/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" Mutect2 \
  --tmp-dir /datapool/yanzeqin/project/ZDW/mutect2/temp1 \
  -R /share/database/openData/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa \
  -I /datapool/yanzeqin/project/ZDW/MI_MGH/step3_MarkDuplicates/MI_MGH_BQSR.bam \
  -I /datapool/yanzeqin/project/ZDW/CTR_MGH/step3_MarkDuplicates/CTR_MGH_BQSR.bam \
  --normal-sample CTR_MGH \
  --germline-resource /share/database/openData/GATK_BestPracticesResource/somatic-hg38/af-only-gnomad.hg38.subset.vcf.gz \
  --f1r2-tar-gz /datapool/yanzeqin/project/ZDW/mutect2/MI_MGH_F1R2.tar.gz \
  -O /datapool/yanzeqin/project/ZDW/mutect2/MI_MGH_Unfiltered.vcf
```
后续需 `FilterMutectCalls`、`bcftools view`（筛选PASS）、`LearnReadOrientationModel` 等步骤，可在工作流中统一执行。

## ANNOVAR 注释示例
- 转换：`convert2annovar.pl -format vcf4old MI_MGH_Filtered.vcf > MI_MGH.annovar`
- 注释：
```
perl /datapool/yanzeqin/software/annovar/table_annovar.pl MI_MGH.annovar \
  /datapool/yanzeqin/software/annovar/humandbhg38_annovar/ \
  -buildver hg38 -out MI_MGH -remove -protocol refGene -operation g -nastring . -polish -otherinfo
```

## 常见问题
- `Wildcards object has no attribute output`：见 `run_snakemake_log` 末尾；说明某条规则的通配符引用不正确，需修正 `new.smk` 第64行附近。
- `python3.10: No such file or directory`：集群环境Python路径缺失，调整 snakemake 可执行的环境或模块。
- BWA/fastp 日志显示完成但后续中断：检查 `config.yaml` 输出路径与文件命名是否匹配规则。

## 产出
- `step1_fastp/`：清洗后的FASTQ与质控报告
- `step2_bwa/`：比对BAM与排序BAM
- `step3_MarkDuplicates/`：去重复BAM、metrics与BQSR（如配置）
- `step5_mutect2withpon/` 或 `mutect2/`：突变VCF与F1R2中间文件
- `annovar/`：注释结果（`*.hg38_multianno.txt` 等）

