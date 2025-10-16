#!/usr/bin/env bash
set -euo pipefail

# 根据 sample1 与 sample2 列表，生成四列制表符分隔的 samples_info：
# 样本ID \t 样本ID \t FQ1绝对路径 \t FQ2绝对路径
# 可按需修改 BASE1 与 BASE2 指向不同批次的清洗后数据目录

BASE1=/datapool/yanzeqin/project/ZDW/RNA/MI-hIGFBP5/00.CleanData
BASE2=/datapool/yanzeqin/project/ZDW/RNA/MI-2nd/00.CleanData

OUT=samples_info
TMP=.samples_info.tmp
> "$TMP"

if [[ -f sample1 ]]; then
  while read -r l1; do
    [[ -z "$l1" ]] && continue
    fq1="${BASE1}/${l1}/${l1}_1.clean.fq.gz"
    fq2="${BASE1}/${l1}/${l1}_2.clean.fq.gz"
    echo -e "${l1}\t${l1}\t${fq1}\t${fq2}" >> "$TMP"
  done < sample1
fi

if [[ -f sample2 ]]; then
  while read -r l1; do
    [[ -z "$l1" ]] && continue
    fq1="${BASE2}/${l1}/${l1}_1.clean.fq.gz"
    fq2="${BASE2}/${l1}/${l1}_2.clean.fq.gz"
    echo -e "${l1}\t${l1}\t${fq1}\t${fq2}" >> "$TMP"
  done < sample2
fi

mv "$TMP" "$OUT"
echo "Wrote $(wc -l < "$OUT") lines to $OUT"
