#!/usr/bin/env bash
set -euo pipefail

# 封装调用 get_salmon.R，避免直接修改 R 脚本
# 用法：bash run_get_salmon.sh -d /path/to/03.salmon

DIR=""
while getopts ":d:" opt; do
  case $opt in
    d) DIR="$OPTARG" ;;
    *) echo "Usage: $0 -d <salmon_output_dir>"; exit 1 ;;
  esac
done

if [[ -z "$DIR" ]]; then
  echo "Usage: $0 -d <salmon_output_dir>"
  exit 1
fi

Rscript sc_multiomics/scRNA/Salmon/get_salmon.R "$DIR"
