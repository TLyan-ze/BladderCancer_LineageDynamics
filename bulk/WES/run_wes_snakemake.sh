#!/usr/bin/env bash
set -euo pipefail

# 运行WES Snakemake工作流（集群环境示例）
# 用法：bash run_wes_snakemake.sh [-j 20]

JOBS=20
if [[ $# -gt 0 ]]; then
  JOBS="$1"
fi

WORKFLOW="bulk/WES/new.smk"
CONF="bulk/WES/config/config.yaml"
CLUSTER_CONF="bulk/WES/config/cluster.yaml"
LOGDIR="bulk/WES/logs"
mkdir -p "$LOGDIR"

snakemake -s "$WORKFLOW" \
  --configfile "$CONF" \
  --cluster-config "$CLUSTER_CONF" \
  --cluster "qsub -clear -cwd -q all_el6.q -l vf={cluster.mem_mb}G,num_proc={threads}" \
  --jobs "$JOBS" --latency-wait 120 --rerun-incomplete --printshellcmds \
  | tee "$LOGDIR/snakemake_run_$(date +%F_%H%M%S).log"
