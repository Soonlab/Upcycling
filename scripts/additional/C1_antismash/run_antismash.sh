#!/usr/bin/env bash
# C1 antiSMASH BGC scan across 111 MAGs
# Inputs: Bakta .gbff per MAG
# Outputs: per-MAG antismash result dir + aggregated counts

set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate antismash_env

BAKTA_DIR="/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
OUT_BASE="/data/data/Upcycling/research/additional/C1_antismash"
PERMAG_DIR="$OUT_BASE/per_mag"
LOG_DIR="$OUT_BASE/logs"
mkdir -p "$PERMAG_DIR" "$LOG_DIR"

THREADS_PER_JOB=4
PARALLEL_JOBS=8   # 4*8 = 32 cores total

run_one() {
  local mag="$1"
  local gbff="$BAKTA_DIR/$mag/$mag.gbff"
  local out="$PERMAG_DIR/$mag"
  if [[ ! -f "$gbff" ]]; then
    echo "[skip] $mag: no gbff" >&2
    return
  fi
  if [[ -f "$out/regions.js" || -f "$out/index.html" ]]; then
    echo "[done] $mag: cached"
    return
  fi
  rm -rf "$out"
  mkdir -p "$out"
  antismash \
    --cb-general --cb-knownclusters --cb-subclusters \
    --asf --pfam2go --cc-mibig \
    --genefinding-tool none \
    --output-dir "$out" \
    -c $THREADS_PER_JOB \
    "$gbff" \
    > "$LOG_DIR/${mag}.log" 2>&1 || echo "[fail] $mag" >&2
}
export -f run_one
export BAKTA_DIR PERMAG_DIR LOG_DIR THREADS_PER_JOB

mags=$(ls "$BAKTA_DIR" | sort)
echo "[$(date)] starting antiSMASH on $(echo "$mags" | wc -l) MAGs, $PARALLEL_JOBS jobs x $THREADS_PER_JOB threads"
echo "$mags" | xargs -n1 -P $PARALLEL_JOBS -I{} bash -c 'run_one "$@"' _ {}
echo "[$(date)] all antiSMASH jobs done"
