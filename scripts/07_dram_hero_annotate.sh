#!/bin/bash
set -u
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate dram_env

OUT=/data/pangenome_work/dram_output
mkdir -p $OUT
CFG=/home/soon/.dram_config

run_one() {
    s=$1
    in=/data/pangenome_work/input_fna_all/${s}.fna
    od=$OUT/$s
    [ -d "$od" ] && rm -rf "$od"
    echo "[start] $s  $(date)"
    DRAM.py annotate -i "$in" -o "$od" --config_loc "$CFG" --threads 8 \
        > $OUT/${s}.log 2>&1 && echo "[done] $s $(date)" || echo "[fail] $s"
}
export -f run_one
export OUT CFG

# 6 hero samples, 3 in parallel
printf "%s\n" C22 M1 S13 S16 S23 S26 | xargs -n1 -P3 -I{} bash -c 'run_one "$@"' _ {}

echo "[all annotate done] $(date)"

# Distill
DRAM.py distill -i "$OUT/*/annotations.tsv" -o $OUT/distillate \
    --rrna_path "$OUT/*/rrnas.tsv" --trna_path "$OUT/*/trnas.tsv" \
    > $OUT/distill.log 2>&1 && echo "[distill done]" || echo "[distill fail]"
