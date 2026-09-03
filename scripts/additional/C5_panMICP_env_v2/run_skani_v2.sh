#!/usr/bin/env bash
# C5 re-run: skani ANI of the six MICP-complete MAGs against the accession-verified
# reference panel.  Same parameters as the original run (run_panMICP_env.sh):
#   skani triangle --full-matrix   and   skani dist -q HERO_* -r <all> --min-af 30
# The only change is the reference panel itself, which is now organism-verified.
set -euo pipefail
OUT=/data/data/Upcycling/research/additional/C5_panMICP_env_v2
REFS="$OUT/refs"
SKANI=/home/soon/miniconda3/envs/dram_env/bin/skani

for h in S13 S16 S23 C22 M1 S26; do
    [[ -f "$REFS/HERO_${h}.fna" ]] || zcat "/data/data/Upcycling/MAGs_FASTA_files/${h}.fasta.gz" > "$REFS/HERO_${h}.fna"
done

ls "$REFS"/*.fna > "$OUT/genome_list.txt"
"$SKANI" triangle -l "$OUT/genome_list.txt" --full-matrix \
    -o "$OUT/skani_panMICP.matrix" 2> "$OUT/logs/skani_v2.log"
"$SKANI" dist -q "$REFS"/HERO_*.fna -r "$REFS"/*.fna --min-af 30 \
    -o "$OUT/skani_hero_vs_refs.tsv" 2>> "$OUT/logs/skani_v2.log"
wc -l "$OUT/skani_hero_vs_refs.tsv"
