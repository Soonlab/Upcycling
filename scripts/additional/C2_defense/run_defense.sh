#!/usr/bin/env bash
# C2 Defense systems (DefenseFinder) + CRISPR arrays (minced)
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate defense_env

BAKTA_DIR="/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
MAG_FA_DIR="/data/data/Upcycling/MAGs_FASTA_files"
OUT_BASE="/data/data/Upcycling/research/additional/C2_defense"
DF_DIR="$OUT_BASE/defensefinder"
CR_DIR="$OUT_BASE/crispr_minced"
LOG_DIR="$OUT_BASE/logs"
mkdir -p "$DF_DIR" "$CR_DIR" "$LOG_DIR"

THREADS=4
PJOBS=8

run_one() {
  local mag="$1"
  local faa="$BAKTA_DIR/$mag/$mag.faa"
  local fna_gz="$MAG_FA_DIR/$mag.fasta.gz"
  local df_out="$DF_DIR/$mag"
  local cr_out="$CR_DIR/$mag"
  mkdir -p "$df_out" "$cr_out"

  # DefenseFinder on protein FAA
  if [[ ! -f "$df_out/defense_finder_systems.tsv" && -f "$faa" ]]; then
    defense-finder run -o "$df_out" -w $THREADS "$faa" \
      > "$LOG_DIR/${mag}_df.log" 2>&1 || echo "[df fail] $mag" >&2
  fi
  # minced CRISPR on nucleotide
  if [[ ! -f "$cr_out/${mag}.minced.gff" && -f "$fna_gz" ]]; then
    local tmp_fna="$cr_out/${mag}.fna"
    zcat "$fna_gz" > "$tmp_fna"
    minced -gff "$tmp_fna" "$cr_out/${mag}.minced.txt" "$cr_out/${mag}.minced.gff" \
      > "$LOG_DIR/${mag}_minced.log" 2>&1 || echo "[minced fail] $mag" >&2
    rm -f "$tmp_fna"
  fi
}
export -f run_one
export BAKTA_DIR MAG_FA_DIR DF_DIR CR_DIR LOG_DIR THREADS

mags=$(ls "$BAKTA_DIR" | sort)
echo "[$(date)] starting defense+CRISPR on $(echo "$mags" | wc -l) MAGs"
echo "$mags" | xargs -n1 -P $PJOBS -I{} bash -c 'run_one "$@"' _ {}
echo "[$(date)] done"
