#!/bin/bash
# Re-run just plasmidfinder (+ bacmet2, argannot) — abricate
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate biosafety

BAKTA=/data/data/Upcycling/MAGs_FASTA_files/bakta_results
OUT=/data/data/Upcycling/research/additional/A1_biosafety
mkdir -p $OUT/plasmidfinder

for FNA in $BAKTA/*/*.fna; do
    MAG=$(basename $FNA .fna)
    abricate --db plasmidfinder --threads 4 --minid 80 --mincov 70 "$FNA" \
        > $OUT/plasmidfinder/${MAG}.tsv 2>/dev/null || true
done

# concat
{
    head -1 $(ls $OUT/plasmidfinder/*.tsv | head -1)
    for f in $OUT/plasmidfinder/*.tsv; do tail -n +2 "$f"; done
} > $OUT/combined/plasmidfinder_all.tsv
echo "plasmidfinder hits: $(($(wc -l < $OUT/combined/plasmidfinder_all.tsv) - 1))"
