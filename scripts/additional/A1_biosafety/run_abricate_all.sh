#!/bin/bash
# A1: biosafety panel — abricate CARD/VFDB/ResFinder/PlasmidFinder on 111 MAGs
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate biosafety

BAKTA=/data/data/Upcycling/MAGs_FASTA_files/bakta_results
OUT=/data/data/Upcycling/research/additional/A1_biosafety
MAGS_DIR=/data/data/Upcycling/MAGs_FASTA_files
mkdir -p $OUT/{card,vfdb,resfinder,plasmidfinder,ncbi,combined}

# collect MAG fna paths (Bakta produces .fna per MAG)
FNAS=$(ls ${BAKTA}/*/*.fna)
N=$(echo "$FNAS" | wc -l)
echo "[A1] scanning $N MAGs with 4 databases"

for DB in card vfdb resfinder plasmidfinder; do
    echo "[A1] === $DB ==="
    for FNA in $FNAS; do
        MAG=$(basename $FNA .fna)
        abricate --db $DB --threads 4 --minid 80 --mincov 70 "$FNA" > $OUT/$DB/${MAG}.tsv 2>/dev/null || true
    done
    # concat
    {
        head -1 $(ls $OUT/$DB/*.tsv | head -1)
        for f in $OUT/$DB/*.tsv; do tail -n +2 "$f"; done
    } > $OUT/combined/${DB}_all.tsv
    echo "[A1] $DB hits: $(($(wc -l < $OUT/combined/${DB}_all.tsv) - 1))"
done

echo "[A1] DONE"
