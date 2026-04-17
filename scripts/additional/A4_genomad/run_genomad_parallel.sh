#!/bin/bash
# A4 speedup: run a second parallel worker targeting S/V MAGs (where heroes are)
# so we don't have to wait for alphabetical progression through C/M.
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mge_tools

WD=/data/data/Upcycling/research/additional/A4_genomad
BAKTA=/data/data/Upcycling/MAGs_FASTA_files/bakta_results
DB=$WD/db/genomad_db

mkdir -p $WD/results
cd $WD

# priority MAG list — hero MAGs first, then S/V series (most likely not yet started by main worker)
# Main worker started from C* so M* and later are untouched
PRIORITY="S13 S16 S23 M1 S26"

# other S/V MAGs as padding
SVS=$(ls $BAKTA | grep -E "^(S|V)" | grep -v -E "^($(echo $PRIORITY | tr ' ' '|'))$")

ORDER="$PRIORITY $SVS"
echo "[A4p] processing ${ORDER}"

for mag in $ORDER; do
    fna=$BAKTA/$mag/$mag.fna
    outdir=$WD/results/$mag
    # skip if main worker already picked it up (even partial)
    if [ -d $outdir ] && [ -f $outdir/${mag}_summary/${mag}_plasmid_summary.tsv ]; then
        echo "[A4p] $mag already DONE, skip"
        continue
    fi
    if [ -d $outdir ] && [ ! -f $outdir/${mag}_summary/${mag}_plasmid_summary.tsv ]; then
        echo "[A4p] $mag started but not done — assume other worker, skip"
        continue
    fi
    echo "[A4p] $mag ..."
    genomad end-to-end --cleanup -t 16 --min-score 0.7 \
        $fna $outdir $DB > $outdir.log 2>&1 || echo "[A4p] $mag FAILED"
done

echo "[A4p] DONE"
