#!/bin/bash
# A3: Download 149 Pseudomonas_E reference genomes via NCBI datasets
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dram_env

WD=/data/data/Upcycling/research/additional/A3_pseudomonas_ani
cd $WD
mkdir -p ref_genomes zip
ACC=pseudomonas_e_accessions.txt
N=$(wc -l < $ACC)
echo "[A3] downloading $N genomes"

# batch download in groups of 20
split -l 20 $ACC $WD/zip/batch_
for batch in $WD/zip/batch_*; do
    name=$(basename $batch)
    datasets download genome accession --inputfile $batch \
        --filename $WD/zip/${name}.zip --no-progressbar 2>&1 | tail -2 || echo "batch $name failed"
done

# extract all
cd $WD/zip
for z in batch_*.zip; do
    unzip -q -o $z -d unzipped 2>/dev/null || true
done
# move all genomic.fna.gz -> ref_genomes/ACCESSION.fna
find unzipped -name "*_genomic.fna*" -not -name "*from*" | while read f; do
    # extract accession
    base=$(basename $f)
    acc=$(echo $base | grep -oE "GC[AF]_[0-9]+\.[0-9]+")
    if [[ -n "$acc" ]]; then
        if [[ $f == *.gz ]]; then
            gunzip -c "$f" > $WD/ref_genomes/${acc}.fna
        else
            cp "$f" $WD/ref_genomes/${acc}.fna
        fi
    fi
done

echo "[A3] downloaded: $(ls $WD/ref_genomes/*.fna 2>/dev/null | wc -l) genomes"
rm -rf $WD/zip/unzipped
echo "[A3] DONE"
