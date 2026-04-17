#!/bin/bash
# B: Download MGnify livestock MAG catalog metadata + KEGG completeness files
set -e
BASE=/data/data/Upcycling/research/additional/B_rarity_screen
cd $BASE
mkdir -p mgnify
FTP=https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes

# (catalog_dir, version) pairs matching the actual layout
declare -A CATS
CATS[cow-rumen]="v1.0.1"
CATS[sheep-rumen]="v1.0"
CATS[pig-gut]="v1.0"
CATS[chicken-gut]="v1.0.1"

for cat in "${!CATS[@]}"; do
    ver=${CATS[$cat]}
    outdir=$BASE/mgnify/$cat
    mkdir -p $outdir
    echo "[B] === $cat $ver ==="
    # metadata
    wget -q -c -O $outdir/genomes-all_metadata.tsv \
        $FTP/$cat/$ver/genomes-all_metadata.tsv
    # kegg completeness
    wget -q -c -O $outdir/kegg_completeness.tar.gz \
        $FTP/$cat/$ver/pangenome_functions/kegg_completeness.tar.gz
    # functional profiles (contains Pfam etc.)
    wget -q -c -O $outdir/functional_profiles.tar.gz \
        $FTP/$cat/$ver/pangenome_functions/functional_profiles.tar.gz
    ls -la $outdir/
done

echo "[B] download DONE"
