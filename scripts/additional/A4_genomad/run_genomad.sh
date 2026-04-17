#!/bin/bash
# A4: geNomad plasmid/prophage detection on 111 MAGs.
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mge_tools

WD=/data/data/Upcycling/research/additional/A4_genomad
BAKTA=/data/data/Upcycling/MAGs_FASTA_files/bakta_results
DB=$WD/db/genomad_db

mkdir -p $WD/results
cd $WD

N=$(ls $BAKTA/*/*.fna | wc -l)
echo "[A4] running geNomad on $N MAGs (CPU, threads=16 per MAG)"

i=0
for fna in $BAKTA/*/*.fna; do
    mag=$(basename $(dirname $fna))
    outdir=$WD/results/$mag
    if [ -f $outdir/${mag}_summary/${mag}_summary.tsv ]; then
        echo "[A4] $mag already done, skip"
        continue
    fi
    i=$((i+1))
    echo "[A4] [$i/$N] running $mag ..."
    genomad end-to-end --cleanup -t 16 --min-score 0.7 \
            $fna $outdir $DB > $outdir.log 2>&1 || \
            echo "[A4] $mag FAILED (continuing)"
done

echo "[A4] aggregating results"
python - <<'PY'
import os, glob
import pandas as pd
WD = "/data/data/Upcycling/research/additional/A4_genomad/results"
HERO = {"S13","S16","S23","C22","M1","S26"}
rows = []
for d in sorted(os.listdir(WD)):
    sdir = os.path.join(WD, d, f"{d}_summary")
    if not os.path.isdir(sdir): continue
    # plasmid summary
    pfn = f"{sdir}/{d}_plasmid_summary.tsv"
    vfn = f"{sdir}/{d}_virus_summary.tsv"
    np = pd.read_csv(pfn, sep="\t").shape[0] if os.path.exists(pfn) else 0
    nv = pd.read_csv(vfn, sep="\t").shape[0] if os.path.exists(vfn) else 0
    rows.append({"MAG": d,
                 "group": "MICP_complete" if d in HERO else "rest",
                 "n_plasmid_contigs": np, "n_virus_contigs": nv})
df = pd.DataFrame(rows)
df.to_csv("/data/data/Upcycling/research/additional/A4_genomad/genomad_summary.csv", index=False)
print("[A4] summary:")
print(df.groupby("group").agg(
    mean_plas=("n_plasmid_contigs","mean"),
    mean_vir=("n_virus_contigs","mean"),
    n=("MAG","count")))
print("\n[A4] MICP-complete individual:")
print(df[df.group=="MICP_complete"].to_string(index=False))
PY
echo "[A4] DONE"
