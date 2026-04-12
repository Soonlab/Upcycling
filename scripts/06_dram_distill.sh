#!/bin/bash
set -u
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate dram_env
D=/data/pangenome_work/dram_output
O=$D/distillate
mkdir -p $O

python - <<'PY'
import os, pandas as pd, glob
D="/data/pangenome_work/dram_output"
anns=[]
for f in sorted(glob.glob(f"{D}/*/annotations.tsv")):
    try:
        df=pd.read_csv(f, sep="\t", index_col=0)
        anns.append(df)
    except Exception as e: print("skip",f,e)
big=pd.concat(anns)
big.to_csv(f"{D}/all_annotations.tsv", sep="\t")
print("merged", big.shape)

# merge trnas/rrnas
for kind in ["trnas","rrnas"]:
    rows=[]
    for f in sorted(glob.glob(f"{D}/*/{kind}.tsv")):
        try: rows.append(pd.read_csv(f, sep="\t"))
        except Exception: pass
    if rows:
        pd.concat(rows).to_csv(f"{D}/all_{kind}.tsv", sep="\t", index=False)
        print(f"merged {kind}:", sum(len(r) for r in rows))
PY

DRAM.py distill -i $D/all_annotations.tsv -o $O \
    --rrna_path $D/all_rrnas.tsv --trna_path $D/all_trnas.tsv 2>&1 | tail -30
