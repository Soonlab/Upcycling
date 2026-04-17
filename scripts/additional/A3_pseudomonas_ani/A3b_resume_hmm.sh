#!/bin/bash
# A3b resume: finish prodigal + run hmmsearch with shared Pfam HMMs
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dram_env

WD=/data/data/Upcycling/research/additional/A3_pseudomonas_ani
HMM=/data/data/Upcycling/research/additional/hmms_shared/urease_CA.hmm
cd $WD
mkdir -p proteins hmm_out

echo "[A3b] completing prodigal"
for fna in ref_genomes/*.fna; do
    acc=$(basename $fna .fna)
    faa=proteins/${acc}.faa
    if [ ! -s $faa ]; then
        prodigal -i $fna -a $faa -p single -q -o /dev/null 2>/dev/null || \
        prodigal -i $fna -a $faa -p meta -q -o /dev/null 2>/dev/null || \
        echo "fail $acc"
    fi
done
echo "[A3b] proteins: $(ls proteins/*.faa 2>/dev/null | wc -l)"

# Prepare HMM
hmmpress -f $HMM >/dev/null 2>&1 || true

echo "[A3b] running hmmscan against urease+CA HMMs"
i=0
total=$(ls proteins/*.faa 2>/dev/null | wc -l)
for faa in proteins/*.faa; do
    acc=$(basename $faa .faa)
    i=$((i+1))
    if [ -f hmm_out/${acc}.dom ] && [ -s hmm_out/${acc}.dom ]; then continue; fi
    hmmscan --cpu 4 --domtblout hmm_out/${acc}.dom --noali \
            -E 1e-10 $HMM $faa > /dev/null 2>&1 || true
    if [ $((i % 20)) -eq 0 ]; then echo "[A3b] $i/$total done"; fi
done
echo "[A3b] hmm_out files: $(ls hmm_out/*.dom 2>/dev/null | wc -l)"

# Aggregate
python - <<'PY'
import os, glob, re
import pandas as pd
WD = "/data/data/Upcycling/research/additional/A3_pseudomonas_ani"
PF_MAP = {
    "PF00449": "UreC_alpha",
    "PF00699": "UreB_beta_gamma",
    "PF01979": "Amidohydro_superfamily",
    "PF00484": "CA_beta",
    "PF00194": "CA_alpha",
    "PF01682": "Urease_acc_DFG",
    "PF00988": "CA_gamma",
}
rows = []
for dom in sorted(glob.glob(f"{WD}/hmm_out/*.dom")):
    acc = os.path.basename(dom).replace(".dom","")
    present = {v:0 for v in PF_MAP.values()}
    try:
        with open(dom) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split()
                if len(parts) < 22: continue
                t_acc = parts[1]
                pf = t_acc.split(".")[0]
                if pf in PF_MAP:
                    present[PF_MAP[pf]] = 1
    except Exception:
        pass
    present["accession"] = acc
    rows.append(present)

df = pd.DataFrame(rows)
cols = ["accession"] + list(PF_MAP.values())
df = df[cols]
df["urease_core"] = ((df["UreC_alpha"]==1) & (df["UreB_beta_gamma"]==1)).astype(int)
df["urease_acc"] = df["Urease_acc_DFG"]
df["any_CA"] = ((df["CA_beta"]==1) | (df["CA_alpha"]==1) | (df["CA_gamma"]==1)).astype(int)
df["MICP_complete_like"] = (df.urease_core & df.any_CA).astype(int)
df.to_csv(f"{WD}/pseudomonas_e_MICP_rarity_screen.csv", index=False)

n = len(df)
print(f"\n[A3b] MICP screen on {n} Pseudomonas_E reference genomes:")
print(f"  UreC (Pfam PF00449):          {df.UreC_alpha.sum()}/{n}  ({100*df.UreC_alpha.mean():.1f}%)")
print(f"  UreB_beta_gamma (PF00699):     {df.UreB_beta_gamma.sum()}/{n}")
print(f"  urease_core (UreC+UreB):      {df.urease_core.sum()}/{n}  ({100*df.urease_core.mean():.1f}%)")
print(f"  urease_accessory DFG (PF01682):{df.urease_acc.sum()}/{n}")
print(f"  any_CA (alpha+beta+gamma):    {df.any_CA.sum()}/{n}  ({100*df.any_CA.mean():.1f}%)")
print(f"  urease_core + CA (MICP-like): {df.MICP_complete_like.sum()}/{n}  ({100*df.MICP_complete_like.mean():.1f}%)")
PY

echo "[A3b] DONE"
