#!/bin/bash
# A3b: Screen MICP operon presence in 146 Pseudomonas_E reference genomes
# Method: prodigal predict genes -> hmmsearch against urease subunit + cah HMMs
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dram_env

WD=/data/data/Upcycling/research/additional/A3_pseudomonas_ani
cd $WD
mkdir -p proteins hmm_out
REFDIR=ref_genomes

echo "[A3b] predicting genes on 146 Pseudomonas_E refs with prodigal"
for fna in $REFDIR/*.fna; do
    acc=$(basename $fna .fna)
    faa=proteins/${acc}.faa
    if [ ! -s $faa ]; then
        prodigal -i $fna -a $faa -p single -q -o /dev/null 2>/dev/null || \
        prodigal -i $fna -a $faa -p meta -q -o /dev/null
    fi
done
echo "[A3b] proteomes: $(ls proteins/*.faa 2>/dev/null | wc -l)"

# Fetch KEGG/Pfam HMMs for urease operon if not present
HMM_DIR=$WD/hmms
mkdir -p $HMM_DIR

# Use DRAM's existing HMMs? Alternative: use PFAM
# KEGG orthologs:
#   K01428 ureC (urease alpha)
#   K01429 ureB (urease beta)
#   K01430 ureA (urease gamma)
#   K14048 ureAB fusion (some bacteria)
#   K03187 ureE (accessory)
#   K03188 ureF
#   K03189 ureG
#   K03190 ureD (=ureH)
#   K01672 cah (carbonic anhydrase, alpha class)
#   K18246 cynT (beta-CA)
#   K01673 canA (gamma-CA)

# Use Pfam profiles for robust screening:
#   PF00449.18 = Urease alpha
#   PF00699.18 = Urease beta/gamma
#   PF01979.22 = Amidohydro_1 (urease superfamily)
#   PF00484.21 = Carbonic anhydrase (beta)
#   PF08445.14 = Fe2+/Zn2+ uptake regulation (not for cah)
#   PF00194.23 = Alpha carbonic anhydrase
#   PF00870.22 = p53  (not used here)

if [ ! -s $HMM_DIR/urease.hmm ]; then
    cd $HMM_DIR
    # PFAM direct download
    for pf in PF00449 PF00699 PF01979 PF00484 PF00194 PF08009; do
        # PF08009 = ureE-like carbohydrate kinase -- not urease. skip.
        wget -q -O ${pf}.hmm.gz "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${pf}?annotation=hmm" 2>/dev/null || true
        if [ -s ${pf}.hmm.gz ]; then
            gunzip -f ${pf}.hmm.gz
        fi
    done
    cat *.hmm > urease_all.hmm 2>/dev/null || true
    ls -la $HMM_DIR/
    cd $WD
fi

echo "[A3b] HMM file size: $(ls -la $HMM_DIR/urease_all.hmm 2>/dev/null | awk '{print $5}')"

# Run hmmsearch
if [ -s $HMM_DIR/urease_all.hmm ]; then
    hmmpress -f $HMM_DIR/urease_all.hmm >/dev/null 2>&1 || true
    echo "[A3b] running hmmsearch against Pfam urease domains"
    for faa in proteins/*.faa; do
        acc=$(basename $faa .faa)
        hmmscan --cpu 4 --domtblout hmm_out/${acc}.dom --noali \
                -E 1e-10 $HMM_DIR/urease_all.hmm $faa > /dev/null 2>&1 || true
    done
    echo "[A3b] hmmsearch results: $(ls hmm_out/*.dom 2>/dev/null | wc -l)"
fi

# Summarize: is ureC (PF00449) present? cah/beta-CA (PF00484 or PF00194)? etc.
python - <<'PY'
import os, glob, re
import pandas as pd

def parse_dom(path):
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 22: continue
            target_name, target_acc = parts[0], parts[1]
            query_name = parts[3]
            evalue = float(parts[6])
            hits.append((target_name, target_acc, query_name, evalue))
    return hits

# column names for Pfam -> readable names
PF_MAP = {
    "PF00449": "UreC_alpha",
    "PF00699": "UreB_beta_gamma",
    "PF01979": "Amidohydro_urease_superfamily",
    "PF00484": "CA_beta_class",
    "PF00194": "CA_alpha_class",
}

rows = []
for dom in sorted(glob.glob("hmm_out/*.dom")):
    acc = os.path.basename(dom).replace(".dom","")
    hits = parse_dom(dom)
    present = {PF_MAP[k]:0 for k in PF_MAP}
    for t_name, t_acc, q_name, ev in hits:
        # t_acc like "PF00449.18"
        pf_id = t_acc.split(".")[0]
        if pf_id in PF_MAP:
            present[PF_MAP[pf_id]] = 1
    present["accession"] = acc
    rows.append(present)

df = pd.DataFrame(rows)
cols = ["accession"] + list(PF_MAP.values())
df = df[cols]
df["urease_operon_all3"] = ((df["UreC_alpha"]==1) & (df["UreB_beta_gamma"]==1)).astype(int)
df["has_CA"] = ((df["CA_beta_class"]==1) | (df["CA_alpha_class"]==1)).astype(int)
df["MICP_complete_like"] = (df["urease_operon_all3"] & df["has_CA"]).astype(int)

df.to_csv("pseudomonas_e_MICP_rarity_screen.csv", index=False)

print("\n[A3b] MICP operon screen across 146 Pseudomonas_E refs:")
print(f"  UreC (α) present:           {df['UreC_alpha'].sum()}/{len(df)}")
print(f"  UreB (β/γ) present:         {df['UreB_beta_gamma'].sum()}/{len(df)}")
print(f"  UreC + UreB (urease core):  {df['urease_operon_all3'].sum()}/{len(df)}")
print(f"  CA (α or β) present:        {df['has_CA'].sum()}/{len(df)}")
print(f"  urease + CA (MICP-like):    {df['MICP_complete_like'].sum()}/{len(df)}  = {100*df['MICP_complete_like'].mean():.1f}%")
print(f"\n[A3b] output: pseudomonas_e_MICP_rarity_screen.csv")
PY

echo "[A3b] DONE"
