#!/bin/bash
# A3: skani ANI of M1, S26 vs 146 Pseudomonas_E references
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dram_env

WD=/data/data/Upcycling/research/additional/A3_pseudomonas_ani
BAKTA=/data/data/Upcycling/MAGs_FASTA_files/bakta_results
cd $WD

# queries: M1, S26 (MICP-complete Pseudomonas_E)
ls $BAKTA/M1/M1.fna > queries.txt
ls $BAKTA/S26/S26.fna >> queries.txt
ls $BAKTA/C22/C22.fna >> queries.txt  # Sphingobacterium but include as outgroup negative control
ls ref_genomes/*.fna > refs.txt
echo "[A3] queries: $(wc -l < queries.txt)"
echo "[A3] refs:    $(wc -l < refs.txt)"

# skani search -- triangle ANI
skani dist --ql queries.txt --rl refs.txt -s 90 -t 16 \
    > skani_pseudomonas_e.tsv 2> skani.log

wc -l skani_pseudomonas_e.tsv
head -5 skani_pseudomonas_e.tsv

# pick closest per query
python - <<'PY'
import pandas as pd
df = pd.read_csv("skani_pseudomonas_e.tsv", sep="\t")
# cols: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query ...
df["query"] = df["Query_file"].str.replace(".*/","", regex=True).str.replace(".fna","", regex=False)
df["ref"]   = df["Ref_file"].str.replace(".*/","", regex=True).str.replace(".fna","", regex=False)
out = df.sort_values(["query","ANI"], ascending=[True, False]).groupby("query").head(5)
out[["query","ref","ANI","Align_fraction_ref","Align_fraction_query"]].to_csv(
    "pseudomonas_e_top5_per_query.csv", index=False)
print("\nTop 5 hits per query:")
print(out[["query","ref","ANI","Align_fraction_query"]].to_string(index=False))

# novelty threshold
print("\nNovelty check (<95% ANI → candidate novel species):")
for q, grp in df.groupby("query"):
    best = grp.sort_values("ANI", ascending=False).head(1).iloc[0]
    ani = best.ANI; ref = best.ref
    novel = "** NOVEL **" if ani < 95 else "assigned"
    print(f"  {q}: best ANI = {ani:.2f}% vs {ref}  → {novel}")
PY

echo "[A3] DONE"
