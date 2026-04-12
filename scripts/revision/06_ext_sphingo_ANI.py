#!/usr/bin/env python3
"""
Revision add-on #3 — External Sphingobacterium comparison.

Uses 63 RefSeq reference Sphingobacterium genomes vs. our 6 Sphingobacterium
MAGs (C13, C22, S13, S16, S23, V3) for a definitive novelty test.

Steps:
  1. Unpack NCBI FNA files, gather assembly metadata
  2. skani triangle: pairwise ANI across (our 6) + (63 refs)
  3. For each Sphingobacterium MAG: nearest reference + ANI
  4. Reproduce S13 / S16 novelty claim against the full public panel
  5. Build NJ-like phylogeny from ANI distance and generate a figure
"""
import os, re, json, shutil, subprocess, gzip, tempfile
from glob import glob
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE="/data/data/Upcycling"
EXT =f"{BASE}/research/revision/ext_sphingo"
OUT =f"{BASE}/research/revision"
NCBI_DATA=f"{EXT}/ncbi_dataset/data"

# -----------------------------------------------------------------------------
# 1. Parse assembly report
# -----------------------------------------------------------------------------
meta_rows=[]
with open(f"{NCBI_DATA}/assembly_data_report.jsonl") as fh:
    for line in fh:
        j=json.loads(line)
        acc = j["accession"]
        org = j.get("organism",{}).get("organismName","")
        category = j.get("assemblyInfo",{}).get("refseqCategory","")
        level = j.get("assemblyInfo",{}).get("assemblyLevel","")
        meta_rows.append({"accession":acc,"organism":org,
                          "category":category,"level":level})
meta = pd.DataFrame(meta_rows).set_index("accession")
print(f"[info] reference genomes: {len(meta)}")

# collect fna paths
ref_fnas = {}
for acc in meta.index:
    fn = glob(f"{NCBI_DATA}/{acc}/*.fna")
    if fn: ref_fnas[acc] = fn[0]
print(f"[info] found {len(ref_fnas)} reference fna files")

# our Sphingobacterium MAGs
gtdb = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
gtdb["Genus"] = gtdb["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "")
OUR = gtdb.loc[gtdb.Genus=="Sphingobacterium","user_genome"].tolist()
print(f"[info] our Sphingobacterium MAGs: {OUR}")

# -----------------------------------------------------------------------------
# 2. Build skani listfile mixing our MAGs (uncompressed) + ref genomes
# -----------------------------------------------------------------------------
tmp = tempfile.mkdtemp(prefix="sk_")
listfile = f"{tmp}/list.txt"
with open(listfile,"w") as fh:
    for m in OUR:
        src = f"{BASE}/MAGs_FASTA_files/{m}.fasta.gz"
        dst = f"{tmp}/OUR_{m}.fa"
        with gzip.open(src,"rt") as g, open(dst,"w") as o: o.write(g.read())
        fh.write(dst+"\n")
    for acc, path in ref_fnas.items():
        # symlink ref with a friendly name
        dst = f"{tmp}/REF_{acc}.fa"
        os.symlink(path, dst)
        fh.write(dst+"\n")

skani = "/home/soon/miniconda3/envs/dram_env/bin/skani"
tri = f"{OUT}/skani_ext_sphingo.tsv"
print("[info] running skani triangle (this may take 1-2 min)...")
subprocess.run([skani,"triangle","-l",listfile,"--full-matrix","-o",tri,"-t","8"],
               check=True, capture_output=True)

# -----------------------------------------------------------------------------
# 3. Parse matrix
# -----------------------------------------------------------------------------
with open(tri) as fh:
    lines=[l.rstrip("\n") for l in fh if l.strip()]
n=int(lines[0]); names=[]; rows=[]
for l in lines[1:]:
    p=l.split("\t")
    nm = os.path.basename(p[0]).replace(".fa","")
    names.append(nm)
    rows.append([float(x) if x not in ("","NA") else np.nan for x in p[1:]])
ani = pd.DataFrame(rows, index=names, columns=names)
ani.to_csv(f"{OUT}/ANI_ext_sphingo_matrix.csv")

# -----------------------------------------------------------------------------
# 4. For each of our MAG: nearest ref + ANI
# -----------------------------------------------------------------------------
our_rows = [n for n in names if n.startswith("OUR_")]
ref_rows = [n for n in names if n.startswith("REF_")]
summary=[]
for ourn in our_rows:
    mag = ourn.replace("OUR_","")
    ours = ani.loc[ourn, ref_rows].dropna()
    if ours.empty:
        summary.append({"MAG":mag,"Nearest_accession":"","Nearest_organism":"",
                        "Nearest_ANI":np.nan,"Max_AF":np.nan,
                        "Novel_species_candidate":True})
        continue
    best_ref = ours.idxmax()
    acc = best_ref.replace("REF_","")
    summary.append({
        "MAG": mag,
        "Nearest_accession": acc,
        "Nearest_organism": meta.loc[acc,"organism"] if acc in meta.index else "",
        "Nearest_ANI": round(ours.max(),3),
        "Novel_species_candidate": ours.max() < 95.0
    })
res = pd.DataFrame(summary).sort_values("Nearest_ANI")
res.to_csv(f"{OUT}/ANI_ext_sphingo_novelty.csv", index=False)
print("\n=== Our 6 Sphingobacterium MAGs vs 63 RefSeq references ===")
print(res.to_string(index=False))

# -----------------------------------------------------------------------------
# 5. Figure — ANI distance of our MAGs to every reference
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(len(our_rows), 1, figsize=(10, max(8, len(our_rows)*1.4)),
                         sharex=True)
HERO = {"S13","S16"}
for ax, ourn in zip(axes, our_rows):
    mag = ourn.replace("OUR_","")
    ours = ani.loc[ourn, ref_rows].dropna().sort_values(ascending=False)
    x = np.arange(len(ours))
    colors = ["#c0392b" if ours.iloc[i]>=95 else "#3498db" for i in range(len(ours))]
    ax.bar(x, ours.values, color=colors, edgecolor="black", lw=0.2)
    ax.axhline(95, ls="--", color="red", lw=0.8)
    ax.set_ylabel("ANI (%)")
    ax.set_ylim(70, 100)
    title_col = "#c0392b" if mag in HERO else "black"
    ax.set_title(f"{mag} vs 63 RefSeq Sphingobacterium  |  max ANI = "
                 f"{ours.max():.2f}% ({'novel species candidate' if ours.max()<95 else 'assigned species'})",
                 fontsize=9, loc="left", color=title_col,
                 fontweight="bold" if mag in HERO else "normal")
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.set_xticks([])
axes[-1].set_xlabel("RefSeq genomes, ranked by ANI")
fig.suptitle("External novelty screen — 6 study MAGs against 63 RefSeq Sphingobacterium",
             fontsize=11)
fig.tight_layout(rect=[0,0,1,0.98])
fig.savefig(f"{OUT}/Fig_ext_Sphingo_ANI.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_ext_Sphingo_ANI.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_ext_Sphingo_ANI")

# cleanup tmp
shutil.rmtree(tmp, ignore_errors=True)
