#!/usr/bin/env python3
"""ANI-based novelty screen + hero-clade pairwise ANI with skani."""
import os, re, subprocess, tempfile, gzip, shutil
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE="/data/data/Upcycling"
OUT=f"{BASE}/research/extra"
os.makedirs(OUT, exist_ok=True)

# 1) novelty from GTDB-Tk (closest_genome_ani vs 95% species cutoff)
g = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
g["Genus"]   = g["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "Unclassified")
g["Species"] = g["classification"].apply(lambda s: (re.search(r"s__([^;]*)",str(s)) or [None,""])[1])
def to_float(x):
    try: return float(x)
    except: return np.nan
g["ANI"] = g["closest_genome_ani"].apply(to_float)
g["Novel_sp_candidate"] = (g["ANI"] < 95.0) | (g["Species"]=="")
HERO=["S13","S16","S23","S26","C22","M1"]
g["Hero"] = g["user_genome"].isin(HERO)
nov = g[["user_genome","Genus","Species","ANI","Novel_sp_candidate","Hero"]]
nov.sort_values(["Hero","Novel_sp_candidate","ANI"],ascending=[False,False,True], inplace=True)
nov.to_csv(f"{OUT}/novelty_ANI_screen.csv", index=False)
print(f"[ok] novelty_ANI_screen.csv  ({nov.Novel_sp_candidate.sum()} candidates)")
print(nov[nov.Hero].to_string(index=False))

# 2) skani pairwise among all 111 MAGs (fast)
src = f"{BASE}/MAGs_FASTA_files"
tmpdir = tempfile.mkdtemp(prefix="skani_")
fnas=[]
for gz in sorted([x for x in os.listdir(src) if x.endswith(".fasta.gz")]):
    samp = gz.replace(".fasta.gz","")
    out = os.path.join(tmpdir, f"{samp}.fa")
    with gzip.open(os.path.join(src,gz),"rt") as fh, open(out,"w") as fo:
        fo.write(fh.read())
    fnas.append(out)
listfile = os.path.join(tmpdir, "list.txt")
with open(listfile,"w") as fh:
    for f in fnas: fh.write(f+"\n")

skani = "/home/soon/miniconda3/envs/dram_env/bin/skani"
tri_out = f"{OUT}/skani_triangle.tsv"
print("[info] running skani triangle ...")
subprocess.run([skani,"triangle","-l",listfile,"--full-matrix","-o",tri_out,"-t","8"], check=True)

shutil.rmtree(tmpdir, ignore_errors=True)

# 3) parse the full matrix
with open(tri_out) as fh:
    lines = [l.rstrip("\n") for l in fh if l.strip()]
n = int(lines[0])
names=[]; rows=[]
for l in lines[1:]:
    parts = l.split("\t")
    names.append(os.path.basename(parts[0]).replace(".fa",""))
    vals = [float(x) if x not in ("","NA") else np.nan for x in parts[1:]]
    rows.append(vals)
ani = pd.DataFrame(rows, index=names, columns=names)
ani.to_csv(f"{OUT}/skani_full_matrix.csv")

# 4) hero-clade submatrix heatmap
sub = ani.loc[HERO, HERO]
fig, ax = plt.subplots(figsize=(6,5))
sns.heatmap(sub, annot=True, fmt=".1f", cmap="viridis", vmin=80, vmax=100,
            cbar_kws={"label":"ANI (%)"}, ax=ax, linewidths=0.4, linecolor="white")
ax.set_title("Pairwise ANI within Sphingobacterium hero clade", fontsize=10)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_T2d_HeroANI.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_T2d_HeroANI.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_T2d_HeroANI")

# 5) novelty overview figure
fig, ax = plt.subplots(figsize=(8,5))
data = nov.dropna(subset=["ANI"])
colors = ["#c0392b" if h else "#7f8c8d" for h in data["Hero"]]
ax.scatter(range(len(data)), data["ANI"], c=colors, s=[60 if h else 16 for h in data["Hero"]],
           edgecolor="black", lw=0.3, zorder=3)
ax.axhline(95, ls="--", color="red", lw=0.8, label="species cutoff (95%)")
ax.set_xlabel("MAG (sorted)"); ax.set_ylabel("ANI to closest GTDB genome (%)")
ax.set_title(f"Novelty screen — {nov.Novel_sp_candidate.sum()} candidates below 95% ANI",
             fontsize=10)
ax.legend(frameon=False)
# label hero points
sdata = data.sort_values("ANI").reset_index(drop=True)
for i,r in sdata.iterrows():
    if r.Hero:
        ax.annotate(r.user_genome, (i, r.ANI), xytext=(3,4), textcoords="offset points",
                    fontsize=8, color="#c0392b", fontweight="bold")
ax.set_xticks([])
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_T2d_Novelty_overview.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_T2d_Novelty_overview.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_T2d_Novelty_overview")
