#!/usr/bin/env python3
"""
Revision add-on #2 — PCoA of the 111 MAGs by waste source (cattle / swine /
sheep / poultry) using Panaroo gene presence-absence + PERMANOVA test.
Also computes DRAM-module PCoA as complementary view.
"""
import os, re, itertools
import numpy as np, pandas as pd
from scipy.spatial.distance import pdist, squareform
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from skbio.stats.distance import DistanceMatrix, permanova
from skbio.stats.ordination import pcoa

BASE="/data/data/Upcycling"
OUT =f"{BASE}/research/revision"
os.makedirs(OUT, exist_ok=True)

# taxonomy + hero
g = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
g["Genus"] = g["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "Unclassified")
tax = g.set_index("user_genome")["Genus"]
HERO = {"C22","M1","S13","S16","S23","S26"}

# waste source from first letter of MAG ID
def source(m):
    return {"C":"Cattle","M":"Swine","S":"Sheep","V":"Poultry"}.get(m[0],"Other")

# -----------------------------------------------------------------------------
# (1) Panaroo Rtab — gene presence/absence
# -----------------------------------------------------------------------------
rtab = pd.read_csv(f"{BASE}/pangenome_work/panaroo_results/gene_presence_absence.Rtab",
                   sep="\t", index_col=0)
print(f"Panaroo: {rtab.shape[0]} genes × {rtab.shape[1]} MAGs")
# Filter genes found in 5–95 % of MAGs (drop core + singletons for informative PCoA)
prevalence = rtab.sum(axis=1)/rtab.shape[1]
keep = (prevalence>=0.05) & (prevalence<=0.95)
rtab_f = rtab.loc[keep]
print(f"Informative genes (5-95% prevalence): {len(rtab_f)}")

mat = rtab_f.T.values.astype(int)  # rows = MAGs
mags = rtab_f.columns.tolist()
# Jaccard distance on binary
dist = pdist(mat, metric="jaccard")
dmat = squareform(dist)
dm = DistanceMatrix(dmat, ids=mags)

# metadata
meta = pd.DataFrame({"MAG":mags,
                     "Source":[source(m) for m in mags],
                     "Genus":[tax.get(m,"Unknown") for m in mags],
                     "Hero":[m in HERO for m in mags]}).set_index("MAG")

# PERMANOVA by source
perm = permanova(dm, grouping=meta.loc[mags,"Source"], permutations=999)
print(f"\nPERMANOVA (Source): pseudo-F = {perm['test statistic']:.3f}, p = {perm['p-value']:.4f}")
perm_genus = permanova(dm, grouping=meta.loc[mags,"Genus"], permutations=999)
print(f"PERMANOVA (Genus):  pseudo-F = {perm_genus['test statistic']:.3f}, p = {perm_genus['p-value']:.4f}")

# PCoA
ord_ = pcoa(dm, number_of_dimensions=3)
coords = ord_.samples.copy()
coords.index = mags
var_exp = ord_.proportion_explained * 100
meta = meta.join(coords[["PC1","PC2","PC3"]])
meta.to_csv(f"{OUT}/PCoA_panaroo_coords.csv")
pd.DataFrame({"pseudo_F_source":perm["test statistic"],
              "p_source":perm["p-value"],
              "pseudo_F_genus":perm_genus["test statistic"],
              "p_genus":perm_genus["p-value"],
              "PC1_var":var_exp.iloc[0],
              "PC2_var":var_exp.iloc[1]}, index=[0]).to_csv(f"{OUT}/PCoA_PERMANOVA.csv", index=False)

# -----------------------------------------------------------------------------
# Figures — two panels: colour by source and by genus
# -----------------------------------------------------------------------------
SOURCE_COLORS = {"Cattle":"#e74c3c","Swine":"#3498db","Sheep":"#2ecc71","Poultry":"#f39c12"}
genera = sorted(meta["Genus"].unique())
gcmap = plt.get_cmap("tab20", len(genera))
GENUS_COLORS = {g: gcmap(i) for i,g in enumerate(genera)}

fig, axes = plt.subplots(1,2, figsize=(14,6))

# panel A — source
ax = axes[0]
for s, col in SOURCE_COLORS.items():
    sub = meta[meta.Source==s]
    ax.scatter(sub.PC1, sub.PC2, s=35, color=col, alpha=0.8, edgecolor="black", lw=0.4,
               label=f"{s} (n={len(sub)})")
# highlight hero
hero = meta[meta.Hero]
ax.scatter(hero.PC1, hero.PC2, s=130, facecolor="none", edgecolor="red", lw=1.8, zorder=5)
for _,r in hero.iterrows():
    ax.annotate(r.name, (r.PC1, r.PC2), xytext=(5,5), textcoords="offset points",
                fontsize=8, color="red", fontweight="bold")
ax.set_xlabel(f"PC1 ({var_exp.iloc[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_exp.iloc[1]:.1f}%)")
ax.set_title(f"a  PCoA coloured by waste source\nPERMANOVA pseudo-F={perm['test statistic']:.2f}, p={perm['p-value']:.3f}",
             fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="best")
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax.grid(ls=":", alpha=0.3)

# panel B — genus
ax = axes[1]
counts = meta["Genus"].value_counts()
for gname, col in GENUS_COLORS.items():
    if counts.get(gname,0) < 2: continue
    sub = meta[meta.Genus==gname]
    ax.scatter(sub.PC1, sub.PC2, s=35, color=col, alpha=0.85, edgecolor="black", lw=0.4,
               label=f"{gname} (n={len(sub)})")
hero = meta[meta.Hero]
ax.scatter(hero.PC1, hero.PC2, s=130, facecolor="none", edgecolor="red", lw=1.8, zorder=5)
ax.set_xlabel(f"PC1 ({var_exp.iloc[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_exp.iloc[1]:.1f}%)")
ax.set_title(f"b  PCoA coloured by GTDB genus\nPERMANOVA pseudo-F={perm_genus['test statistic']:.2f}, p={perm_genus['p-value']:.3f}",
             fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=7, loc="best", ncol=1)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax.grid(ls=":", alpha=0.3)

fig.suptitle("Pan-genome PCoA (Jaccard on Panaroo presence/absence, 5-95 % prevalence filter)",
             fontsize=11)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_PCoA_source_genus.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_PCoA_source_genus.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_PCoA_source_genus")

# -----------------------------------------------------------------------------
# (2) Pairwise PERMANOVA between sources (post-hoc)
# -----------------------------------------------------------------------------
sources = sorted(meta["Source"].unique())
rows=[]
for a,b in itertools.combinations(sources,2):
    ids = meta.index[meta.Source.isin([a,b])].tolist()
    sub_dm = dm.filter(ids)
    grp = meta.loc[ids,"Source"]
    res = permanova(sub_dm, grouping=grp, permutations=999)
    rows.append({"A":a,"B":b,"n_A":(grp==a).sum(),"n_B":(grp==b).sum(),
                 "pseudo_F":round(res["test statistic"],3),"p":round(res["p-value"],4)})
ph = pd.DataFrame(rows)
# BH FDR
from scipy.stats import false_discovery_control
ph["q_BH"] = false_discovery_control(ph["p"].values).round(4)
ph.to_csv(f"{OUT}/PCoA_pairwise_PERMANOVA.csv", index=False)
print("\n=== Pairwise PERMANOVA between waste sources ===")
print(ph.to_string(index=False))

# -----------------------------------------------------------------------------
# (3) Trait-module based PCoA (orthogonal view)
# -----------------------------------------------------------------------------
counts = pd.read_csv(f"{BASE}/research/extra/gene_category_counts.csv")
subcols = [c for c in counts.columns if "::" in c]
cds = counts.set_index("Sample")["CDS_total"]
norm = counts.set_index("Sample")[subcols].div(cds, axis=0)*1000.0
norm = norm.loc[[m for m in mags if m in norm.index]]
dist_t = pdist(norm.values, metric="euclidean")
dm_t = DistanceMatrix(squareform(dist_t), ids=norm.index.tolist())

perm_t = permanova(dm_t, grouping=meta.loc[norm.index,"Source"], permutations=999)
ord_t = pcoa(dm_t, number_of_dimensions=3)
coords_t = ord_t.samples.copy(); coords_t.index = norm.index
meta_t = meta.loc[norm.index].join(coords_t[["PC1","PC2"]], rsuffix="_trait")
var_t = ord_t.proportion_explained*100

fig, ax = plt.subplots(figsize=(7,6))
for s,col in SOURCE_COLORS.items():
    sub = meta_t[meta_t.Source==s]
    ax.scatter(sub["PC1_trait"], sub["PC2_trait"], s=35, color=col, alpha=0.85,
               edgecolor="black", lw=0.4, label=f"{s} (n={len(sub)})")
hero_t = meta_t[meta_t.Hero]
ax.scatter(hero_t["PC1_trait"], hero_t["PC2_trait"], s=140, facecolor="none",
           edgecolor="red", lw=1.8, zorder=5)
for _,r in hero_t.iterrows():
    ax.annotate(r.name, (r["PC1_trait"], r["PC2_trait"]), xytext=(5,5),
                textcoords="offset points", fontsize=8, color="red", fontweight="bold")
ax.set_xlabel(f"PC1 ({var_t.iloc[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_t.iloc[1]:.1f}%)")
ax.set_title(f"Trait-module PCoA (Euclidean, 38 modules × per-1k CDS)\n"
             f"Source PERMANOVA pseudo-F={perm_t['test statistic']:.2f}, p={perm_t['p-value']:.3f}",
             fontsize=10)
ax.legend(frameon=False, fontsize=8)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax.grid(ls=":", alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_PCoA_trait_source.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_PCoA_trait_source.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_PCoA_trait_source")
print(f"Trait PERMANOVA Source: pseudo-F={perm_t['test statistic']:.3f}, p={perm_t['p-value']:.4f}")
