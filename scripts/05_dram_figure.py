#!/usr/bin/env python3
"""Fig 4 — DRAM metabolic module heatmap, tree-ordered, hero-highlighted."""
import os, re
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import Phylo

BASE="/data/data/Upcycling"
OUT=f"{BASE}/research/extra"
prod = pd.read_csv("/data/pangenome_work/dram_output/distillate/product.tsv", sep="\t", index_col=0)
print("product shape:", prod.shape)

gtdb = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
gtdb["Genus"] = gtdb["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "Unclassified")
tax = gtdb.set_index("user_genome")["Genus"]

HERO={"S13","S16","S23","S26","C22","M1"}

# pick informative modules: keep columns with std>0 and at least 3 non-zero entries
keep = [c for c in prod.columns if prod[c].std()>0 and (prod[c]>0).sum()>=3]
prod = prod[keep]

# choose a curated MICP-relevant subset + top variable modules
curated = [c for c in prod.columns if re.search(
    r"urease|carbonic anhyd|Nitrogen|Na.*H.*antiport|glutamate|glutamine|"
    r"Cobalamin|B12|Urea cycle|Arginine|Carbohydrate|CAZy|biofilm|"
    r"Flagell|Methionine|Serine|Fatty acid|Pyruvate|Glycolysis|TCA", c, re.I)]
curated = list(dict.fromkeys(curated))[:35]
prod_sel = prod[curated].copy()
prod_sel["Genus"] = [tax.get(i,"Unclassified") for i in prod_sel.index]
prod_sel["Hero"]  = [i in HERO for i in prod_sel.index]

# order: hero first (highlight), then by genus
prod_sel = prod_sel.sort_values(["Hero","Genus"], ascending=[False, True])

fig, ax = plt.subplots(figsize=(max(10, len(curated)*0.35), max(10, len(prod_sel)*0.11)))
mat = prod_sel[curated].values.astype(float)
im = ax.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
ax.set_xticks(range(len(curated)))
ax.set_xticklabels(curated, rotation=60, ha="right", fontsize=6.5)
ax.set_yticks(range(len(prod_sel)))
labels = [f"{i}  [{g}]" for i,g in zip(prod_sel.index, prod_sel["Genus"])]
ax.set_yticklabels(labels, fontsize=5.2)
for i,(_,r) in enumerate(prod_sel.iterrows()):
    if r.Hero:
        ax.get_yticklabels()[i].set_color("#c0392b")
        ax.get_yticklabels()[i].set_fontweight("bold")
cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.01)
cbar.set_label("Module completeness", fontsize=8)
ax.set_title("Fig. 4 | DRAM metabolic module completeness across 111 MAGs "
             "(hero clade highlighted)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig4_DRAM_metabolism_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig4_DRAM_metabolism_heatmap.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig4_DRAM_metabolism_heatmap")

# focused comparison: hero vs rest for MICP-critical modules
crit = [c for c in curated if re.search(r"urease|carbonic anhyd|Nitrogen|antiport|Cobalamin|B12|CAZy|Carbohydrate", c, re.I)]
if crit:
    hero_mean = prod_sel.loc[prod_sel.Hero, crit].mean()
    rest_mean = prod_sel.loc[~prod_sel.Hero, crit].mean()
    fig, ax = plt.subplots(figsize=(9, max(4, len(crit)*0.35)))
    y=np.arange(len(crit)); w=0.4
    ax.barh(y-w/2, hero_mean.values, w, label=f"Hero (n={prod_sel.Hero.sum()})",
            color="#c0392b", edgecolor="black", lw=0.4)
    ax.barh(y+w/2, rest_mean.values, w, label=f"Rest (n={(~prod_sel.Hero).sum()})",
            color="#7f8c8d", edgecolor="black", lw=0.4)
    ax.set_yticks(y); ax.set_yticklabels(crit, fontsize=7)
    ax.set_xlabel("Mean module completeness")
    ax.set_title("MICP-critical modules: hero vs rest", fontsize=10)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/Fig4b_DRAM_HeroVsRest.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Fig4b_DRAM_HeroVsRest.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[ok] Fig4b")
