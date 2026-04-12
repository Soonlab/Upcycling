#!/usr/bin/env python3
"""
Regenerate supplementary trait-module heatmaps (S1, S2, S3, S5) with
genus-aggregated rows for readability, matching Figure 4a style.
"""
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE="/data/data/Upcycling"
OUT =f"{BASE}/SUBMISSION/Figures_supplementary"
os.makedirs(OUT, exist_ok=True)

gtdb = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
gtdb["Genus"] = gtdb["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "Unclassified")
tax = gtdb.set_index("user_genome")["Genus"]
LINEAGE={"C22","M1","S13","S16","S23","S26"}
COLOR_HIGHLIGHT="#c0392b"

counts=pd.read_csv(f"{BASE}/research/extra/gene_category_counts.csv").set_index("Sample")
subcols=[c for c in counts.columns if "::" in c]
norm=counts[subcols].div(counts["CDS_total"], axis=0)*1000.0
norm["Genus"]=[tax.get(m,"Unclassified") for m in norm.index]

CAT_CONF = {
    "S1": ("Biofilm_EPS",    "Biofilm / EPS gene modules"),
    "S2": ("Ammonia_N",      "Ammonia & nitrogen-assimilation modules"),
    "S3": ("Alkaline_Osmo",  "Alkaline / osmotic-stress tolerance modules"),
    "S5": ("MetalResist_AMR","Heavy-metal & antibiotic-resistance modules"),
}

def plot_heatmap(cat, title, fig_id):
    cols=[c for c in norm.columns if c.startswith(f"{cat}::")]
    sub=norm[cols+["Genus"]].copy()
    agg=sub.groupby("Genus")[cols].mean()
    nmag=sub.groupby("Genus").size()
    agg=agg.loc[nmag[nmag>=2].index]
    # flag genera containing MICP-complete members
    contains_mc={g: any(tax.get(i)==g for i in LINEAGE) for g in agg.index}
    # priority sort
    agg["_pri"]=[0 if contains_mc[g] else 1 for g in agg.index]
    agg=agg.sort_values("_pri").drop(columns="_pri")

    ylbls=[f"{g}  (n = {nmag[g]})" for g in agg.index]
    clean_xcols=[c.replace(f"{cat}::","") for c in cols]
    fig, ax=plt.subplots(figsize=(max(7, len(cols)*0.8), max(4, len(agg)*0.5+1.5)))
    # log10(x+0.1) so low-abundance modules remain visible
    mat = np.log10(agg.values.astype(float) + 0.1)
    im=ax.imshow(mat, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(clean_xcols, rotation=40, ha="right", fontsize=9)
    ax.set_yticks(range(len(agg))); ax.set_yticklabels(ylbls, fontsize=10)
    for i,g in enumerate(agg.index):
        if contains_mc[g]:
            ax.get_yticklabels()[i].set_color(COLOR_HIGHLIGHT)
            ax.get_yticklabels()[i].set_fontweight("bold")
    ax.set_xticks(np.arange(len(cols))-0.5, minor=True)
    ax.set_yticks(np.arange(len(agg))-0.5, minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar=fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("log₁₀( gene hits per 10³ CDS + 0.1 )", fontsize=8)
    ax.set_title(f"Figure {fig_id}. {title} — genus-aggregated ({len(agg)} genera, n ≥ 2).",
                 fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_{fig_id}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_{fig_id}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Figure {fig_id} done")

for fig_id,(cat,title) in CAT_CONF.items():
    plot_heatmap(cat, title, fig_id)
