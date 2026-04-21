#!/usr/bin/env python3
"""C1 antiSMASH figure: (a) n_regions hero vs rest, (b) top BGC class enrichments."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT = Path("/data/data/Upcycling/research/additional/C1_antismash")
FIG = Path("/home/soon/Upcycling_repo/figures/additional")
FIG.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(ROOT/"antismash_per_MAG.csv")
stats = pd.read_csv(ROOT/"antismash_hero_vs_rest.csv")
stats = stats.sort_values("MWU_p").reset_index(drop=True)
HEROES = {"S13","S16","S23","C22","M1","S26"}

GOOD = "#0a7d6e"; GREY = "#bdbdbd"; BAD = "#c53e1f"

fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))

# (a) n_regions boxplot + hero scatter
ax = axes[0]
h = df[df.is_hero]["n_regions"].values
r = df[~df.is_hero]["n_regions"].values
bp = ax.boxplot([r, h], tick_labels=[f"Rest (n={len(r)})", f"MICP-complete (n={len(h)})"],
                patch_artist=True, widths=0.55)
for patch, c in zip(bp["boxes"], [GREY, GOOD]):
    patch.set_facecolor(c); patch.set_alpha(0.7)
# scatter with hero labels
for i, (lbl, vals) in enumerate(zip(["rest","hero"], [r, h]), start=1):
    x = np.random.normal(i, 0.06, len(vals))
    ax.scatter(x, vals, s=20 if lbl=="rest" else 55,
               c="k" if lbl=="rest" else BAD,
               alpha=0.5 if lbl=="rest" else 0.95,
               edgecolors="none" if lbl=="rest" else "k", zorder=3)
# Label heroes
hero_rows = df[df.is_hero].sort_values("n_regions").reset_index(drop=True)
for idx, row in hero_rows.iterrows():
    ax.annotate(row["MAG"], (2 + 0.08, row["n_regions"]),
                fontsize=8, va="center")
ax.set_ylabel("BGC regions per MAG (antiSMASH 7)")
ax.set_title("(a) Total BGC count — heroes not BGC-rich or BGC-poor\n"
             "(hero mean 6.2, rest mean 6.9; MWU n.s.)")

# (b) top enriched BGC classes
ax = axes[1]
top = stats.head(8).copy()
top["class"] = top["metric"].str.replace("BGC_", "")
y = np.arange(len(top))
w_h = top["hero_mean"].values
w_r = top["rest_mean"].values
ax.barh(y - 0.18, w_h, height=0.36, color=GOOD, alpha=0.85, label="MICP-complete")
ax.barh(y + 0.18, w_r, height=0.36, color=GREY, alpha=0.85, label="Rest")
ax.set_yticks(y)
ax.set_yticklabels(top["class"], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("Mean BGC copies per MAG")
# annotate p-values
for yi, p in zip(y, top["MWU_p"].values):
    ax.text(max(top[["hero_mean","rest_mean"]].max().max(), 0.1) * 1.02, yi,
            f"p = {p:.2g}", va="center", fontsize=8)
ax.set_title("(b) BGC-class enrichment in MICP-complete lineages\n"
             "T3PKS + RRE signature characterises the 4 Sphingobacterium heroes")
ax.legend(loc="lower right", fontsize=9)

fig.suptitle("C1 antiSMASH secondary-metabolism profile (111 MAGs)", y=1.02, fontsize=11)
fig.tight_layout()
fig.savefig(FIG/"C1_antismash.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG/"C1_antismash.pdf", bbox_inches="tight")
print("saved C1_antismash")
