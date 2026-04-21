#!/usr/bin/env python3
"""C4 composite figure: bar of TM-score + RMSD per hero vs PDB 4CEU."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT = Path("/data/data/Upcycling/research/additional/C4_esmfold")
FIG = Path("/home/soon/Upcycling_repo/figures/additional")
FIG.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(ROOT/"ureC_vs_4CEU_tm.csv")
df["hero"] = df["MAG"].str.replace("_UreC","")
df = df.sort_values("hero").reset_index(drop=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

colors_sph = ["#0a7d6e"]*4 + ["#b8860b"]*2  # Sphingobacterium vs Pseudomonas_E
bars1 = ax1.bar(df["hero"], df["tm_norm_ref"], color=colors_sph, alpha=0.85, edgecolor="k")
ax1.axhline(0.5, ls="--", color="k", lw=1, label="TM = 0.5 (same fold)")
ax1.set_ylabel("TM-score (normalised to 4CEU-C, 569 aa)")
ax1.set_title("C4 ESMFold UreC backbone similarity to PDB 4CEU")
ax1.set_ylim(0, 1.0)
for b, v in zip(bars1, df["tm_norm_ref"]):
    ax1.text(b.get_x() + b.get_width()/2, v + 0.01, f"{v:.2f}",
             ha="center", va="bottom", fontsize=9)
ax1.legend(loc="upper right", fontsize=8)

bars2 = ax2.bar(df["hero"], df["rmsd"], color=colors_sph, alpha=0.85, edgecolor="k")
ax2.set_ylabel(r"Backbone RMSD to 4CEU-C (Å)")
ax2.set_title("C4 ESMFold vs 4CEU-C RMSD")
for b, v in zip(bars2, df["rmsd"]):
    ax2.text(b.get_x() + b.get_width()/2, v + 0.1, f"{v:.1f}",
             ha="center", va="bottom", fontsize=9)

fig.suptitle("C4 UreC structural comparison vs Sporosarcina pasteurii urease α (PDB 4CEU, chain C)",
             y=1.03, fontsize=11)
fig.tight_layout()
fig.savefig(FIG/"C4_esmfold_superposition.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG/"C4_esmfold_superposition.pdf", bbox_inches="tight")
print(f"saved {FIG}/C4_esmfold_superposition")
print(df[["hero","pred_len","tm_norm_ref","rmsd"]].to_string(index=False))
