#!/usr/bin/env python3
"""Generate extra C3 figures (codeml M0 omega + yn00 hero/rest pairwise distribution)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT = Path("/data/data/Upcycling/research/additional/C3_dnds_codon")
FIG = Path("/home/soon/Upcycling_repo/figures/additional")
FIG.mkdir(parents=True, exist_ok=True)

GOOD = "#0a7d6e"
BAD = "#c53e1f"
GREY = "#bdbdbd"

def save(fig, name):
    fig.savefig(FIG/f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG/f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("saved", name)

# codeml M0
df_m = pd.read_csv(ROOT/"codeml_M0_summary.csv")
if not df_m.empty:
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    genes = df_m["gene"].tolist()
    omegas = df_m["omega_M0"].astype(float).values
    bars = ax.bar(genes, omegas, color=GOOD, alpha=0.85, edgecolor="k")
    ax.axhline(1.0, ls="--", color="k", lw=1, label="ω = 1 (neutral)")
    for b, v in zip(bars, omegas):
        ax.text(b.get_x() + b.get_width()/2, v + 0.005, f"{v:.3f}",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("codeml M0 ω (dN/dS)")
    ax.set_title("C3 codeml M0: urease gene family selection pressure\n(6 MICP-complete + 12 non-hero reps, FastTree topology)")
    ax.set_ylim(0, max(0.15, max(omegas)*1.3))
    ax.legend(loc="upper right", fontsize=8)
    save(fig, "C3_codeml_omega")

# yn00 hero vs rest boxplot
pw = pd.read_csv(ROOT/"yn00_pairwise.csv")
if not pw.empty:
    HEROES = {"S13","S16","S23","C22","M1","S26"}
    pw = pw[(pw["omega"] < 99) & (pw["omega"] > 0) & (pw["dS"] > 0.01)].copy()
    def bucket(r):
        if r.a in HEROES and r.b in HEROES: return "hero-hero"
        if r.a in HEROES or r.b in HEROES: return "hero-rest"
        return "rest-rest"
    pw["bucket"] = pw.apply(bucket, axis=1)
    genes = ["ureA","ureB","ureC","ureG"]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), sharey=False)
    for ax, g in zip(axes, genes):
        gdf = pw[pw.gene == g]
        if gdf.empty: continue
        buckets = ["hero-hero","hero-rest","rest-rest"]
        data = [gdf[gdf.bucket == b]["omega"].values for b in buckets]
        labels = [f"{b}\n(n={len(d)})" for b, d in zip(buckets, data)]
        bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.55, showfliers=False)
        for patch, c in zip(bp["boxes"], [BAD, GREY, GOOD]):
            patch.set_facecolor(c); patch.set_alpha(0.7)
        for i, d in enumerate(data, start=1):
            if len(d) == 0: continue
            x = np.random.normal(i, 0.05, len(d))
            ax.scatter(x, d, s=10, alpha=0.35, c="k")
        ax.axhline(1.0, ls="--", color="k", lw=0.8)
        ax.set_title(g)
        ax.set_ylabel("yn00 ω" if g == "ureA" else "")
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle("C3 yn00 pairwise ω across hero/rest partitions (ureG hero-hero elevation, MWU p = 7.7e-8)", y=1.02)
    fig.tight_layout()
    save(fig, "C3_yn00_hero_vs_rest")
