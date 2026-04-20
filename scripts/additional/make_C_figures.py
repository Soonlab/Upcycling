#!/usr/bin/env python3
"""Generate publication figures for C1-C6 analyses.
Tolerant to missing inputs (skips silently)."""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/data/data/Upcycling/research/additional")
FIG  = Path("/home/soon/Upcycling_repo/figures/additional")
FIG.mkdir(parents=True, exist_ok=True)

def safe_read_csv(p):
    try:
        if Path(p).exists() and Path(p).stat().st_size > 0:
            return pd.read_csv(p)
    except Exception:
        pass
    return None

def save(fig, name):
    fig.savefig(FIG/f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG/f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"saved {name}")

def boxplot_hero_vs_rest(df, value_col, title, fname, ylabel=None):
    if df is None or value_col not in df.columns: return
    h = df[df.is_hero][value_col].astype(float).dropna().values
    r = df[~df.is_hero][value_col].astype(float).dropna().values
    if len(h)<1 or len(r)<1: return
    fig, ax = plt.subplots(figsize=(4.4,4.0))
    bp = ax.boxplot([r, h], labels=["Rest (n=%d)"%len(r), "MICP-complete (n=%d)"%len(h)],
                    patch_artist=True, widths=0.55)
    colors = ["#bdbdbd","#0a7d6e"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    # scatter
    ax.scatter(np.random.normal(1, 0.05, len(r)), r, s=10, c="#444", alpha=0.5)
    ax.scatter(np.random.normal(2, 0.05, len(h)), h, s=24, c="#c53e1f", edgecolors="k", lw=0.5, zorder=3)
    ax.set_ylabel(ylabel or value_col)
    ax.set_title(title)
    save(fig, fname)

# --- C1 antiSMASH ---
df = safe_read_csv(ROOT/"C1_antismash"/"antismash_per_MAG.csv")
boxplot_hero_vs_rest(df, "n_regions",
                     "C1 antiSMASH BGC regions per MAG",
                     "C1_antismash_regions",
                     "BGC regions per MAG")
if df is not None:
    bgc_cols = [c for c in df.columns if c.startswith("BGC_")]
    if bgc_cols:
        # stacked bar of BGC class composition for hero vs rest
        h_means = df[df.is_hero][bgc_cols].mean()
        r_means = df[~df.is_hero][bgc_cols].mean()
        top = (h_means+r_means).sort_values(ascending=False).head(10).index
        fig, ax = plt.subplots(figsize=(7,4))
        x = np.arange(len(top)); w=0.4
        ax.bar(x-w/2, r_means[top], w, label="Rest", color="#bdbdbd")
        ax.bar(x+w/2, h_means[top], w, label="MICP-complete", color="#0a7d6e")
        ax.set_xticks(x); ax.set_xticklabels([t.replace("BGC_","") for t in top], rotation=35, ha="right")
        ax.set_ylabel("Mean regions / MAG"); ax.legend()
        ax.set_title("C1 Top BGC classes")
        save(fig, "C1_antismash_topclasses")

# --- C2 Defense ---
df = safe_read_csv(ROOT/"C2_defense"/"defense_per_MAG.csv")
boxplot_hero_vs_rest(df, "n_defense_systems",
                     "C2 DefenseFinder systems per MAG",
                     "C2_defense_systems",
                     "Defense systems per MAG")
boxplot_hero_vs_rest(df, "n_crispr_arrays",
                     "C2 CRISPR arrays (minced) per MAG",
                     "C2_crispr_arrays",
                     "CRISPR arrays per MAG")

# --- C3 Codon usage ---
df = safe_read_csv(ROOT/"C3_dnds_codon"/"codon_usage_per_MAG.csv")
boxplot_hero_vs_rest(df, "ENC", "C3 Effective Number of Codons", "C3_ENC", "ENC")
boxplot_hero_vs_rest(df, "GC3_pct", "C3 GC content at 3rd codon position", "C3_GC3", "GC3 (%)")

# C3 dN/dS bar
codeml = safe_read_csv(ROOT/"C3_dnds_codon"/"codeml_M0_summary.csv")
if codeml is not None and "omega_M0" in codeml.columns:
    cdf = codeml.dropna(subset=["omega_M0"])
    if len(cdf):
        fig, ax = plt.subplots(figsize=(4.4,3.4))
        ax.bar(cdf.gene, cdf.omega_M0, color="#0a7d6e", alpha=0.85, edgecolor="k")
        ax.axhline(1.0, ls="--", c="#c53e1f", lw=1.2, label="ω = 1 (neutral)")
        ax.set_ylabel("ω (dN/dS)")
        ax.set_title("C3 codeml M0 single-omega per urease gene")
        ax.legend()
        save(fig, "C3_dnds_M0")

# --- C4 ESMFold vs 4CEU ---
df = safe_read_csv(ROOT/"C4_esmfold"/"esmfold_vs_4CEU.csv")
if df is not None and len(df):
    fig, axes = plt.subplots(1,2, figsize=(8,3.5))
    axes[0].bar(df.MAG, df.RMSD_to_4CEU, color="#0a7d6e", edgecolor="k")
    axes[0].set_ylabel("RMSD vs 4CEU (Å)"); axes[0].set_title("C4 RMSD")
    axes[0].tick_params(axis="x", rotation=30)
    axes[1].bar(df.MAG, df.TM_score_norm_4CEU, color="#0a7d6e", edgecolor="k")
    axes[1].axhline(0.5, ls="--", c="#c53e1f", lw=1.0, label="TM=0.5 same fold")
    axes[1].set_ylabel("TM-score (norm by 4CEU)"); axes[1].set_title("C4 TM-score")
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].legend()
    fig.suptitle("ESMFold UreC vs S. pasteurii UreC (PDB 4CEU)")
    fig.tight_layout()
    save(fig, "C4_esmfold_vs_4CEU")

# --- C5 panMICP env ---
sk = safe_read_csv(ROOT/"C5_panMICP_env"/"skani_hero_vs_refs.tsv")
if sk is None:
    # try TSV
    p = ROOT/"C5_panMICP_env"/"skani_hero_vs_refs.tsv"
    if p.exists():
        try: sk = pd.read_csv(p, sep="\t")
        except Exception: sk = None
if sk is not None and len(sk):
    fig, ax = plt.subplots(figsize=(7,4.5))
    cols = sk.columns.str.lower().tolist()
    qcol = "Query_file" if "Query_file" in sk.columns else sk.columns[0]
    rcol = "Ref_file" if "Ref_file" in sk.columns else sk.columns[1]
    acol = "ANI" if "ANI" in sk.columns else [c for c in sk.columns if "ani" in c.lower()][0]
    sk["query"] = sk[qcol].astype(str).apply(lambda s: Path(s).stem)
    sk["ref"] = sk[rcol].astype(str).apply(lambda s: Path(s).stem)
    pivot = sk.pivot_table(index="ref", columns="query", values=acol, aggfunc="max").fillna(0)
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis", vmin=70, vmax=100)
    ax.set_xticks(range(len(pivot.columns))); ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(pivot.index))); ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_title("C5 ANI: hero MAGs vs curated MICP references (per environment)")
    fig.colorbar(im, ax=ax, label="ANI (%)")
    fig.tight_layout()
    save(fig, "C5_panMICP_ANI_heatmap")

# --- C6 abundance proxy ---
df = safe_read_csv(ROOT/"C6_abundance_proxy"/"abundance_proxy_per_MAG.csv")
if df is not None:
    boxplot_hero_vs_rest(df, "length_weighted_cov",
                         "C6 In-situ depth proxy (length-weighted SPAdes cov)",
                         "C6_abundance_proxy",
                         "Length-weighted contig coverage")
    # Per-source breakdown
    fig, ax = plt.subplots(figsize=(5,3.6))
    grp = df.groupby("source")["length_weighted_cov"].apply(list).to_dict()
    keys = sorted(grp.keys())
    bp = ax.boxplot([grp[k] for k in keys], labels=keys, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#bdbdbd"); patch.set_alpha(0.7)
    # Overlay heroes
    for i,k in enumerate(keys, start=1):
        sub = df[(df.source==k) & df.is_hero]
        ax.scatter(np.random.normal(i, 0.05, len(sub)), sub.length_weighted_cov,
                   c="#c53e1f", s=40, edgecolors="k", zorder=3, label="hero" if i==1 else None)
    ax.set_ylabel("length-weighted cov")
    ax.set_title("C6 Per-source depth proxy")
    ax.legend()
    save(fig, "C6_abundance_proxy_by_source")

print("done figures")
