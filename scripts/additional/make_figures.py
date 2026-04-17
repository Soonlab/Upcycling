#!/usr/bin/env python3
"""Generate figures for all additional analyses.

Outputs to /data/data/Upcycling/research/additional/figures/
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from scipy.stats import mannwhitneyu

ROOT = "/data/data/Upcycling/research/additional"
OUT = f"{ROOT}/figures"
os.makedirs(OUT, exist_ok=True)
mpl.rcParams["figure.dpi"] = 150
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False

HERO = ["S13", "S16", "S23", "C22", "M1", "S26"]
COLOR_HERO = "#c0392b"
COLOR_REST = "#7f8c8d"

def fig_A1_biosafety():
    df = pd.read_csv(f"{ROOT}/A1_biosafety/biosafety_counts_per_MAG.csv", index_col=0)
    df = df.reindex(columns=["card","vfdb","resfinder","plasmidfinder"]).fillna(0).astype(int)
    df["group"] = ["MICP-complete" if m in HERO else "rest" for m in df.index]

    fig, axes = plt.subplots(1, 4, figsize=(12, 4), sharey=False)
    for ax, db in zip(axes, ["card","vfdb","resfinder","plasmidfinder"]):
        hero = df[df.group=="MICP-complete"][db].values
        rest = df[df.group=="rest"][db].values
        ax.boxplot([hero, rest], labels=["MICP\ncomplete\n(n=6)","rest\n(n=105)"],
                   widths=0.5, patch_artist=True,
                   boxprops=dict(facecolor="#ecf0f1"), medianprops=dict(color="black"))
        # scatter
        for i, vals, c in [(0, hero, COLOR_HERO), (1, rest, COLOR_REST)]:
            xs = np.random.normal(i+1, 0.06, size=len(vals))
            ax.scatter(xs, vals, color=c, alpha=0.6, s=18, edgecolor="black", linewidth=0.3)
        ax.set_title(db.upper(), fontsize=11, fontweight="bold")
        ax.set_ylabel("n hits per MAG")
    plt.suptitle("A1. Biosafety panel — hits per MAG by database",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(f"{OUT}/A1_biosafety.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/A1_biosafety.pdf", bbox_inches="tight")
    plt.close()
    print(f"[fig] A1 saved")

def fig_A2_active_site():
    df = pd.read_csv(f"{ROOT}/A2_structure/UreC_active_site_residues.csv")
    sites = df["site"].tolist()
    mags = HERO
    M = np.zeros((len(mags), len(sites)), dtype=int)
    for i, m in enumerate(mags):
        for j, s in enumerate(sites):
            v = df.iloc[j][m]
            exp = df.iloc[j]["expected"]
            M[i, j] = 1 if v == exp else 0
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.imshow(M, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    for i in range(len(mags)):
        for j in range(len(sites)):
            v = df.iloc[j][mags[i]]
            ax.text(j, i, v, ha="center", va="center", fontsize=11, fontweight="bold")
    ax.set_xticks(range(len(sites)))
    ax.set_xticklabels([f"{s}\n({r})" for s, r in zip(sites, df["expected"])], fontsize=10)
    ax.set_yticks(range(len(mags)))
    ax.set_yticklabels(mags)
    ax.set_title("A2. UreC active-site residue conservation vs S. pasteurii P41020\n(42/42 = 100%)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig.savefig(f"{OUT}/A2_ureC_active_site.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/A2_ureC_active_site.pdf", bbox_inches="tight")
    plt.close()
    print(f"[fig] A2 saved")

def fig_A5_alkaliphile():
    df = pd.read_csv(f"{ROOT}/A5_alkaliphile/alkaliphile_signature_per_MAG.csv")
    metrics = [("Mrp_count","Mrp antiporter copies\n(11.7× enrichment, p=5.3e-4)"),
               ("Nha_count","Nha antiporter copies"),
               ("pI_median","proteome pI (median)"),
               ("pI_acidic_frac","acidic protein fraction (pI<5)")]
    fig, axes = plt.subplots(1, 4, figsize=(12, 4))
    for ax, (m, title) in zip(axes, metrics):
        hero_v = df[df.group=="MICP_complete"][m].dropna().values
        rest_v = df[df.group=="rest"][m].dropna().values
        try:
            u, p = mannwhitneyu(hero_v, rest_v, alternative="two-sided")
        except Exception:
            p = np.nan
        bp = ax.boxplot([hero_v, rest_v], labels=["MICP\ncomplete","rest"],
                        widths=0.5, patch_artist=True,
                        boxprops=dict(facecolor="#ecf0f1"),
                        medianprops=dict(color="black"))
        for i, vals, c in [(0, hero_v, COLOR_HERO), (1, rest_v, COLOR_REST)]:
            xs = np.random.normal(i+1, 0.06, size=len(vals))
            ax.scatter(xs, vals, color=c, alpha=0.7, s=18, edgecolor="black", linewidth=0.3)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(m.replace("_", " "))
        if not np.isnan(p):
            ax.text(0.5, 0.97, f"MWU p = {p:.2e}", transform=ax.transAxes,
                    ha="center", va="top", fontsize=9, style="italic")
    plt.suptitle("A5. Alkaliphile signatures: MICP-complete vs rest",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(f"{OUT}/A5_alkaliphile.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/A5_alkaliphile.pdf", bbox_inches="tight")
    plt.close()
    print(f"[fig] A5 saved")

def fig_A7_gRodon():
    if not os.path.exists(f"{ROOT}/A7_grodon/gRodon_growth_rates_per_MAG.csv"):
        print("[fig] A7 missing, skip"); return
    df = pd.read_csv(f"{ROOT}/A7_grodon/gRodon_growth_rates_per_MAG.csv")
    fig, ax = plt.subplots(figsize=(5, 5))
    hero_v = df[df.group=="MICP_complete"]["d_hours"].values
    rest_v = df[df.group=="rest"]["d_hours"].values
    bp = ax.boxplot([hero_v, rest_v],
                    labels=[f"MICP-complete\n(n={len(hero_v)})", f"rest\n(n={len(rest_v)})"],
                    widths=0.5, patch_artist=True,
                    boxprops=dict(facecolor="#ecf0f1"),
                    medianprops=dict(color="black"))
    for i, vals, c in [(0, hero_v, COLOR_HERO), (1, rest_v, COLOR_REST)]:
        xs = np.random.normal(i+1, 0.06, size=len(vals))
        ax.scatter(xs, vals, color=c, alpha=0.7, s=30, edgecolor="black", linewidth=0.3)
    ax.set_ylabel("predicted doubling time (h)")
    ax.set_title("A7. Minimum doubling time (gRodon)\nNo growth-rate penalty for MICP trait",
                 fontsize=11)
    try:
        u, p = mannwhitneyu(hero_v, rest_v, alternative="two-sided")
        ax.text(0.5, 0.97, f"MWU p = {p:.2f}", transform=ax.transAxes,
                ha="center", va="top", fontsize=10, style="italic")
    except Exception:
        pass
    plt.tight_layout()
    fig.savefig(f"{OUT}/A7_growth_rate.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/A7_growth_rate.pdf", bbox_inches="tight")
    plt.close()
    print(f"[fig] A7 saved")

def fig_B_rarity():
    df = pd.read_csv(f"{ROOT}/B_rarity_screen/MICP_rarity_summary_mgnify.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel A: bar of % MICP-complete per biome + pooled
    ax = axes[0]
    total = df.n_species_clusters.sum()
    pooled_micp = df.n_MICP_gene_complete.sum()
    cats = df.catalog.tolist() + ["pooled\n(4 biomes)"]
    pcts = df.pct_MICP_gene_complete.tolist() + [100*pooled_micp/total]
    ns   = df.n_species_clusters.tolist() + [total]
    colors_ = ["#3498db","#9b59b6","#e67e22","#16a085","#2c3e50"]
    bars = ax.bar(cats, pcts, color=colors_, edgecolor="black")
    ax.set_ylabel("% species clusters with MICP gene-complete profile")
    ax.set_title("B. MICP-complete rarity across MGnify livestock MAGs",
                 fontsize=11, fontweight="bold")
    for bar, pct, n in zip(bars, pcts, ns):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
                f"{pct:.2f}%\n(n={n})", ha="center", va="bottom", fontsize=9)
    # our MAGs reference line
    ax.axhline(100, color=COLOR_HERO, ls="--", lw=1.2)
    ax.text(len(cats)-0.5, 95, "this study:\n6/6 = 100%",
            ha="right", va="top", fontsize=9, color=COLOR_HERO, fontweight="bold")
    ax.set_ylim(0, 12)

    # Panel B: single-contig ureC+CA %
    ax = axes[1]
    pcts2 = df.pct_MICP_single_contig.tolist() + [100*df.n_MICP_single_contig_ureC_CA.sum()/total]
    bars = ax.bar(cats, pcts2, color=colors_, edgecolor="black")
    ax.set_ylabel("% with single-contig ureC+CA (operon-like)")
    ax.set_title("Single-contig operon architecture — even rarer",
                 fontsize=11, fontweight="bold")
    for bar, pct in zip(bars, pcts2):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
                f"{pct:.2f}%", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, 8)

    plt.tight_layout()
    fig.savefig(f"{OUT}/B_rarity_mgnify.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/B_rarity_mgnify.pdf", bbox_inches="tight")
    plt.close()
    print(f"[fig] B saved")

def fig_A6_stoich():
    df = pd.read_csv(f"{ROOT}/A6_metabolic/stoichiometry_per_MAG.csv")
    # select per-gene columns for hero MAGs
    cols = ["ureA","ureB","ureC","ureD_H","ureE","ureF","ureG",
            "CA_generic","Ca_transporter","Ca_ATPase","Na_H_antiporter_Mrp"]
    hero_df = df[df.group=="MICP_complete"].set_index("MAG")[cols].reindex(HERO)
    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(hero_df.values, cmap="YlGnBu", aspect="auto", vmin=0, vmax=hero_df.values.max())
    for i in range(len(HERO)):
        for j in range(len(cols)):
            v = hero_df.values[i, j]
            ax.text(j, i, str(v), ha="center", va="center", fontsize=10,
                    color="black" if v < hero_df.values.max()/2 else "white")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(len(HERO)))
    ax.set_yticklabels(HERO, fontsize=11)
    ax.set_title("A6. Per-gene dosage for MICP pathway in 6 MICP-complete MAGs\n"
                 "(urea + 2H₂O + Ca²⁺ → 2NH₄⁺ + CaCO₃ + OH⁻)",
                 fontsize=11, fontweight="bold")
    plt.colorbar(im, ax=ax, label="gene copies")
    plt.tight_layout()
    fig.savefig(f"{OUT}/A6_stoichiometry.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/A6_stoichiometry.pdf", bbox_inches="tight")
    plt.close()
    print(f"[fig] A6 saved")

def fig_A3b_pseudomonas_rarity():
    fp = f"{ROOT}/A3_pseudomonas_ani/pseudomonas_e_single_contig.csv"
    if not os.path.exists(fp):
        print("[fig] A3b missing, skip"); return
    df = pd.read_csv(fp)
    n = len(df)
    cats = ["UreC\npresent", "any CA\npresent", "UreC+UreB\nsame contig", "UreC+CA\nsame contig\n(M1/S26 architecture)"]
    vals = [df.has_UreC.sum(), df.has_CA_any.sum(),
            df.ureC_and_ureB_single_contig.sum(),
            df.ureC_and_CA_single_contig.sum()]
    pcts = [100*v/n for v in vals]
    fig, ax = plt.subplots(figsize=(8, 5))
    colors_ = ["#3498db","#3498db","#27ae60","#e74c3c"]
    bars = ax.bar(cats, pcts, color=colors_, edgecolor="black")
    for bar, pct, v in zip(bars, pcts, vals):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1.5,
                f"{pct:.1f}%\n(n={v})", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("% of 146 Pseudomonas_E reference genomes")
    ax.set_ylim(0, 115)
    ax.set_title(f"A3b. MICP feature prevalence in Pseudomonas_E genus (n={n})\n"
                 "Gene presence is near-universal; co-localized architecture rarer",
                 fontsize=11, fontweight="bold")
    ax.axhline(100, color="grey", ls="--", lw=0.5)
    plt.tight_layout()
    fig.savefig(f"{OUT}/A3b_pseudomonas_rarity.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/A3b_pseudomonas_rarity.pdf", bbox_inches="tight")
    plt.close()
    print("[fig] A3b saved")

def fig_A4_genomad():
    # only run if summary csv exists
    fp = f"{ROOT}/A4_genomad/genomad_summary.csv"
    if not os.path.exists(fp):
        print("[fig] A4 not ready yet, skip")
        return
    df = pd.read_csv(fp)
    fig, ax = plt.subplots(figsize=(6, 5))
    h = df[df.group=="MICP_complete"]
    r = df[df.group=="rest"]
    bp = ax.boxplot([h.n_plasmid_contigs, r.n_plasmid_contigs,
                     h.n_virus_contigs, r.n_virus_contigs],
                    tick_labels=["MICP\nplasmid","rest\nplasmid","MICP\nvirus","rest\nvirus"],
                    widths=0.5, patch_artist=True,
                    boxprops=dict(facecolor="#ecf0f1"),
                    medianprops=dict(color="black"))
    for i, vals, c in [(0, h.n_plasmid_contigs, COLOR_HERO),
                       (1, r.n_plasmid_contigs, COLOR_REST),
                       (2, h.n_virus_contigs, COLOR_HERO),
                       (3, r.n_virus_contigs, COLOR_REST)]:
        xs = np.random.normal(i+1, 0.06, size=len(vals))
        ax.scatter(xs, vals, color=c, alpha=0.6, s=18, edgecolor="black", linewidth=0.3)
    ax.set_ylabel("n MGE-predicted contigs per MAG")
    ax.set_title("A4. geNomad plasmid/virus calls per MAG", fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig.savefig(f"{OUT}/A4_genomad.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{OUT}/A4_genomad.pdf", bbox_inches="tight")
    plt.close()
    print("[fig] A4 saved")

if __name__ == "__main__":
    fig_A1_biosafety()
    fig_A2_active_site()
    fig_A5_alkaliphile()
    fig_A7_gRodon()
    fig_B_rarity()
    fig_A6_stoich()
    fig_A3b_pseudomonas_rarity()
    fig_A4_genomad()
    print("\nAll figures in", OUT)
