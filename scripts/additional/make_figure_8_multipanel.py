#!/usr/bin/env python3
"""
Figure 8 — Main figure combining structural and external validation of the
two MICP-complete lineages.

Panels:
(a) A2 : UreC active-site residue conservation heatmap
         (6 MICP-complete MAGs x 7 canonical catalytic residues,
          colour = match/mismatch vs S. pasteurii P41020 / PDB 4CEU).
(b) C4 : ESMFold UreC backbone TM-score (normalised to 4CEU chain C, 569 aa)
         and all-residue RMSD (A) for the same 6 MAGs, dashed line at
         TM = 0.5 same-fold threshold.
(c) B  : MGnify 7,599-cluster external rarity: % MICP gene-complete and
         % single-contig ureC+CA per biome + pooled, horizontal line at
         the 6/6 (100%) detection rate of the present study.
(d) C1 : antiSMASH 7 class-level BGC enrichment - hero mean vs rest mean
         per class for the significantly enriched / depleted classes,
         T3PKS and RRE-containing highlighted.

Inputs are the committed supplementary tables; no raw FASTA / PDB required.
Output: Figure_8.png (300 dpi), Figure_8.pdf (vector).
"""
from __future__ import annotations

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
SUPP = "/data/data/Upcycling/SUBMISSION/Supplementary_tables"
OUT_DIR = "/data/data/Upcycling/research/figures_final"
os.makedirs(OUT_DIR, exist_ok=True)

A2_PATH = os.path.join(SUPP, "Table_S12_UreC_active_site_residues.csv")
C4_PATH = os.path.join(SUPP, "Table_S22_ureC_vs_4CEU_tm.csv")
B_PATH = os.path.join(SUPP, "Table_S14a_mgnify_catalog_summary.csv")
C1_PATH = os.path.join(SUPP, "Table_S23b_antismash_hero_vs_rest.csv")

HEROES = ["S13", "S16", "S23", "C22", "M1", "S26"]
SPHINGO = {"S13", "S16", "S23", "C22"}
PSEUDO = {"M1", "S26"}

SPHINGO_COLOR = "#2b6cb0"  # blue
PSEUDO_COLOR = "#c05621"   # orange
MATCH_COLOR = "#2f855a"    # green
MISMATCH_COLOR = "#c53030"  # red
ACCENT = "#6b46c1"         # purple for key enrichments
NEUTRAL_GREY = "#4a5568"

# -----------------------------------------------------------------------------
# Panel (a) - A2 UreC active-site heatmap
# -----------------------------------------------------------------------------

def panel_a_active_site(ax: plt.Axes) -> None:
    df = pd.read_csv(A2_PATH)
    residues = df["site"].tolist()
    expected = df["expected"].tolist()
    mat = np.zeros((len(HEROES), len(residues)), dtype=int)
    for i, mag in enumerate(HEROES):
        for j, exp in enumerate(expected):
            observed = df.iloc[j][mag]
            mat[i, j] = 1 if observed == exp else 0

    cmap = ListedColormap([MISMATCH_COLOR, MATCH_COLOR])
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)

    for i, mag in enumerate(HEROES):
        for j, exp in enumerate(expected):
            observed = df.iloc[j][mag]
            ax.text(j, i, observed, ha="center", va="center",
                    color="white", fontsize=8, fontweight="bold")

    ax.set_xticks(range(len(residues)))
    ax.set_xticklabels([f"{r}\n({e})" for r, e in zip(residues, expected)],
                       fontsize=8)
    ax.set_yticks(range(len(HEROES)))
    # colour tick labels by lineage
    ax.set_yticklabels(HEROES, fontsize=9)
    for tick, mag in zip(ax.get_yticklabels(), HEROES):
        tick.set_color(SPHINGO_COLOR if mag in SPHINGO else PSEUDO_COLOR)
        tick.set_fontweight("bold")

    ax.set_xlabel(
        "UreC active-site residue (P41020 / 4CEU numbering)",
        fontsize=9)
    ax.set_title(
        "(a) Active-site conservation — 42/42 matches",
        fontsize=11, fontweight="bold", loc="left")


# -----------------------------------------------------------------------------
# Panel (b) - C4 ESMFold TM-score + RMSD
# -----------------------------------------------------------------------------

def panel_b_esmfold(ax: plt.Axes) -> None:
    df = pd.read_csv(C4_PATH)
    df["MAG_short"] = df["MAG"].str.replace("_UreC", "", regex=False)
    # keep the hero order
    df = df.set_index("MAG_short").loc[HEROES].reset_index()

    x = np.arange(len(HEROES))
    bar_colors = [SPHINGO_COLOR if m in SPHINGO else PSEUDO_COLOR
                  for m in HEROES]

    width = 0.38
    tm_bars = ax.bar(x - width / 2, df["tm_norm_ref"], width,
                     color=bar_colors, edgecolor="black",
                     label="TM-score (norm. 4CEU chain C)")
    ax.axhline(0.5, color="black", linestyle="--", linewidth=1,
               label="Same-fold threshold (TM = 0.5)")
    ax.set_ylabel("Backbone TM-score", fontsize=9)
    ax.set_ylim(0.0, 0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(HEROES, fontsize=9)
    for tick, mag in zip(ax.get_xticklabels(), HEROES):
        tick.set_color(SPHINGO_COLOR if mag in SPHINGO else PSEUDO_COLOR)
        tick.set_fontweight("bold")

    for bar, val in zip(tm_bars, df["tm_norm_ref"]):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.01,
                f"{val:.3f}", ha="center", va="bottom", fontsize=7.5)

    ax2 = ax.twinx()
    rmsd_bars = ax2.bar(x + width / 2, df["rmsd"], width,
                        color="white", edgecolor="black", hatch="////",
                        label="Backbone RMSD (Å)")
    ax2.set_ylabel("Backbone RMSD (Å)", fontsize=9)
    ax2.set_ylim(0, 6)
    for bar, val in zip(rmsd_bars, df["rmsd"]):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 0.05,
                 f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)

    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2,
              loc="upper right", fontsize=7.5, frameon=False)

    ax.set_title(
        "(b) ESMFold UreC vs 4CEU — 6/6 TM > 0.5",
        fontsize=11, fontweight="bold", loc="left")


# -----------------------------------------------------------------------------
# Panel (c) - B MGnify 7,599-cluster rarity
# -----------------------------------------------------------------------------

def panel_c_mgnify(ax: plt.Axes) -> None:
    df = pd.read_csv(B_PATH)
    rename = {
        "cow-rumen": "cow\nrumen",
        "sheep-rumen": "sheep\nrumen",
        "pig-gut": "pig\ngut",
        "chicken-gut": "chicken\ngut",
    }
    df["label"] = df["catalog"].map(rename)
    pooled_n = int(df["n_species_clusters"].sum())
    pooled_mgc = int(df["n_MICP_gene_complete"].sum())
    pooled_sc = int(df["n_MICP_single_contig_ureC_CA"].sum())
    pooled_pct_mgc = 100 * pooled_mgc / pooled_n
    pooled_pct_sc = 100 * pooled_sc / pooled_n
    pooled_row = {
        "catalog": "pooled",
        "n_species_clusters": pooled_n,
        "n_MICP_gene_complete": pooled_mgc,
        "n_MICP_single_contig_ureC_CA": pooled_sc,
        "pct_MICP_gene_complete": pooled_pct_mgc,
        "pct_MICP_single_contig": pooled_pct_sc,
        "label": f"pooled\n(n={pooled_n:,})",
    }
    df = pd.concat([df, pd.DataFrame([pooled_row])], ignore_index=True)

    x = np.arange(len(df))
    width = 0.38
    bars_gc = ax.bar(x - width / 2, df["pct_MICP_gene_complete"], width,
                     color="#3182ce", edgecolor="black",
                     label="% MICP gene-complete")
    bars_sc = ax.bar(x + width / 2, df["pct_MICP_single_contig"], width,
                     color="#63b3ed", edgecolor="black",
                     label="% single-contig ureC + CA")

    for bar, val in zip(bars_gc, df["pct_MICP_gene_complete"]):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.3,
                f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)
    for bar, val in zip(bars_sc, df["pct_MICP_single_contig"]):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.3,
                f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)

    ax.axhline(100, color=ACCENT, linestyle="--", linewidth=1.2)
    ax.text(len(df) - 0.5, 100, "  present study 6/6 = 100%",
            color=ACCENT, fontsize=8, va="center", ha="right")

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], fontsize=8)
    ax.set_ylabel("Species-cluster representatives (%)", fontsize=9)
    ax.set_ylim(0, 110)
    ax.set_title(
        "(c) MGnify rarity — 7,599 clusters, 3.07% MICP-complete",
        fontsize=11, fontweight="bold", loc="left")
    ax.legend(loc="upper left", fontsize=7.5, frameon=False)


# -----------------------------------------------------------------------------
# Panel (d) - C1 antiSMASH class-level enrichment
# -----------------------------------------------------------------------------

def panel_d_antismash(ax: plt.Axes) -> None:
    df = pd.read_csv(C1_PATH)
    # Focus on: T3PKS, RRE, arylpolyene, terpene, NRP-metallophore, betalactone,
    #           HCN, NAGGN, RiPP-like, NRPS, T1PKS (if present)
    focus = [
        "BGC_T3PKS",
        "BGC_RRE-containing",
        "BGC_arylpolyene",
        "BGC_terpene",
        "BGC_NAGGN",
        "BGC_RiPP-like",
        "BGC_NRPS",
        "BGC_betalactone",
        "BGC_hydrogen-cyanide",
        "BGC_NRP-metallophore",
    ]
    sub = df[df["metric"].isin(focus)].copy()
    # Safety: if a focus class is absent, just skip it
    sub["class_label"] = sub["metric"].str.replace("BGC_", "", regex=False)
    sub = sub.set_index("class_label").reindex(
        [c.replace("BGC_", "") for c in focus if
         c in sub["metric"].values or c.replace("BGC_", "") in sub.index]
    ).dropna(how="all").reset_index()

    # Sort by hero_mean descending
    sub = sub.sort_values("hero_mean", ascending=True).reset_index(drop=True)

    y = np.arange(len(sub))
    height = 0.4

    def _star(p):
        try:
            pv = float(p)
        except Exception:
            return ""
        if pv < 1e-9:
            return "***"
        if pv < 1e-4:
            return "**"
        if pv < 0.05:
            return "*"
        return ""

    bars_h = ax.barh(y - height / 2, sub["hero_mean"], height,
                     color=ACCENT, edgecolor="black",
                     label="MICP-complete mean (n=6)")
    bars_r = ax.barh(y + height / 2, sub["rest_mean"], height,
                     color="#a0aec0", edgecolor="black",
                     label="Rest mean (n=105)")

    for yi, (cl, row) in enumerate(sub.iterrows()):
        star = _star(row["MWU_p"])
        if star:
            ax.text(max(row["hero_mean"], row["rest_mean"]) + 0.03, yi,
                    f"{star}  p={row['MWU_p']:.1e}",
                    va="center", fontsize=7.5, color=ACCENT)

    ax.set_yticks(y)
    ax.set_yticklabels(sub["class_label"], fontsize=8)
    # Bold the two key classes
    for tick_lbl in ax.get_yticklabels():
        if tick_lbl.get_text() in {"T3PKS", "RRE-containing"}:
            tick_lbl.set_color(ACCENT)
            tick_lbl.set_fontweight("bold")

    ax.set_xlabel("BGC regions per MAG (mean)", fontsize=9)
    ax.set_xlim(0, max(sub["hero_mean"].max(), sub["rest_mean"].max()) * 1.7)
    ax.set_title(
        "(d) antiSMASH 7 — T3PKS 23× + RRE 8× Sphingobacterium fingerprint",
        fontsize=11, fontweight="bold", loc="left")
    ax.legend(loc="lower right", fontsize=7.5, frameon=False)


# -----------------------------------------------------------------------------
# Main figure assembly
# -----------------------------------------------------------------------------

def main() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(16, 11), constrained_layout=True)

    panel_a_active_site(axes[0, 0])
    panel_b_esmfold(axes[0, 1])
    panel_c_mgnify(axes[1, 0])
    panel_d_antismash(axes[1, 1])

    # Shared lineage legend (bottom)
    sphingo_patch = mpatches.Patch(color=SPHINGO_COLOR,
                                   label="Sphingobacterium (S13/S16/S23/C22)")
    pseudo_patch = mpatches.Patch(color=PSEUDO_COLOR,
                                  label="Pseudomonas_E (M1/S26)")
    fig.legend(handles=[sphingo_patch, pseudo_patch],
               loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=10)

    png_path = os.path.join(OUT_DIR, "Figure_8.png")
    pdf_path = os.path.join(OUT_DIR, "Figure_8.pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
