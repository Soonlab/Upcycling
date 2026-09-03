"""Suppl Fig S8 - biosafety screen across the 111 MAGs (analysis A1).

One block per abricate database (CARD, VFDB, ResFinder, PlasmidFinder); within a block the
per-MAG hit count is shown for the MICP-complete group and for the remaining 105 MAGs, as a
box with every MAG overlaid as a point.  MICP-complete points are coloured by lineage.

Provenance.  CARD / VFDB / ResFinder counts come from Table S11.  Table S11 carries NO
PlasmidFinder column, although the A1 run produced one: the counts are recomputed here from
the raw abricate output `research/additional/A1_biosafety/combined/plasmidfinder_all.tsv`
(one row per hit, the MAG taken from the #FILE path), which is the file the aggregation
script `aggregate_biosafety.py` reads.  The recomputed CARD/VFDB/ResFinder counts from the
same raw files are asserted against Table S11 before PlasmidFinder is added, so the fourth
block rests on a source that reproduces the three published ones.
Mann-Whitney P per database is recomputed (two-sided) and shown as a stat annotation.

Colour meanings: blue = Sphingobacterium MICP-complete MAG, orange = Pseudomonas_E
MICP-complete MAG, grey = one of the remaining 105 MAGs.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import (REST, SPHINGO, PSEUDO, GREY, TEXT, AXIS, LIGHT,
                    FS_BODY, FS_STAT, HEROES, hero_col)

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
A1 = Path("/data/data/Upcycling/research/additional/A1_biosafety/combined")

DBS = ["card", "vfdb", "resfinder", "plasmidfinder"]
DB_LAB = {"card": "CARD", "vfdb": "VFDB", "resfinder": "ResFinder",
          "plasmidfinder": "PlasmidFinder"}
RNG = np.random.default_rng(0)     # reproducible jitter (style constant)

# ------------------------------------------------------------------ data
s11 = pd.read_csv(SUPP / "Table_S11_biosafety_counts_per_MAG.csv").set_index("MAG")
mags = list(s11.index)


def raw_counts(db):
    """Per-MAG abricate hit count, recomputed from the combined raw output."""
    df = pd.read_csv(A1 / f"{db}_all.tsv", sep="\t", low_memory=False)
    mag = df["#FILE"].str.extract(r"bakta_results/([^/]+)/")[0]
    return mag.value_counts().reindex(mags).fillna(0).astype(int)


counts = pd.DataFrame({db: raw_counts(db) for db in DBS}, index=mags)
for db in ["card", "vfdb", "resfinder"]:            # the three published columns
    assert (counts[db] == s11[db]).all(), db
counts["group"] = s11["group"]
is_hero = counts.group == "MICP_complete"
assert sorted(counts.index[is_hero]) == sorted(HEROES)
n_hero, n_rest = int(is_hero.sum()), int((~is_hero).sum())

# ------------------------------------------------------------------ page
H = 76.0
fig, ax_mm, text_mm, letter = st.page(H)
letter(11, 5, "A")

W, GAP, LEFT, TOP, HGT = 33.0, 10.0, 22.0, 14.0, 44.0
for k, db in enumerate(DBS):
    ax = ax_mm(LEFT + k * (W + GAP), TOP, W, HGT)
    hero_v = counts.loc[is_hero, db].values
    rest_v = counts.loc[~is_hero, db].values
    bp = ax.boxplot([hero_v, rest_v], positions=[0, 1], widths=0.55,
                    showfliers=False, patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor="white", edgecolor=AXIS, linewidth=0.7)
    for part in ("whiskers", "caps", "medians"):
        for ln in bp[part]:
            ln.set(color=AXIS, linewidth=0.7)
    for xi, vals, cols in ((0, hero_v, [hero_col(m) for m in counts.index[is_hero]]),
                           (1, rest_v, [REST] * n_rest)):
        jit = RNG.uniform(-0.16, 0.16, len(vals))
        ax.scatter(xi + jit, vals, s=5, c=cols, alpha=0.85, linewidths=0, zorder=3)
    p = mannwhitneyu(hero_v, rest_v, alternative="two-sided").pvalue
    top = max(float(np.max(hero_v)), float(np.max(rest_v)))
    ax.set_ylim(-0.06 * max(top, 1), max(top, 1) * 1.30)
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"MICP-complete\nn = {n_hero}", f"rest\nn = {n_rest}"])
    ax.text(0.5, max(top, 1) * 1.20, f"P = {p:.2f}", ha="center", va="center",
            fontsize=FS_STAT, color=TEXT)
    ax.text(0.5, max(top, 1) * 1.30, DB_LAB[db], ha="center", va="bottom",
            fontsize=FS_BODY, color=TEXT)
    if k == 0:
        ax.set_ylabel("abricate hits per MAG")
    st.style_axis(ax)

fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, color=SPHINGO,
                           label="Sphingobacterium"),
                    Line2D([], [], marker="o", ls="", ms=3, color=PSEUDO,
                           label="Pseudomonas_E"),
                    Line2D([], [], marker="o", ls="", ms=3, color=REST,
                           label="other MAG")],
           loc="lower center", ncol=3, frameon=False, fontsize=FS_BODY,
           bbox_to_anchor=(0.5, 0.005), handletextpad=0.4, columnspacing=1.6)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S8")
