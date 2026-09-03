"""Suppl Fig S15 - geNomad mobile-element calls and the urease cross-check (analysis A4).

A  plasmid- and virus-flagged contigs per MAG, MICP-complete group vs the rest
B  per-MICP-complete-MAG cross-check: how many contigs carry the urease core and the
   carbonic anhydrase, how many contigs geNomad flagged, and how many of the urease / CA
   contigs are among the flagged ones

Provenance: Tables S17a and S17b.  Group means and the two-sided Mann-Whitney P in A are
recomputed from S17a.  The contig counts in B are recomputed by splitting the comma-joined
contig lists of S17b; the two contamination columns are read as stored.  The flag columns
urease_on_plasmid / urease_on_virus / CA_on_virus are empty for every MAG in the source and
carry no information, so they are not drawn; CA_on_plasmid is non-empty for S26 alone and is
represented by that MAG's CA contamination count.

Colour meanings: blue = Sphingobacterium MICP-complete MAG, orange = Pseudomonas_E
MICP-complete MAG, grey = one of the remaining 105 MAGs (A); coral = a non-zero
mobile-element contamination count (B).  Nothing else is coloured in B.
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
from _style import (HERO, REST, SPHINGO, PSEUDO, TEXT, AXIS, LIGHT,
                    FS_BODY, FS_STAT, HEROES, hero_col)

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
RNG = np.random.default_rng(0)

METRICS = [("n_plasmid_contigs", "plasmid-flagged contigs"),
           ("n_virus_contigs", "virus-flagged contigs")]

# ------------------------------------------------------------------ data
a = pd.read_csv(SUPP / "Table_S17a_genomad_summary_per_MAG.csv")
is_hero = a.group == "MICP_complete"
assert sorted(a.MAG[is_hero]) == sorted(HEROES)
n_hero, n_rest = int(is_hero.sum()), int((~is_hero).sum())

b = pd.read_csv(SUPP / "Table_S17b_ureCah_vs_MGE_overlap.csv").set_index("MAG").loc[HEROES]


def n_contigs(cell):
    return 0 if pd.isna(cell) else len(str(cell).split(","))


tab_cols = [("urease contigs", [n_contigs(v) for v in b.urease_core_contigs]),
            ("CA contigs", [n_contigs(v) for v in b.CA_contigs]),
            ("plasmid contigs", list(b.n_plasmid_contigs)),
            ("virus contigs", list(b.n_virus_contigs)),
            ("urease on MGE", list(b.urease_core_MGE_contamination)),
            ("CA on MGE", list(b.CA_MGE_contamination))]
# the per-MAG plasmid/virus counts of the two tables must agree
assert (b.n_plasmid_contigs.values ==
        a.set_index("MAG").loc[HEROES].n_plasmid_contigs.values).all()

# ------------------------------------------------------------------ page
H = 80.0
fig, ax_mm, text_mm, letter = st.page(H)

# ---- A -------------------------------------------------------------------
letter(11, 5, "A")
for k, (col, ylab) in enumerate(METRICS):
    ax = ax_mm(22 + k * 42.0, 13, 30, 45)
    hero_v = a.loc[is_hero, col].values
    rest_v = a.loc[~is_hero, col].values
    bp = ax.boxplot([hero_v, rest_v], positions=[0, 1], widths=0.55, showfliers=False,
                    patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor="white", edgecolor=AXIS, linewidth=0.7)
    for part in ("whiskers", "caps", "medians"):
        for ln in bp[part]:
            ln.set(color=AXIS, linewidth=0.7)
    for xi, vals, cols in ((0, hero_v, [hero_col(m) for m in a.MAG[is_hero]]),
                           (1, rest_v, [REST] * n_rest)):
        jit = RNG.uniform(-0.16, 0.16, len(vals))
        ax.scatter(xi + jit, vals, s=5, c=cols, alpha=0.85, linewidths=0, zorder=3)
    p = mannwhitneyu(hero_v, rest_v, alternative="two-sided").pvalue
    hi = float(max(hero_v.max(), rest_v.max()))
    ax.set_ylim(-0.04 * hi, hi * 1.22)
    ax.set_yticks([t for t in ax.get_yticks() if t >= 0])   # drop out-of-view ticks
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"MICP-complete\nn = {n_hero}", f"rest\nn = {n_rest}"])
    ax.text(0.5, hi * 1.15, f"P = {p:.2f}", ha="center", va="center", fontsize=FS_STAT,
            color=TEXT)
    ax.set_ylabel(ylab)
    st.style_axis(ax)

fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, color=SPHINGO,
                           label="Sphingobacterium"),
                    Line2D([], [], marker="o", ls="", ms=3, color=PSEUDO,
                           label="Pseudomonas_E"),
                    Line2D([], [], marker="o", ls="", ms=3, color=REST,
                           label="other MAG")],
           loc="lower left", ncol=1, frameon=False, fontsize=FS_BODY,
           bbox_to_anchor=(0.10, 0.02), handletextpad=0.4)

# ---- B: borderless mini table --------------------------------------------
letter(103, 5, "B")
axB = ax_mm(116, 13, 58, 45)
axB.set_xlim(0, len(tab_cols))
axB.set_ylim(len(HEROES), -1.4)
axB.axis("off")
for j, (head, vals) in enumerate(tab_cols):
    axB.text(j + 0.5, -1.0, head.replace(" ", "\n", 1), ha="center", va="center",
             fontsize=FS_STAT, color=TEXT)
    for i, v in enumerate(vals):
        axB.text(j + 0.5, i, str(int(v)), ha="center", va="center", fontsize=FS_BODY,
                 color=HERO if (int(v) > 0 and "on MGE" in head) else TEXT)
for i, mag in enumerate(HEROES):
    axB.text(-0.15, i, mag, ha="right", va="center", fontsize=FS_BODY,
             color=hero_col(mag))
    axB.axhline(i + 0.5, color=LIGHT, lw=0.5, xmin=0.0, xmax=1.0)
axB.axhline(-0.5, color=AXIS, lw=0.7)
for t in axB.texts:
    t.set_clip_on(False)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S15")
