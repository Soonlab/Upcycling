"""Suppl Fig S11 - alkaliphile signature per MAG (analysis A5).

Four blocks, one per feature: multi-subunit Mrp antiporter copy count, Nha copy count,
proteome pI median and proteome acidic-pI fraction.  Within a block the MICP-complete group
and the remaining 105 MAGs are shown as a box with every MAG overlaid as a point, and the
two-sided Mann-Whitney P is a stat annotation bound to the block.

Provenance: Table S15a.  Group means and the Mann-Whitney P are recomputed here; the Mrp
fold difference the manuscript quotes (11.7x) is recomputed and asserted against the ratio
of the recomputed group means, and the two-sided P is the value the manuscript reports.

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
from _style import REST, SPHINGO, PSEUDO, TEXT, AXIS, FS_BODY, FS_STAT, HEROES, hero_col

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

FEATURES = [("Mrp_count", "Mrp copies"),
            ("Nha_count", "Nha copies"),
            ("pI_median", "proteome pI (median)"),
            ("pI_acidic_frac", "acidic-pI fraction")]
RNG = np.random.default_rng(0)

# ------------------------------------------------------------------ data
df = pd.read_csv(SUPP / "Table_S15a_alkaliphile_signature_per_MAG.csv")
is_hero = df.group == "MICP_complete"
assert sorted(df.MAG[is_hero]) == sorted(HEROES)
n_hero, n_rest = int(is_hero.sum()), int((~is_hero).sum())
mrp_fold = df.loc[is_hero, "Mrp_count"].mean() / df.loc[~is_hero, "Mrp_count"].mean()
assert abs(mrp_fold - 11.7) < 0.05, mrp_fold          # the value quoted in the manuscript

# ------------------------------------------------------------------ page
H = 76.0
fig, ax_mm, text_mm, letter = st.page(H)
letter(11, 5, "A")

W, GAP, LEFT, TOP, HGT = 28.0, 10.0, 20.0, 14.0, 44.0
for k, (col, ylab) in enumerate(FEATURES):
    ax = ax_mm(LEFT + k * (W + GAP), TOP, W, HGT)
    hero_v = df.loc[is_hero, col].values
    rest_v = df.loc[~is_hero, col].values
    bp = ax.boxplot([hero_v, rest_v], positions=[0, 1], widths=0.55,
                    showfliers=False, patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor="white", edgecolor=AXIS, linewidth=0.7)
    for part in ("whiskers", "caps", "medians"):
        for ln in bp[part]:
            ln.set(color=AXIS, linewidth=0.7)
    for xi, vals, cols in ((0, hero_v, [hero_col(m) for m in df.MAG[is_hero]]),
                           (1, rest_v, [REST] * n_rest)):
        jit = RNG.uniform(-0.16, 0.16, len(vals))
        ax.scatter(xi + jit, vals, s=5, c=cols, alpha=0.85, linewidths=0, zorder=3)
    p = mannwhitneyu(hero_v, rest_v, alternative="two-sided").pvalue
    lo = min(float(np.min(hero_v)), float(np.min(rest_v)))
    hi = max(float(np.max(hero_v)), float(np.max(rest_v)))
    span = hi - lo
    ax.set_ylim(lo - 0.08 * span, hi + 0.26 * span)
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"MICP-complete\nn = {n_hero}", f"rest\nn = {n_rest}"])
    txt = f"P = {p:.1e}" if p < 0.01 else f"P = {p:.2f}"
    ax.text(0.5, hi + 0.18 * span, txt, ha="center", va="center", fontsize=FS_STAT,
            color=TEXT)
    ax.set_ylabel(ylab)
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
st.save(fig, OUT, "Fig_S11")
