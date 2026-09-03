"""Suppl Fig S12 - predicted minimum doubling time by group (gRodon2, analysis A7).

One block: predicted minimum doubling time for the MICP-complete group and for the rest,
box with every MAG overlaid as a point, group n in the tick label and the two-sided
Mann-Whitney P as a stat annotation.  The y axis is logarithmic: the rest group spans
0.4-17.4 h, and on a linear axis its upper outliers compress every MICP-complete point onto
the baseline.

Provenance: Table S16, which holds the 85 MAGs that passed the >= 10 ribosomal-protein
anchor filter.  Medians and P are recomputed; they are asserted against the values quoted in
the figure legend (median 1.06 h vs 1.10 h, P = 0.58).  S13 and S16 are absent from the
table - the filter removed them - and this figure therefore shows 4 of the 6 MICP-complete
MAGs; that is a property of the source, not a drop made here.

Colour meanings: blue = Sphingobacterium MICP-complete MAG, orange = Pseudomonas_E
MICP-complete MAG, grey = one of the remaining MAGs.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import REST, SPHINGO, PSEUDO, TEXT, AXIS, FS_BODY, FS_STAT, HEROES, hero_col

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
RNG = np.random.default_rng(0)

# ------------------------------------------------------------------ data
df = pd.read_csv(SUPP / "Table_S16_gRodon_growth_rates_per_MAG.csv")
is_hero = df.group == "MICP_complete"
hero = df[is_hero]
rest = df[~is_hero]
assert set(hero.MAG) <= set(HEROES)
p = mannwhitneyu(hero.d_hours, rest.d_hours, alternative="two-sided").pvalue
assert abs(hero.d_hours.median() - 1.06) < 0.005
assert abs(rest.d_hours.median() - 1.10) < 0.005
assert abs(p - 0.58) < 0.005
dropped = [m for m in HEROES if m not in set(hero.MAG)]
assert dropped, "the legend records which MICP-complete MAGs the filter removed"

# ------------------------------------------------------------------ page
H = 64.0
fig, ax_mm, text_mm, letter = st.page(H)
letter(11, 5, "A")

ax = ax_mm(30, 13, 104, 38)
bp = ax.boxplot([hero.d_hours.values, rest.d_hours.values], positions=[0, 1], widths=0.5,
                showfliers=False, patch_artist=True)
for box in bp["boxes"]:
    box.set(facecolor="white", edgecolor=AXIS, linewidth=0.7)
for part in ("whiskers", "caps", "medians"):
    for ln in bp[part]:
        ln.set(color=AXIS, linewidth=0.7)

jit = RNG.uniform(-0.13, 0.13, len(rest))
ax.scatter(1 + jit, rest.d_hours, s=5, c=REST, alpha=0.85, linewidths=0, zorder=3)
hx = np.linspace(-0.18, 0.18, len(hero))
ax.scatter(hx, hero.d_hours, s=12, c=[hero_col(m) for m in hero.MAG], zorder=4,
           linewidths=0)

ax.set_xticks([0, 1])
ax.set_xticklabels([f"MICP-complete\nn = {len(hero)}", f"rest\nn = {len(rest)}"])
ax.set_ylabel("predicted minimum doubling time (h)")
ax.set_xlim(-0.6, 1.6)
ax.set_yscale("log")
ax.set_yticks([0.5, 1, 2, 5, 10, 20])
ax.get_yaxis().set_major_formatter(ScalarFormatter())
lo = float(min(hero.d_hours.min(), rest.d_hours.min()))
hi = float(max(hero.d_hours.max(), rest.d_hours.max()))
ax.set_ylim(lo * 0.75, hi * 2.2)
ax.text(0.5, 0.97, f"P = {p:.2f}", ha="center", va="top", fontsize=FS_STAT,
        color=TEXT, transform=ax.transAxes)
st.style_axis(ax)

fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, color=SPHINGO,
                           label="Sphingobacterium"),
                    Line2D([], [], marker="o", ls="", ms=3, color=PSEUDO,
                           label="Pseudomonas_E"),
                    Line2D([], [], marker="o", ls="", ms=3, color=REST,
                           label="other MAG")],
           loc="lower center", ncol=3, frameon=False, fontsize=FS_BODY,
           bbox_to_anchor=(0.5, 0.015), handletextpad=0.4, columnspacing=1.6)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S12")
