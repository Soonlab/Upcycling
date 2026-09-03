"""Main Fig 7 - pan-genome and trait-module ordination by waste source.

One composed page, 180 mm wide, two panels:
  A  Jaccard PCoA of the Panaroo gene presence/absence matrix, stored coordinates
  B  z-scored Euclidean PCoA of the 38 per-10^3-CDS trait modules, RECOMPUTED here

Sources
  Table_S9a_PCoA_coordinates.csv    PC1-PC3, waste source, genus, MICP flag (panel A)
  Table_S9b_PERMANOVA_global.csv    pan-genome pseudo-F and p for source and genus, and the
                                    PC1 / PC2 variance explained (panel A axis titles)
  Table_S2b_trait_module_per1kCDS.csv   the 38 trait modules per MAG (panel B input)

Panel B recomputation
  Table S9a stores the PAN-GENOME ordination only; the trait ordination of the shipped
  panel b has no stored coordinates anywhere in the repository.  It is recomputed here with
  the recipe of `research/revision/05_PCoA_by_source.py` and the same z-scoring the shipped
  builder used (`scripts/revision/regenerate_panelized_figures.py`): per-module z-score,
  Euclidean distance, classical PCoA, PERMANOVA on waste source with 999 permutations.
  The recomputation reproduces the published statistic exactly - pseudo-F = 2.7056 against
  the 2.71 of the figure legend, p = 0.001, PC1 16.47 % and PC2 10.38 % against the 16.5 %
  and 10.4 % of the shipped axis titles - and the script asserts each of those.

Colour meanings on this page:
  cattle brown / swine pink / sheep green / poultry purple   waste source, both panels
  coral ring + coral label                                   a MICP-complete MAG
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore
from skbio.stats.distance import DistanceMatrix, permanova
from skbio.stats.ordination import pcoa

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, TEXT, SOURCE, FS_BODY, FS_STAT, HEROES

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

SEED = 0            # PERMANOVA permutation seed, fixed so the p value is reproducible
N_PERM = 999
CLUSTER_PX = 28     # marker separation below which two MAG labels would collide
LINE_PT = 7         # vertical step used to fan out a cluster of labels

# ------------------------------------------------------------------ data
coords = pd.read_csv(SUPP / "Table_S9a_PCoA_coordinates.csv").set_index("MAG")
glob = pd.read_csv(SUPP / "Table_S9b_PERMANOVA_global.csv").iloc[0]
mags = list(coords.index)
assert set(coords.index[coords.Hero]) == set(HEROES)

# panel B - recompute the trait ordination
traits = pd.read_csv(SUPP / "Table_S2b_trait_module_per1kCDS.csv").set_index("Sample")
mod_cols = [c for c in traits.columns if "::" in c]
X = traits.loc[mags, mod_cols]
Z = np.nan_to_num(zscore(X.values, axis=0, ddof=0, nan_policy="omit"), nan=0.0)
dm_t = DistanceMatrix(squareform(pdist(Z, metric="euclidean")), ids=mags)
np.random.seed(SEED)
perm_t = permanova(dm_t, grouping=coords.loc[mags, "Source"], permutations=N_PERM)
ord_t = pcoa(dm_t, number_of_dimensions=2)
var_t = ord_t.proportion_explained * 100
coords_t = ord_t.samples.copy()
coords_t.index = mags

F_T, P_T = float(perm_t["test statistic"]), float(perm_t["p-value"])
# the values the shipped figure and the legend report, reproduced to their printed precision
assert round(F_T, 2) == 2.71, F_T
assert P_T <= 0.001, P_T
assert round(var_t.iloc[0], 1) == 16.5, var_t.iloc[0]
assert round(var_t.iloc[1], 1) == 10.4, var_t.iloc[1]

SRC_ORDER = ["Cattle", "Swine", "Sheep", "Poultry"]
assert set(coords.Source) == set(SRC_ORDER)
COL = {s: SOURCE[s.lower()] for s in SRC_ORDER}
n_src = coords.Source.value_counts()

# ------------------------------------------------------------------ layout
TOP = 10.0
W_P, H_P = 66.0, 62.0
X1, X2 = 14.0, 104.0
H = TOP + H_P + 24.0

fig, ax_mm, text_mm, letter = st.page(H)


def scatter(ax, xy, xlab, ylab):
    for s in SRC_ORDER:
        idx = coords.index[coords.Source == s]
        ax.scatter(xy.loc[idx, "PC1"], xy.loc[idx, "PC2"], s=9, color=COL[s],
                   linewidth=0, alpha=0.9, zorder=2, label=f"{s} n = {n_src[s]}")
    h = list(HEROES)
    ax.scatter(xy.loc[h, "PC1"], xy.loc[h, "PC2"], s=34, facecolor="none",
               edgecolor=HERO, linewidth=0.8, zorder=4)
    # several MICP-complete MAGs project onto nearly the same point, so their labels are
    # fanned out: MAGs whose markers fall within CLUSTER_PX of one another get successive
    # vertical offsets instead of printing on top of each other
    ax.figure.canvas.draw()
    pts = {m: ax.transData.transform((xy.loc[m, "PC1"], xy.loc[m, "PC2"])) for m in h}
    placed = []
    for m in sorted(h, key=lambda k: (-pts[k][1], pts[k][0])):
        k = sum(1 for q in placed
                if abs(pts[q][0] - pts[m][0]) < CLUSTER_PX
                and abs(pts[q][1] - pts[m][1]) < CLUSTER_PX)
        placed.append(m)
        ax.annotate(m, (xy.loc[m, "PC1"], xy.loc[m, "PC2"]),
                    xytext=(4, 2 - k * LINE_PT), textcoords="offset points",
                    fontsize=FS_STAT, color=HERO, zorder=5)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    # prune the extreme ticks: the first x tick and the first y tick otherwise print on
    # top of each other in the corner where the two spines meet
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    st.style_axis(ax)


axA = ax_mm(X1, TOP, W_P, H_P)
scatter(axA, coords, f"PC1  {glob.PC1_var:.1f} %", f"PC2  {glob.PC2_var:.1f} %")

axB = ax_mm(X2, TOP, W_P, H_P)
scatter(axB, coords_t, f"PC1  {var_t.iloc[0]:.1f} %", f"PC2  {var_t.iloc[1]:.1f} %")

# PERMANOVA results as a stat block bound to each panel
STAT_T = TOP + 1.0
text_mm(X1 + 3.0, STAT_T, "PERMANOVA", fontsize=FS_STAT, color=TEXT)
for k, (lab, f, p) in enumerate((("source", glob.pseudo_F_source, glob.p_source),
                                 ("genus", glob.pseudo_F_genus, glob.p_genus))):
    text_mm(X1 + 3.0, STAT_T + 3.4 + k * 3.4,
            f"{lab}  F = {f:.2f}  p = {p:.3f}", fontsize=FS_STAT, color=TEXT)
text_mm(X2 + 3.0, STAT_T, "PERMANOVA", fontsize=FS_STAT, color=TEXT)
text_mm(X2 + 3.0, STAT_T + 3.4, f"source  F = {F_T:.2f}  p = {P_T:.3f}",
        fontsize=FS_STAT, color=TEXT)

handles, labels = axA.get_legend_handles_labels()
handles.append(Line2D([0], [0], marker="o", ls="", ms=5, mfc="none", mec=HERO, mew=0.8))
labels.append(f"MICP-complete n = {len(HEROES)}")
fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 6 / H), ncol=5,
           fontsize=FS_BODY, handletextpad=0.3, columnspacing=1.4, borderpad=0.2)

letter(4, 4, "A")
letter(X2 - 12.0, 4, "B")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig7")
