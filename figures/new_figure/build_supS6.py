"""Suppl Fig S6 - pairwise skani ANI within the six MICP-complete MAGs.

Source: SUBMISSION/Supplementary_tables/Table_S4a_skani_ANI_matrix_111MAGs.csv,
        subset to the six MICP-complete MAGs.

skani reports no value for a pair below its alignment threshold; the shipped matrix
stores those as 0.0.  A stored 0.0 is therefore rendered as an empty cell with an em
dash, not as "0 % identity" - the pair is simply too divergent for skani to align.
The matrix is checked for symmetry before drawing.

Colour meanings on this page:
  white -> green   ANI (%) between two MICP-complete MAGs (the only encoding)
  blue             Sphingobacterium MAG label (S13, S16, S23, C22)
  orange           Pseudomonas_E MAG label (M1, S26)
  light grey       a pair skani could not align
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
NO_ALIGNMENT = 0.0

ani = pd.read_csv(SUPP / "Table_S4a_skani_ANI_matrix_111MAGs.csv", index_col=0)
mat = ani.loc[st.HEROES, st.HEROES]
assert np.allclose(mat.values, mat.values.T)
assert np.allclose(np.diag(mat.values), 100.0)

vals = mat.values.astype(float)
shown = np.where(vals > NO_ALIGNMENT, vals, np.nan)
off = shown[~np.eye(len(st.HEROES), dtype=bool)]
vmin = float(np.nanmin(off))

CELL = 14.0
TOP, LEFT = 22.0, 30.0
SIDE = CELL * len(st.HEROES)
H = TOP + SIDE + 20.0

fig, ax_mm, text_mm, letter = st.page(H)
ax = ax_mm(LEFT, TOP, SIDE, SIDE)
cax = ax_mm(LEFT + SIDE + 8.0, TOP, 3.0, SIDE * 0.6)
letter(4, 4, "A")

cmap = st.seq_cmap()
cmap.set_bad(color="#F4F4F4")
im = ax.imshow(shown, cmap=cmap, vmin=vmin, vmax=100.0, aspect="equal")

ax.set_xticks(range(len(st.HEROES)))
ax.set_xticklabels(st.HEROES, fontsize=st.FS_BODY)
ax.xaxis.set_ticks_position("top")
ax.set_yticks(range(len(st.HEROES)))
ax.set_yticklabels(st.HEROES, fontsize=st.FS_BODY)
ax.tick_params(length=0, pad=3)
for labs in (ax.get_xticklabels(), ax.get_yticklabels()):
    for lab, mag in zip(labs, st.HEROES):
        lab.set_color(st.hero_col(mag))
        lab.set_fontweight("bold")
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xticks(np.arange(-0.5, len(st.HEROES), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(st.HEROES), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.2)
ax.tick_params(which="minor", length=0)

for i in range(len(st.HEROES)):
    for j in range(len(st.HEROES)):
        v = vals[i, j]
        if v > NO_ALIGNMENT:
            frac = (v - vmin) / (100.0 - vmin)
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=st.FS_STAT,
                    color="white" if frac > 0.55 else st.TEXT)
        else:
            ax.text(j, i, "—", ha="center", va="center", fontsize=st.FS_STAT,
                    color=st.GREY)

grad = np.linspace(100.0, vmin, 128).reshape(-1, 1)
cax.imshow(grad, cmap=st.seq_cmap(), vmin=vmin, vmax=100.0, aspect="auto",
           extent=[0, 1, vmin, 100.0])
cax.set_xticks([])
cax.yaxis.set_ticks_position("right")
cax.yaxis.set_major_locator(MaxNLocator(nbins=4))
cax.tick_params(labelsize=st.FS_STAT, length=2, pad=1)
for s in cax.spines.values():
    s.set_linewidth(0.6)
    s.set_color(st.AXIS)
text_mm(LEFT + SIDE + 18.0, TOP + SIDE * 0.3, "ANI (%)", rotation=90,
        ha="center", va="center", fontsize=st.FS_BODY)

handles = [Line2D([], [], color=st.SPHINGO, lw=4, label="Sphingobacterium"),
           Line2D([], [], color=st.PSEUDO, lw=4, label="Pseudomonas_E"),
           Rectangle((0, 0), 1, 1, facecolor="#F4F4F4", edgecolor=st.LIGHT,
                     label="— no alignment")]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.06),
          fontsize=st.FS_BODY, handlelength=1.2, ncol=3, columnspacing=1.4)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S6")
