"""Suppl Fig S20 - ESMFold UreC backbone agreement with PDB 4CEU (C4).

Panels (old file Figure_S20.png -> new page):
  A  backbone TM-score per MICP-complete UreC, normalised to 4CEU chain C
  B  all-residue backbone RMSD per MICP-complete UreC

Bars are ordered by descending TM-score in A, and B keeps that same MAG order so the two
panels read as one comparison.  The TM = 0.5 same-fold threshold is drawn as a reference
line in A; it is a published constant, not a measurement from this study, so it is labelled
as a threshold and coloured grey.

Sources
  Table_S22_ureC_vs_4CEU_tm.csv   MAG, pred_len, ref_len, tm_norm_pred, tm_norm_ref, rmsd
The `_UreC` suffix is stripped from the MAG identifier so the labels match every other
figure in the set.  `ref_len` is asserted constant across rows (a single reference chain).

Colour meanings on this page
  blue   = Sphingobacterium lineage (S13, S16, S23, C22)
  orange = Pseudomonas_E lineage (M1, S26)
  grey   = the TM = 0.5 same-fold threshold, which is not a lineage
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import GREY, TEXT, FS_STAT
from matplotlib.patches import Patch

st.setup()
OUT = HERE / "figures"
SUPP = Path(gh.SUPP)

TM_SAME_FOLD = 0.5  # Xu & Zhang 2010 threshold, a published constant

# ------------------------------------------------------------------ data
d = pd.read_csv(SUPP / "Table_S22_ureC_vs_4CEU_tm.csv")
d["MAG"] = d.MAG.str.replace("_UreC", "", regex=False)
assert sorted(d.MAG) == sorted(st.HEROES), sorted(d.MAG)
assert d.ref_len.nunique() == 1, d.ref_len.unique()
d = d.sort_values("tm_norm_ref", ascending=False).reset_index(drop=True)
cols = [st.hero_col(m) for m in d.MAG]

# ------------------------------------------------------------------ page
H = 66.0
fig, ax_mm, text_mm, letter = st.page(H)

LEFT = [17.0, 105.0]
TOP = 13.0
W, PH = 62.0, 40.0
x = np.arange(len(d), dtype=float)

# ---- A: TM-score
axA = ax_mm(LEFT[0], TOP, W, PH)
axA.bar(x, d.tm_norm_ref, width=0.66, color=cols, edgecolor="none", zorder=2)
gh.value_labels(axA, x, d.tm_norm_ref.to_numpy(), fmt="{v:.3f}", dy=0.01)
axA.axhline(TM_SAME_FOLD, color=GREY, lw=0.8, ls="--", zorder=3)
axA.text(1.02, TM_SAME_FOLD, "TM = 0.5", transform=axA.get_yaxis_transform(),
         ha="left", va="center", fontsize=FS_STAT, color=GREY)
axA.set_xticks(x)
axA.set_xticklabels(d.MAG)
axA.set_ylabel("TM-score vs 4CEU chain C")
axA.set_ylim(0, 0.78)
st.style_axis(axA)
axA.set_xlim(-0.7, len(d) - 0.3)
letter(LEFT[0] - 13.0, TOP - 6.0, "A")

# ---- B: RMSD
axB = ax_mm(LEFT[1], TOP, W, PH)
axB.bar(x, d.rmsd, width=0.66, color=cols, edgecolor="none", zorder=2)
gh.value_labels(axB, x, d.rmsd.to_numpy(), fmt="{v:.2f}", dy=0.07)
axB.set_xticks(x)
axB.set_xticklabels(d.MAG)
axB.set_ylabel("Backbone RMSD (Å)")
axB.set_ylim(0, 5.0)
st.style_axis(axB)
axB.set_xlim(-0.7, len(d) - 0.3)
letter(LEFT[1] - 13.0, TOP - 6.0, "B")

handles = [Patch(facecolor=st.LINEAGE_COL[g], edgecolor="none", label=lab)
           for g, lab in [("Sphingobacterium", "Sphingobacterium"),
                          ("Pseudomonas_E", "Pseudomonas_E")]]
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, 0.015), handlelength=1.2, columnspacing=1.6)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S20")
