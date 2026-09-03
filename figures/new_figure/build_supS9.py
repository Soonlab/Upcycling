"""Suppl Fig S9 - UreC active-site residue conservation, extended view (analysis A2).

Rows are the seven canonical catalytic positions of the S. pasteurii reference
(UniProt P41020 / PDB 4CEU chain C), ordered by residue number; columns are the six
MICP-complete MAGs, ordered Sphingobacterium first then Pseudomonas_E.  Cell colour
encodes agreement with the reference and the cell carries the observed residue.
Two stat columns to the right of the map carry the reference residue and the column of the
multiple alignment the comparison was made in.

Provenance: Table S12 in full - no value is recomputed, and none is dropped.  The site
order in the table (H137, H139, K220, H249, H275, D363, C322) is re-sorted by residue
number so that reading order matches the numbering.

Colour meanings: green = observed residue matches the reference, coral = differs (legend
entry only; no cell takes it), blue / orange = the two lineages, on the column labels.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, GREEN, TEXT, FS_BODY, FS_STAT, HEROES, hero_col

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

# ------------------------------------------------------------------ data
df = pd.read_csv(SUPP / "Table_S12_UreC_active_site_residues.csv")
df["pos"] = df.site.str.extract(r"(\d+)").astype(int)
df = df.sort_values("pos").reset_index(drop=True)
assert (df.expected == df.ref_aa).all()        # the table's two reference columns agree
obs = df[HEROES].values                        # rows = site, cols = MAG
match = obs == df.expected.values[:, None]
n_site, n_mag = obs.shape

# ------------------------------------------------------------------ page
H = 64.0
fig, ax_mm, text_mm, letter = st.page(H)
letter(11, 5, "A")

ax = ax_mm(26, 14, 88, 42)
ax.imshow(np.where(match, 1, 0), cmap=st.seq_cmap("match", hi=GREEN), vmin=0, vmax=1,
          aspect="auto")
for i in range(n_site):
    for j in range(n_mag):
        ax.text(j, i, obs[i, j], ha="center", va="center", fontsize=FS_BODY,
                color="white" if match[i, j] else TEXT)
ax.set_xticks(range(n_mag))
ax.set_xticklabels(HEROES)
for tick, mag in zip(ax.get_xticklabels(), HEROES):
    tick.set_color(hero_col(mag))
ax.set_yticks(range(n_site))
ax.set_yticklabels(df.site)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(handles=[Patch(facecolor=GREEN, label="matches reference"),
                   Patch(facecolor=HERO, label="differs")],
          loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.1,
          handleheight=0.9, columnspacing=1.2, fontsize=FS_BODY)

# stat columns bound to the rows
for dx, head, vals in ((n_mag - 0.1, "reference", df.expected),
                       (n_mag + 1.1, "alignment column", df.ref_column)):
    ax.text(dx, -0.75, head, fontsize=FS_STAT, ha="left", va="center", color=TEXT)
    for i, v in enumerate(vals):
        ax.text(dx, i, str(v), fontsize=FS_STAT, ha="left", va="center", color=TEXT)
for t in ax.texts:
    t.set_clip_on(False)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S9")
