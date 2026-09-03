"""Suppl Fig S14 - MICP pathway gene dosage per MICP-complete MAG (analysis A6).

Gene-copy heat map: rows are the six MICP-complete MAGs, columns the urease structural and
accessory subunits, the three carbonic-anhydrase families plus the generic CA hits, and the
calcium / carbonate / cation-transport genes.  Cell colour encodes copy number and the cell
carries the count.  Two stat columns to the right carry the two completeness calls stored in
the table: the urease core and the calcium pathway.

Provenance: Table S15b.  The Ca_pathway call is recomputed here (Ca_transporter or
Ca_ATPase present) and asserted against the stored flag, which is what makes the
"at least one Ca-handling gene in 5 of 6" statement checkable: S26 carries neither, and its
five CO3_transporter copies do not enter the stored call.

Colour meaning: green intensity = gene copy number (one meaning, whole page).  The stat
columns are text, not colour, so nothing competes with the copy-number scale.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import GREEN, TEXT, FS_BODY, FS_STAT, HEROES, hero_col

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

# column groups (category names, not data)
URE = ["ureA", "ureB", "ureC", "ureD_H", "ureE", "ureF", "ureG"]
CA = ["cah_alphaCA", "canA_gammaCA", "cynT_betaCA", "CA_generic"]
ION = ["Ca_transporter", "Ca_ATPase", "CO3_transporter", "Na_H_antiporter_Mrp",
       "K_transport"]
COLS = URE + CA + ION
NICE = {"ureD_H": "ureD/H", "cah_alphaCA": "cah (α-CA)", "canA_gammaCA": "canA (γ-CA)",
        "cynT_betaCA": "cynT (β-CA)", "CA_generic": "CA generic",
        "Ca_transporter": "Ca transporter", "Ca_ATPase": "Ca ATPase",
        "CO3_transporter": "CO₃ transporter", "Na_H_antiporter_Mrp": "Mrp Na⁺/H⁺",
        "K_transport": "K⁺ transport"}

# ------------------------------------------------------------------ data
df = pd.read_csv(SUPP / "Table_S15b_stoichiometry_per_MAG.csv")
df = df[df.group == "MICP_complete"].set_index("MAG").loc[HEROES]
recomputed_ca = ((df.Ca_transporter > 0) | (df.Ca_ATPase > 0)).astype(int)
assert (recomputed_ca == df.Ca_pathway).all()      # the stored call is Ca transporter/ATPase
mat = df[COLS].values.astype(int)

# ------------------------------------------------------------------ page
H = 60.0
fig, ax_mm, text_mm, letter = st.page(H)
letter(11, 5, "A")

ax = ax_mm(24, 12, 104, 28)
ax.imshow(mat, cmap=st.seq_cmap("copies", hi=GREEN), vmin=0, vmax=mat.max(),
          aspect="auto")
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j, i, mat[i, j], ha="center", va="center", fontsize=FS_BODY,
                color="white" if mat[i, j] > mat.max() / 2 else TEXT)
ax.set_xticks(range(len(COLS)))
ax.set_xticklabels([NICE.get(c, c) for c in COLS], fontsize=FS_BODY,
                   rotation=90)
ax.set_yticks(range(len(HEROES)))
ax.set_yticklabels(HEROES)
for tick, mag in zip(ax.get_yticklabels(), HEROES):
    tick.set_color(hero_col(mag))
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

# stat columns bound to the rows
for dx, head, vals in ((len(COLS) + 0.4, "urease\ncore", df.urease_core_complete),
                       (len(COLS) + 2.6, "Ca\npathway", df.Ca_pathway)):
    ax.text(dx, -1.15, head, fontsize=FS_STAT, ha="left", va="center", color=TEXT)
    for i, v in enumerate(vals):
        ax.text(dx, i, str(int(v)), fontsize=FS_STAT, ha="left", va="center", color=TEXT)
for t in ax.texts:
    t.set_clip_on(False)

cb = fig.colorbar(ax.images[0], ax=ax, fraction=0.03, pad=0.22, aspect=12)
cb.set_label("gene copies", fontsize=FS_BODY)
cb.ax.tick_params(labelsize=FS_BODY, length=2)
cb.outline.set_visible(False)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S14")
