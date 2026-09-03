"""Suppl Fig S10 - MICP feature prevalence within the Pseudomonas_E genus (analysis A3b).

Prevalence of each MICP feature across the 146 NCBI Pseudomonas_E reference genomes
screened by Pfam hmmscan, as horizontal bars with the percentage and the hit / total count
bound to each bar.  The two architecture features (ureC + ureB on one contig, ureC + CA on
one contig) come from the single-contig table, the gene-presence features from the rarity
screen; both tables cover the same 146 accessions, which is asserted before merging.

Provenance: Tables S13a and S13b.  Every percentage is recomputed from the 0/1 columns;
nothing is taken from a pre-tabulated summary.

Colour meaning: coral = a MICP feature carried by the reference genomes.  The bar for
"ureC + CA on one contig" - the architecture M1 and S26 carry - is not given a separate
colour; it is identified by its own row label.
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, TEXT, AXIS, FS_BODY, FS_STAT

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

# feature -> (source table, column, row label); category names, not data
FEATURES = [("a", "UreC_alpha", "UreC (PF00449)"),
            ("a", "UreB_beta_gamma", "UreB (PF00699)"),
            ("a", "urease_core", "urease core (UreC + UreB)"),
            ("a", "has_CA", "any carbonic anhydrase"),
            ("b", "ureC_and_ureB_single_contig", "ureC + ureB on one contig"),
            ("b", "ureC_and_CA_single_contig", "ureC + CA on one contig")]

# ------------------------------------------------------------------ data
a = pd.read_csv(SUPP / "Table_S13a_pseudomonas_e_MICP_rarity_screen.csv")
b = pd.read_csv(SUPP / "Table_S13b_pseudomonas_e_single_contig.csv")
assert sorted(a.accession) == sorted(b.accession)
N = len(a)
src = {"a": a, "b": b}

rows = []
for tbl, col, label in FEATURES:
    hits = int(src[tbl][col].sum())
    rows.append((label, hits, 100 * hits / N))
rows = rows[::-1]                       # first feature ends up at the top of the axis
labels = [r[0] for r in rows]
hits = [r[1] for r in rows]
pct = [r[2] for r in rows]

# ------------------------------------------------------------------ page
H = 57.0
fig, ax_mm, text_mm, letter = st.page(H)
letter(11, 5, "A")

ax = ax_mm(52, 11, 86, 35)
y = range(len(rows))
ax.barh(list(y), pct, 0.6, color=HERO, edgecolor=AXIS, linewidth=0.5)
for yi, p, h in zip(y, pct, hits):
    ax.text(p + 1.5, yi, f"{p:.1f}%   {h}/{N}", va="center", ha="left",
            fontsize=FS_STAT, color=TEXT)
ax.set_yticks(list(y))
ax.set_yticklabels(labels)
ax.set_xlabel("Pseudomonas_E reference genomes (%)")
ax.set_xlim(0, 118)
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_xticks([0, 25, 50, 75, 100])
st.style_axis(ax)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S10")
