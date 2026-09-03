"""Suppl Fig S5 - heavy-metal and antibiotic-resistance modules, genus-aggregated.

Source: SUBMISSION/Supplementary_tables/Table_S2b_trait_module_per1kCDS.csv
        (per-MAG hits per 10^3 CDS; the MetalResist_AMR:: subcategory block - the header name
        of the category the legend calls "heavy-metal and antibiotic resistance").

Aggregation and colour transform reproduce the superseded PowerPoint builder
(`scripts/pptx_builder/build_editable_figures.py::_genus_aggregate` and
`_draw_genus_heatmap`): mean per genus over genera with >= 2 MAGs, the two
MICP-complete-containing genera first and the rest by descending row sum, colour on
log10(value + 1).  Cell labels carry the untransformed mean, so nothing is read off
the colour alone.

Colour meanings on this page:
  white -> green   mean hits per 10^3 CDS for that genus x subcategory (the only encoding)
  coral            genus label of a genus containing MICP-complete MAGs
  black            every other genus label
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _supp_traits as tr

st.setup()
OUT = HERE / "figures"
CATEGORY = "MetalResist_AMR::"

df = tr.load_traits()
cols = [c for c in df.columns if c.startswith(CATEGORY)]
table, n_mags = tr.genus_aggregate(df, cols)

# the aggregate must be reproducible from the per-MAG rows
for g in table.index:
    rows = df[df["Genus"] == g]
    assert np.allclose(table.loc[g].values, rows[cols].mean().values)
    assert n_mags[g] == len(rows)

H_ROW, TOP, LEFT = 6.0, 36.0, 36.0
W = 126.0
H = TOP + H_ROW * table.shape[0] + 8.0

fig, ax_mm, text_mm, letter = st.page(H)
ax = ax_mm(LEFT, TOP, W, H_ROW * table.shape[0])
cax = ax_mm(LEFT + W + 6.0, TOP, 3.0, min(30.0, H_ROW * table.shape[0]))
letter(4, 4, "A")

tr.draw_genus_heatmap(fig, ax, cax, table, n_mags)
text_mm(LEFT + W + 15.5, TOP + min(30.0, H_ROW * table.shape[0]) / 2,
        "log₁₀(mean hits per 10³ CDS + 1)", rotation=90, ha="center",
        va="center", fontsize=st.FS_BODY)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S5")
