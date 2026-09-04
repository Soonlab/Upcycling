"""Consolidated Suppl Fig S1 - trait-module landscape across genera.

One composed page, 180 mm wide, five genus-aggregated trait heat maps
(reading order A -> E, left to right then top to bottom):

  A  biofilm / EPS modules                          (old build_supS1.py)
  B  ammonia handling and nitrogen assimilation     (old build_supS2.py)
  C  alkaline and osmotic stress tolerance          (old build_supS3.py)
  D  CAZy keyword-proxy subcategories               (old build_supS4.py panel A)
  E  heavy-metal and antibiotic-resistance modules  (old build_supS5.py)

The dbCAN panels of old build_supS4.py (its panels B and C) are not carried here; they
belong to the consolidated Suppl Fig S2 (DESIGN.md).

Source: SUBMISSION/Supplementary_tables/Table_S2b_trait_module_per1kCDS.csv
        (per-MAG hits per 10^3 CDS), one subcategory block per panel.

Construction, unchanged from the five superseded scripts and shared through
`_supp_traits.py`: mean per genus over genera with >= 2 MAGs, the two MICP-complete-
containing genera first and the rest by descending row sum, colour on log10(value + 1).
Cell labels carry the untransformed mean, so nothing is read off the colour alone.  Each
panel is aggregated and colour-scaled independently, exactly as before; because the log
range differs from block to block (the biofilm block reaches two orders of magnitude
above the ammonia block) the panels do NOT share one colour scale, and each therefore
carries its own colour bar.  The row order likewise differs between panels because the
descending-row-sum rank is computed within each block.

Provenance: no value is hard-coded; every panel re-derives its genus means from the
per-MAG rows of Table S2b and asserts the aggregate against those rows, as the old
scripts did.

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
OUT = HERE / "figures_v2"

# subcategory blocks of Table S2b, in reading order (category names, not data)
PANELS = [("A", "Biofilm_EPS::"),
          ("B", "Ammonia_N::"),
          ("C", "Alkaline_Osmo::"),
          ("D", "CAZyme_proxy::"),
          ("E", "MetalResist_AMR::")]
CBAR_LABEL = "log₁₀(mean hits per 10³ CDS + 1)"

# Row pitch: the old one-panel-per-page scripts used 6.0 mm, which put five stacked
# panels 43 mm over the one-page ceiling of 235 mm.  4.8 mm still leaves a 7 pt genus
# label (~3.0 mm line height) and a 6.5 pt cell value clear of their neighbours; the
# audit is the check.
H_ROW = 4.8          # heat-map row height, mm
WRAP_CHARS = 15      # a column label longer than this is wrapped at a space ...
WRAP_MIN_PITCH = 8.0 # ... but only where the cell pitch can hold a two-line block
CBAR_W = 3.0
CBAR_GAP = 6.0       # heat map -> colour bar, mm (unchanged)
CBAR_LAB_DX = 15.5   # heat map -> colour-bar caption, mm (unchanged)
CBAR_H_MAX = 30.0

# ------------------------------------------------------------------ data
df = tr.load_traits()
tables = {}
for key, cat in PANELS:
    cols = [c for c in df.columns if c.startswith(cat)]
    assert cols, cat
    table, n_mags = tr.genus_aggregate(df, cols)
    # the aggregate must be reproducible from the per-MAG rows
    for g in table.index:
        rows = df[df["Genus"] == g]
        assert np.allclose(table.loc[g].values, rows[cols].mean().values)
        assert n_mags[g] == len(rows)
    tables[key] = (table, n_mags)

n_rows = {k: tables[k][0].shape[0] for k, _ in PANELS}
assert len(set(n_rows.values())) == 1          # every block spans the same genera

# ------------------------------------------------------------------ page
BLOCK_H = H_ROW * n_rows["A"]

# per-row header space for the rotated column labels, mm.  Band A is short because its
# one long label ("polysaccharide intercellular") is wrapped onto two lines - panel A has
# the widest cells on the page, so the two-line block fits between its columns.
HDR_A, HDR_BC, HDR_DE = 18.0, 25.0, 25.0
TOP_A = 4.0 + HDR_A
TOP_BC = TOP_A + BLOCK_H + 8.0 + HDR_BC
TOP_DE = TOP_BC + BLOCK_H + 8.0 + HDR_DE
H = TOP_DE + BLOCK_H + 10.0

# x geometry: A spans the page, B/C and D/E sit in two columns
X_A, W_A = 34.0, 112.0
X_L, W_HALF = 33.0, 36.0
X_R = 122.0
POS = {"A": (X_A, TOP_A, W_A, 4.0, 4.0),
       "B": (X_L, TOP_BC, W_HALF, 4.0, TOP_BC - HDR_BC),
       "C": (X_R, TOP_BC, W_HALF, 88.0, TOP_BC - HDR_BC),
       "D": (X_L, TOP_DE, W_HALF, 4.0, TOP_DE - HDR_DE),
       "E": (X_R, TOP_DE, W_HALF, 88.0, TOP_DE - HDR_DE)}

fig, ax_mm, text_mm, letter = st.page(H)


def wrap_ticklabels(ax, table):
    """Wrap over-long column labels at a space; layout only, the text is unchanged."""
    out = []
    for t in ax.get_xticklabels():
        s = t.get_text()
        if len(s) > WRAP_CHARS and " " in s:
            i = s.rfind(" ", 0, WRAP_CHARS + 1)
            i = i if i > 0 else s.find(" ")
            s = s[:i] + "\n" + s[i + 1:]
        out.append(s)
    ax.set_xticklabels(out, rotation=90, ha="center", va="bottom",
                       fontsize=st.FS_BODY)


for key, _cat in PANELS:
    table, n_mags = tables[key]
    x, top, w, lx, ly = POS[key]
    cbar_h = min(CBAR_H_MAX, BLOCK_H)
    ax = ax_mm(x, top, w, BLOCK_H)
    cax = ax_mm(x + w + CBAR_GAP, top, CBAR_W, cbar_h)
    tr.draw_genus_heatmap(fig, ax, cax, table, n_mags)
    if w / table.shape[1] >= WRAP_MIN_PITCH:
        wrap_ticklabels(ax, table)
    text_mm(x + w + CBAR_LAB_DX, top + cbar_h / 2, CBAR_LABEL, rotation=90,
            ha="center", va="center", fontsize=st.FS_BODY)
    letter(lx, ly, key)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S1")
