"""Suppl Fig S13 - external MICP rarity across the MGnify livestock catalogs (analysis B).

A  percentage of species-cluster representatives carrying the MICP gene-complete profile
B  percentage carrying the single-contig ureC + CA architecture
Both are shown per catalog and pooled, with the catalog size and the hit count bound to
each bar.

Provenance: Table S14a.  Every percentage is recomputed from the count columns and asserted
against the stored pct columns; the pooled bar is recomputed from the summed counts.  The
"present study 6/6 = 100%" reference line carried by the shipped figure is dropped: it is a
statement about a different denominator (this panel's 111-MAG study set, not an MGnify
species cluster), and it belongs in the legend.

Colour meanings: dark green = MICP gene-complete prevalence (A), light green =
single-contig ureC + CA prevalence (B).  The two greens are the same pair as main Fig 8c.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import GREEN, TEXT, AXIS, FS_BODY, FS_STAT

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
GREEN_LT = tuple(0.45 + 0.55 * c for c in to_rgb(GREEN))

# ------------------------------------------------------------------ data
mg = pd.read_csv(SUPP / "Table_S14a_mgnify_catalog_summary.csv")
pct_complete = 100 * mg.n_MICP_gene_complete / mg.n_species_clusters
pct_single = 100 * mg.n_MICP_single_contig_ureC_CA / mg.n_species_clusters
assert np.allclose(pct_complete, mg.pct_MICP_gene_complete, atol=1e-3)
assert np.allclose(pct_single, mg.pct_MICP_single_contig, atol=1e-3)

n_pool = int(mg.n_species_clusters.sum())
labels = [c.replace("-", "\n") for c in mg.catalog] + ["pooled"]
n_clusters = list(mg.n_species_clusters) + [n_pool]
PANELS = [("A", "n_MICP_gene_complete", list(pct_complete), GREEN,
           "MICP gene-complete (%)"),
          ("B", "n_MICP_single_contig_ureC_CA", list(pct_single), GREEN_LT,
           "single-contig ureC + CA (%)")]

# ------------------------------------------------------------------ page
H = 62.0
fig, ax_mm, text_mm, letter = st.page(H)

for k, (lab, count_col, vals, colour, ylab) in enumerate(PANELS):
    left = 22.0 + k * 88.0
    letter(left - 11, 5, lab)
    ax = ax_mm(left, 13, 62, 34)
    counts = list(mg[count_col]) + [int(mg[count_col].sum())]
    v = vals + [100 * counts[-1] / n_pool]
    x = np.arange(len(v))
    ax.bar(x, v, 0.6, color=colour, edgecolor=AXIS, linewidth=0.5)
    for xi, vi, ci in zip(x, v, counts):
        ax.text(xi, vi, f"{vi:.2f}\n{ci:,}", ha="center", va="bottom",
                fontsize=FS_STAT, color=TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{l}\nn = {n:,}" for l, n in zip(labels, n_clusters)])
    ax.set_ylabel(ylab)
    ax.set_ylim(0, max(v) * 1.35)
    ax.set_xlim(-0.7, len(v) - 0.3)
    st.style_axis(ax)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S13")
