"""Layout helpers shared by the genus-aggregated trait heat maps (Suppl S1-S5).

Layout only: this module contains no data value.  The aggregation reproduces
`scripts/pptx_builder/build_editable_figures.py::_genus_aggregate` (mean per genus,
genera with >= 2 MAGs, the two MICP-complete-containing genera first, the rest by
descending row sum) and the log10(v + 1) colour transform of `_draw_genus_heatmap`.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

import _style as st

SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
TRAIT_TABLE = SUPP / "Table_S2b_trait_module_per1kCDS.csv"
PSEUDOCOUNT = 1.0            # colour transform log10(v + PSEUDOCOUNT), from the old builder
MIN_MAGS_PER_GENUS = 2       # readability filter, from the old builder
HERO_GENERA = ["Sphingobacterium", "Pseudomonas_E"]


def load_traits():
    return pd.read_csv(TRAIT_TABLE)


def genus_aggregate(df, cols):
    """Genus-mean table for `cols`; returns (table, n_MAGs per genus)."""
    counts = df["Genus"].value_counts()
    keep = counts[counts >= MIN_MAGS_PER_GENUS].index.tolist()
    sub = df[df["Genus"].isin(keep)]
    agg = sub.groupby("Genus")[cols].mean()
    first = [g for g in HERO_GENERA if g in agg.index]
    rest = sorted([g for g in agg.index if g not in first],
                  key=lambda g: -agg.loc[g].sum())
    order = first + rest
    return agg.loc[order], counts.loc[order]


def draw_genus_heatmap(fig, ax, cax, table, n_mags, cell_text_fmt="{:.1f}"):
    """Genus x subcategory heat map with value labels and a colour bar.

    Returns the log10-transformed matrix so the caller can assert against it.
    """
    mat = np.log10(table.values + PSEUDOCOUNT)
    vmax = float(mat.max()) if mat.size else 1.0
    vmax = vmax if vmax > 0 else 1.0
    cmap = st.seq_cmap()
    ax.imshow(mat, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")

    ax.set_xticks(range(table.shape[1]))
    ax.set_xticklabels([c.split("::", 1)[-1].replace("_", " ")
                        for c in table.columns], rotation=90,
                       ha="center", va="bottom", fontsize=st.FS_BODY)
    ax.xaxis.set_ticks_position("top")
    ax.tick_params(axis="x", length=0, pad=2)

    ax.set_yticks(range(table.shape[0]))
    ax.set_yticklabels([f"{g}  (n = {n_mags[g]})" for g in table.index],
                       fontsize=st.FS_BODY)
    ax.tick_params(axis="y", length=0, pad=2)
    for lab, g in zip(ax.get_yticklabels(), table.index):
        if g in HERO_GENERA:
            lab.set_color(st.HERO)
            lab.set_fontweight("bold")

    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, table.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, table.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", length=0)

    for ri in range(table.shape[0]):
        for ci in range(table.shape[1]):
            v = table.values[ri, ci]
            frac = mat[ri, ci] / vmax
            ax.text(ci, ri, cell_text_fmt.format(v), ha="center", va="center",
                    fontsize=st.FS_STAT,
                    color="white" if frac > 0.55 else st.TEXT)

    grad = np.linspace(vmax, 0, 128).reshape(-1, 1)
    cax.imshow(grad, cmap=cmap, vmin=0, vmax=vmax, aspect="auto",
               extent=[0, 1, 0, vmax])
    cax.set_xticks([])
    cax.yaxis.set_ticks_position("right")
    cax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=False))
    cax.tick_params(labelsize=st.FS_STAT, length=2, pad=1)
    for s in cax.spines.values():
        s.set_linewidth(0.6)
        s.set_color(st.AXIS)
    return mat, vmax
