"""Suppl Fig S19 - in-situ relative-abundance proxy from SPAdes contig coverage (C6).

Panels (old files -> new page):
  Figure_S19.png           A  length-weighted coverage, MICP-complete vs rest
  Figure_S19_by_source.png B  length-weighted coverage by waste source, MICP-complete MAGs
                              overlaid on the source distribution

The old figure had three panels: coverage by group, coverage by source, and the heroes drawn
over the source distribution.  The last two are the same plot with and without the hero
points, so they are merged into one panel B in which the six MICP-complete MAGs are marked
inside their own source column.  The merge is recorded in the journal.

Coverage spans two orders of magnitude, so both panels use a log10 y-axis; the whiskers
therefore show the full range without compressing the bulk of the distribution.

SOURCE LABELLING DEFECT (corrected here)
The `source` column of Table S21a attaches "sheep" to the 15 M-prefixed MAGs and "swine" to
the 32 S-prefixed MAGs.  That is the opposite of the sample-code key used throughout the
manuscript (C cattle, M swine, S sheep, V poultry) and of the independent per-MAG source
column in Table S9a, which both make M swine and S sheep.  This page therefore derives the
waste source from the MAG identifier - a rule that can be checked against the identifier
itself - and asserts it against Table S9a.  Table S21b inherits the same swap, so its stored
per-source statistics are asserted against the recomputed groups under the swapped names,
which both reproduces the stored numbers and demonstrates the swap.  The consequence for the
manuscript is recorded in the journal: the claim that sheep MAGs are the most deeply covered
(mean 57.4x) is really about the swine MAGs.

Sources
  Table_S21a_abundance_proxy_per_MAG.csv    111 MAGs, per-MAG coverage statistics
  Table_S21b_abundance_proxy_per_source.csv stored per-source count/mean/median/std
                                            (recomputed and asserted against, under the swap)
  Table_S9a_PCoA_coordinates.csv            independent per-MAG waste source (asserted against)

Colour meanings on this page
  A: coral = MICP-complete group, grey = the remaining 105 MAGs
  B: four waste-source colours (cattle brown, swine pink, sheep green, poultry purple) for
     the distributions, and black-edged points for the six MICP-complete MAGs.  Coral is not
     reused in B, so no colour carries two meanings on the page; the MICP-complete MAGs are
     marked by outline rather than by fill because their fill states which source they came from.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import HERO, REST, TEXT, GREY, FS_STAT

st.setup()
OUT = HERE / "figures"
SUPP = Path(gh.SUPP)

COL = "length_weighted_cov"
SOURCE_ORDER = ["cattle", "swine", "sheep", "poultry"]
PREFIX_SOURCE = {"C": "cattle", "M": "swine", "S": "sheep", "V": "poultry"}
# how the labels in Table S21a/b map onto the corrected groups
STORED_LABEL = {"cattle": "cattle", "swine": "sheep", "sheep": "swine", "poultry": "poultry"}

# ------------------------------------------------------------------ data
d = pd.read_csv(SUPP / "Table_S21a_abundance_proxy_per_MAG.csv")
per_source = pd.read_csv(SUPP / "Table_S21b_abundance_proxy_per_source.csv").set_index("source")

hero = d.is_hero.astype(bool)
assert sorted(d.loc[hero, "MAG"]) == sorted(st.HEROES)
assert set(d.source) == set(SOURCE_ORDER), sorted(set(d.source))

# waste source taken from the MAG identifier, checked against the independent Table S9a
d["waste_source"] = d.MAG.str[0].map(PREFIX_SOURCE)
assert d.waste_source.notna().all()
s9 = pd.read_csv(SUPP / "Table_S9a_PCoA_coordinates.csv")
s9_src = dict(zip(s9.MAG, s9.Source.str.lower()))
assert all(s9_src[m] == w for m, w in zip(d.MAG, d.waste_source)), "S9a disagrees with the prefix rule"
n_swapped = int((d.source != d.waste_source).sum())
print(f"  Table S21a source column disagrees with the MAG identifier for {n_swapped} of {len(d)} MAGs "
      f"(swine and sheep are exchanged)")

h = d.loc[hero, COL].to_numpy(float)
r = d.loc[~hero, COL].to_numpy(float)
p_group = mannwhitneyu(h, r, alternative="two-sided").pvalue

# the stored per-source table is recomputed from the per-MAG rows, under the label swap
for src in SOURCE_ORDER:
    v = d.loc[d.waste_source == src, COL]
    s = per_source.loc[STORED_LABEL[src]]
    assert len(v) == s["count"], (src, len(v), s["count"])
    assert abs(v.mean() - s["mean"]) < 1e-6, (src, v.mean(), s["mean"])
    assert abs(v.median() - s["median"]) < 1e-6, (src, v.median(), s["median"])
    assert abs(v.std(ddof=1) - s["std"]) < 1e-6, (src, v.std(ddof=1), s["std"])
print(f"  group means: MICP-complete {h.mean():.1f}x  rest {r.mean():.1f}x  "
      f"MWU P = {p_group:.4f}")

# ------------------------------------------------------------------ page
H = 76.0
fig, ax_mm, text_mm, letter = st.page(H)

TOP, PH = 15.0, 47.0

# ---- A: group comparison
axA = ax_mm(17.0, TOP, 46.0, PH)
groups, colours = [("MICP-complete", h), ("Rest", r)], [HERO, REST]
xs, _ = gh.strip_box(axA, groups, colours, log=True)
axA.set_ylabel("Length-weighted contig coverage (×)")
axA.set_ylim(1.0, 3000.0)
axA.set_yticks([1, 10, 100, 1000])
gh.group_counts(axA, xs, groups, -0.12)
gh.stat_bracket(axA, xs[0], xs[1], 1100.0, gh.fmt_p(p_group), drop=250.0)
letter(5.0, TOP - 6.0, "A")

# ---- B: by waste source, MICP-complete MAGs marked in their own column
axB = ax_mm(84.0, TOP, 82.0, PH)
src_groups = [(s, d.loc[d.waste_source == s, COL].to_numpy(float)) for s in SOURCE_ORDER]
src_cols = [st.SOURCE[s] for s in SOURCE_ORDER]
xsB, jitB = gh.strip_box(axB, src_groups, src_cols, log=True)
for si, s in enumerate(SOURCE_ORDER):
    in_src = (d.waste_source == s).to_numpy()
    is_hero_in_src = hero.to_numpy()[in_src]
    if not is_hero_in_src.any():
        continue
    axB.scatter(jitB[si][is_hero_in_src], d.loc[in_src, COL].to_numpy()[is_hero_in_src],
                s=26, facecolor="none", edgecolor=TEXT, linewidth=0.9, zorder=4)
axB.set_ylabel("Length-weighted contig coverage (×)")
axB.set_ylim(1.0, 3000.0)
axB.set_yticks([1, 10, 100, 1000])
gh.group_counts(axB, xsB, src_groups, -0.12)
letter(72.0, TOP - 6.0, "B")

handles = [Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
                  markeredgecolor=TEXT, markersize=5, label="MICP-complete MAG")]
axB.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.06), frameon=False,
           handletextpad=0.4)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S19")
