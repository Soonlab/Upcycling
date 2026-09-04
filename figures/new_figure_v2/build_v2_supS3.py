"""Consolidated Suppl Fig S3 - biosafety, mobile elements and defence systems.

Consolidation of 2026-09-04 (see consolidation_260904/DESIGN.md).  Four panels, one
180 mm page, two rows of four box-and-point columns:

  A  abricate hits per MAG across CARD, VFDB, ResFinder and PlasmidFinder, MICP-complete
     group vs the remaining 105 MAGs, with the two-sided Mann-Whitney P per database
     (old Suppl Fig S8)
  B  geNomad plasmid-flagged and virus-flagged contigs per MAG, same two groups
     (old Suppl Fig S15A)
  C  DefenseFinder anti-phage systems per MAG (old Suppl Fig S16A)
  D  minced CRISPR arrays per MAG (old Suppl Fig S16B)

Old Suppl Fig S15B (the per-MAG urease / CA vs geNomad cross-check table) is NOT on this
page; it moves to a main figure.  Its underlying Table S17b is still read here for the one
cross-table consistency check that involves a quantity panel B draws.

C and D are a true null: every one of the 111 MAGs carries zero detected defence systems
and zero CRISPR arrays.  Drawn honestly - the y-axis is held at 0-1 and the point cloud
sits on the zero line - rather than rescaled to manufacture spread.  The count of MAGs
above zero is printed under each column so a reader sees the null is a real 0 / n, not a
plotting failure.

Provenance
  A  CARD / VFDB / ResFinder counts are in Table S11.  Table S11 carries NO PlasmidFinder
     column, although the A1 run produced one: the counts are recomputed here from the raw
     abricate output `research/additional/A1_biosafety/combined/plasmidfinder_all.tsv` (one
     row per hit, the MAG taken from the #FILE path), which is the file the aggregation
     script `aggregate_biosafety.py` reads.  The recomputed CARD / VFDB / ResFinder counts
     from the same raw files are asserted against Table S11 before PlasmidFinder is added,
     so the fourth block rests on a source that reproduces the three published ones.
  B  Table S17a, with the per-MAG plasmid counts cross-checked against Table S17b.
  C, D  Table S18a per-MAG counts, with the group means and the Mann-Whitney P asserted
     against the stored Table S18b.
Every Mann-Whitney P on the page is recomputed here, two-sided.

Colour meanings on this page
  blue    a Sphingobacterium MICP-complete MAG (A, B, per-point lineage colour)
  orange  a Pseudomonas_E MICP-complete MAG (A, B, per-point lineage colour)
  coral   the MICP-complete group taken as a whole (C, D)
  grey    one of the remaining 105 MAGs, in every panel
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
from _style import (HERO, REST, SPHINGO, PSEUDO, TEXT, AXIS, FS_BODY, FS_STAT,
                    HEROES, hero_col)

st.setup()
OUT = HERE / "figures_v2"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
A1 = Path("/data/data/Upcycling/research/additional/A1_biosafety/combined")
RNG = np.random.default_rng(0)     # reproducible jitter (style constant)

DBS = ["card", "vfdb", "resfinder", "plasmidfinder"]
DB_LAB = {"card": "CARD", "vfdb": "VFDB", "resfinder": "ResFinder",
          "plasmidfinder": "PlasmidFinder"}
MGE_METRICS = [("n_plasmid_contigs", "plasmid-flagged contigs"),
               ("n_virus_contigs", "virus-flagged contigs")]
NULL_METRICS = [("n_defense_systems", "Defense systems per MAG"),
                ("n_crispr_arrays", "CRISPR arrays per MAG")]

# ------------------------------------------------------------------ data: A
s11 = pd.read_csv(SUPP / "Table_S11_biosafety_counts_per_MAG.csv").set_index("MAG")
mags = list(s11.index)


def raw_counts(db):
    """Per-MAG abricate hit count, recomputed from the combined raw output."""
    df = pd.read_csv(A1 / f"{db}_all.tsv", sep="\t", low_memory=False)
    mag = df["#FILE"].str.extract(r"bakta_results/([^/]+)/")[0]
    return mag.value_counts().reindex(mags).fillna(0).astype(int)


counts = pd.DataFrame({db: raw_counts(db) for db in DBS}, index=mags)
for db in ["card", "vfdb", "resfinder"]:            # the three published columns
    assert (counts[db] == s11[db]).all(), db
counts["group"] = s11["group"]
is_hero_a = counts.group == "MICP_complete"
assert sorted(counts.index[is_hero_a]) == sorted(HEROES)
n_hero, n_rest = int(is_hero_a.sum()), int((~is_hero_a).sum())

# ------------------------------------------------------------------ data: B
gnm = pd.read_csv(SUPP / "Table_S17a_genomad_summary_per_MAG.csv")
is_hero_b = gnm.group == "MICP_complete"
assert sorted(gnm.MAG[is_hero_b]) == sorted(HEROES)
assert int(is_hero_b.sum()) == n_hero and int((~is_hero_b).sum()) == n_rest

overlap = pd.read_csv(SUPP / "Table_S17b_ureCah_vs_MGE_overlap.csv").set_index("MAG").loc[HEROES]
# the per-MAG plasmid counts of the two tables must agree
assert (overlap.n_plasmid_contigs.values ==
        gnm.set_index("MAG").loc[HEROES].n_plasmid_contigs.values).all()

# ------------------------------------------------------------------ data: C and D
dfn = pd.read_csv(SUPP / "Table_S18a_defense_per_MAG.csv")
stored = pd.read_csv(SUPP / "Table_S18b_defense_hero_vs_rest.csv").set_index("metric")
is_hero_d = dfn.is_hero.astype(bool)
assert int(is_hero_d.sum()) == len(HEROES), is_hero_d.sum()
assert sorted(dfn.loc[is_hero_d, "MAG"]) == sorted(HEROES)

null_panels = []
for col, ylab in NULL_METRICS:
    h = dfn.loc[is_hero_d, col].to_numpy(float)
    r = dfn.loc[~is_hero_d, col].to_numpy(float)
    p = mannwhitneyu(h, r, alternative="two-sided").pvalue
    s = stored.loc[col]
    assert abs(h.mean() - s.hero_mean) < 1e-9, (col, h.mean(), s.hero_mean)
    assert abs(r.mean() - s.rest_mean) < 1e-9, (col, r.mean(), s.rest_mean)
    assert abs(p - s.MWU_p) < 1e-9, (col, p, s.MWU_p)
    null_panels.append(dict(col=col, ylab=ylab, h=h, r=r, p=p,
                            n_pos_h=int((h > 0).sum()), n_pos_r=int((r > 0).sum())))

# ------------------------------------------------------------------ page
# one four-column grid per row, so that the eight box-and-point columns of the page share
# a width and a baseline; the old S8 grid ran to 184 mm and clipped its fourth block
LEFT, W, GAP = 22.0, 30.0, 10.0
COL_X = [LEFT + k * (W + GAP) for k in range(4)]
assert COL_X[-1] + W <= st.PAGE_W_MM
ROW1_T, ROW2_T, PH = 15.0, 84.0, 44.0
H = 148.0

fig, ax_mm, text_mm, letter = st.page(H)


def lineage_columns(ax, hero_v, rest_v, hero_mags, block_label, ylab=None):
    """Two columns - MICP-complete and rest - with every MAG overlaid as a point, the
    MICP-complete points coloured by lineage.  Returns nothing; draws in place."""
    bp = ax.boxplot([hero_v, rest_v], positions=[0, 1], widths=0.55,
                    showfliers=False, patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor="white", edgecolor=AXIS, linewidth=0.7)
    for part in ("whiskers", "caps", "medians"):
        for ln in bp[part]:
            ln.set(color=AXIS, linewidth=0.7)
    for xi, vals, cols in ((0, hero_v, [hero_col(m) for m in hero_mags]),
                           (1, rest_v, [REST] * len(rest_v))):
        jit = RNG.uniform(-0.16, 0.16, len(vals))
        ax.scatter(xi + jit, vals, s=5, c=cols, alpha=0.85, linewidths=0, zorder=3)
    p = mannwhitneyu(hero_v, rest_v, alternative="two-sided").pvalue
    top = max(float(np.max(hero_v)), float(np.max(rest_v)), 1.0)
    ax.set_ylim(-0.06 * top, top * 1.34)
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"MICP-complete\nn = {n_hero}", f"rest\nn = {n_rest}"])
    ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= top * 1.34])
    ax.text(0.5, top * 1.16, f"P = {p:.2f}", ha="center", va="center", fontsize=FS_STAT,
            color=TEXT)
    ax.text(0.5, top * 1.34, block_label, ha="center", va="bottom", fontsize=FS_BODY,
            color=TEXT)
    if ylab:
        ax.set_ylabel(ylab)
    st.style_axis(ax)


# ---- A: abricate, one block per database ---------------------------------
for k, db in enumerate(DBS):
    ax = ax_mm(COL_X[k], ROW1_T, W, PH)
    lineage_columns(ax, counts.loc[is_hero_a, db].values,
                    counts.loc[~is_hero_a, db].values,
                    list(counts.index[is_hero_a]), DB_LAB[db],
                    ylab="abricate hits per MAG" if k == 0 else None)

# ---- B: geNomad plasmid- and virus-flagged contigs ------------------------
for k, (col, lab) in enumerate(MGE_METRICS):
    ax = ax_mm(COL_X[k], ROW2_T, W, PH)
    lineage_columns(ax, gnm.loc[is_hero_b, col].values,
                    gnm.loc[~is_hero_b, col].values,
                    list(gnm.MAG[is_hero_b]), lab,
                    ylab="geNomad contigs per MAG" if k == 0 else None)

# ---- C and D: the two null assays ----------------------------------------
for k, pan in enumerate(null_panels):
    ax = ax_mm(COL_X[2 + k], ROW2_T, W, PH)
    groups = [("MICP-complete", pan["h"]), ("Rest", pan["r"])]
    xs, _ = gh.strip_box(ax, groups, [HERO, REST])
    ax.set_ylabel(pan["ylab"])
    ax.set_ylim(-0.08, 1.0)
    ax.set_yticks([0, 0.5, 1.0])
    gh.group_counts(ax, xs, groups, -0.15)
    # count above zero, bound to each column: shows the null is 0 / n, not a missing series
    for x, n_pos, n_tot in zip(xs, [pan["n_pos_h"], pan["n_pos_r"]],
                               [len(pan["h"]), len(pan["r"])]):
        ax.text(x, 0.10, f"{n_pos} / {n_tot} > 0", ha="center", va="bottom",
                fontsize=FS_STAT, color=TEXT)
    gh.stat_bracket(ax, xs[0], xs[1], 0.80, gh.fmt_p(pan["p"]), drop=0.05)

fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, color=SPHINGO,
                           label="Sphingobacterium"),
                    Line2D([], [], marker="o", ls="", ms=3, color=PSEUDO,
                           label="Pseudomonas_E"),
                    Line2D([], [], marker="s", ls="", ms=3, color=HERO,
                           label="MICP-complete group"),
                    Line2D([], [], marker="o", ls="", ms=3, color=REST,
                           label="other MAG")],
           loc="lower center", ncol=4, frameon=False, fontsize=FS_BODY,
           bbox_to_anchor=(0.5, 0.004), handletextpad=0.4, columnspacing=2.2)

letter(COL_X[0] - 14.0, ROW1_T - 8.0, "A")
letter(COL_X[0] - 14.0, ROW2_T - 8.0, "B")
letter(COL_X[2] - 14.0, ROW2_T - 8.0, "C")
letter(COL_X[3] - 14.0, ROW2_T - 8.0, "D")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S3")
