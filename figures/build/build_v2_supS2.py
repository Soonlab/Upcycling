"""Consolidated Suppl Fig S2 - DRAM metabolic context and dbCAN CAZyme validation.

Consolidation of 2026-09-04 (see consolidation_260904/DESIGN.md).  Four panels, one
180 mm page:

  A  genus-aggregated mean completeness of the DRAM modules that are scored as a fraction
     (old main Fig 4A)
  B  genus-aggregated prevalence of the CAZy and nitrogen-metabolism modules, which DRAM
     scores as present / absent rather than as a fraction (old main Fig 4B)
  C  dbCAN class abundance, MICP-complete vs rest, with the permutation statistics
     (old Suppl Fig S4B)
  D  dbCAN family abundance, MICP-complete vs rest, top families by MICP-complete mean
     (old Suppl Fig S4C)

Old main Fig 4C (MICP-critical trait modules) is NOT on this page; it moves to a main
figure.  Old Suppl Fig S4A (keyword CAZy proxy) moves to consolidated Suppl Fig S1.

SOURCE-FIDELITY NOTE for A and B (carried over from build_fig4.py, _job/journal_main_4_7.md)
  The shipped supplementary table `Table_S5a_DRAM_product.tsv` carries ALL-ZERO rows for the
  six MICP-complete MAGs (C22, M1, S13, S16, S23, S26): the hero-only DRAM re-run
  (`scripts/07_dram_hero_annotate.sh`) overwrote the distillate of the full panel and was
  then interrupted, leaving `/data/pangenome_work/dram_output/distillate/product.tsv` - the
  file that was copied into SUBMISSION - with empty hero rows.  The intact full-panel
  distillate survives at `/data/pangenome_work/dram_distilled/product.tsv`; the two files
  are identical for the other 105 MAGs and differed only in those six rows.  This page is
  built from the INTACT file, exactly as build_fig4.py was.  The re-export asked for in that
  journal HAS since landed: as of this build the shipped Table S5a reproduces the intact
  distillate cell for cell, and the assertion below now checks that equality instead of the
  old "differs on exactly the six hero rows" check, which no longer holds.

PROVENANCE NOTE for D (carried over from build_supS4.py)
  The shipped Table_S6c (family counts) contains only the 105 non-MICP-complete MAGs, which
  is why the superseded Figure_S4c rendered as an empty axis - its hero mean was taken over
  an empty set.  The six MICP-complete MAGs were annotated by direct hmmsearch against dbCAN
  (research/revision/dbcan_direct/<MAG>.tbl, best HMM per protein), the same route that
  produced the published hero CLASS counts.  This script re-parses those tables with the
  recipe of research/revision/04b_dbcan_final.py and asserts that the per-class sums
  reproduce Table_S6a exactly for all six MAGs, so the family panel rests on the same
  annotation as the published class panel.  The hero/rest method asymmetry (hmmsearch vs
  DRAM cazy_best_hit) is inherited from that published analysis and is unchanged here.

Sources
  /data/pangenome_work/dram_distilled/product.tsv   DRAM v1.5 distillate, 111 MAGs x 98 modules
  Table_S5a_DRAM_product.tsv                        shipped copy, used only as a check
  Table_S1d_GTDB_Tk_classification.tsv              genus per MAG
  Table_S6a_dbCAN_class_counts.csv                  hero class counts, asserted against
  Table_S6b_dbCAN_class_per1kCDS.csv                per-MAG class hits per 10^3 CDS (C)
  Table_S6c_dbCAN_family_counts.csv                 rest family counts (D)
  Table_S6d_dbCAN_hero_vs_rest.csv                  stored fold change / q / P (C)
  research/revision/dbcan_direct/<MAG>.tbl          hero hmmsearch tables (D)
  research/extra/gene_category_counts.csv           CDS_total per MAG (D)

Colour meanings on this page:
  white -> green   A: mean module completeness (0-1); B: fraction of the genus carrying the
                   module (0-1).  Both ramps encode "more of the module", and each panel
                   names its own quantity on its own colour bar.
  coral            the MICP-complete group (n = 6): its genus row labels in A and B, its
                   bars in C and its dots and mean column in D
  grey             the remaining 105 MAGs
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, REST, TEXT, AXIS, LIGHT, FS_BODY, FS_STAT, HEROES

st.setup()
OUT = HERE / "figures_v2"
BASE = Path("/data/data/Upcycling")
SUPP = BASE / "SUBMISSION/Supplementary_tables"
DRAM_INTACT = Path("/data/pangenome_work/dram_distilled/product.tsv")
DRAM_SHIPPED = SUPP / "Table_S5a_DRAM_product.tsv"
TBL_DIR = BASE / "research/revision/dbcan_direct"
CDS_SRC = BASE / "research/extra/gene_category_counts.csv"

MIN_PER_GENUS = 2          # genera drawn in A and B
MIN_MAGS_WITH_MODULE = 3   # a module is drawn if at least this many MAGs carry it
CLASSES = ["GH", "GT", "PL", "CE", "AA", "CBM"]
N_FAMILIES = 14

# ------------------------------------------------------------------ data: A and B
prod = pd.read_csv(DRAM_INTACT, sep="\t", index_col=0)
shipped = pd.read_csv(DRAM_SHIPPED, sep="\t", index_col=0)

# the re-export has since landed: the shipped Table S5a must now reproduce the intact
# distillate cell for cell, including the six MICP-complete rows that used to be empty
assert list(prod.columns) == list(shipped.columns)
assert list(prod.index) == list(shipped.index)
num_cols = list(prod.select_dtypes("number").columns)
bool_cols = list(prod.select_dtypes(bool).columns)
diff = (prod[num_cols] - shipped.loc[prod.index, num_cols]).abs().sum(axis=1) \
    + (prod[bool_cols].astype(int) - shipped.loc[prod.index, bool_cols].astype(int)).abs().sum(axis=1)
assert (diff == 0).all(), sorted(diff[diff > 0].index)
assert (shipped.loc[HEROES, num_cols].to_numpy().sum() > 0)
assert (shipped.loc[HEROES, bool_cols].to_numpy().sum() > 0)

gtdb = pd.read_csv(SUPP / "Table_S1d_GTDB_Tk_classification.tsv", sep="\t")
gtdb["Genus"] = gtdb.classification.str.extract(r"g__([^;]*)")[0].replace("", np.nan)
genus = gtdb.set_index("user_genome").Genus.fillna("Unclassified")
genus = genus.reindex(prod.index)
assert genus.notna().all()

hero_flag = pd.Series([m in HEROES for m in prod.index], index=prod.index)
HERO_GENERA = sorted(set(genus[hero_flag]))

# modules with signal
frac_mods = [c for c in num_cols
             if prod[c].std() > 0 and (prod[c] > 0).sum() >= MIN_MAGS_WITH_MODULE]
pres_all = [c for c in bool_cols if 0 < prod[c].sum() < len(prod)]
# panel B keeps the two categories the figure is about; the remaining present/absent
# modules (methanogenesis, SCFA, sulfur, other reductases) are in Table S5 but out of scope
pres_mods = [c for c in pres_all if c.startswith(("CAZy:", "Nitrogen metabolism:"))]

n_per_genus = genus.value_counts()
keep_genera = [g for g in n_per_genus.index if n_per_genus[g] >= MIN_PER_GENUS]
# hero-bearing genera first, then the rest alphabetically
order = sorted(keep_genera, key=lambda g: (g not in HERO_GENERA, g))

A = prod.groupby(genus)[frac_mods].mean().loc[order]
B = prod.groupby(genus)[pres_mods].mean().loc[order]      # bool mean = fraction present
row_lbl = [f"{g}  n = {n_per_genus[g]}" for g in order]

# ------------------------------------------------------------------ data: C
per1k = pd.read_csv(SUPP / "Table_S6b_dbCAN_class_per1kCDS.csv").set_index("fasta")
stats_c = pd.read_csv(SUPP / "Table_S6d_dbCAN_hero_vs_rest.csv").set_index("Class")
is_hero = per1k.index.isin(HEROES)
hero_mean = per1k[is_hero][CLASSES].mean()
rest_mean = per1k[~is_hero][CLASSES].mean()
assert is_hero.sum() == len(HEROES)
for c in CLASSES:                      # published means and fold changes must reproduce
    assert abs(hero_mean[c] - stats_c.loc[c, "Hero_mean_per1kCDS"]) < 5e-3
    assert abs(rest_mean[c] - stats_c.loc[c, "Rest_mean_per1kCDS"]) < 5e-3
    assert abs(hero_mean[c] / rest_mean[c] - stats_c.loc[c, "Fold_change"]) < 5e-3
N_HERO, N_REST = int(is_hero.sum()), int((~is_hero).sum())

# ------------------------------------------------------------------ data: D
def best_hmm_per_protein(path):
    best = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        target, hmm, evalue = f[0], f[2], float(f[4])
        if target not in best or evalue < best[target][1]:
            best[target] = (hmm, evalue)
    return best


FAMILY_RE = re.compile(r"(GH|GT|PL|CE|AA|CBM)(\d+.*)?$")
# dbCAN HMM names carry a subfamily suffix (GT2_Glycos_transf_2) while the DRAM-derived
# family table for the other 105 MAGs carries the bare family (GT2).  Both sides are
# reduced to the bare family so the two groups are counted in the same unit.
FAMILY_KEY = re.compile(r"^(GH|GT|PL|CE|AA|CBM)\d+")


def family_key(name):
    m = FAMILY_KEY.match(name)
    return m.group(0) if m else None


hero_fam = defaultdict(lambda: defaultdict(int))
hero_cls = defaultdict(lambda: defaultdict(int))
for mag in HEROES:
    for _, (hmm, _ev) in best_hmm_per_protein(TBL_DIR / f"{mag}.tbl").items():
        name = hmm.replace(".hmm", "")
        m = FAMILY_RE.match(name)
        if m:
            key = family_key(name)
            if key:
                hero_fam[mag][key] += 1
            hero_cls[mag][m.group(1)] += 1

class_counts = pd.read_csv(SUPP / "Table_S6a_dbCAN_class_counts.csv").set_index("fasta")
hero_fam_df = pd.DataFrame(hero_fam).T.fillna(0).astype(int)
for mag in HEROES:                     # the re-parse must reproduce the published counts
    for c in CLASSES:
        assert hero_cls[mag].get(c, 0) == class_counts.loc[mag, c], (mag, c)
    assert hero_fam_df.loc[mag].sum() == class_counts.loc[mag, "Total"]
assert all(family_key(f) == f for f in hero_fam_df.columns)

fam_tbl = pd.read_csv(SUPP / "Table_S6c_dbCAN_family_counts.csv").set_index("fasta")
# the DRAM-derived columns are already bare families; keep only those, so the two groups
# are compared family-for-family
fam_tbl = fam_tbl[[c for c in fam_tbl.columns if family_key(c) == c]]
# Table S6c has since been re-exported and now carries the six MICP-complete MAGs as well
# (it held only the 105 others when build_supS4.py was written).  The hmmsearch re-parse
# above must reproduce those published hero rows family for family; the rest group is then
# the 105 non-MICP-complete rows of the same table.
assert set(HEROES) <= set(fam_tbl.index)
_cols = sorted(set(hero_fam_df.columns) | set(fam_tbl.columns))
assert (hero_fam_df.reindex(columns=_cols, fill_value=0).loc[HEROES].to_numpy()
        == fam_tbl.reindex(columns=_cols, fill_value=0).loc[HEROES].to_numpy()).all()
rest_fam_df = fam_tbl.drop(index=HEROES)
assert len(rest_fam_df) == len(fam_tbl) - len(HEROES)
cds = pd.read_csv(CDS_SRC).set_index("Sample")["CDS_total"]
hero_n1k = hero_fam_df.div(cds.loc[hero_fam_df.index], axis=0) * 1000.0
rest_n1k = rest_fam_df.div(cds.loc[rest_fam_df.index], axis=0) * 1000.0
# per-class sums of the normalised family table must reproduce the panel-C class means
for c in CLASSES:
    fams = [f for f in hero_n1k.columns if FAMILY_KEY.match(f).group(1) == c]
    assert abs(hero_n1k[fams].sum(axis=1).mean() - hero_mean[c]) < 5e-3

fam_all = sorted(set(hero_n1k.columns) | set(rest_n1k.columns))
hero_fam_mean = hero_n1k.reindex(columns=fam_all, fill_value=0).mean()
rest_fam_mean = rest_n1k.reindex(columns=fam_all, fill_value=0).mean()
n_hero_with = (hero_n1k.reindex(columns=fam_all, fill_value=0) > 0).sum()
top = hero_fam_mean.sort_values(ascending=False).head(N_FAMILIES).index.tolist()


# ------------------------------------------------------------------ heat-map helpers
MIN_BLOCK = 4   # a category becomes a block header only if it spans this many columns


def blocks_of(names):
    """Consecutive runs of the same category prefix, as (prefix, first, last)."""
    out, start = [], 0
    pref = [n.split(": ", 1)[0] if ": " in n else "" for n in names]
    for i in range(1, len(names) + 1):
        if i == len(names) or pref[i] != pref[start]:
            out.append((pref[start], start, i - 1))
            start = i
    return out


def header_blocks(names):
    """Only wide runs get a header; a narrow run keeps its prefix on the tick label, so
    that two short headers cannot collide above adjacent one-column groups."""
    return [b for b in blocks_of(names) if b[0] and (b[2] - b[1] + 1) >= MIN_BLOCK]


def tick_labels(names):
    headed = {i for _, i0, i1 in header_blocks(names) for i in range(i0, i1 + 1)}
    return [n.split(": ", 1)[1] if (i in headed and ": " in n) else n
            for i, n in enumerate(names)]


def heat(ax, mat, col_names, cmap, ylabels=True):
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(tick_labels(col_names), rotation=90,
                       ha="center", va="top", fontsize=FS_BODY)
    ax.set_yticks(range(mat.shape[0]))
    if ylabels:
        ax.set_yticklabels(row_lbl, fontsize=FS_BODY)
        for i, g in enumerate(order):
            if g in HERO_GENERA:
                ax.get_yticklabels()[i].set_color(HERO)
    else:
        ax.set_yticklabels([])
    ax.set_xticks(np.arange(mat.shape[1]) - 0.5, minor=True)
    ax.set_yticks(np.arange(mat.shape[0]) - 0.5, minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


# ------------------------------------------------------------------ layout constants
CELL_W = 2.9
LAB_W = 31.0
XA = LAB_W + 3.0
WA = CELL_W * len(frac_mods)
XB = XA + WA + 9.0
WB = CELL_W * len(pres_mods)
# Row pitch: the one-panel-per-page build_fig4.py used 5.0 mm; 4.8 mm - the pitch the
# consolidated Fig S1 settled on - still leaves a 7 pt genus label (~3.0 mm line height)
# clear of its neighbours and brings this page under the 235 mm one-page ceiling.
ROW_H = 4.8
H_HEAT = ROW_H * len(order)
T_HEAT = 10.0
CMAP = st.seq_cmap()


def draw_heatmaps(fig, ax_mm, text_mm):
    axA = ax_mm(XA, T_HEAT, WA, H_HEAT)
    imA = heat(axA, A.values, frac_mods, CMAP)
    axB = ax_mm(XB, T_HEAT, WB, H_HEAT)
    imB = heat(axB, B.values, pres_mods, CMAP, ylabels=False)
    for ax, names, x0 in ((axA, frac_mods, XA), (axB, pres_mods, XB)):
        for pref, i0, i1 in header_blocks(names):
            xm = x0 + (i0 + i1 + 1) / 2 * CELL_W
            text_mm(xm, T_HEAT - 1.2, pref, fontsize=FS_STAT, ha="center", va="bottom",
                    color=TEXT)
    return axA, axB, imA, imB


def label_depth_mm(fig, axes):
    """Millimetres from the top of the page to the bottom of the deepest rotated tick
    label.  Measured with the real renderer rather than estimated, so the colour bars and
    the lower half of the page are placed under the labels as drawn, not under a guess."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    deepest = fig.bbox.height
    for ax in axes:
        for t in ax.get_xticklabels():
            bb = t.get_window_extent(renderer=rend)
            deepest = min(deepest, bb.y0)   # display y grows upward, so y0 is the bottom
    return (fig.bbox.height - deepest) / fig.dpi * 25.4


# pass 1: measure how deep the rotated module names run, on a throwaway page
probe, probe_ax_mm, probe_text_mm, _ = st.page(300.0)
pa, pb, _, _ = draw_heatmaps(probe, probe_ax_mm, probe_text_mm)
LAB_BOTTOM = label_depth_mm(probe, (pa, pb))
st.plt.close(probe)

# ------------------------------------------------------------------ page geometry
CB_Y = LAB_BOTTOM + 4.0            # colour bars sit under the measured label bottom
CB_H = 2.2
T_CD = CB_Y + CB_H + 18.0          # top of the two lower panels

C_LEFT, C_W, C_H = 18.0, 34.0, 44.0
C_STAT_X = [60.0, 71.5, 82.0]
D_LEFT, D_W = 104.0, 28.0
D_H = 4.6 * N_FAMILIES
D_STAT_X = [141.0, 154.0, 167.0]
LEG_Y = T_CD + C_H + 16.0          # shared key, in the free band under C
H = T_CD + max(C_H, D_H) + 12.0
assert H <= 235.0, H             # single-page ceiling

fig, ax_mm, text_mm, letter = st.page(H)
axA, axB, imA, imB = draw_heatmaps(fig, ax_mm, text_mm)
assert abs(label_depth_mm(fig, (axA, axB)) - LAB_BOTTOM) < 0.5


def colourbar(im, l, t, w, h, label):
    cax = ax_mm(l, t, w, h)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.4)
    cb.set_ticks([0, 0.5, 1])
    cax.tick_params(labelsize=FS_STAT, length=1.5, pad=1.2)
    cax.set_xlabel(label, fontsize=FS_STAT, labelpad=1.5)


colourbar(imA, XA, CB_Y, 26, CB_H, "Mean module completeness")
colourbar(imB, XB, CB_Y, 26, CB_H, "Fraction of genus with module")

# ------------------------------------------------------------------ C: dbCAN classes
axC = ax_mm(C_LEFT, T_CD, C_W, C_H)
y = np.arange(len(CLASSES))
bh = 0.36
axC.barh(y - bh / 2, hero_mean[CLASSES].values, height=bh, color=HERO, zorder=3)
axC.barh(y + bh / 2, rest_mean[CLASSES].values, height=bh, color=REST, zorder=3)
for i, c in enumerate(CLASSES):
    axC.text(hero_mean[c] + 0.4, i - bh / 2, f"{hero_mean[c]:.2f}", va="center",
             ha="left", fontsize=FS_STAT, color=HERO)
    axC.text(rest_mean[c] + 0.4, i + bh / 2, f"{rest_mean[c]:.2f}", va="center",
             ha="left", fontsize=FS_STAT, color=TEXT)
axC.set_yticks(y)
axC.set_yticklabels(CLASSES, fontsize=FS_BODY)
axC.invert_yaxis()
axC.set_xlabel("mean hits per 10³ CDS")
axC.set_xlim(0, float(hero_mean.max()) * 1.45)
st.style_axis(axC, left=False)
axC.tick_params(axis="y", length=0)
for x, head in zip(C_STAT_X, ["fold\nchange", "permut.\nq (BH)", "MWU\nP"]):
    text_mm(x, T_CD - 3.0, head, ha="center", va="bottom", fontsize=FS_STAT,
            color=AXIS, linespacing=1.25)
for i, c in enumerate(CLASSES):
    yy = T_CD + (i + 0.5) * (C_H / len(CLASSES))
    q = stats_c.loc[c, "Permutation_q_BH"]
    text_mm(C_STAT_X[0], yy, f"{stats_c.loc[c, 'Fold_change']:.2f}", ha="center",
            va="center", fontsize=FS_STAT)
    text_mm(C_STAT_X[1], yy, f"{q:.4f}", ha="center", va="center", fontsize=FS_STAT,
            color=HERO if q < 0.05 else TEXT)
    text_mm(C_STAT_X[2], yy, f"{stats_c.loc[c, 'MannWhitney_greater_p']:.4f}",
            ha="center", va="center", fontsize=FS_STAT)

# ------------------------------------------------------------------ D: dbCAN families
axD = ax_mm(D_LEFT, T_CD, D_W, D_H)
yd = np.arange(len(top))
for i, f in enumerate(top):
    axD.plot([rest_fam_mean[f], hero_fam_mean[f]], [i, i], color=LIGHT, lw=2.4,
             solid_capstyle="round", zorder=1)
axD.scatter(rest_fam_mean[top].values, yd, s=16, color=REST, zorder=3)
axD.scatter(hero_fam_mean[top].values, yd, s=16, color=HERO, zorder=4)
axD.set_yticks(yd)
axD.set_yticklabels(top, fontsize=FS_BODY)
axD.invert_yaxis()
axD.set_xlabel("mean hits per 10³ CDS")
axD.set_xlim(-0.05, float(hero_fam_mean[top].max()) * 1.12)
st.style_axis(axD, left=False)
axD.tick_params(axis="y", length=0)
for x, head in zip(D_STAT_X, ["MICP-compl.\nmean", "rest\nmean",
                              "MICP-compl.\nMAGs with\nfamily"]):
    text_mm(x, T_CD - 3.0, head, ha="center", va="bottom", fontsize=FS_STAT,
            color=AXIS, linespacing=1.25)
for i, f in enumerate(top):
    yy = T_CD + (i + 0.5) * (D_H / len(top))
    text_mm(D_STAT_X[0], yy, f"{hero_fam_mean[f]:.2f}", ha="center", va="center",
            fontsize=FS_STAT, color=HERO)
    text_mm(D_STAT_X[1], yy, f"{rest_fam_mean[f]:.2f}", ha="center", va="center",
            fontsize=FS_STAT)
    text_mm(D_STAT_X[2], yy, f"{int(n_hero_with[f])} / {len(HEROES)}", ha="center",
            va="center", fontsize=FS_STAT)

# one key for both lower panels; an in-axes legend would sit on top of the shortest
# bars in C and the tied dots in D
fig.legend(handles=[Line2D([], [], color=HERO, lw=4, label=f"MICP-complete n = {N_HERO}"),
                    Line2D([], [], color=REST, lw=4, label=f"rest n = {N_REST}")],
           loc="upper left", bbox_to_anchor=(18.0 / st.PAGE_W_MM, 1 - LEG_Y / H),
           ncol=2, frameon=False, fontsize=FS_BODY, handlelength=1.4,
           handletextpad=0.8, columnspacing=2.4)

letter(4, 4, "A")
letter(XB - 5.0, 4, "B")
letter(4, T_CD - 6.0, "C")
letter(92, T_CD - 6.0, "D")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S2")
