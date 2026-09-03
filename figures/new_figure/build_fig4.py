"""Main Fig 4 - DRAM metabolic reconstruction across the 111-MAG panel.

One composed page, 180 mm wide, three panels:
  A  genus-aggregated mean completeness of the DRAM modules that are scored as a fraction
  B  genus-aggregated prevalence of the CAZy and nitrogen-metabolism modules, which DRAM
     scores as present / absent rather than as a fraction
  C  MICP-critical trait modules, mean hits per 10^3 CDS, MICP-complete vs rest

SOURCE-FIDELITY NOTE (see _job/journal_main_4_7.md)
  The shipped supplementary table `Table_S5a_DRAM_product.tsv` carries ALL-ZERO rows for the
  six MICP-complete MAGs (C22, M1, S13, S16, S23, S26): the hero-only DRAM re-run
  (`scripts/07_dram_hero_annotate.sh`) overwrote the distillate of the full panel and was
  then interrupted, leaving `/data/pangenome_work/dram_output/distillate/product.tsv` - the
  file that was copied into SUBMISSION - with empty hero rows.  The intact full-panel
  distillate survives at `/data/pangenome_work/dram_distilled/product.tsv`; the two files
  are identical for the other 105 MAGs and differ only in those six rows.  This page is
  built from the INTACT file.  Table S5a must be re-exported from it before submission.

Sources
  /data/pangenome_work/dram_distilled/product.tsv   DRAM v1.5 distillate, 111 MAGs x 98 modules
  Table_S1d_GTDB_Tk_classification.tsv              genus per MAG
  Table_S2a_trait_module_counts.csv                 Bakta keyword hits and CDS_total (panel C)
  Table_S2b_trait_module_per1kCDS.csv               stored normalisation, used as a check
  Table_S2c_permutation_statistics.csv              stored group means, used as a check

Colour meanings on this page:
  white -> green   panel A: mean module completeness (0-1); panel B: fraction of the genus
                   carrying the module (0-1).  Both ramps encode "more of the module", and
                   each panel names its own quantity on its own colour bar.
  coral            the MICP-complete group (n = 6): its genus row labels in A and B, and its
                   bars in C
  grey             the remaining MAGs
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, REST, TEXT, GREY, FS_BODY, FS_STAT, HEROES

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
DRAM_INTACT = Path("/data/pangenome_work/dram_distilled/product.tsv")
DRAM_SHIPPED = SUPP / "Table_S5a_DRAM_product.tsv"

MIN_PER_GENUS = 2          # genera drawn in A and B
MIN_MAGS_WITH_MODULE = 3   # a module is drawn if at least this many MAGs carry it

# ------------------------------------------------------------------ data
prod = pd.read_csv(DRAM_INTACT, sep="\t", index_col=0)
shipped = pd.read_csv(DRAM_SHIPPED, sep="\t", index_col=0)

# prove the clobbering claim rather than asserting it in prose: the two files agree
# everywhere except on the six MICP-complete rows, which are empty in the shipped copy
assert list(prod.columns) == list(shipped.columns)
num_cols = list(prod.select_dtypes("number").columns)
bool_cols = list(prod.select_dtypes(bool).columns)
diff = (prod[num_cols] - shipped.loc[prod.index, num_cols]).abs().sum(axis=1) \
    + (prod[bool_cols].astype(int) - shipped.loc[prod.index, bool_cols].astype(int)).abs().sum(axis=1)
assert set(diff[diff > 0].index) == set(HEROES), sorted(diff[diff > 0].index)
assert (shipped.loc[HEROES, num_cols].sum(axis=1) == 0).all()
assert (shipped.loc[HEROES, bool_cols].sum(axis=1) == 0).all()

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

# panel C - MICP-critical trait modules, recomputed from raw counts
counts = pd.read_csv(SUPP / "Table_S2a_trait_module_counts.csv").set_index("Sample")
per1k_stored = pd.read_csv(SUPP / "Table_S2b_trait_module_per1kCDS.csv").set_index("Sample")
sub_cols = [c for c in counts.columns if "::" in c]
per1k = counts[sub_cols].div(counts["CDS_total"], axis=0) * 1000.0
assert np.allclose(per1k.values, per1k_stored.loc[per1k.index, sub_cols].values)
is_hero_c = counts["Hero"].astype(bool)
assert set(counts.index[is_hero_c]) == set(HEROES)

CRIT = [("Alkaline_Osmo::Mrp_complex", "Mrp Na+/H+ antiporter"),
        ("Alkaline_Osmo::Na_H_antiporter", "nhaA-C antiporter"),
        ("Alkaline_Osmo::oxidative", "Oxidative-stress defence"),
        ("Alkaline_Osmo::compatible_solute", "Compatible-solute synthesis"),
        ("CAZyme_proxy::glycoside_hydrolase", "Glycoside hydrolase"),
        ("CAZyme_proxy::carb_binding", "Carbohydrate-binding module"),
        ("Ammonia_N::urea_transport", "Urea transport, urtABC"),
        ("Ammonia_N::GS_GOGAT", "GS-GOGAT, glnA gltBD"),
        ("Biofilm_EPS::quorum", "Quorum sensing"),
        ("Biofilm_EPS::cellulose", "Cellulose synthesis, bcs")]
keys = [k for k, _ in CRIT]
crit_lbl = [l for _, l in CRIT]
hero_mean = per1k.loc[is_hero_c, keys].mean()
rest_mean = per1k.loc[~is_hero_c, keys].mean()

# assert the two group means against the permutation table, which stores them rounded
perm = pd.read_csv(SUPP / "Table_S2c_permutation_statistics.csv")
perm["key"] = perm.Category + "::" + perm.Subcategory
perm = perm.set_index("key")
for k in keys:
    assert abs(round(hero_mean[k], 3) - perm.loc[k, "Hero_mean_per1kCDS"]) < 1e-9, k
    assert abs(round(rest_mean[k], 3) - perm.loc[k, "Rest_mean_per1kCDS"]) < 1e-9, k

N_HERO, N_REST = int(is_hero_c.sum()), int((~is_hero_c).sum())


# ------------------------------------------------------------------ helpers
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


def heat(ax, mat, col_names, cmap):
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(tick_labels(col_names), rotation=90,
                       ha="center", va="top", fontsize=FS_BODY)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(row_lbl, fontsize=FS_BODY)
    for i, g in enumerate(order):
        if g in HERO_GENERA:
            ax.get_yticklabels()[i].set_color(HERO)
    ax.set_xticks(np.arange(mat.shape[1]) - 0.5, minor=True)
    ax.set_yticks(np.arange(mat.shape[0]) - 0.5, minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def colourbar(fig, im, l, t, w, h, label):
    cax = ax_mm(l, t, w, h)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.4)
    cb.set_ticks([0, 0.5, 1])
    cax.tick_params(labelsize=FS_STAT, length=1.5, pad=1.2)
    cax.set_xlabel(label, fontsize=FS_STAT, labelpad=1.5)


# ------------------------------------------------------------------ layout
CELL_W = 2.9
LAB_W = 31.0
XA = LAB_W + 3.0
WA = CELL_W * len(frac_mods)
XB = XA + WA + 9.0
WB = CELL_W * len(pres_mods)
ROW_H = 5.0
H_HEAT = ROW_H * len(order)
T_HEAT = 10.0
LAB_H = 96.0                       # rotated module names below the heat maps; the
                                   # colour bars are then placed under the MEASURED bottom
                                   # of those labels, not under this estimate
T_C = T_HEAT + H_HEAT + LAB_H + 12.0
H_C = len(CRIT) * 5.6
H = T_C + H_C + 22.0

fig, ax_mm, text_mm, letter = st.page(H)
cmap = st.seq_cmap()

axA = ax_mm(XA, T_HEAT, WA, H_HEAT)
imA = heat(axA, A.values, frac_mods, cmap)
axB = ax_mm(XB, T_HEAT, WB, H_HEAT)
imB = heat(axB, B.values, pres_mods, cmap)
axB.set_yticklabels([])

# category block headers above the column groups
for ax, names, x0, w in ((axA, frac_mods, XA, WA), (axB, pres_mods, XB, WB)):
    for pref, i0, i1 in header_blocks(names):
        xm = x0 + (i0 + i1 + 1) / 2 * CELL_W
        text_mm(xm, T_HEAT - 1.2, pref, fontsize=FS_STAT, ha="center", va="bottom",
                color=TEXT)

# place the colour bars below the deepest rotated tick label.  st.audit compares text
# against text only, so a bar drawn on top of a long label would pass unnoticed.
fig.canvas.draw()
rend = fig.canvas.get_renderer()
deepest = fig.bbox.height
for ax in (axA, axB):
    for t in ax.get_xticklabels():
        bb = t.get_window_extent(renderer=rend)
        deepest = min(deepest, bb.y0)          # display y grows upward, so y0 is the bottom
CB_Y = (fig.bbox.height - deepest) / fig.dpi * 25.4 + 4.0
assert CB_Y + 8.0 < H, (CB_Y, H)
colourbar(fig, imA, XA, CB_Y, 26, 2.2, "Mean module completeness")
colourbar(fig, imB, XB, CB_Y, 26, 2.2, "Fraction of genus with module")

# panel C
axC = ax_mm(XA + 22.0, T_C, 96.0, H_C)
y = np.arange(len(keys))
bh = 0.38
axC.barh(y - bh / 2, hero_mean[keys].values, bh, color=HERO, label=f"MICP-complete n = {N_HERO}")
axC.barh(y + bh / 2, rest_mean[keys].values, bh, color=REST, label=f"Rest n = {N_REST}")
axC.set_yticks(y)
axC.set_yticklabels(crit_lbl, fontsize=FS_BODY)
axC.invert_yaxis()
axC.set_xscale("log")
# a Unicode superscript keeps the label in the page font; mathtext would emit a
# second font-family into the SVG and break the Arial-first rule
axC.set_xlabel("Gene hits per 10\u00b3 CDS, Bakta keyword scan (log scale)")
axC.tick_params(axis="y", length=0)
st.style_axis(axC, left=False, bottom=True)
axC.xaxis.set_major_locator(FixedLocator([0.001, 0.01, 0.1, 1, 10, 100]))
axC.xaxis.set_minor_locator(NullLocator())
axC.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
axC.set_xlim(min(rest_mean.min(), hero_mean.min()) * 0.5, 100)

# the two group means sit in a stat column rather than beside the bars, where the two
# values of a near-tied row would collide
X_STAT = XA + 22.0 + 96.0 + 4.0
text_mm(X_STAT, T_C - 1.5, "MICP-compl.", fontsize=FS_STAT, ha="left", va="bottom")
text_mm(X_STAT + 15.0, T_C - 1.5, "rest", fontsize=FS_STAT, ha="left", va="bottom")
for yi, k in enumerate(keys):
    yy = T_C + (yi + 0.5) * (H_C / len(keys))
    text_mm(X_STAT, yy, f"{hero_mean[k]:.2f}", fontsize=FS_STAT, ha="left", va="center",
            color=HERO)
    text_mm(X_STAT + 15.0, yy, f"{rest_mean[k]:.2f}", fontsize=FS_STAT, ha="left",
            va="center", color=REST)
axC.legend(handles=[Patch(color=HERO, label=f"MICP-complete n = {N_HERO}"),
                    Patch(color=REST, label=f"Rest n = {N_REST}")],
           loc="lower right", fontsize=FS_BODY, handlelength=1.2, handletextpad=0.4,
           borderpad=0.2, labelspacing=0.25)

letter(4, 4, "A")
letter(XB - 5.0, 4, "B")
letter(4, T_C - 6.0, "C")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig4")
