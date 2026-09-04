"""Consolidated main Fig 1 - phylogenomic distribution and completeness of the MICP module.

Consolidation of 2026-09-04 (consolidation_260904/DESIGN.md): old main Fig 1 and old main
Fig 2 become one page, and the page must fit inside a single 180 x 235 mm journal page.

Panels (reading order, left to right then top to bottom):
  A  maximum-likelihood bac120 tree (midpoint-rooted for display), a genus colour strip,
     and the MAG identifier with its GTDB species where one was assigned   [old Fig 1A]
  B  presence of ureA-G and cah for the same rows                          [old Fig 1B]
  C  MICP module score (ureA-G + cah, 0-8) per GTDB-Tk genus, box + jittered points,
     the six MICP-complete MAGs overplotted as coral rings                 [old Fig 2A]
  D  per-gene prevalence, MICP-complete group (n = 6) vs the rest (n = 105)[old Fig 2B]

A and B share one row order and one row height and stay side by side exactly as the old
Fig 1 drew them.  All 111 tips are kept.  To reach the 235 mm page ceiling the tip rows
are carried in TWO adjacent tree columns instead of one: the left column holds the top
part of the tip order and the right column continues it, at an identical row pitch and an
identical branch-length scale (one scale bar, under the left tree, applies to both).  The
cut is not arbitrary - CUT_LO..CUT_HI is searched for the tip boundary crossed by the
fewest clades, so the split falls at a near-natural break in the topology and only a
handful of backbone edges run off a column edge.  Panel letters A and B sit over the left
column and label the same two element types in the right column.

Sources
  pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile   IQ-TREE topology
                                    and branch lengths; tip names carry the MAG id and
                                    the GTDB species
  Table_S1d_GTDB_Tk_classification.tsv    genus and species per MAG
  MAGs_FASTA_files/bakta_results/*/*.tsv  gene presence, via _micp_presence.presence():
                                    a CDS-only keyword scan of all 111 annotations
  Table_S15a_alkaliphile_signature_per_MAG.csv   the MICP-complete / rest group flag

Provenance note (see _job/JOURNAL.md).  Two earlier sources for panels B-D are unusable.
pangenome_work/MICP_Pangenome_Final_Summary.csv, used by the shipped figure, holds only
100 of the 111 MAGs - C1 and C10-C19 are missing and were drawn as all-absent.
Table_S1a_ace_samples_list.csv, used by the first version of this rebuild, lists only 45
MAGs (13 of the 66 it omits carry a ure CDS) and counts the Bakta 5_ureB_sRNA non-coding
feature as a copy of the beta subunit, which is why it scores S26 8/8 where the Bakta
annotation has no protein-coding ureB.  Panels B-D are therefore built from the
annotations directly; _micp_presence documents both defects and asserts the relationship
to S1a.  presence() is called once here, so A/B and C/D cannot disagree.

Colour meanings on this page (one meaning per colour):
  genus strip = one colour per GTDB genus; Sphingobacterium blue and Pseudomonas_E orange
                are the two MICP-complete lineages, every other genus takes a muted hue
                and no genus is given green or coral
  green       = the gene is present in that MAG (panel B); white = absent
  coral       = the MICP-complete group (tip labels in A, rings in C, bars in D)
  grey        = the remaining 105 MAGs (bars in D)
  black       = tree branches, medians and individual MAG points
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.patches import Patch, Rectangle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, REST, SPHINGO, PSEUDO, GREEN, TEXT, GREY, AXIS, LIGHT, \
    FS_BODY, FS_STAT
from _micp_presence import presence

st.setup()
OUT = HERE / "figures_v2"
BASE = Path("/data/data/Upcycling")
SUPP = BASE / "SUBMISSION/Supplementary_tables"
TREE = BASE / "pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile"

GENES = ["ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG", "cah"]
# genus strip: the two MICP-complete lineages keep the lineage colours; every other genus
# takes a muted hue chosen to avoid green (gene present) and coral (MICP-complete)
GENUS_COL = {"Sphingobacterium": SPHINGO, "Pseudomonas_E": PSEUDO,
             "Stenotrophomonas": "#9B8AB8", "Acinetobacter": "#C2A878",
             "Achromobacter": "#B77BA8", "Comamonas": "#5FA8A0",
             "Chryseobacterium": "#8C8C8C", "Alcaligenes": "#C9C24D",
             "Paraburkholderia": "#7C5C3E"}
OTHER_COL = "#D9D9D9"
OTHER_LABEL = "Other genera"

# ------------------------------------------------------------------ data
tax = pd.read_csv(SUPP / "Table_S1d_GTDB_Tk_classification.tsv", sep="\t",
                  index_col="user_genome")
genus = (tax["classification"].str.extract(r"g__([^;]*)")[0]
         .replace("", np.nan).fillna("Unclassified"))
species = tax["classification"].str.extract(r"s__([^;]*)")[0].fillna("")
PANEL = sorted(tax.index)

grp = pd.read_csv(SUPP / "Table_S15a_alkaliphile_signature_per_MAG.csv", index_col=0)
heroes_list = sorted(grp.index[grp["group"] == "MICP_complete"])
heroes = set(heroes_list)
assert heroes_list == sorted(st.HEROES), heroes_list
assert len(genus) == len(grp) == 111, (len(genus), len(grp))
n_hero, n_rest = len(heroes_list), len(PANEL) - len(heroes_list)

tree = Phylo.read(TREE, "newick")
tree.root_at_midpoint()
tips = [t.name for t in tree.get_terminals()]
mag_of = {t: t.split("_s__")[0] for t in tips}
assert len(tips) == 111 == len(tax), (len(tips), len(tax))
assert set(mag_of.values()) == set(tax.index)

pres = presence().reindex(PANEL)[GENES]
assert pres.notna().all().all()
score = pres.sum(axis=1)
# five of the six MICP-complete MAGs carry the full module; S26 has no protein-coding
# urease beta subunit in its Bakta annotation (see _micp_presence)
assert (score[heroes_list] == len(GENES)).sum() == 5, score[heroes_list].to_dict()
assert pres.loc["S26", "ureB"] == 0 and score["S26"] == len(GENES) - 1

is_hero = pd.Series(pres.index.isin(heroes_list), index=pres.index)
prev_hero = pres[is_hero.values].mean() * 100
prev_rest = pres[~is_hero.values].mean() * 100
assert len(pres[is_hero.values]) == n_hero
assert len(pres[~is_hero.values]) == n_rest

# tip order as drawn, top to bottom
ordered = []


def collect(clade):
    if clade.is_terminal():
        ordered.append(clade.name)
    else:
        for c in clade.clades:
            collect(c)


collect(tree.root)
assert len(ordered) == len(tips)
y_of = {name: i for i, name in enumerate(ordered)}

depths = tree.depths()
xmax = max(depths.values())

gcount_all = genus.value_counts()
strip_genera = [g for g in GENUS_COL if g in set(genus)]
assert set(strip_genera) == set(gcount_all.head(len(GENUS_COL)).index), \
    (sorted(strip_genera), sorted(gcount_all.head(len(GENUS_COL)).index))
n_other = int((~genus.isin(strip_genera)).sum())

# ---- where to break the tip order into two columns.  Every clade covers a contiguous
# block of the drawn tip order; a boundary crossed by few clades cuts few branches.
CUT_LO, CUT_HI = 50, 62      # layout constants: keeps both columns near half height
clade_span = []
for cl in tree.get_nonterminals():
    ii = [y_of[x.name] for x in cl.get_terminals()]
    clade_span.append((min(ii), max(ii)))


def n_crossing(i):
    return sum(1 for lo, hi in clade_span if lo <= i - 1 and hi >= i)


SPLIT = min(range(CUT_LO, CUT_HI + 1),
            key=lambda i: (n_crossing(i), abs(i - len(ordered) / 2)))
BLOCKS = [(0, SPLIT), (SPLIT, len(ordered))]
assert sum(hi - lo for lo, hi in BLOCKS) == len(ordered) == 111
n_rows_max = max(hi - lo for lo, hi in BLOCKS)

# panel C grouping
df = pd.DataFrame({"genus": genus.reindex(PANEL), "score": score})
gcount = df["genus"].value_counts()
TOP_N = 9  # style constant: genera drawn individually; the remainder are pooled
top = list(gcount.head(TOP_N).index)
df["grp"] = np.where(df["genus"].isin(top), df["genus"], "Other genera")
order_c = df.groupby("grp")["score"].mean().sort_values(ascending=False).index.tolist()

# ------------------------------------------------------------------ page geometry
FS_TIP = 5.5                     # tip labels only; 111 rows on one page need a tighter
                                 # pitch than the 7 pt body size allows
ROWS_H = 119.0                   # height of the taller of the two tree columns
ROW = ROWS_H / n_rows_max        # mm per tip, identical in both columns
TOP = 11.0
COL_X = [4.0, 92.0]              # left edge of each tree column block
TREE_W = 19.0
STRIP_DX, STRIP_W = 20.0, 2.2
LAB_DX, LAB_W = 23.4, 40.0
HEAT_DX, HEAT_W = 64.0, 20.0

AB_END = TOP + ROWS_H            # bottom of the taller column
LEG_Y = AB_END + 10.0            # top of the A/B key (clears the tree scale bar)
CD_LET_Y = LEG_Y + 16.0          # panel letters of the second row
CD_TOP = CD_LET_Y + 7.0          # top of the C / D axes
CD_H = 50.0
H = CD_TOP + CD_H + 14.0
assert H <= 235.0, H             # single-page ceiling

fig, ax_mm, text_mm, letter = st.page(H)


def fx(x_mm):
    return x_mm / st.PAGE_W_MM


def fy(y_mm):
    return 1.0 - y_mm / H


letter(COL_X[0], 6.0, "A")
letter(COL_X[0] + HEAT_DX - 5.0, 6.0, "B")
letter(4.0, CD_LET_Y, "C")
letter(102.0, CD_LET_Y, "D")


# ---- A: tree, drawn in full in each column and clipped to that column's tip range
def draw(ax, clade, x_parent):
    x = depths[clade]
    if clade.is_terminal():
        y = y_of[clade.name]
        ax.plot([x_parent, x], [y, y], color=TEXT, lw=0.4, solid_capstyle="butt")
        return y
    ys = [draw(ax, c, x) for c in clade.clades]
    ax.plot([x, x], [min(ys), max(ys)], color=TEXT, lw=0.4, solid_capstyle="butt")
    ymid = (min(ys) + max(ys)) / 2
    ax.plot([x_parent, x], [ymid, ymid], color=TEXT, lw=0.4, solid_capstyle="butt")
    return ymid


tip_labels = []
for k, (lo, hi) in enumerate(BLOCKS):
    x0 = COL_X[k]
    h = (hi - lo) * ROW
    ax_t = ax_mm(x0, TOP, TREE_W, h)
    ax_s = ax_mm(x0 + STRIP_DX, TOP, STRIP_W, h)
    ax_l = ax_mm(x0 + LAB_DX, TOP, LAB_W, h)
    ax_h = ax_mm(x0 + HEAT_DX, TOP, HEAT_W, h)

    draw(ax_t, tree.root, 0.0)
    ax_t.set_xlim(-xmax * 0.02, xmax * 1.02)
    ax_t.set_ylim(hi - 0.5, lo - 0.5)
    ax_t.set_xticks([])
    ax_t.set_yticks([])
    for s in ax_t.spines.values():
        s.set_visible(False)

    if k == 0:
        # scale bar: a round substitutions-per-site distance, under the left tree only;
        # both columns share one branch-length scale
        step = 10 ** np.floor(np.log10(xmax / 4))
        bar = float(max(kk * step for kk in (1, 2, 5) if kk * step <= xmax / 3))
        y_bar = hi + 1.6 / ROW
        ax_t.plot([0, bar], [y_bar, y_bar], color=TEXT, lw=0.9, clip_on=False)
        ax_t.text(0, y_bar + 1.2 / ROW, f"{bar:g} substitutions/site",
                  ha="left", va="top", fontsize=FS_STAT, color=TEXT, clip_on=False)

    for name in ordered[lo:hi]:
        mag = mag_of[name]
        y = y_of[name]
        g = genus[mag]
        ax_s.add_patch(Rectangle((0, y - 0.5), 1, 1,
                                 facecolor=GENUS_COL.get(g, OTHER_COL), edgecolor="none"))
        sp = species[mag]
        label = f"{mag}  {sp}" if sp else mag
        tip_labels.append(
            ax_l.text(0, y, label, va="center", ha="left", fontsize=FS_TIP,
                      color=HERO if mag in heroes else TEXT,
                      fontweight="bold" if mag in heroes else "normal"))

    for ax in (ax_s, ax_l):
        ax.set_xlim(0, 1)
        ax.set_ylim(hi - 0.5, lo - 0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    # ---- B: gene presence for the same rows
    mat = np.array([pres.loc[mag_of[n], GENES].values for n in ordered[lo:hi]],
                   dtype=float)
    ax_h.imshow(mat, aspect="auto", cmap=st.seq_cmap(), vmin=0, vmax=1,
                extent=[0, len(GENES), hi - 0.5, lo - 0.5], interpolation="nearest")
    for kk in range(len(GENES) + 1):
        ax_h.plot([kk, kk], [lo - 0.5, hi - 0.5], color=LIGHT, lw=0.3)
    ax_h.set_xticks(np.arange(len(GENES)) + 0.5)
    ax_h.set_xticklabels(GENES, rotation=90, fontstyle="italic", fontsize=FS_STAT)
    ax_h.xaxis.set_ticks_position("top")
    ax_h.tick_params(axis="x", length=0, pad=1.5)
    ax_h.set_yticks([])
    ax_h.set_ylim(hi - 0.5, lo - 0.5)
    for s in ax_h.spines.values():
        s.set_visible(False)

# ---- key for A and B, between the two rows
handles = [Patch(facecolor=GENUS_COL[g], label=f"{g} ({int(gcount_all[g])})")
           for g in sorted(strip_genera, key=lambda g: -gcount_all[g])]
handles.append(Patch(facecolor=OTHER_COL, label=f"{OTHER_LABEL} ({n_other})"))
handles.append(Patch(facecolor=GREEN, label="gene present"))
handles.append(Patch(facecolor="white", edgecolor=LIGHT, lw=0.5, label="gene absent"))
handles.append(Patch(facecolor="none", edgecolor="none", label="MICP-complete MAG"))
leg = fig.legend(handles=handles, loc="upper left", ncol=4, fontsize=FS_STAT,
                 frameon=False, bbox_to_anchor=(fx(8.0), fy(LEG_Y)), handlelength=1.1,
                 columnspacing=1.2, handletextpad=0.45, labelspacing=0.4)
for txt in leg.get_texts():
    if txt.get_text() == "MICP-complete MAG":
        txt.set_color(HERO)
        txt.set_fontweight("bold")

# ---- C: module score by genus (horizontal; genus names would collide when rotated)
axC = ax_mm(34.0, CD_TOP, 62.0, CD_H)
order_c = order_c[::-1]  # highest mean at the top of a horizontal axis
data = [df.loc[df["grp"] == g, "score"].values for g in order_c]
bp = axC.boxplot(data, positions=range(len(order_c)), widths=0.62, patch_artist=True,
                 vert=False,
                 medianprops=dict(color=TEXT, lw=0.9),
                 whiskerprops=dict(color=GREY, lw=0.7),
                 capprops=dict(color=GREY, lw=0.7),
                 flierprops=dict(marker="", markersize=0))
for patch in bp["boxes"]:
    patch.set_facecolor(LIGHT)
    patch.set_edgecolor(GREY)
    patch.set_linewidth(0.7)

rng = np.random.default_rng(0)  # jitter only; carries no data
for i, g in enumerate(order_c):
    sub = df[df["grp"] == g]
    y = i + rng.normal(0, 0.13, size=len(sub))
    hero_mask = sub.index.isin(heroes_list)
    axC.scatter(sub["score"].values[~hero_mask], y[~hero_mask], s=5,
                color=TEXT, alpha=0.55, lw=0, zorder=3)
    if hero_mask.any():
        axC.scatter(sub["score"].values[hero_mask], y[hero_mask], s=26,
                    facecolor="none", edgecolor=HERO, lw=1.1, zorder=4)

axC.set_yticks(range(len(order_c)))
axC.set_yticklabels([f"{g} ({int((df['grp'] == g).sum())})" for g in order_c],
                    fontsize=FS_BODY)
axC.set_xlabel("MICP module score (max 8)")
axC.set_xlim(-0.4, 8.6)
axC.set_xticks([0, 2, 4, 6, 8])
axC.set_ylim(-0.7, len(order_c) - 0.3)
st.style_axis(axC)
ring = Patch(facecolor="none", edgecolor=HERO, lw=1.1,
             label=f"MICP-complete (n = {n_hero})")
axC.legend(handles=[ring], loc="upper left", bbox_to_anchor=(0.0, -0.16),
           fontsize=FS_STAT, handlelength=1.0, borderpad=0.2)

# ---- D: per-gene prevalence
axD = ax_mm(114.0, CD_TOP, 58.0, CD_H)
x = np.arange(len(GENES))
w = 0.38
axD.bar(x - w / 2, prev_hero[GENES].values, w, color=HERO,
        label=f"MICP-complete (n = {n_hero})")
axD.bar(x + w / 2, prev_rest[GENES].values, w, color=REST,
        label=f"Rest (n = {n_rest})")
for xi, v in zip(x - w / 2, prev_hero[GENES].values):
    axD.text(xi, v + 1.5, f"{v:.0f}", ha="center", va="bottom", fontsize=FS_STAT,
             color=TEXT)
for xi, v in zip(x + w / 2, prev_rest[GENES].values):
    axD.text(xi, v + 1.5, f"{v:.0f}", ha="center", va="bottom", fontsize=FS_STAT,
             color=TEXT)
axD.set_xticks(x)
axD.set_xticklabels(GENES, rotation=45, ha="right", fontstyle="italic")
axD.set_ylabel("MAGs with the gene (%)")
axD.set_ylim(0, 112)
axD.set_yticks([0, 25, 50, 75, 100])
st.style_axis(axD)
axD.legend(loc="upper right", bbox_to_anchor=(1.02, 1.12), fontsize=FS_STAT,
           handlelength=1.0, borderpad=0.2)

# the tip label column is the one slot st.audit cannot police (a label that runs past it
# would sit over the presence matrix, where there is no other text to collide with), so
# its rendered width is measured against the slot it was given
fig.canvas.draw()
w_mm = max(t.get_window_extent(renderer=fig.canvas.get_renderer()).width
           for t in tip_labels) / fig.dpi * 25.4
assert w_mm <= LAB_W, (w_mm, LAB_W)

print(f"  page height {H:.1f} mm | tip columns {[hi - lo for lo, hi in BLOCKS]} "
      f"split at tip {SPLIT} ({n_crossing(SPLIT)} clades crossed) | row pitch "
      f"{ROW:.2f} mm | widest tip label {w_mm:.1f} of {LAB_W:.1f} mm")
st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig1")
