"""Consolidated main Fig 1 - phylogenomic distribution and completeness of the MICP module.

Consolidation of 2026-09-04 (consolidation_260904/DESIGN.md): old main Fig 1 and old main
Fig 2 become one page, and the page must fit inside a single 180 x 235 mm journal page.

Panels (reading order, left to right then top to bottom):
  A  maximum-likelihood bac120 tree (midpoint-rooted for display), a genus colour strip,
     the MAG identifier with its GTDB species where one was assigned, and the presence
     of ureA-G and cah as concentric rings                        [old Fig 1A + 1B]
  B  MICP module score (ureA-G + cah, 0-8) per GTDB-Tk genus, box + jittered points,
     the six MICP-complete MAGs overplotted as coral rings                 [old Fig 2A]
  C  per-gene prevalence, MICP-complete group (n = 6) vs the rest (n = 105)[old Fig 2B]

Revision of 2026-09-09 (second round): the tree and the presence rings are ONE panel
(A); the earlier split into A and B put a second letter inside the circle.  The keys sit
to the right of the circle, titled, and the circle is shifted left to make room.

A and B are drawn as ONE circular layout (revision of 2026-09-09): the tree is a radial
phylogram from the page centre outwards, the genus strip and the eight presence rings sit
concentrically outside the tips, and every tip label radiates outwards beyond the rings.
The 111 tips occupy a 360 - GAP_DEG arc; the gap at twelve o'clock carries the ring
names.  This
replaced the earlier two-column rectangular layout, whose letter B could only sit over
the left column while the right column carried the same two element types unlabelled.
Tip labels carry the MAG identifier and the GTDB species with the genus abbreviated to
its initial (the genus is encoded by the colour strip and named in full in the key).

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
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Wedge

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

# panel C grouping
df = pd.DataFrame({"genus": genus.reindex(PANEL), "score": score})
gcount = df["genus"].value_counts()
TOP_N = 9  # style constant: genera drawn individually; the remainder are pooled
top = list(gcount.head(TOP_N).index)
df["grp"] = np.where(df["genus"].isin(top), df["genus"], "Other genera")
order_c = df.groupby("grp")["score"].mean().sort_values(ascending=False).index.tolist()

# ------------------------------------------------------------------ page geometry
FS_TIP = 5.5                     # tip labels only; 111 tips on one page need a tighter
                                 # pitch than the 7 pt body size allows
GAP_DEG = 16.0                   # arc left empty at twelve o'clock for the ring names
R_IN, R_TREE = 3.0, 28.0         # root radius and tip radius of the radial tree (mm)
STRIP_R0, STRIP_W = 28.6, 2.0    # genus strip
RING_R0, RING_W = 31.2, 1.9      # first presence ring and ring pitch
RING_R1 = RING_R0 + RING_W * len(GENES)
R_LAB = RING_R1 + 1.2            # tip labels start here
TOP = 6.0                        # top of the circle's bounding square
CX = 76.0                        # circle centre, x; the keys take the right margin
KEY_X = 150.0                    # left edge of the keys

# the outer radius is the label start plus the widest tip label, measured below; the
# circle's bounding square is sized from that measurement so nothing is clipped
CD_H = 44.0
fig = None


def label_text(mag):
    sp = species[mag]
    if not sp:
        return mag
    g, _, rest = sp.partition(" ")
    return f"{mag}  {g[0]}. {rest}" if rest else f"{mag}  {sp}"


# measure the widest tip label on a throw-away figure at the tip size
_f = plt.figure(figsize=(1, 1))
_r = _f.canvas.get_renderer()
w_lab = max(_f.text(0, 0, label_text(m), fontsize=FS_TIP,
                    fontweight="bold" if m in heroes else "normal")
            .get_window_extent(renderer=_r).width for m in PANEL) / _f.dpi * 25.4
plt.close(_f)
R_OUT = R_LAB + w_lab + 1.0
assert CX - R_OUT >= 3.0 and CX + R_OUT <= KEY_X - 3.0, (CX, R_OUT, KEY_X)
CY = TOP + R_OUT                 # circle centre, y (mm from the top of the page)
CIRC_END = TOP + 2 * R_OUT
CD_LET_Y = CIRC_END + 6.0        # panel letters of the second row
CD_TOP = CD_LET_Y + 7.0          # top of the C / D axes
H = CD_TOP + CD_H + 14.0
assert H <= 235.0, H             # single-page ceiling

fig, ax_mm, text_mm, letter = st.page(H)


def fx(x_mm):
    return x_mm / st.PAGE_W_MM


def fy(y_mm):
    return 1.0 - y_mm / H


letter(4.0, 5.0, "A")
letter(4.0, CD_LET_Y, "B")
letter(102.0, CD_LET_Y, "C")

# ---- A and B: one square Axes in mm units, y up, centre at (0, 0)
axR = ax_mm(CX - R_OUT, TOP, 2 * R_OUT, 2 * R_OUT)
axR.set_xlim(-R_OUT, R_OUT)
axR.set_ylim(-R_OUT, R_OUT)
axR.set_aspect("equal")
axR.axis("off")

n_tip = len(ordered)
pitch = (360.0 - GAP_DEG) / n_tip                     # degrees per tip
theta0 = 90.0 - GAP_DEG / 2.0                          # arc starts right of the gap
# tip i (top to bottom in the rectangular order) runs clockwise from the gap
ang_of = {name: theta0 - (i + 0.5) * pitch for i, name in enumerate(ordered)}
# the last tip must end just left of the gap, so the two flanking half-gaps are equal
assert abs((theta0 - n_tip * pitch) - (90.0 + GAP_DEG / 2.0 - 360.0)) < 1e-9


def r_of(clade):
    return R_IN + depths[clade] / xmax * (R_TREE - R_IN)


def pol(r, a_deg):
    a = np.deg2rad(a_deg)
    return r * np.cos(a), r * np.sin(a)


def draw(clade, r_parent):
    """Radial phylogram: an arc at the parent radius spanning the children's angles,
    then a radial spoke to each child.  Returns the clade's angle."""
    r = r_of(clade)
    if clade.is_terminal():
        a = ang_of[clade.name]
    else:
        angs = [draw(c, r) for c in clade.clades]
        a = (min(angs) + max(angs)) / 2.0
        arc = np.linspace(min(angs), max(angs), max(2, int(abs(max(angs) - min(angs)) * 2)))
        xs, ys = pol(r, arc)
        axR.plot(xs, ys, color=TEXT, lw=0.4, solid_capstyle="butt")
    x0, y0 = pol(r_parent, a)
    x1, y1 = pol(r, a)
    axR.plot([x0, x1], [y0, y1], color=TEXT, lw=0.4, solid_capstyle="butt")
    return a


draw(tree.root, R_IN)

tip_labels = []
for name in ordered:
    mag = mag_of[name]
    a = ang_of[name]
    g = genus[mag]
    axR.add_patch(Wedge((0, 0), STRIP_R0 + STRIP_W, a - pitch / 2, a + pitch / 2,
                        width=STRIP_W, facecolor=GENUS_COL.get(g, OTHER_COL),
                        edgecolor="none"))
    # ---- B: gene presence as concentric rings, one ring per gene, inner to outer
    for j, gene in enumerate(GENES):
        r0 = RING_R0 + j * RING_W
        axR.add_patch(Wedge((0, 0), r0 + RING_W, a - pitch / 2, a + pitch / 2,
                            width=RING_W,
                            facecolor=GREEN if pres.loc[mag, gene] else "white",
                            edgecolor=LIGHT, lw=0.2))
    # tip label, radiating outwards; flipped on the left half so it reads left to right
    x, y = pol(R_LAB, a)
    right = np.cos(np.deg2rad(a)) >= 0
    tip_labels.append(
        axR.text(x, y, label_text(mag), rotation=a if right else a - 180,
                 rotation_mode="anchor", ha="left" if right else "right", va="center",
                 fontsize=FS_TIP, color=HERO if mag in heroes else TEXT,
                 fontweight="bold" if mag in heroes else "normal"))

# ring names in the gap at twelve o'clock, one per ring
for j, gene in enumerate(GENES):
    axR.text(0, RING_R0 + (j + 0.5) * RING_W, gene, ha="center", va="center",
             fontsize=FS_TIP, fontstyle="italic", color=TEXT)

# scale bar: a round substitutions-per-site distance, drawn at the tree's radial scale in
# the free lower-left corner of the circle's bounding square
step = 10 ** np.floor(np.log10(xmax / 4))
bar = float(max(kk * step for kk in (1, 2, 5) if kk * step <= xmax / 3))
bar_mm = bar / xmax * (R_TREE - R_IN)
sx, sy = -R_OUT + 2.0, -R_OUT + 4.0
axR.plot([sx, sx + bar_mm], [sy, sy], color=TEXT, lw=0.9)
axR.text(sx, sy - 1.0, f"{bar:g} substitutions/site", ha="left", va="top",
         fontsize=FS_STAT, color=TEXT)

# ---- keys for A, to the right of the circle: the genus strip, then the rings
handles = [Patch(facecolor=GENUS_COL[g], label=f"{g} ({int(gcount_all[g])})")
           for g in sorted(strip_genera, key=lambda g: -gcount_all[g])]
handles.append(Patch(facecolor=OTHER_COL, label=f"{OTHER_LABEL} ({n_other})"))
leg1 = fig.legend(handles=handles, loc="upper left", ncol=1, fontsize=FS_STAT,
                  frameon=False, bbox_to_anchor=(fx(KEY_X), fy(CY - 34.0)),
                  handlelength=1.1, handletextpad=0.45, labelspacing=0.4,
                  title="GTDB-Tk genus (MAGs)", title_fontsize=FS_STAT, alignment="left")
handles2 = [Patch(facecolor=GREEN, label="present"),
            Patch(facecolor="white", edgecolor=LIGHT, lw=0.5, label="absent"),
            Patch(facecolor="none", edgecolor="none", label="MICP-complete MAG")]
leg2 = fig.legend(handles=handles2, loc="upper left", ncol=1, fontsize=FS_STAT,
                  frameon=False, bbox_to_anchor=(fx(KEY_X), fy(CY + 10.0)),
                  handlelength=1.1, handletextpad=0.45, labelspacing=0.4,
                  title="ureA–G, cah rings", title_fontsize=FS_STAT, alignment="left")
for leg in (leg1, leg2):
    leg.get_title().set_fontweight("bold")
for txt in leg2.get_texts():
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
# the top rows of the box plot carry no data left of score 6, so the key sits there
axC.legend(handles=[ring], loc="upper left", bbox_to_anchor=(0.01, 0.99),
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
hD, lD = axD.get_legend_handles_labels()
fig.legend(hD, lD, loc="upper left", bbox_to_anchor=(fx(112.0), fy(CD_LET_Y - 0.5)),
           ncol=2, fontsize=FS_STAT, frameon=False, handlelength=1.0, borderpad=0.2,
           columnspacing=1.2, handletextpad=0.4)

# the tip labels are the one element st.audit cannot police against a slot (they radiate
# into empty page), so their rendered extent is checked against the bounding square
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_box = axR.get_window_extent(renderer=_r)
for t in tip_labels:
    bb = t.get_window_extent(renderer=_r)
    assert bb.x0 >= _box.x0 - 1 and bb.x1 <= _box.x1 + 1 and \
        bb.y0 >= _box.y0 - 1 and bb.y1 <= _box.y1 + 1, (t.get_text(), bb, _box)

print(f"  page height {H:.1f} mm | circle radius {R_OUT:.1f} mm | tip pitch "
      f"{2 * np.pi * R_LAB * (360 - GAP_DEG) / 360 / n_tip:.2f} mm at the label ring | "
      f"widest tip label {w_lab:.1f} mm")
st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig1")
