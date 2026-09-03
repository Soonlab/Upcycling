"""Main Fig 1 - bac120 phylogeny of the 111 MAGs beside the MICP gene presence matrix.

Panels (reading order left to right):
  A  maximum-likelihood bac120 tree (midpoint-rooted for display), a genus colour strip,
     and the MAG identifier with its GTDB species where one was assigned
  B  presence of ureA-G and cah for the same rows

Sources
  pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile   IQ-TREE topology
                                    and branch lengths; tip names carry the MAG id and
                                    the GTDB species
  Table_S1d_GTDB_Tk_classification.tsv    genus and species per MAG
  MAGs_FASTA_files/bakta_results/*/*.tsv  gene presence, via _micp_presence.presence():
                                    a CDS-only keyword scan of all 111 annotations
  Table_S15a_alkaliphile_signature_per_MAG.csv   the MICP-complete / rest group flag

Provenance note (see _job/JOURNAL.md).  Two earlier sources for panel B are unusable.
pangenome_work/MICP_Pangenome_Final_Summary.csv, used by the shipped figure, holds only
100 of the 111 MAGs - C1 and C10-C19 are missing and were drawn as all-absent.
Table_S1a_ace_samples_list.csv, used by the first version of this rebuild, lists only 45
MAGs (13 of the 66 it omits carry a ure CDS) and counts the Bakta 5_ureB_sRNA non-coding
feature as a copy of the beta subunit, which is why it scores S26 8/8 where the Bakta
annotation has no protein-coding ureB.  Panel B is therefore built from the annotations
directly; _micp_presence documents both defects and asserts the relationship to S1a.

Colour meanings on this page (one meaning per colour):
  genus strip = one colour per GTDB genus; Sphingobacterium blue and Pseudomonas_E orange
                are the two MICP-complete lineages, every other genus takes a muted hue
                and no genus is given green or coral
  green       = the gene is present in that MAG (panel B); white = absent
  coral       = a MICP-complete MAG (tip label only)
  black       = tree branches
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
from _style import HERO, SPHINGO, PSEUDO, GREEN, TEXT, GREY, AXIS, LIGHT, FS_BODY, FS_STAT
from _micp_presence import presence

st.setup()
OUT = HERE / "figures"
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

grp = pd.read_csv(SUPP / "Table_S15a_alkaliphile_signature_per_MAG.csv", index_col=0)
heroes = set(grp.index[grp["group"] == "MICP_complete"])
assert sorted(heroes) == sorted(st.HEROES), sorted(heroes)

tree = Phylo.read(TREE, "newick")
tree.root_at_midpoint()
tips = [t.name for t in tree.get_terminals()]
mag_of = {t: t.split("_s__")[0] for t in tips}
assert len(tips) == 111 == len(tax), (len(tips), len(tax))
assert set(mag_of.values()) == set(tax.index)

pres = presence().reindex(sorted(tax.index))[GENES]
assert pres.notna().all().all()
# five of the six MICP-complete MAGs carry the full module; S26 has no protein-coding
# urease beta subunit in its Bakta annotation (see _micp_presence)
hero_score = pres.loc[sorted(heroes)].sum(axis=1)
assert (hero_score == len(GENES)).sum() == 5, hero_score.to_dict()
assert pres.loc["S26", "ureB"] == 0 and hero_score["S26"] == len(GENES) - 1

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

gcount = genus.value_counts()
strip_genera = [g for g in GENUS_COL if g in set(genus)]
assert set(strip_genera) == set(gcount.head(len(GENUS_COL)).index), \
    (sorted(strip_genera), sorted(gcount.head(len(GENUS_COL)).index))
n_other = int((~genus.isin(strip_genera)).sum())

# ------------------------------------------------------------------ page
ROW = 2.25                       # mm per tip
TOP = 13.0
PLOT_H = len(ordered) * ROW
TREE_L, TREE_W = 9.0, 50.0
STRIP_L, STRIP_W = TREE_L + TREE_W + 1.5, 2.6
LAB_L = STRIP_L + STRIP_W + 1.5
HEAT_W = 30.0
HEAT_L = 180.0 - 8.0 - HEAT_W
H = TOP + PLOT_H + 26.0

fig, ax_mm, text_mm, letter = st.page(H)
letter(3.0, 6.0, "A")
letter(HEAT_L - 2.0, 6.0, "B")

ax_t = ax_mm(TREE_L, TOP, TREE_W, PLOT_H)
ax_s = ax_mm(STRIP_L, TOP, STRIP_W, PLOT_H)
ax_l = ax_mm(LAB_L, TOP, HEAT_L - LAB_L - 1.0, PLOT_H)
ax_h = ax_mm(HEAT_L, TOP, HEAT_W, PLOT_H)

# ---- A: tree
def draw(clade, x_parent):
    x = depths[clade]
    if clade.is_terminal():
        y = y_of[clade.name]
        ax_t.plot([x_parent, x], [y, y], color=TEXT, lw=0.45, solid_capstyle="butt")
        return y
    ys = [draw(c, x) for c in clade.clades]
    ax_t.plot([x, x], [min(ys), max(ys)], color=TEXT, lw=0.45, solid_capstyle="butt")
    ymid = (min(ys) + max(ys)) / 2
    ax_t.plot([x_parent, x], [ymid, ymid], color=TEXT, lw=0.45, solid_capstyle="butt")
    return ymid


draw(tree.root, 0.0)
ax_t.set_xlim(-xmax * 0.02, xmax * 1.02)
ax_t.set_ylim(len(ordered) - 0.5, -0.5)
ax_t.set_yticks([])
ax_t.set_xticks([])
for s in ax_t.spines.values():
    s.set_visible(False)

# scale bar: a round substitutions-per-site distance, drawn under the tree
step = 10 ** np.floor(np.log10(xmax / 4))
bar = float(max(k * step for k in (1, 2, 5) if k * step <= xmax / 3))
y_bar = len(ordered) + 1.5
ax_t.plot([0, bar], [y_bar, y_bar], color=TEXT, lw=0.9, clip_on=False)
ax_t.text(bar / 2, y_bar + 1.6, f"{bar:g} substitutions/site", ha="center", va="top",
          fontsize=FS_STAT, color=TEXT, clip_on=False)

# ---- A: genus strip and tip labels
for name in ordered:
    mag = mag_of[name]
    y = y_of[name]
    g = genus[mag]
    ax_s.add_patch(Rectangle((0, y - 0.5), 1, 1,
                             facecolor=GENUS_COL.get(g, OTHER_COL), edgecolor="none"))
    sp = species[mag]
    label = f"{mag}  {sp}" if sp else mag
    ax_l.text(0, y, label, va="center", ha="left", fontsize=FS_STAT,
              color=HERO if mag in heroes else TEXT,
              fontweight="bold" if mag in heroes else "normal")

for ax in (ax_s, ax_l):
    ax.set_xlim(0, 1)
    ax.set_ylim(len(ordered) - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

# ---- B: gene presence
mat = np.array([pres.loc[mag_of[n], GENES].values for n in ordered], dtype=float)
ax_h.imshow(mat, aspect="auto", cmap=st.seq_cmap(), vmin=0, vmax=1,
            extent=[0, len(GENES), len(ordered) - 0.5, -0.5], interpolation="nearest")
for k in range(len(GENES) + 1):
    ax_h.plot([k, k], [-0.5, len(ordered) - 0.5], color=LIGHT, lw=0.3)
ax_h.set_xticks(np.arange(len(GENES)) + 0.5)
ax_h.set_xticklabels(GENES, rotation=90, fontstyle="italic",
                     fontsize=FS_STAT)
ax_h.xaxis.set_ticks_position("top")
ax_h.tick_params(axis="x", length=0, pad=1.5)
ax_h.set_yticks([])
ax_h.set_ylim(len(ordered) - 0.5, -0.5)
for s in ax_h.spines.values():
    s.set_visible(False)

# ---- legends, below the plot
handles = [Patch(facecolor=GENUS_COL[g], label=f"{g} ({int(gcount[g])})")
           for g in sorted(strip_genera, key=lambda g: -gcount[g])]
handles.append(Patch(facecolor=OTHER_COL, label=f"{OTHER_LABEL} ({n_other})"))
handles.append(Patch(facecolor=GREEN, label="gene present"))
handles.append(Patch(facecolor="white", edgecolor=LIGHT, lw=0.5, label="gene absent"))
handles.append(Patch(facecolor="none", edgecolor="none", label="MICP-complete MAG"))
leg = fig.legend(handles=handles, loc="lower left", ncol=4, fontsize=FS_STAT,
                 frameon=False, bbox_to_anchor=(0.045, 0.004), handlelength=1.1,
                 columnspacing=1.2, handletextpad=0.45, labelspacing=0.5)
for txt, h in zip(leg.get_texts(), handles):
    if txt.get_text() == "MICP-complete MAG":
        txt.set_color(HERO)
        txt.set_fontweight("bold")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig1")
