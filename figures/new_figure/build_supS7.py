"""Suppl Fig S7 - ureC gene tree and its congruence with the species tree.

Panels
  A  maximum-likelihood ureC gene tree, 46 ureC-encoding MAGs, drawn from
     Table_S7c_ureC_gene_tree.newick; node labels are the IQ-TREE ultrafast bootstrap
     values, shown where >= 70 (a display threshold, not a datum)
  B  topology test and Robinson-Foulds distance, as stat columns:
     Table_S7f_SH_AU_test.iqtree (the USER TREES block) and Table_S7a_RF_distance.txt

The two candidate trees of the topology test are, in the order IQ-TREE received them
(research/revision/02_ureC_gene_tree.py writes ureC_pruned.tre then species_pruned.tre
into candidate_trees.tre): tree 1 = the ureC gene tree, tree 2 = the bac120 species-tree
topology.  Both are pruned to the same 46 taxa.

The phylogram is drawn from the clade branch lengths directly rather than with
Phylo.draw, so that the tip label colours, the bootstrap placement and the millimetre
layout are under the same control as every other page in this set.  Tips are ladderized
for readability; ladderization reorders display only and does not change the topology.

Colour meanings on this page:
  blue    Sphingobacterium MICP-complete tip (S13, S16, S23, C22)
  orange  Pseudomonas_E MICP-complete tip (M1, S26)
  black   every other tip, and all branches
  grey    bootstrap values and the scale bar
  coral   a topology-test column that excludes the tree at the 95 % level
"""

import re
import sys
from pathlib import Path

import numpy as np
from Bio import Phylo
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
BOOTSTRAP_SHOWN_AT = 70          # display threshold for node labels
SCALE_BAR = 0.1                  # substitutions per site, drawn as a ruler

# ------------------------------------------------------------------ tree
tree = Phylo.read(SUPP / "Table_S7c_ureC_gene_tree.newick", "newick")
tree.ladderize()
tips = [t.name for t in tree.get_terminals()]
assert len(tips) == len(set(tips))
assert set(st.HEROES).issubset(set(tips))

depths = tree.depths()                      # root-to-node branch-length sums
y_of = {}
for i, t in enumerate(tree.get_terminals()):
    y_of[id(t)] = float(i)


def node_y(clade):
    if id(clade) in y_of:
        return y_of[id(clade)]
    ys = [node_y(c) for c in clade.clades]
    y = (min(ys) + max(ys)) / 2.0
    y_of[id(clade)] = y
    return y


node_y(tree.root)
max_depth = max(depths.values())

# ------------------------------------------------------------------ stats
rf_txt = (SUPP / "Table_S7a_RF_distance.txt").read_text()
rf = float(re.search(r"^RF = ([\d.]+)", rf_txt, re.M).group(1))
rf_max = float(re.search(r"^max_RF = ([\d.]+)", rf_txt, re.M).group(1))
rf_norm = float(re.search(r"^normalized_RF = ([\d.]+)", rf_txt, re.M).group(1))
n_taxa = int(re.search(r"n=(\d+) taxa", rf_txt).group(1))
assert n_taxa == len(tips)
assert abs(rf / rf_max - rf_norm) < 1e-9

TEST_ROW = re.compile(r"^\s+([12])\s+(-[\d.]+)\s+(\S+)\s+(\S+)\s+([+-])\s+(\S+)\s+([+-])"
                      r"\s+(\S+)\s+([+-])\s+\S+\s+[+-]\s+\S+\s+[+-]\s+\S+\s+[+-]"
                      r"\s+(\S+)\s+([+-])")
rows = []
for line in (SUPP / "Table_S7f_SH_AU_test.iqtree").read_text().splitlines():
    m = TEST_ROW.match(line)
    if m:
        rows.append(dict(tree=int(m.group(1)), logL=float(m.group(2)),
                         dL=float(m.group(3)), p_kh=float(m.group(6)),
                         p_sh=float(m.group(8)), p_au=float(m.group(10)),
                         excluded=(m.group(11) == "-")))
assert len(rows) == 2 and rows[0]["dL"] == 0.0
assert rows[0]["logL"] > rows[1]["logL"]
assert abs((rows[0]["logL"] - rows[1]["logL"]) - rows[1]["dL"]) < 0.01
TREE_NAMES = {1: "ureC gene tree", 2: "species tree (bac120)"}

# ------------------------------------------------------------------ page
ROW = 3.9
A_TOP, A_LEFT, A_W = 16.0, 6.0, 86.0
A_H = ROW * (len(tips) - 1)
B_LEFT = 108.0
H = A_TOP + A_H + 20.0

fig, ax_mm, text_mm, letter = st.page(H)
ax = ax_mm(A_LEFT, A_TOP, A_W, A_H)
letter(4, 4, "A")
letter(B_LEFT - 4, 4, "B")

# bootstrap labels are nudged apart deterministically when two nodes land on the same
# spot; the nudge moves the label, never the node
placed = []


def label_y(x, y, x_tol, y_tol=0.55, step=0.5):
    while any(abs(x - px) < x_tol and abs(y - py) < y_tol for px, py in placed):
        y -= step
    placed.append((x, y))
    return y


for clade in tree.find_clades(order="level"):
    x1 = depths[clade]
    if clade is not tree.root:
        parent = tree.get_path(clade)[-2] if len(tree.get_path(clade)) > 1 else tree.root
        x0 = depths[parent]
        ax.plot([x0, x1], [y_of[id(clade)]] * 2, color=st.TEXT, lw=0.6, zorder=2)
    if clade.clades:
        ys = [y_of[id(c)] for c in clade.clades]
        ax.plot([x1, x1], [min(ys), max(ys)], color=st.TEXT, lw=0.6, zorder=2)
        bs = clade.confidence
        if bs is not None and bs >= BOOTSTRAP_SHOWN_AT and clade is not tree.root:
            # centred on the branch that leads to the node: branch midpoints are further
            # apart than the node x positions, which keeps neighbouring labels clear
            bx_ = (x0 + x1) / 2
            by_ = label_y(bx_, y_of[id(clade)] - 0.22, max_depth * 0.035)
            ax.text(bx_, by_, f"{int(bs)}", ha="center", va="bottom", fontsize=4.6,
                    color=st.GREY, zorder=3)

for t in tree.get_terminals():
    hero = t.name in st.HEROES
    ax.text(depths[t] + max_depth * 0.012, y_of[id(t)], t.name, va="center",
            ha="left", fontsize=st.FS_BODY,
            color=st.hero_col(t.name) if hero else st.TEXT,
            fontweight="bold" if hero else "normal")

ax.set_xlim(-max_depth * 0.01, max_depth * 1.18)
ax.set_ylim(len(tips) - 0.4, -1.6)
ax.axis("off")

# scale bar, bound to the same x units as the branches
sb_y = -1.1
ax.plot([0, SCALE_BAR], [sb_y, sb_y], color=st.GREY, lw=1.0)
for x in (0, SCALE_BAR):
    ax.plot([x, x], [sb_y - 0.3, sb_y + 0.3], color=st.GREY, lw=1.0)
ax.text(SCALE_BAR * 1.15, sb_y, f"{SCALE_BAR} substitutions per site", ha="left",
        va="center", fontsize=st.FS_STAT, color=st.GREY)

def fmt_p(v):
    """IQ-TREE prints an exhausted RELL p-value as a bare 0; report it at the
    resolution the test actually has (1,000 replicates) rather than as an exact zero."""
    return "< 0.001" if v == 0 else f"{v:.3g}"


# --- B: topology test and RF, as stat columns
bx = [B_LEFT, B_LEFT + 30.0, B_LEFT + 42.0, B_LEFT + 52.0, B_LEFT + 62.0]
by = 16.0
text_mm(bx[0], by, "topology test", ha="left", va="bottom", fontsize=st.FS_BODY,
        fontweight="bold", color=st.AXIS)
by += 5.0
for x, head, al in zip(bx, ["candidate tree", "logL", "ΔlogL", "p-SH", "p-AU"],
                       ["left", "right", "right", "right", "right"]):
    text_mm(x if al == "left" else x + 8.0, by, head, ha=al, va="bottom",
            fontsize=st.FS_STAT, color=st.AXIS)
by += 1.5
for r in rows:
    by += 5.0
    col = st.HERO if r["excluded"] else st.TEXT
    text_mm(bx[0], by, TREE_NAMES[r["tree"]], ha="left", va="center",
            fontsize=st.FS_STAT)
    for x, v in zip(bx[1:], [f"{r['logL']:.2f}", f"{r['dL']:.2f}",
                             fmt_p(r['p_sh']), fmt_p(r['p_au'])]):
        text_mm(x + 8.0, by, v, ha="right", va="center", fontsize=st.FS_STAT, color=col)

by += 12.0
text_mm(bx[0], by, "Robinson–Foulds, ureC vs species tree", ha="left", va="bottom",
        fontsize=st.FS_BODY, fontweight="bold", color=st.AXIS)
for label, value in [("taxa compared", f"{n_taxa}"), ("RF", f"{rf:.0f}"),
                     ("maximum RF", f"{rf_max:.0f}"), ("normalised RF", f"{rf_norm:.3f}")]:
    by += 5.0
    text_mm(bx[0], by, label, ha="left", va="center", fontsize=st.FS_STAT)
    text_mm(bx[1] + 8.0, by, value, ha="right", va="center", fontsize=st.FS_STAT)

# tip-colour legend, placed in the right-hand column below the stat blocks
by += 12.0
fig.legend(handles=[Line2D([], [], color=st.SPHINGO, lw=4, label="Sphingobacterium"),
                    Line2D([], [], color=st.PSEUDO, lw=4, label="Pseudomonas_E"),
                    Line2D([], [], color=st.GREY, lw=1.2,
                           label=f"bootstrap ≥ {BOOTSTRAP_SHOWN_AT} shown")],
           loc="upper left", fontsize=st.FS_BODY, handlelength=1.2,
           bbox_to_anchor=(B_LEFT * st.MM / (st.PAGE_W_MM * st.MM),
                           1 - by * st.MM / (H * st.MM)))

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S7")
