#!/usr/bin/env python3
"""
Figure 1 redesign — circular (fan) phylogenetic tree with concentric
MICP-gene rings, genus colour arc, and labels placed radially.

Input : GTDB-Tk bac120 renamed tree + MICP_Pangenome_Final_Summary.csv
Output: /data/data/Upcycling/research/figures_final/Figure_1.(png|pdf)
"""
import os, re, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Patch, Circle
from matplotlib.colors import LinearSegmentedColormap
from Bio import Phylo

BASE="/data/data/Upcycling"
OUT =f"{BASE}/research/figures_final"
os.makedirs(OUT, exist_ok=True)

# ---------------- load data ----------------
tree_path=f"{BASE}/pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile"
tree=Phylo.read(tree_path,"newick")
try: tree.root_at_midpoint()
except: pass
sid=lambda n: n.split("_s__")[0]

gtdb=pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
gtdb["Genus"]=gtdb["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "Unclassified")
tax=gtdb.set_index("user_genome")["Genus"]

micp=pd.read_csv(f"{BASE}/pangenome_work/MICP_Pangenome_Final_Summary.csv", index_col=0)
gene_cols=["ureA","ureB","ureC","ureD","ureE","ureF","ureG","cah"]

MICP_COMPLETE = {"C22","M1","S13","S16","S23","S26"}

# --------------------------------------------------------------------------
# Radial layout
# --------------------------------------------------------------------------
leaves_full=[t.name for t in tree.get_terminals()]
leaf_ids   =[sid(x) for x in leaves_full]

# order leaves as tree is drawn (ladderize first)
tree.ladderize(reverse=False)
leaf_order=[t.name for t in tree.get_terminals()]
n = len(leaf_order)

# equally spaced angles (leave a small gap at theta=0 for orientation)
GAP_DEG = 20
total_deg = 360 - GAP_DEG
ang = {leaf: np.deg2rad(GAP_DEG/2 + total_deg * i/(n-1)) for i,leaf in enumerate(leaf_order)}

# compute depths (cumulative branch length from root)
depths = tree.depths()
max_depth = max(depths.values()) or 1.0

# tree radius mapping
R_TREE_MIN  = 0.30
R_TREE_MAX  = 0.72
def depth2r(d):  return R_TREE_MIN + (d/max_depth)*(R_TREE_MAX-R_TREE_MIN)

# compute angle for each internal node = mean of descendant leaf angles
def node_angle(clade):
    if clade.is_terminal(): return ang[clade.name]
    a=[node_angle(c) for c in clade.clades]
    return (min(a)+max(a))/2
for cl in tree.find_clades():
    cl._ang = node_angle(cl)
    cl._r   = depth2r(depths[cl])

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(15, 17))
# main tree takes top ~90 % so there is room for a legend strip below
ax = fig.add_axes([0.03, 0.16, 0.94, 0.80], projection="polar")
ax.set_ylim(0, 1.02)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_axis_off()

# ---- edges ----
for cl in tree.get_nonterminals():
    r_p = cl._r
    a_p = cl._ang
    # collect child angles to draw the horizontal arc
    child_angles=[c._ang for c in cl.clades]
    a_lo, a_hi = min(child_angles), max(child_angles)
    arc_steps = max(8, int(np.rad2deg(a_hi-a_lo)*2))
    arc_theta = np.linspace(a_lo, a_hi, arc_steps)
    ax.plot(arc_theta, [r_p]*len(arc_theta), color="black", lw=0.7, solid_capstyle="round")
    # radial lines to each child
    for c in cl.clades:
        ax.plot([c._ang, c._ang], [r_p, c._r], color="black", lw=0.7)

# ---- genus colour arc (between R_TREE_MAX and R_GENUS) ----
R_GENUS_IN  = R_TREE_MAX + 0.005
R_GENUS_OUT = R_TREE_MAX + 0.025

genera = sorted(tax.unique())
cmap   = plt.get_cmap("tab20", len(genera))
gcol   = {g: cmap(i) for i,g in enumerate(genera)}

# angular width per leaf (half step each side)
d_ang = total_deg/(n-1) * np.pi/180  # radians between leaves
half  = d_ang/2
for leaf in leaf_order:
    g = tax.get(sid(leaf),"Unclassified")
    theta0 = ang[leaf]-half; theta1=ang[leaf]+half
    # polar bar
    ax.bar(x=(theta0+theta1)/2, height=R_GENUS_OUT-R_GENUS_IN, width=theta1-theta0,
           bottom=R_GENUS_IN, color=gcol.get(g,"#cccccc"), edgecolor="none", align="center")

# ---- MICP heatmap rings ----
R_RINGS_IN = R_GENUS_OUT + 0.004
R_RING_W   = 0.020
heat_cmap = LinearSegmentedColormap.from_list("greens", ["#f0f8f0","#1e6b1e"])
for gi,gene in enumerate(gene_cols):
    r0 = R_RINGS_IN + gi*R_RING_W
    r1 = r0 + R_RING_W*0.9
    for leaf in leaf_order:
        sid_=sid(leaf)
        v = micp.loc[sid_, gene] if sid_ in micp.index else 0
        color = heat_cmap(0.95 if v else 0.05)
        theta0=ang[leaf]-half; theta1=ang[leaf]+half
        ax.bar(x=(theta0+theta1)/2, height=r1-r0, width=theta1-theta0,
               bottom=r0, color=color, edgecolor="none", align="center")
    # gene label at gap position (top of figure, inside the GAP_DEG wedge)
    lbl_theta = np.deg2rad(GAP_DEG/2 - 2)   # just inside the gap, left of first leaf
    ax.text(lbl_theta, (r0+r1)/2, gene, ha="right", va="center",
            fontsize=10, fontweight="bold", style="italic", color="#1e6b1e")

# ---- leaf labels (radial, aligned outward) ----
R_LABEL = R_RINGS_IN + len(gene_cols)*R_RING_W + 0.01
for leaf in leaf_order:
    sid_ = sid(leaf)
    a = ang[leaf]
    # rotation: text aligned with radius
    deg = np.rad2deg(a)
    rot = 90 - deg if deg <= 180 else 270 - deg
    ha  = "left" if deg <= 180 else "right"
    is_complete = sid_ in MICP_COMPLETE
    color = "#c0392b" if is_complete else "#2c3e50"
    fw    = "bold"    if is_complete else "normal"
    fs    = 10       if is_complete else 6.2
    ax.text(a, R_LABEL, sid_, rotation=rot, rotation_mode="anchor",
            ha=ha, va="center", fontsize=fs, color=color, fontweight=fw)

# ---- add big star + connecting line to MICP-complete leaves ----
for leaf in leaf_order:
    if sid(leaf) in MICP_COMPLETE:
        a=ang[leaf]
        ax.plot(a, R_TREE_MAX+0.002, marker="*", markersize=14,
                color="#c0392b", markeredgecolor="black", markeredgewidth=0.6, zorder=5)

# ---- title + legend ----
# genus legend — only genera with >=2 members shown
counts = tax.value_counts()
legend_genera = [g for g in genera if counts.get(g,0) >= 2]
handles = [Patch(facecolor=gcol[g], edgecolor="none",
                 label=f"{g} (n={counts[g]})") for g in legend_genera]
handles.append(Patch(facecolor="#1e6b1e", edgecolor="black",
                     label="MICP gene present"))
handles.append(Patch(facecolor="white", edgecolor="#c0392b", lw=1.8,
                     label="★  MICP-complete lineage member (n=6)"))
leg = fig.legend(handles=handles, loc="lower center",
                 bbox_to_anchor=(0.5, 0.01), fontsize=10, frameon=False,
                 ncol=5, title="Genus (GTDB-Tk)  ·  MICP gene ring legend",
                 title_fontsize=11)
leg._legend_box.align = "center"

fig.suptitle("Figure 1. Phylogenomic distribution of MICP-related genes across 111 livestock-waste MAGs.",
             fontsize=13, y=0.97, fontweight="bold")

# small scale bar in the centre
ax.plot([0, np.deg2rad(1)], [0.02, 0.02 + 0.10], color="black", lw=2)
ax.text(np.deg2rad(0.5), 0.085, f"{max_depth*0.10/(R_TREE_MAX-R_TREE_MIN):.2f}\nsubst./site",
        ha="center", va="center", fontsize=7, color="black")

fig.savefig(f"{OUT}/Figure_1.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Figure_1.pdf", bbox_inches="tight")
plt.close(fig)
print(f"Figure 1 (circular) written to {OUT}")
