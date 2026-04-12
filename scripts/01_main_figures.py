#!/usr/bin/env python3
"""
Figures for the livestock-waste MICP study (111 MAGs).
Targets: Bioresource Technology / Microbial Biotechnology.

Outputs (300 dpi PNG + PDF):
  Fig1_Phylogeny_MICP_heatmap.(png|pdf)
  Fig2_MICP_completeness_by_genus.(png|pdf)
  Fig3_ureCah_cluster_synteny.(png|pdf)
"""
import os, re, gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch, Rectangle
from matplotlib.lines import Line2D
from Bio import Phylo
from io import StringIO

BASE = "/data/data/Upcycling"
OUT  = f"{BASE}/research"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
micp = pd.read_csv(f"{BASE}/pangenome_work/MICP_Pangenome_Final_Summary.csv",
                   index_col=0)
gene_cols = ["ureA","ureB","ureC","ureD","ureE","ureF","ureG","cah"]

gtdb = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv",
                   sep="\t")
def parse_genus(s):
    m = re.search(r"g__([^;]*)", str(s)); return m.group(1) if m and m.group(1) else "Unclassified"
def parse_species(s):
    m = re.search(r"s__([^;]*)", str(s)); return m.group(1) if m and m.group(1) else ""
gtdb["Genus"]   = gtdb["classification"].apply(parse_genus)
gtdb["Species"] = gtdb["classification"].apply(parse_species)
tax = gtdb.set_index("user_genome")[["Genus","Species"]]

tree_path = f"{BASE}/pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile"
tree = Phylo.read(tree_path, "newick")
# map stripped leaf name (everything before _s__) -> full leaf
leaf_full = [t.name for t in tree.get_terminals()]
def sid(full): return full.split("_s__")[0]
full_by_id = {sid(f): f for f in leaf_full}

# midpoint root for nicer layout
try: tree.root_at_midpoint()
except Exception: pass

# order leaves as drawn (top->bottom)
ordered_full = []
def collect(clade):
    if clade.is_terminal(): ordered_full.append(clade.name)
    else:
        for c in clade.clades: collect(c)
collect(tree.root)
ordered_ids = [sid(x) for x in ordered_full]

# genus palette
genera = sorted(tax["Genus"].unique())
cmap = plt.get_cmap("tab20", len(genera))
genus_color = {g: cmap(i) for i, g in enumerate(genera)}

# highlight "hero" Sphingobacterium clade
HERO = {"S13","S16","S23","S26","C22","M1"}

# ---------------------------------------------------------------------------
# Figure 1 — Phylogeny + MICP heatmap
# ---------------------------------------------------------------------------
def draw_tree(ax, tree, ordered_full):
    # compute x coords via cumulative branch length; y by leaf order
    y_of = {name: i for i, name in enumerate(ordered_full[::-1])}  # bottom->top
    depths = tree.depths()
    xmax = max(depths.values()) or 1.0
    def _draw(clade, x0):
        if clade.is_terminal():
            y = y_of[clade.name]
            ax.plot([x0, depths[clade]], [y, y], color="black", lw=0.6)
            return y
        ys = [_draw(c, depths[clade]) for c in clade.clades]
        ymid = (min(ys)+max(ys))/2
        ax.plot([depths[clade], depths[clade]], [min(ys), max(ys)], color="black", lw=0.6)
        if clade is not tree.root:
            ax.plot([x0, depths[clade]], [ymid, ymid], color="black", lw=0.6)
        return ymid
    _draw(tree.root, 0.0)
    ax.set_xlim(0, xmax*1.02)
    ax.set_ylim(-0.5, len(ordered_full)-0.5)
    ax.invert_yaxis()
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
    ax.set_yticks([]); ax.set_xlabel("Substitutions / site", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    return y_of

n = len(ordered_full)
fig = plt.figure(figsize=(12, max(10, n*0.14)))
gs  = fig.add_gridspec(1, 4, width_ratios=[2.2, 0.12, 0.12, 1.4], wspace=0.02)
ax_t = fig.add_subplot(gs[0,0])
ax_g = fig.add_subplot(gs[0,1], sharey=ax_t)
ax_l = fig.add_subplot(gs[0,2], sharey=ax_t)
ax_h = fig.add_subplot(gs[0,3], sharey=ax_t)

y_of = draw_tree(ax_t, tree, ordered_full)

# Genus strip
for full in ordered_full:
    sid_ = sid(full); y = y_of[full]
    g = tax.loc[sid_, "Genus"] if sid_ in tax.index else "Unclassified"
    ax_g.add_patch(Rectangle((0, y-0.5), 1, 1, color=genus_color.get(g,"#cccccc")))
ax_g.set_xlim(0,1); ax_g.set_xticks([]); ax_g.set_yticks([])
for s in ax_g.spines.values(): s.set_visible(False)
ax_g.set_title("Genus", fontsize=8, rotation=90, loc="left", pad=4)

# Leaf labels
ax_l.set_xlim(0,1); ax_l.set_xticks([]); ax_l.set_yticks([])
for s in ax_l.spines.values(): s.set_visible(False)
for full in ordered_full:
    sid_ = sid(full); y = y_of[full]
    sp = tax.loc[sid_,"Species"] if sid_ in tax.index else ""
    color = "#c0392b" if sid_ in HERO else "black"
    weight = "bold" if sid_ in HERO else "normal"
    label = f"{sid_}  {sp}" if sp else sid_
    ax_l.text(0, y, label, va="center", ha="left", fontsize=6.2,
              color=color, fontweight=weight)

# MICP heatmap
mat = np.zeros((n, len(gene_cols)))
for i, full in enumerate(ordered_full):
    sid_ = sid(full)
    if sid_ in micp.index:
        mat[y_of[full], :] = micp.loc[sid_, gene_cols].values
ax_h.imshow(mat, aspect="auto", cmap="Greens", vmin=0, vmax=1,
            extent=[-0.5, len(gene_cols)-0.5, n-0.5, -0.5])
ax_h.set_xticks(range(len(gene_cols)))
ax_h.set_xticklabels(gene_cols, rotation=45, ha="right", fontsize=8, style="italic")
ax_h.set_yticks([])
for s in ax_h.spines.values(): s.set_visible(False)
ax_h.set_title("MICP gene presence", fontsize=9)

# Legend (top genera)
counts = tax["Genus"].value_counts()
top = counts.head(10).index.tolist()
handles = [Patch(facecolor=genus_color[g], label=f"{g} (n={counts[g]})") for g in top]
fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7, frameon=False,
           bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Fig. 1 | Phylogeny of 111 MAGs and distribution of MICP-related genes",
             fontsize=11, y=0.995)
fig.tight_layout(rect=[0,0.02,1,0.98])
fig.savefig(f"{OUT}/Fig1_Phylogeny_MICP_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig1_Phylogeny_MICP_heatmap.pdf", bbox_inches="tight")
plt.close(fig)
print("Fig1 done")

# ---------------------------------------------------------------------------
# Figure 2 — MICP completeness by genus
# ---------------------------------------------------------------------------
df = micp.copy()
df["Genus"] = [tax.loc[i,"Genus"] if i in tax.index else "Unclassified" for i in df.index]
df["Score"] = df[gene_cols].sum(axis=1)
order = df.groupby("Genus")["Score"].mean().sort_values(ascending=False).index.tolist()

fig, axes = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios":[1.2,1]})
# 2a: boxplot of total score by genus
ax = axes[0]
data = [df.loc[df.Genus==g,"Score"].values for g in order]
bp = ax.boxplot(data, positions=range(len(order)), widths=0.6, patch_artist=True,
                medianprops=dict(color="black"))
for patch, g in zip(bp["boxes"], order):
    patch.set_facecolor(genus_color.get(g,"#bbbbbb")); patch.set_alpha(0.85)
# overlay individual points
for i, g in enumerate(order):
    y = df.loc[df.Genus==g,"Score"].values
    x = np.random.normal(i, 0.07, size=len(y))
    hero_mask = np.array([idx in HERO for idx in df.loc[df.Genus==g].index])
    ax.scatter(x[~hero_mask], y[~hero_mask], s=14, color="black", alpha=0.5, zorder=3)
    if hero_mask.any():
        ax.scatter(x[hero_mask], y[hero_mask], s=55, color="#c0392b",
                   edgecolor="black", lw=0.6, zorder=4, label="_nolegend_")
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, rotation=40, ha="right", fontsize=9)
ax.set_ylabel("MICP module score (ureA–G + cah, max 8)", fontsize=10)
ax.set_title("a  MICP module completeness per genus", fontsize=10, loc="left")
ax.axhline(8, ls="--", lw=0.6, color="gray"); ax.set_ylim(-0.3, 8.8)
ax.grid(axis="y", ls=":", alpha=0.5)

# 2b: per-gene prevalence in hero clade vs rest
ax = axes[1]
hero_df = df.loc[df.index.isin(HERO), gene_cols].mean()*100
rest_df = df.loc[~df.index.isin(HERO), gene_cols].mean()*100
x = np.arange(len(gene_cols)); w = 0.38
ax.bar(x-w/2, hero_df.values, w, label=f"Sphingobacterium hero clade (n={len(HERO)})",
       color="#c0392b", edgecolor="black", lw=0.5)
ax.bar(x+w/2, rest_df.values, w, label=f"Other MAGs (n={len(df)-len(HERO)})",
       color="#7f8c8d", edgecolor="black", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(gene_cols, style="italic")
ax.set_ylabel("Prevalence (% of genomes)")
ax.set_ylim(0,105); ax.set_title("b  Per-gene prevalence: hero clade vs. others",
                                  fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.grid(axis="y", ls=":", alpha=0.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)

fig.suptitle("Fig. 2 | MICP functional completeness across the 111-MAG panel",
             fontsize=11)
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(f"{OUT}/Fig2_MICP_completeness_by_genus.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig2_MICP_completeness_by_genus.pdf", bbox_inches="tight")
plt.close(fig)
print("Fig2 done")

# ---------------------------------------------------------------------------
# Figure 3 — ure-cah conserved cluster synteny (Sphingobacterium hero samples)
# ---------------------------------------------------------------------------
HERO_LIST = ["C22","M1","S13","S16","S23","S26"]
GENE_COLORS = {
    "ureA":"#1f77b4","ureB":"#2ca02c","ureC":"#d62728",
    "ureD":"#9467bd","ureE":"#8c564b","ureF":"#e377c2",
    "ureG":"#17becf","cah":"#ff7f0e","other":"#bbbbbb"
}
def classify(prod, gene):
    p = (prod or "").lower(); g = (gene or "").lower()
    if "urease subunit alpha" in p: return "ureC"
    if "urease subunit beta"  in p: return "ureB"
    if "urease subunit gamma" in p: return "ureA"
    if g in ["ured","uree","uref","ureg"]: return g[0:3]+g[3].upper()  # ureD etc
    if "ured" in p: return "ureD"
    if "uree" in p: return "ureE"
    if "uref" in p: return "ureF"
    if "ureg" in p: return "ureG"
    if "carbonic anhyd" in p: return "cah"
    return "other"

def parse_gff(path):
    rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS": continue
            attrs = dict(kv.split("=",1) for kv in p[8].split(";") if "=" in kv)
            rows.append({"contig":p[0],"start":int(p[3]),"end":int(p[4]),
                         "strand":p[6],"gene":attrs.get("gene",""),
                         "product":attrs.get("product","")})
    return pd.DataFrame(rows)

def find_cluster(df):
    # find the contig with the densest ure cluster (>=3 ure genes within 20 kb)
    df = df.copy()
    df["cls"] = [classify(p,g) for p,g in zip(df["product"], df["gene"])]
    ure = df[df.cls.str.startswith("ure")]
    best = None
    for ctg, sub in ure.groupby("contig"):
        sub = sub.sort_values("start")
        if len(sub) < 3: continue
        span = sub["end"].max() - sub["start"].min()
        if span > 30000: continue
        score = len(set(sub.cls))
        if best is None or score > best[0]:
            best = (score, ctg, sub["start"].min(), sub["end"].max())
    if not best: return None
    _, ctg, s, e = best
    win_s, win_e = max(0, s-2000), e+2000
    region = df[(df.contig==ctg)&(df.start>=win_s)&(df.end<=win_e)].sort_values("start")
    return ctg, win_s, win_e, region

fig, axes = plt.subplots(len(HERO_LIST), 1, figsize=(12, 1.2*len(HERO_LIST)+0.8),
                         sharex=False)
if len(HERO_LIST)==1: axes=[axes]

for ax, samp in zip(axes, HERO_LIST):
    gff = f"{BASE}/MAGs_FASTA_files/bakta_results/{samp}/{samp}.gff3"
    if not os.path.exists(gff):
        ax.text(0.5,0.5,f"{samp}: GFF not found", transform=ax.transAxes, ha="center")
        ax.axis("off"); continue
    df = parse_gff(gff)
    res = find_cluster(df)
    if res is None:
        ax.text(0.5,0.5,f"{samp}: cluster not found", transform=ax.transAxes, ha="center")
        ax.axis("off"); continue
    ctg, ws, we, region = res
    # normalize coordinates to kb from cluster start
    off = region["start"].min()
    L   = region["end"].max() - off
    for _, r in region.iterrows():
        cls = classify(r["product"], r["gene"])
        col = GENE_COLORS.get(cls, GENE_COLORS["other"])
        x0 = (r["start"]-off)/1000; x1 = (r["end"]-off)/1000
        if r["strand"] == "-":
            x_start, dx = x1, -(x1-x0)
        else:
            x_start, dx = x0, (x1-x0)
        head = min(0.35, abs(dx)*0.4)
        ax.add_patch(FancyArrow(x_start, 0, dx, 0,
                                width=0.55, head_width=0.8, head_length=head,
                                length_includes_head=True,
                                facecolor=col, edgecolor="black", lw=0.4))
        if cls != "other":
            ax.text((x0+x1)/2, 0.9, cls, ha="center", va="bottom",
                    fontsize=7, style="italic")
    ax.set_xlim(-0.5, L/1000+0.5)
    ax.set_ylim(-1.2, 1.6)
    ax.set_yticks([])
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("kb" if ax is axes[-1] else "", fontsize=8)
    ax.text(-0.02, 0.5, samp, transform=ax.transAxes, ha="right", va="center",
            fontsize=9, fontweight="bold", color="#c0392b")

# legend
handles = [Patch(facecolor=GENE_COLORS[g], edgecolor="black", label=g)
           for g in ["ureA","ureB","ureC","ureD","ureE","ureF","ureG","cah","other"]]
fig.legend(handles=handles, loc="lower center", ncol=9, fontsize=8, frameon=False,
           bbox_to_anchor=(0.5,-0.02))
fig.suptitle("Fig. 3 | Conserved ure–cah gene cluster in Sphingobacterium MAGs",
             fontsize=11, y=0.995)
fig.tight_layout(rect=[0.02, 0.03, 1, 0.97])
fig.savefig(f"{OUT}/Fig3_ureCah_cluster_synteny.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig3_ureCah_cluster_synteny.pdf", bbox_inches="tight")
plt.close(fig)
print("Fig3 done")
print("All figures written to", OUT)
