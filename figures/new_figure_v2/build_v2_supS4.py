"""Suppl Fig S4 (consolidated) - genome relatedness, gene-tree topology and reference context.

Consolidation of four old supplementary pages onto one page, per
consolidation_260904/DESIGN.md.

Panels (old page -> new panel)
  A  old S6   pairwise skani ANI among the six MICP-complete MAGs, values printed, an em
              dash where skani reported no alignment
  B  old S7A  maximum-likelihood ureC gene tree, 46 ureC-encoding MAGs, ultrafast bootstrap
              shown at >= 70, MICP-complete tips coloured by lineage
  C  old S18A ANI of the six against the organism-verified MICP reference panel
  D  old S18B the reported MAG-reference pairs, ranked, against the 95 % species boundary
  E  old S10  MICP feature prevalence within the 146 Pseudomonas_E reference genomes

The old S18 carried two panels.  They are kept as two lettered panels (C and D) rather than
merged: the reference matrix needs a column of long organism names, and at the type size of
this figure the ranked-pair strip cannot sit beside it inside 180 mm.  The last panel is
therefore lettered E, the alternative the task allowed.

The old S7 panel B (SH/AU topology test and Robinson-Foulds distance against the species
tree) is not drawn on this page; the task specifies the tree alone for panel B.  Its files
are still read and every assertion of the old script is still evaluated, and the values are
printed to stdout, so the numbers remain checked and are available for the legend.

A stored ANI of 0.0 means skani reported no alignment for that pair, not 0 % identity, and
is drawn as an em dash (A) or a dash (C) rather than as a value.  That rendering is a
deliberate honesty choice and is kept from the old pages.

Colour meanings on this page
  white -> green   ANI (%) between two MICP-complete MAGs (panel A only)
  blue             Sphingobacterium MAG (S13, S16, S23, C22): its label, cells and points
  orange           Pseudomonas_E MAG (M1, S26): its label, cells and points
  coral            a MICP feature carried by the Pseudomonas_E reference genomes (E)
  grey             no-alignment dashes, bootstrap values, the scale bar and the 95 %
                   species boundary - never a lineage

Sources
  Table_S4a_skani_ANI_matrix_111MAGs.csv                        A
  Table_S7c_ureC_gene_tree.newick                               B
  Table_S7a_RF_distance.txt, Table_S7f_SH_AU_test.iqtree        B (checked, printed)
  research/additional/C5_panMICP_env_v2/reference_manifest.tsv  C, D (organism-verified panel)
  research/additional/C5_panMICP_env_v2/skani_panMICP.matrix    C, D
  Table_S20_skani_hero_vs_refs.tsv                              C, D (asserted against)
  Table_S13a / Table_S13b pseudomonas_e screens                 E
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import HERO, GREY, TEXT, LIGHT, AXIS, FS_BODY, FS_STAT

st.setup()
OUT = HERE / "figures_v2"
SUPP = Path(gh.SUPP)
C5 = Path(gh.ADDITIONAL) / "C5_panMICP_env_v2"

NO_ALIGNMENT = 0.0
SPECIES_ANI = 95.0        # prokaryotic species boundary (Konstantinidis & Tiedje 2005)
BOOTSTRAP_SHOWN_AT = 70   # display threshold for node labels
SCALE_BAR = 0.1           # substitutions per site, drawn as a ruler

# ================================================================== A: ANI among the six
ani6 = pd.read_csv(SUPP / "Table_S4a_skani_ANI_matrix_111MAGs.csv", index_col=0)
mat = ani6.loc[st.HEROES, st.HEROES]
assert np.allclose(mat.values, mat.values.T)
assert np.allclose(np.diag(mat.values), 100.0)

vals6 = mat.values.astype(float)
shown6 = np.where(vals6 > NO_ALIGNMENT, vals6, np.nan)
off6 = shown6[~np.eye(len(st.HEROES), dtype=bool)]
vmin6 = float(np.nanmin(off6))

# ================================================================== B: ureC gene tree
tree = Phylo.read(SUPP / "Table_S7c_ureC_gene_tree.newick", "newick")
tree.ladderize()
tips = [t.name for t in tree.get_terminals()]
assert len(tips) == len(set(tips))
assert set(st.HEROES).issubset(set(tips))

depths = tree.depths()
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

# congruence statistics: checked and reported, not drawn on this page
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
rows_top = []
for line in (SUPP / "Table_S7f_SH_AU_test.iqtree").read_text().splitlines():
    m = TEST_ROW.match(line)
    if m:
        rows_top.append(dict(tree=int(m.group(1)), logL=float(m.group(2)),
                             dL=float(m.group(3)), p_kh=float(m.group(6)),
                             p_sh=float(m.group(8)), p_au=float(m.group(10)),
                             excluded=(m.group(11) == "-")))
assert len(rows_top) == 2 and rows_top[0]["dL"] == 0.0
assert rows_top[0]["logL"] > rows_top[1]["logL"]
assert abs((rows_top[0]["logL"] - rows_top[1]["logL"]) - rows_top[1]["dL"]) < 0.01
print(f"  ureC vs species tree: RF {rf:.0f}/{rf_max:.0f} (normalised {rf_norm:.3f}), "
      f"{n_taxa} taxa | p-AU tree 1 {rows_top[0]['p_au']:.3g}, tree 2 {rows_top[1]['p_au']:.3g}"
      f" | excluded at 95 %: {[r['tree'] for r in rows_top if r['excluded']]}")

# ================================================================== C, D: reference panel
man = pd.read_csv(C5 / "reference_manifest.tsv", sep="\t")
assert (man.status == "verified").all(), man.loc[man.status != "verified", "accession"].tolist()
n_refs = len(man)


def binomial(organism):
    """Genus + species of the organism NCBI reports, for use as a row label."""
    return " ".join(str(organism).split()[:2])


man["label"] = [binomial(o) for o in man.reported_organism]
# where two accessions share a binomial, the strain designation distinguishes them; it is
# read from the assembly report, not typed in, so the label always names the genome drawn
dupe = man.label.duplicated(keep=False)
man.loc[dupe, "label"] = [f"{lab} {str(strain).replace('type strain: ', '')}"
                          for lab, strain in zip(man.loc[dupe, "label"],
                                                 man.loc[dupe, "strain"])]
assert man.label.is_unique

lines = (C5 / "skani_panMICP.matrix").read_text().splitlines()
n_declared = int(lines[0])
names, rows = [], []
for ln in lines[1:]:
    parts = ln.split("\t")
    names.append(Path(parts[0]).stem)
    rows.append([float(v) for v in parts[1:]])
M = np.array(rows)
assert M.shape == (n_declared, len(names)) and M.shape[0] == M.shape[1], M.shape
assert np.allclose(M, M.T), "skani matrix is not symmetric"
idx = {nm: i for i, nm in enumerate(names)}
assert sorted(n.replace("HERO_", "") for n in names if n.startswith("HERO_")) == sorted(st.HEROES)

tab = pd.read_csv(SUPP / "Table_S20_skani_hero_vs_refs.tsv", sep="\t")
# pandas reads the TRUE/FALSE column as bool; normalise so the check below is unambiguous
tab["Alignment_reported"] = tab.Alignment_reported.astype(str).str.upper() == "TRUE"
assert len(tab) == len(st.HEROES) * n_refs, (len(tab), n_refs)
assert set(tab.Ref_accession) == set(man.accession)

ani = pd.DataFrame(np.nan, index=st.HEROES, columns=list(man.label))
acc2label = dict(zip(man.accession, man.label))
for r in tab.itertuples():
    m = M[idx[f"HERO_{r.MAG}"], idx[r.Ref_accession]]
    reported = bool(r.Alignment_reported)
    assert reported == (m > 0), (r.MAG, r.Ref_accession, m, r.Alignment_reported)
    if reported:
        assert abs(m - float(r.ANI)) < 0.01, (r.MAG, r.Ref_accession, m, r.ANI)
        ani.loc[r.MAG, acc2label[r.Ref_accession]] = m

reported = ani.stack().sort_values(ascending=False)
assert len(reported) == int(tab.Alignment_reported.sum())
print(f"  verified references {n_refs} | reported MAG-reference pairs {len(reported)}")

# every reported pair is with a reference of the MAG's own genus - stated as a check,
# derived from the manifest rather than asserted from a list
for (mag, lab), v in reported.items():
    genus_of_ref = lab.split()[0]
    expected = "Sphingobacterium" if st.LINEAGE[mag] == "Sphingobacterium" else "Pseudomonas"
    assert genus_of_ref == expected, (mag, lab)

# ================================================================== E: Pseudomonas_E screen
# feature -> (source table, column, row label); category names, not data
FEATURES = [("a", "UreC_alpha", "UreC (PF00449)"),
            ("a", "UreB_beta_gamma", "UreB (PF00699)"),
            ("a", "urease_core", "urease core (UreC + UreB)"),
            ("a", "has_CA", "any carbonic anhydrase"),
            ("b", "ureC_and_ureB_single_contig", "ureC + ureB on one contig"),
            ("b", "ureC_and_CA_single_contig", "ureC + CA on one contig")]

pa = pd.read_csv(SUPP / "Table_S13a_pseudomonas_e_MICP_rarity_screen.csv")
pb = pd.read_csv(SUPP / "Table_S13b_pseudomonas_e_single_contig.csv")
assert sorted(pa.accession) == sorted(pb.accession)
N_PSE = len(pa)
src = {"a": pa, "b": pb}

feat_rows = []
for tbl, col, label in FEATURES:
    hits = int(src[tbl][col].sum())
    feat_rows.append((label, hits, 100 * hits / N_PSE))
feat_rows = feat_rows[::-1]                 # first feature ends up at the top of the axis
feat_labels = [r[0] for r in feat_rows]
feat_hits = [r[1] for r in feat_rows]
feat_pct = [r[2] for r in feat_rows]

# ================================================================== page
LCOL, RCOL = 6.0, 92.0
H = 230.0
fig, ax_mm, text_mm, letter = st.page(H)

# ------------------------------------------------------------------ A
axA = ax_mm(15.0, 20.0, 57.0, 57.0)
caxA = ax_mm(77.0, 20.0, 3.0, 34.0)
letter(LCOL, 7.0, "A")

cmap = st.seq_cmap()
cmap.set_bad(color="#F4F4F4")
axA.imshow(shown6, cmap=cmap, vmin=vmin6, vmax=100.0, aspect="equal")
axA.set_xticks(range(len(st.HEROES)))
axA.set_xticklabels(st.HEROES, fontsize=FS_BODY)
axA.xaxis.set_ticks_position("top")
axA.set_yticks(range(len(st.HEROES)))
axA.set_yticklabels(st.HEROES, fontsize=FS_BODY)
axA.tick_params(length=0, pad=3)
for labs in (axA.get_xticklabels(), axA.get_yticklabels()):
    for lab, mag in zip(labs, st.HEROES):
        lab.set_color(st.hero_col(mag))
        lab.set_fontweight("bold")
for s in axA.spines.values():
    s.set_visible(False)
axA.set_xticks(np.arange(-0.5, len(st.HEROES), 1), minor=True)
axA.set_yticks(np.arange(-0.5, len(st.HEROES), 1), minor=True)
axA.grid(which="minor", color="white", linewidth=1.2)
axA.tick_params(which="minor", length=0)

for i in range(len(st.HEROES)):
    for j in range(len(st.HEROES)):
        v = vals6[i, j]
        if v > NO_ALIGNMENT:
            frac = (v - vmin6) / (100.0 - vmin6)
            axA.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=FS_STAT,
                     color="white" if frac > 0.55 else TEXT)
        else:
            axA.text(j, i, "—", ha="center", va="center", fontsize=FS_STAT, color=GREY)

grad = np.linspace(100.0, vmin6, 128).reshape(-1, 1)
caxA.imshow(grad, cmap=st.seq_cmap(), vmin=vmin6, vmax=100.0, aspect="auto",
            extent=[0, 1, vmin6, 100.0])
caxA.set_xticks([])
caxA.yaxis.set_ticks_position("right")
caxA.yaxis.set_major_locator(MaxNLocator(nbins=4))
caxA.tick_params(labelsize=FS_STAT, length=2, pad=1)
for s in caxA.spines.values():
    s.set_linewidth(0.6)
    s.set_color(AXIS)
text_mm(88.5, 37.0, "ANI (%)", rotation=90, ha="center", va="center", fontsize=FS_BODY)

axA.legend(handles=[Line2D([], [], color=st.SPHINGO, lw=4, label="Sphingobacterium"),
                    Line2D([], [], color=st.PSEUDO, lw=4, label="Pseudomonas_E"),
                    Rectangle((0, 0), 1, 1, facecolor="#F4F4F4", edgecolor=LIGHT,
                              label="— no alignment")],
           loc="upper left", bbox_to_anchor=(-0.14, -0.045), fontsize=FS_STAT,
           handlelength=1.2, ncol=2, columnspacing=1.2, handletextpad=0.5)

# ------------------------------------------------------------------ B: ureC gene tree
axB = ax_mm(RCOL, 14.0, 86.0, 155.0)
letter(RCOL, 7.0, "B")

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
        axB.plot([x0, x1], [y_of[id(clade)]] * 2, color=TEXT, lw=0.6, zorder=2)
    if clade.clades:
        ys = [y_of[id(c)] for c in clade.clades]
        axB.plot([x1, x1], [min(ys), max(ys)], color=TEXT, lw=0.6, zorder=2)
        bs = clade.confidence
        if bs is not None and bs >= BOOTSTRAP_SHOWN_AT and clade is not tree.root:
            # centred on the branch that leads to the node: branch midpoints are further
            # apart than the node x positions, which keeps neighbouring labels clear
            bx_ = (x0 + x1) / 2
            by_ = label_y(bx_, y_of[id(clade)] - 0.22, max_depth * 0.035)
            axB.text(bx_, by_, f"{int(bs)}", ha="center", va="bottom", fontsize=4.6,
                     color=GREY, zorder=3)

for t in tree.get_terminals():
    hero = t.name in st.HEROES
    axB.text(depths[t] + max_depth * 0.012, y_of[id(t)], t.name, va="center", ha="left",
             fontsize=FS_BODY, color=st.hero_col(t.name) if hero else TEXT,
             fontweight="bold" if hero else "normal")

axB.set_xlim(-max_depth * 0.01, max_depth * 1.18)
axB.set_ylim(len(tips) - 0.4, -1.6)
axB.axis("off")

sb_y = -1.1
axB.plot([0, SCALE_BAR], [sb_y, sb_y], color=GREY, lw=1.0)
for x in (0, SCALE_BAR):
    axB.plot([x, x], [sb_y - 0.3, sb_y + 0.3], color=GREY, lw=1.0)
axB.text(SCALE_BAR * 1.30, sb_y, f"{SCALE_BAR} substitutions per site", ha="left",
         va="center", fontsize=FS_STAT, color=GREY)
axB.text(max_depth * 1.18, sb_y, f"bootstrap ≥ {BOOTSTRAP_SHOWN_AT} shown",
         ha="right", va="center", fontsize=FS_STAT, color=GREY)

# ------------------------------------------------------------------ C: six x reference panel
axC = ax_mm(45.0, 102.0, 42.0, 66.0)
letter(LCOL, 92.0, "C")
nrC, ncC = n_refs, len(st.HEROES)
for i, lab in enumerate(man.label):
    for j, mag in enumerate(st.HEROES):
        v = ani.loc[mag, lab]
        if np.isnan(v):
            axC.text(j, i, "–", ha="center", va="center", fontsize=FS_STAT, color=GREY)
        else:
            axC.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                    facecolor=st.hero_col(mag), edgecolor="none"))
            axC.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=FS_STAT,
                     color="white")
axC.set_xlim(-0.5, ncC - 0.5)
axC.set_ylim(nrC - 0.5, -0.5)
axC.xaxis.set_ticks_position("top")
axC.set_xticks(range(ncC))
axC.set_xticklabels(st.HEROES, fontsize=FS_BODY)
for tick, mag in zip(axC.get_xticklabels(), st.HEROES):
    tick.set_color(st.hero_col(mag))
    tick.set_fontweight("bold")
axC.set_yticks(range(nrC))
axC.set_yticklabels(man.label, fontsize=FS_STAT, fontstyle="italic")
axC.tick_params(length=0, pad=2)
for sp in axC.spines.values():
    sp.set_visible(False)
axC.set_xticks(np.arange(-0.5, ncC, 1), minor=True)
axC.set_yticks(np.arange(-0.5, nrC, 1), minor=True)
axC.grid(which="minor", color=LIGHT, lw=0.5)
axC.tick_params(which="minor", length=0)
text_mm(45.0, 171.0, "–  no alignment reported by skani", fontsize=FS_STAT, color=GREY)

# ------------------------------------------------------------------ D: reported pairs, ranked
axD = ax_mm(48.0, 187.0, 38.0, 32.0)
letter(LCOL, 180.0, "D")
yD = np.arange(len(reported), dtype=float)
XMIN = 88.0
for yy, ((mag, lab), v) in zip(yD, reported.items()):
    col = st.hero_col(mag)
    axD.plot([XMIN, v], [yy, yy], lw=1.0, color=col, zorder=2)
    axD.scatter([v], [yy], s=20, color=col, zorder=3)
    axD.text(v + 0.3, yy, f"{v:.2f}", ha="left", va="center", fontsize=FS_STAT, color=TEXT)
axD.axvline(SPECIES_ANI, color=GREY, lw=0.8, ls="--", zorder=1)
axD.text(SPECIES_ANI, 1.05, "95 % species boundary", transform=axD.get_xaxis_transform(),
         ha="center", va="bottom", fontsize=FS_STAT, color=GREY)
axD.set_yticks(yD)
axD.set_yticklabels([f"{mag}  vs  {lab}" for (mag, lab), _ in reported.items()],
                    fontsize=FS_STAT)
for tick, ((mag, _lab), _v) in zip(axD.get_yticklabels(), reported.items()):
    tick.set_color(st.hero_col(mag))
axD.set_ylim(len(reported) - 0.5, -0.5)
axD.set_xlim(XMIN, 101.0)
axD.set_xlabel("skani ANI (%)")
axD.grid(axis="x", color=LIGHT, lw=0.5, zorder=0)
axD.set_axisbelow(True)
st.style_axis(axD, left=False)
axD.tick_params(left=False)

# ------------------------------------------------------------------ E: Pseudomonas_E screen
axE = ax_mm(126.0, 187.0, 34.0, 32.0)
letter(RCOL, 180.0, "E")
yE = range(len(feat_rows))
axE.barh(list(yE), feat_pct, 0.6, color=HERO, edgecolor=AXIS, linewidth=0.5)
for yi, p, h in zip(yE, feat_pct, feat_hits):
    axE.text(p + 1.5, yi, f"{p:.1f}%   {h}/{N_PSE}", va="center", ha="left",
             fontsize=FS_STAT, color=TEXT)
axE.set_yticks(list(yE))
axE.set_yticklabels(feat_labels)
axE.set_xlabel("Pseudomonas_E reference genomes (%)")
axE.set_xlim(0, 118)
axE.set_ylim(-0.6, len(feat_rows) - 0.4)
axE.set_xticks([0, 25, 50, 75, 100])
st.style_axis(axE)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S4")
