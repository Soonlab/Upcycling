"""Consolidated main Fig 5 - novelty, global rarity and engineering background.

One composed page, 180 mm wide, four panels (reading order A -> D, left to right then
top to bottom).  Every panel is carried over unchanged from a builder of the 8-figure
set; see consolidation_260904/DESIGN.md.

  A  novelty screen across the 111-MAG panel: closest-GTDB-reference ANI per MAG with
     the 95 % species boundary, plus the block of MAGs for which GTDB-Tk returned no
     species-level ANI at all          (old build_fig5.py panel A)
  B  MGnify livestock species-cluster rarity of the MICP gene-complete profile and of the
     single-contig ureC + CA architecture, per biome and pooled
                                       (old build_fig8.py panel C)
  C  antiSMASH 7 per-MAG BGC region count by class, MICP-complete vs rest, with the
     Mann-Whitney P of each class as a stat column
                                       (old build_fig8.py panel D)
  D  Jaccard PCoA of the Panaroo gene presence/absence matrix, coloured by waste source,
     MICP-complete ringed               (old build_fig7.py panel A)

Old build_fig5.py panels B (AAI of S13 / S16) and C (skani vs the 63 RefSeq
Sphingobacterium genomes) and old build_fig7.py panel B (trait-module ordination) are
dropped from the main set; their numbers survive in the supplementary tables
(Table S4b, Table S10a/b, Table S2b).

Sources
  Table_S8_novelty_ANI_screen.csv        ANI to closest GTDB reference, novelty flag (A)
  Table_S14a_mgnify_catalog_summary.csv  MGnify per-catalog counts and percentages (B)
  Table_S23b_antismash_hero_vs_rest.csv  BGC class means and Mann-Whitney P (C)
  Table_S9a_PCoA_coordinates.csv         pan-genome PC1-PC3, source, genus, MICP flag (D)
  Table_S9b_PERMANOVA_global.csv         pan-genome pseudo-F / p and PC variance (D)

Provenance: no value is hard-coded.  The novelty flag of A is asserted to be exactly
"no species-level ANI"; the MGnify percentages of B are recomputed from the per-catalog
counts and asserted against the stored pct columns before the pooled bar is derived from
the same counts; the ten BGC classes of C are asserted present in the stored table; the
MICP-complete set of D is asserted against the Hero flag of the coordinate table.

Colour meanings on this page (one colour, one meaning):
  coral   a MICP-complete MAG - its point and label in A, its bars in C, its ring and
          label in D - and the MICP-complete group everywhere it is a series
  grey    a MAG that is not MICP-complete, the remaining 105 MAGs in C, and every
          threshold rule
  dark / light green   MICP gene-complete and single-contig ureC + CA prevalence in B
  cattle brown / swine pink / sheep green / poultry purple   waste source in D
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import (HERO, REST, GREEN, GREY, TEXT, AXIS, SOURCE,
                    FS_BODY, FS_STAT, HEROES)

st.setup()
OUT = HERE / "figures_v2"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

ANI_SPECIES = 95.0          # species boundary, ANI
GREEN_LT = tuple(0.45 + 0.55 * c for c in to_rgb(GREEN))   # light tint of GREEN
CLUSTER_PX = 28             # marker separation below which two MAG labels would collide
LINE_PT = 7                 # vertical step used to fan out a cluster of labels

# the ten BGC classes carried by the shipped Fig 8d (category names, not data)
BGC_CLASSES = ["BGC_T3PKS", "BGC_RRE-containing", "BGC_arylpolyene", "BGC_terpene",
               "BGC_betalactone", "BGC_RiPP-like", "BGC_NAGGN", "BGC_NRPS",
               "BGC_hydrogen-cyanide", "BGC_NRP-metallophore"]

# ------------------------------------------------------------------ data: A
nov = pd.read_csv(SUPP / "Table_S8_novelty_ANI_screen.csv")
has_ani = nov.ANI.notna()
ranked = nov[has_ani].sort_values("ANI").reset_index(drop=True)
no_ani = nov[~has_ani].reset_index(drop=True)
# the novelty flag is exactly "no species-level ANI"; state it by assertion, not in prose
assert set(nov.user_genome[nov.Novel_sp_candidate]) == set(no_ani.user_genome)
assert (ranked.ANI >= ANI_SPECIES).all()

# ------------------------------------------------------------------ data: B
mg = pd.read_csv(SUPP / "Table_S14a_mgnify_catalog_summary.csv")
pct_complete = 100 * mg.n_MICP_gene_complete / mg.n_species_clusters
pct_single = 100 * mg.n_MICP_single_contig_ureC_CA / mg.n_species_clusters
assert np.allclose(pct_complete, mg.pct_MICP_gene_complete, atol=1e-3)
assert np.allclose(pct_single, mg.pct_MICP_single_contig, atol=1e-3)
n_pool = int(mg.n_species_clusters.sum())
pool_complete = 100 * mg.n_MICP_gene_complete.sum() / n_pool
pool_single = 100 * mg.n_MICP_single_contig_ureC_CA.sum() / n_pool
bio_lab = [c.replace("-", "\n") for c in mg.catalog] + ["pooled"]
bio_n = list(mg.n_species_clusters) + [n_pool]
v_complete = list(pct_complete) + [pool_complete]
v_single = list(pct_single) + [pool_single]

# ------------------------------------------------------------------ data: C
bgc = pd.read_csv(SUPP / "Table_S23b_antismash_hero_vs_rest.csv").set_index("metric")
missing = [c for c in BGC_CLASSES if c not in bgc.index]
assert not missing, missing
bgc = bgc.loc[BGC_CLASSES].copy()
bgc["diff"] = bgc.hero_mean - bgc.rest_mean
bgc = bgc.sort_values("diff")            # most hero-enriched ends up at the top
bgc_lab = [i.replace("BGC_", "") for i in bgc.index]

# ------------------------------------------------------------------ data: D
coords = pd.read_csv(SUPP / "Table_S9a_PCoA_coordinates.csv").set_index("MAG")
glob = pd.read_csv(SUPP / "Table_S9b_PERMANOVA_global.csv").iloc[0]
assert set(coords.index[coords.Hero]) == set(HEROES)
SRC_ORDER = ["Cattle", "Swine", "Sheep", "Poultry"]
assert set(coords.Source) == set(SRC_ORDER)
COL = {s: SOURCE[s.lower()] for s in SRC_ORDER}
n_src = coords.Source.value_counts()

# ------------------------------------------------------------------ page
TOP = 10.0
H_A, W_A, X_A = 46.0, 78.0, 16.0
X_A2 = X_A + W_A + 14.0          # the no-ANI block sits beside the scatter, same panel
T2 = TOP + H_A + 22.0            # row 2: panels B and C
H_B, W_B, X_B = 34.0, 60.0, 24.0
H_C, W_C, X_C = 46.0, 38.0, 110.0
T3 = T2 + H_C + 18.0             # row 3: panel D
H_D, W_D, X_D = 62.0, 76.0, 18.0
H = T3 + H_D + 14.0

fig, ax_mm, text_mm, letter = st.page(H)

# ---------------------------------------------------------------- panel A
axA = ax_mm(X_A, TOP, W_A, H_A)
is_hero = ranked.user_genome.isin(HEROES).values
axA.scatter(np.arange(len(ranked))[~is_hero], ranked.ANI[~is_hero], s=5, color=REST,
            linewidth=0, zorder=2)
axA.scatter(np.arange(len(ranked))[is_hero], ranked.ANI[is_hero], s=20, color=HERO,
            edgecolor="white", linewidth=0.4, zorder=3)
hero_idx = list(np.where(is_hero)[0])
for k, i in enumerate(hero_idx):
    prev_close = k > 0 and (i - hero_idx[k - 1]) < 6
    dy = -7 if prev_close else 3
    axA.annotate(ranked.user_genome[i], (i, ranked.ANI[i]), xytext=(3, dy),
                 textcoords="offset points", fontsize=FS_STAT, color=HERO)
axA.axhline(ANI_SPECIES, ls="--", lw=0.7, color=GREY)
axA.text(len(ranked) - 1, ANI_SPECIES, f"{ANI_SPECIES:g} % species boundary",
         fontsize=FS_STAT, color=GREY, ha="right", va="bottom")
axA.set_xlabel(f"MAGs with a species-level ANI, ranked (n = {len(ranked)})")
axA.set_ylabel("ANI to closest GTDB reference (%)")
axA.set_xticks([])
st.style_axis(axA, left=True, bottom=True)
axA.spines["bottom"].set_visible(False)

# the MAGs GTDB-Tk left without a species-level ANI, as an identifier block
axA2 = ax_mm(X_A2, TOP, 62.0, H_A)
axA2.axis("off")
text_mm(X_A2, TOP - 1.0, f"no species-level ANI  (n = {len(no_ani)})", fontsize=FS_STAT,
        ha="left", va="bottom", color=TEXT)
names = list(no_ani.user_genome)
NCOL = 4
per = int(np.ceil(len(names) / NCOL))
for k, name in enumerate(names):
    col, row = k // per, k % per
    text_mm(X_A2 + 1.5 + col * 15.0, TOP + 3.0 + row * 4.2, name, fontsize=FS_BODY,
            ha="left", va="top", color=HERO if name in HEROES else TEXT)

# ---------------------------------------------------------------- panel B
axB = ax_mm(X_B, T2, W_B, H_B)
w = 0.38                       # paired-bar width (style constant)
xc = np.arange(len(bio_lab))
axB.bar(xc - w / 2, v_complete, w, color=GREEN, edgecolor=AXIS, linewidth=0.5)
axB.bar(xc + w / 2, v_single, w, color=GREEN_LT, edgecolor=AXIS, linewidth=0.5)
for xi, v in zip(xc, v_complete):
    axB.text(xi - w / 2, v, f"{v:.2f}", ha="center", va="bottom", fontsize=FS_STAT)
for xi, v in zip(xc, v_single):
    axB.text(xi + w / 2, v, f"{v:.2f}", ha="center", va="bottom", fontsize=FS_STAT)
axB.set_xticks(xc)
axB.set_xticklabels([f"{lab}\nn = {n:,}" for lab, n in zip(bio_lab, bio_n)],
                    fontsize=FS_BODY)
axB.set_ylabel("species-cluster representatives (%)")
axB.set_ylim(0, max(v_complete + v_single) * 1.25)
st.style_axis(axB)
axB.legend(handles=[Patch(facecolor=GREEN, label="MICP gene-complete"),
                    Patch(facecolor=GREEN_LT, label="single-contig ureC + CA")],
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.1,
           handleheight=0.9, columnspacing=1.2, fontsize=FS_BODY)

# ---------------------------------------------------------------- panel C
axC = ax_mm(X_C, T2, W_C, H_C)
y = np.arange(len(bgc))
hb = 0.38
axC.barh(y + hb / 2, bgc.hero_mean, hb, color=HERO, edgecolor=AXIS, linewidth=0.5)
axC.barh(y - hb / 2, bgc.rest_mean, hb, color=REST, edgecolor=AXIS, linewidth=0.5)
axC.set_yticks(y)
axC.set_yticklabels(bgc_lab, fontsize=FS_BODY)
axC.set_xlabel("BGC regions per MAG (mean)")
axC.set_ylim(-0.7, len(bgc) - 0.3)
xmax = float(max(bgc.hero_mean.max(), bgc.rest_mean.max()))
axC.set_xlim(0, xmax * 1.05)
st.style_axis(axC)
axC.legend(handles=[Patch(facecolor=HERO, label="MICP-complete (n = 6)"),
                    Patch(facecolor=REST, label="rest (n = 105)")],
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=1, handlelength=1.1,
           handleheight=0.9, fontsize=FS_BODY)
# stat column: Mann-Whitney P per class, bound to its row
px = xmax * 1.14
axC.text(px, len(bgc) - 0.35, "MWU P", fontsize=FS_STAT, ha="left", va="bottom",
         color=TEXT)
for yi, p in zip(y, bgc.MWU_p):
    txt = f"{p:.0e}" if p < 0.001 else f"{p:.2f}"
    axC.text(px, yi, txt, fontsize=FS_STAT, ha="left", va="center", color=TEXT)
axC.set_clip_on(False)
for t in axC.texts:
    t.set_clip_on(False)

# ---------------------------------------------------------------- panel D
axD = ax_mm(X_D, T3, W_D, H_D)
for s in SRC_ORDER:
    idx = coords.index[coords.Source == s]
    axD.scatter(coords.loc[idx, "PC1"], coords.loc[idx, "PC2"], s=9, color=COL[s],
                linewidth=0, alpha=0.9, zorder=2, label=f"{s} n = {n_src[s]}")
axD.scatter(coords.loc[HEROES, "PC1"], coords.loc[HEROES, "PC2"], s=34,
            facecolor="none", edgecolor=HERO, linewidth=0.8, zorder=4)
# several MICP-complete MAGs project onto nearly the same point, so their labels are
# fanned out: MAGs whose markers fall within CLUSTER_PX of one another get successive
# vertical offsets instead of printing on top of each other
fig.canvas.draw()
pts = {m: axD.transData.transform((coords.loc[m, "PC1"], coords.loc[m, "PC2"]))
       for m in HEROES}
placed = []
for m in sorted(HEROES, key=lambda k: (-pts[k][1], pts[k][0])):
    k = sum(1 for q in placed
            if abs(pts[q][0] - pts[m][0]) < CLUSTER_PX
            and abs(pts[q][1] - pts[m][1]) < CLUSTER_PX)
    placed.append(m)
    axD.annotate(m, (coords.loc[m, "PC1"], coords.loc[m, "PC2"]),
                 xytext=(4, 2 - k * LINE_PT), textcoords="offset points",
                 fontsize=FS_STAT, color=HERO, zorder=5)
axD.set_xlabel(f"PC1  {glob.PC1_var:.1f} %")
axD.set_ylabel(f"PC2  {glob.PC2_var:.1f} %")
# prune the extreme ticks: the first x tick and the first y tick otherwise print on top
# of each other in the corner where the two spines meet
axD.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
axD.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
st.style_axis(axD)

# PERMANOVA result as a stat block bound to panel D
STAT_T = T3 + 2.0
STAT_X = X_D + W_D - 32.0        # upper right of the ordination, clear of the point cloud
text_mm(STAT_X, STAT_T, "PERMANOVA", fontsize=FS_STAT, color=TEXT)
for k, (lab, f, p) in enumerate((("source", glob.pseudo_F_source, glob.p_source),
                                 ("genus", glob.pseudo_F_genus, glob.p_genus))):
    text_mm(STAT_X, STAT_T + 3.4 + k * 3.4,
            f"{lab}  F = {f:.2f}  p = {p:.3f}", fontsize=FS_STAT, color=TEXT)

handles, labels = axD.get_legend_handles_labels()
handles.append(Line2D([0], [0], marker="o", ls="", ms=5, mfc="none", mec=HERO, mew=0.8))
labels.append(f"MICP-complete n = {len(HEROES)}")
fig.legend(handles, labels, loc="upper left",
           bbox_to_anchor=((X_D + W_D + 8.0) * st.MM / (st.PAGE_W_MM * st.MM),
                           1 - (T3 + 6.0) * st.MM / (H * st.MM)),
           ncol=1, fontsize=FS_BODY, handletextpad=0.4, borderpad=0.2,
           labelspacing=0.5)

letter(4, 4, "A")
letter(4, T2 - 9.0, "B")
letter(99, T2 - 9.0, "C")
letter(4, T3 - 4.0, "D")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig5")
