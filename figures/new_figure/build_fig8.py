"""Main Fig 8 - convergent retention of urease catalytic function, external rarity
and the lineage-specific secondary-metabolism fingerprint.

Panels (reading order A -> D, left to right then top to bottom):
  A  UreC active-site residues, 6 MICP-complete MAGs x 7 canonical catalytic sites
     (Table S12); cell colour encodes agreement with the S. pasteurii reference
  B  ESMFold UreC backbone agreement with PDB 4CEU chain C (Table S22): TM-score
     (filled, left axis) and all-residue backbone RMSD (hatched, right axis)
  C  MGnify livestock species-cluster prevalence of the MICP gene-complete profile and
     of the single-contig ureC + CA architecture, per biome and pooled (Table S14a)
  D  antiSMASH 7 per-MAG BGC region count by class, MICP-complete vs rest, with the
     Mann-Whitney P of each class as a stat column (Table S23b)

Provenance: every number is read from the supplementary tables.  Pooled MGnify
percentages in C are recomputed from the per-catalog counts and asserted against the
stored pct columns; the per-catalog percentages are likewise recomputed and asserted.
Panel B recomputes nothing (TM-score and RMSD are the stored measurements) but the MAG
labels are stripped of the "_UreC" suffix carried in the table.

Colour meanings on this page (one meaning per colour):
  green   = observed residue matches the reference (A only)
  coral   = observed residue differs from the reference (A, legend entry only - no cell
            takes this colour), and the MICP-complete group in D
  blue    = Sphingobacterium lineage (row labels in A, bars in B)
  orange  = Pseudomonas_E lineage (row labels in A, bars in B)
  grey    = the remaining 105 MAGs in D, and the TM = 0.5 reference line in B
  dark/light green = gene-complete / single-contig prevalence in C
Blue is deliberately NOT used for the C prevalence series, so that blue means
"Sphingobacterium lineage" and nothing else on this page.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import (HERO, REST, SPHINGO, PSEUDO, GREEN, GREY, TEXT, AXIS,
                    FS_BODY, FS_STAT, HEROES, hero_col)

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

TM_SAME_FOLD = 0.5          # Xu & Zhang 2010 same-fold threshold (published constant)
GREEN_LT = tuple(0.45 + 0.55 * c for c in to_rgb(GREEN))   # light tint of GREEN

# the ten BGC classes carried by the shipped Fig 8d (category names, not data)
BGC_CLASSES = ["BGC_T3PKS", "BGC_RRE-containing", "BGC_arylpolyene", "BGC_terpene",
               "BGC_betalactone", "BGC_RiPP-like", "BGC_NAGGN", "BGC_NRPS",
               "BGC_hydrogen-cyanide", "BGC_NRP-metallophore"]

# ------------------------------------------------------------------ data: A
sites = pd.read_csv(SUPP / "Table_S12_UreC_active_site_residues.csv")
sites["pos"] = sites.site.str.extract(r"(\d+)").astype(int)
sites = sites.sort_values("pos").reset_index(drop=True)      # reading order by position
obs = sites[HEROES].values.T                                  # rows = MAG, cols = site
exp = sites.expected.values
match = obs == exp[None, :]

# ------------------------------------------------------------------ data: B
tm = pd.read_csv(SUPP / "Table_S22_ureC_vs_4CEU_tm.csv")
tm["MAG"] = tm.MAG.str.replace("_UreC", "", regex=False)
tm = tm.set_index("MAG").loc[HEROES].reset_index()
assert (tm.ref_len == tm.ref_len.iloc[0]).all()               # one reference chain

# ------------------------------------------------------------------ data: C
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

# ------------------------------------------------------------------ data: D
bgc = pd.read_csv(SUPP / "Table_S23b_antismash_hero_vs_rest.csv").set_index("metric")
missing = [c for c in BGC_CLASSES if c not in bgc.index]
assert not missing, missing
bgc = bgc.loc[BGC_CLASSES].copy()
bgc["diff"] = bgc.hero_mean - bgc.rest_mean
bgc = bgc.sort_values("diff")            # most hero-enriched ends up at the top
bgc_lab = [i.replace("BGC_", "") for i in bgc.index]

# ------------------------------------------------------------------ page
H = 122.0
fig, ax_mm, text_mm, letter = st.page(H)
L1, L2 = 13.0, 101.0                     # the two panel-letter x positions
COL1, COL2 = 24.0, 110.0

# ---- A -------------------------------------------------------------------
letter(L1, 4, "A")
axA = ax_mm(COL1, 13, 60, 32)
axA.imshow(np.where(match, 1, 0), cmap=st.seq_cmap("match", hi=GREEN),
           vmin=0, vmax=1, aspect="auto")
for i in range(obs.shape[0]):
    for j in range(obs.shape[1]):
        axA.text(j, i, obs[i, j], ha="center", va="center", fontsize=FS_BODY,
                 color="white" if match[i, j] else TEXT)
axA.set_xticks(range(len(sites)))
axA.set_xticklabels(sites.site, fontsize=FS_BODY)
axA.set_yticks(range(len(HEROES)))
axA.set_yticklabels(HEROES, fontsize=FS_BODY)
for tick, mag in zip(axA.get_yticklabels(), HEROES):
    tick.set_color(hero_col(mag))
axA.set_xlabel("UreC active-site residue (P41020 / 4CEU numbering)")
axA.tick_params(length=0)
for s in axA.spines.values():
    s.set_visible(False)
axA.legend(handles=[Patch(facecolor=GREEN, label="matches reference"),
                    Patch(facecolor=HERO, label="differs")],
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.1,
           handleheight=0.9, columnspacing=1.2, fontsize=FS_BODY)

# ---- B -------------------------------------------------------------------
letter(L2, 4, "B")
axB = ax_mm(COL2, 13, 24, 32)
axB2 = ax_mm(COL2 + 34, 13, 24, 32)
x = np.arange(len(tm))
w = 0.38                       # paired-bar width (style constant)
cols = [hero_col(m) for m in tm.MAG]
axB.bar(x, tm.tm_norm_ref, 0.62, color=cols, edgecolor=AXIS, linewidth=0.5)
axB.axhline(TM_SAME_FOLD, color=GREY, lw=0.8, ls="--")
axB2.bar(x, tm.rmsd, 0.62, facecolor="white", edgecolor=AXIS, linewidth=0.5, hatch="////")
for xi, v in zip(x, tm.tm_norm_ref):
    axB.text(xi, v, f"{v:.3f}", ha="center", va="bottom", fontsize=FS_STAT, rotation=90)
for xi, v in zip(x, tm.rmsd):
    axB2.text(xi, v, f"{v:.2f}", ha="center", va="bottom", fontsize=FS_STAT, rotation=90)
for ax_, ylab, top in ((axB, "TM-score (norm. 4CEU chain C)", float(tm.tm_norm_ref.max())),
                       (axB2, "backbone RMSD (Å)", float(tm.rmsd.max()))):
    ax_.set_xticks(x)
    ax_.set_xticklabels(tm.MAG, rotation=90)
    for tick, mag in zip(ax_.get_xticklabels(), tm.MAG):
        tick.set_color(hero_col(mag))
    ax_.set_ylabel(ylab)
    ax_.set_ylim(0, top * 1.45)
    ax_.set_xlim(-0.7, len(x) - 0.3)
    st.style_axis(ax_)
axB.legend(handles=[Patch(facecolor=SPHINGO, label="Sphingobacterium"),
                    Patch(facecolor=PSEUDO, label="Pseudomonas_E"),
                    Line2D([], [], color=GREY, ls="--", lw=0.8, label="TM = 0.5")],
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, handlelength=1.1,
           handleheight=0.9, columnspacing=1.0, fontsize=FS_BODY)

# ---- C -------------------------------------------------------------------
letter(L1, 54, "C")
axC = ax_mm(COL1, 64, 60, 34)
xc = np.arange(len(bio_lab))
axC.bar(xc - w / 2, v_complete, w, color=GREEN, edgecolor=AXIS, linewidth=0.5)
axC.bar(xc + w / 2, v_single, w, color=GREEN_LT, edgecolor=AXIS, linewidth=0.5)
for xi, v in zip(xc, v_complete):
    axC.text(xi - w / 2, v, f"{v:.2f}", ha="center", va="bottom", fontsize=FS_STAT)
for xi, v in zip(xc, v_single):
    axC.text(xi + w / 2, v, f"{v:.2f}", ha="center", va="bottom", fontsize=FS_STAT)
axC.set_xticks(xc)
axC.set_xticklabels([f"{lab}\nn = {n:,}" for lab, n in zip(bio_lab, bio_n)],
                    fontsize=FS_BODY)
axC.set_ylabel("species-cluster representatives (%)")
axC.set_ylim(0, max(v_complete + v_single) * 1.25)
st.style_axis(axC)
axC.legend(handles=[Patch(facecolor=GREEN, label="MICP gene-complete"),
                    Patch(facecolor=GREEN_LT, label="single-contig ureC + CA")],
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.1,
           handleheight=0.9, columnspacing=1.2, fontsize=FS_BODY)

# ---- D -------------------------------------------------------------------
letter(L2, 54, "D")
axD = ax_mm(COL2, 64, 38, 46)
y = np.arange(len(bgc))
hb = 0.38
axD.barh(y + hb / 2, bgc.hero_mean, hb, color=HERO, edgecolor=AXIS, linewidth=0.5)
axD.barh(y - hb / 2, bgc.rest_mean, hb, color=REST, edgecolor=AXIS, linewidth=0.5)
axD.set_yticks(y)
axD.set_yticklabels(bgc_lab, fontsize=FS_BODY)
axD.set_xlabel("BGC regions per MAG (mean)")
axD.set_ylim(-0.7, len(bgc) - 0.3)
xmax = float(max(bgc.hero_mean.max(), bgc.rest_mean.max()))
axD.set_xlim(0, xmax * 1.05)
st.style_axis(axD)
axD.legend(handles=[Patch(facecolor=HERO, label="MICP-complete (n = 6)"),
                    Patch(facecolor=REST, label="rest (n = 105)")],
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=1, handlelength=1.1,
           handleheight=0.9, fontsize=FS_BODY)
# stat column: Mann-Whitney P per class, bound to its row
px = xmax * 1.14
axD.text(px, len(bgc) - 0.35, "MWU P", fontsize=FS_STAT, ha="left", va="bottom",
         color=TEXT)
for yi, p in zip(y, bgc.MWU_p):
    txt = f"{p:.0e}" if p < 0.001 else f"{p:.2f}"
    axD.text(px, yi, txt, fontsize=FS_STAT, ha="left", va="center", color=TEXT)
axD.set_clip_on(False)
for t in axD.texts:
    t.set_clip_on(False)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig8")
