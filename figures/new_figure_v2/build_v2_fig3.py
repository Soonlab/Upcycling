"""Consolidated Fig 3 - catalytic and structural conservation of urease.

Consolidation of 2026-09-04 (see /data/data/Upcycling/consolidation_260904/DESIGN.md).
Old panels -> new page:

  old Fig 8A  ->  A  UreC active-site residues, 6 MICP-complete MAGs x 7 canonical
                     catalytic sites (Table S12), each shown with its FLANK flanking
                     alignment columns on either side so that the invariant catalytic
                     column stands against its variable neighbours (revision of
                     2026-09-09: the seven-column version was uniformly green and showed
                     nothing).  The reference sequence is the top row; cell colour
                     encodes identity to the reference residue in that column.
  old Fig 8B  ->  B  ESMFold UreC backbone agreement with PDB 4CEU chain C (Table S22):
                     TM-score (left sub-axis) and all-residue backbone RMSD (right
                     sub-axis).  The two sub-axes are one lettered panel.
  old S17C    ->  C  PAML codeml M0 omega per urease gene (Table S19b), a vertical
                     lollipop in the gene order of D and at the height of D
  old S17D    ->  D  yn00 pairwise omega by gene and pair class, log10 omega axis,
                     drawn from the per-pair yn00 output rather than the stored medians

Provenance: every number is read from a repository source.  Panel B recomputes nothing
(TM-score and RMSD are the stored measurements) but the MAG labels are stripped of the
"_UreC" suffix carried in the table.  Panel D reapplies the pair filter of the analysis
script (`run_dnds_v3.py`: 0 < omega < 99 and dS > 0.01) and asserts that it reproduces
every count, median and Mann-Whitney P stored in Table S19c.

Sources
  Table_S12_UreC_active_site_residues.csv        7 sites x 6 MAGs + expected residue and
                                                 the 0-based MSA column of each site
  research/additional/A2_structure/UreC_aligned.faa   the MSA behind Table S12 (A)
  Table_S22_ureC_vs_4CEU_tm.csv                  TM-score and backbone RMSD per MAG
  Table_S19b_codeml_M0_summary.csv               4 genes, codeml M0 omega
  Table_S19c_yn00_hero_vs_rest_summary.csv       stored medians / n / MWU P (asserted)
  research/additional/C3_dnds_codon/yn00_pairwise.csv   the pairwise distribution in D

Colour meanings on this page (one meaning per colour):
  green        residue identical to the reference residue in that column (A)
  light grey   residue differs from the reference (A); white, an alignment gap
  coral        in D a pair whose two members are both MICP-complete; in C the coral
               lollipop marks the urease genes of the MICP-complete lineages under study
  light coral  a MICP-complete x rest pair (D)
  blue         Sphingobacterium lineage (row labels in A, bars in B)
  orange       Pseudomonas_E lineage (row labels in A, bars in B)
  grey         reference lines (TM = 0.5 in B, omega = 1 in D) and rest x rest pairs (D)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from Bio import SeqIO
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import (HERO, HERO_LT, REST, SPHINGO, PSEUDO, GREEN, GREY, TEXT, AXIS,
                    LIGHT, FS_BODY, FS_STAT, HEROES, hero_col)

st.setup()
OUT = HERE / "figures_v2"
SUPP = Path(gh.SUPP)
C3 = Path(gh.ADDITIONAL) / "C3_dnds_codon"

TM_SAME_FOLD = 0.5          # Xu & Zhang 2010 same-fold threshold (published constant)
# the pair filter of run_dnds_v3.py, reproduced so panel D matches Table S19c
OMEGA_MAX, OMEGA_MIN, DS_MIN = 99.0, 0.0, 0.01
PAIR_CLASSES = [("hero-hero", HERO), ("hero-rest", HERO_LT), ("rest-rest", REST)]
PAIR_LABEL = {"hero-hero": "MICP-complete × MICP-complete",
              "hero-rest": "MICP-complete × rest", "rest-rest": "rest × rest"}

# ------------------------------------------------------------------ data: A
FLANK = 2                    # alignment columns drawn either side of each catalytic site
sites = pd.read_csv(SUPP / "Table_S12_UreC_active_site_residues.csv")
sites["pos"] = sites.site.str.extract(r"(\d+)").astype(int)
sites = sites.sort_values("pos").reset_index(drop=True)      # reading order by position
msa = {r.id.split("__")[0]: str(r.seq)
       for r in SeqIO.parse(Path(gh.ADDITIONAL) / "A2_structure/UreC_aligned.faa", "fasta")}
REF = [k for k in msa if "P41020" in k][0]
ROWS = [REF] + HEROES
assert all(m in msa for m in HEROES), sorted(msa)
assert len({len(v) for v in msa.values()}) == 1
# the stored table records, per site, the 0-based MSA column and the residue read there
for _, r in sites.iterrows():
    c = int(r.ref_column)
    assert msa[REF][c] == r.ref_aa == r.expected, (r.site, msa[REF][c], r.ref_aa)
    for m in HEROES:
        assert msa[m][c] == r[m], (r.site, m)
blocks = []                  # per site: list of (residue per row) for the FLANK window
for _, r in sites.iterrows():
    c = int(r.ref_column)
    blocks.append([[msa[m][k] for m in ROWS] for k in range(c - FLANK, c + FLANK + 1)])
n_match = sum(msa[m][int(r.ref_column)] == r.expected
              for _, r in sites.iterrows() for m in HEROES)
assert n_match == len(sites) * len(HEROES)                     # the 42/42 of the text

# ------------------------------------------------------------------ data: B
tm = pd.read_csv(SUPP / "Table_S22_ureC_vs_4CEU_tm.csv")
tm["MAG"] = tm.MAG.str.replace("_UreC", "", regex=False)
tm = tm.set_index("MAG").loc[HEROES].reset_index()
assert (tm.ref_len == tm.ref_len.iloc[0]).all()               # one reference chain

# ------------------------------------------------------------------ data: C
m0 = pd.read_csv(SUPP / "Table_S19b_codeml_M0_summary.csv")
m0 = m0.sort_values("omega_M0").reset_index(drop=True)
assert (m0.omega_M0 < 1).all(), "codeml M0 omega above 1 would not be purifying selection"

# ------------------------------------------------------------------ data: D
yn = pd.read_csv(C3 / "yn00_pairwise.csv")
n_all = len(yn)
yn = yn[(yn.omega < OMEGA_MAX) & (yn.omega > OMEGA_MIN) & (yn.dS > DS_MIN)]
heroes = set(HEROES)
in_h_a = yn.a.isin(heroes)
in_h_b = yn.b.isin(heroes)
yn = yn.assign(pair=np.where(in_h_a & in_h_b, "hero-hero",
                             np.where(~in_h_a & ~in_h_b, "rest-rest", "hero-rest")))
print(f"  yn00 pairs: {n_all} computed | {len(yn)} pass the run_dnds_v3 filter")

stored_yn = pd.read_csv(SUPP / "Table_S19c_yn00_hero_vs_rest_summary.csv").set_index("gene")
GENES = stored_yn.index.tolist()
yn_p = {}
for g in GENES:
    x = yn[yn.gene == g]
    hh = x.loc[x.pair == "hero-hero", "omega"]
    rr = x.loc[x.pair == "rest-rest", "omega"]
    hr = x.loc[x.pair == "hero-rest", "omega"]
    s = stored_yn.loc[g]
    assert len(hh) == s.hero_hero_n and len(rr) == s.rest_rest_n and len(hr) == s.hero_rest_n
    assert abs(hh.median() - s.hero_hero_median) < 1e-9
    assert abs(rr.median() - s.rest_rest_median) < 1e-9
    assert abs(hr.median() - s.hero_rest_median) < 1e-9
    p = mannwhitneyu(hh, rr, alternative="two-sided").pvalue
    assert abs(p - s.MWU_hh_vs_rr_p) < 1e-15, (g, p, s.MWU_hh_vs_rr_p)
    yn_p[g] = p

# ------------------------------------------------------------------ page
T1, PH1 = 15.0, 34.0                     # row 1: A and B
T2, PH2 = 74.0, 46.0                     # row 2: C and D, one height
H = T2 + PH2 + 15.0
fig, ax_mm, text_mm, letter = st.page(H)

L_A, X_A, W_A = 4.0, 24.0, 80.0
L_B, X_B, W_B = 108.0, 120.0, 21.0
GAP_B = 32.0
L_C, X_C, W_C = 4.0, 18.0, 36.0
L_D, X_D, W_D = 60.0, 72.0, 97.0

# ---- A: active-site residues with their flanking alignment columns ---------
letter(L_A, 4, "A")
axA = ax_mm(X_A, T1, W_A, PH1)
n_col = 2 * FLANK + 1
BLOCK_GAP = 1.0                          # empty column between site blocks
xpos = []                                # x of every drawn column
for b in range(len(blocks)):
    x0 = b * (n_col + BLOCK_GAP)
    xpos.append([x0 + k for k in range(n_col)])
ref_seq = msa[REF]
for b, (blk, xs) in enumerate(zip(blocks, xpos)):
    c = int(sites.ref_column.iloc[b])
    for k, (col, x) in enumerate(zip(blk, xs)):
        ref_aa = ref_seq[c - FLANK + k]
        for i, aa in enumerate(col):
            if aa == "-":
                fc, tc = "white", GREY
            elif aa == ref_aa:
                fc, tc = GREEN, "white"
            else:
                fc, tc = LIGHT, TEXT
            axA.add_patch(Rectangle((x - 0.5, i - 0.5), 1, 1, facecolor=fc,
                                    edgecolor="white", lw=0.4))
            axA.text(x, i, aa, ha="center", va="center", fontsize=FS_STAT, color=tc)
    # the catalytic column is outlined
    axA.add_patch(Rectangle((xs[FLANK] - 0.5, -0.5), 1, len(ROWS), facecolor="none",
                            edgecolor=TEXT, lw=0.8, zorder=3))
axA.set_xlim(-0.6, xpos[-1][-1] + 0.6)
axA.set_ylim(len(ROWS) - 0.5, -0.5)
axA.set_xticks([xs[FLANK] for xs in xpos])
axA.set_xticklabels(sites.site, fontsize=FS_BODY)
axA.set_yticks(range(len(ROWS)))
axA.set_yticklabels(["P41020"] + HEROES, fontsize=FS_BODY)
for tick, m in zip(axA.get_yticklabels(), ROWS):
    tick.set_color(TEXT if m == REF else hero_col(m))
axA.set_xlabel(f"UreC active-site residue ± {FLANK} alignment columns "
               "(P41020 / 4CEU numbering)")
axA.tick_params(length=0)
for sp in axA.spines.values():
    sp.set_visible(False)
axA.legend(handles=[Patch(facecolor=GREEN, label="identical to reference"),
                    Patch(facecolor=LIGHT, label="differs"),
                    Patch(facecolor="white", edgecolor=LIGHT, label="gap")],
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, handlelength=1.1,
           handleheight=0.9, columnspacing=1.2, fontsize=FS_BODY)

# ---- B: ESMFold TM-score and backbone RMSD (two sub-axes, one panel) ------
letter(L_B, 4, "B")
axB = ax_mm(X_B, T1, W_B, PH1)
axB2 = ax_mm(X_B + GAP_B, T1, W_B, PH1)
x = np.arange(len(tm))
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
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.1,
           handleheight=0.9, columnspacing=1.0, fontsize=FS_BODY)

# ---- C: codeml M0 omega per urease gene, in the gene order of D -------------
letter(L_C, T2 - 8.0, "C")
axC = ax_mm(X_C, T2, W_C, PH2)
m0 = m0.set_index("gene").loc[GENES].reset_index()
xC = np.arange(len(m0), dtype=float)
axC.vlines(xC, 0, m0.omega_M0, color=HERO, lw=1.2, zorder=2)
axC.scatter(xC, m0.omega_M0, s=26, color=HERO, zorder=3)
for xx, v in zip(xC, m0.omega_M0):
    axC.text(xx, v + m0.omega_M0.max() * 0.06, f"{v:.3f}", ha="center", va="bottom",
             fontsize=FS_STAT, color=TEXT)
axC.set_xticks(xC)
axC.set_xticklabels(m0.gene, style="italic")
axC.set_xlim(-0.6, len(m0) - 0.4)
axC.set_ylim(0, m0.omega_M0.max() * 1.35)
axC.set_ylabel("codeml M0 ω")
st.style_axis(axC)

# ---- D: yn00 pairwise omega, gene x pair class ---------------------------
letter(L_D, T2 - 8.0, "D")
axD = ax_mm(X_D, T2, W_D, PH2)
rng = np.random.default_rng(gh.JITTER_SEED)
step = 0.26
centres = np.arange(len(GENES), dtype=float)
for gi, g in enumerate(GENES):
    for ci, (pc, col) in enumerate(PAIR_CLASSES):
        v = yn.loc[(yn.gene == g) & (yn.pair == pc), "omega"].to_numpy(float)
        xx = centres[gi] + (ci - 1) * step
        bp = axD.boxplot([v], positions=[xx], widths=step * 0.72, showfliers=False,
                         patch_artist=True, whis=(0, 100), zorder=1)
        bp["boxes"][0].set(facecolor="white", edgecolor=col, linewidth=0.8)
        for key in ("whiskers", "caps"):
            for a in bp[key]:
                a.set(color=col, linewidth=0.8)
        bp["medians"][0].set(color=col, linewidth=1.4)
        axD.scatter(xx + rng.uniform(-step * 0.26, step * 0.26, len(v)), v, s=7,
                    facecolor=col, edgecolor="none", alpha=0.6, zorder=2)
        axD.text(xx, -0.05, f"{len(v)}", ha="center", va="top", fontsize=FS_STAT,
                 color=TEXT, transform=axD.get_xaxis_transform())

axD.set_xticks(centres)
axD.set_xticklabels(GENES, style="italic")
axD.tick_params(axis="x", pad=13)
axD.set_ylabel("yn00 pairwise ω")
axD.set_yscale("log")
axD.set_ylim(yn.omega.min() * 0.6, yn.omega.max() * 6.0)
axD.set_yticks([0.001, 0.01, 0.1, 1.0])
axD.axhline(1.0, color=GREY, lw=0.8, ls="--", zorder=1)
axD.text(1.01, 1.0, "ω = 1", transform=axD.get_yaxis_transform(), ha="left", va="center",
         fontsize=FS_STAT, color=GREY)
st.style_axis(axD)
axD.set_xlim(-0.6, len(GENES) - 0.4)
axD.text(-0.62, -0.05, "pairs", ha="right", va="top", fontsize=FS_STAT, color=TEXT,
         transform=axD.get_xaxis_transform())
for gi, g in enumerate(GENES):
    gh.stat_bracket(axD, centres[gi] - step, centres[gi] + step, yn.omega.max() * 2.2,
                    gh.fmt_p(yn_p[g]), drop=yn.omega.max() * 0.7)

axD.legend(handles=[Patch(facecolor=c, edgecolor="none", label=f"{PAIR_LABEL[pc]} pairs")
                    for pc, c in PAIR_CLASSES],
           loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2, frameon=False,
           handlelength=1.1, columnspacing=1.4, fontsize=FS_BODY)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig3")
