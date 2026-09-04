"""Consolidated Fig 3 - catalytic and structural conservation of urease.

Consolidation of 2026-09-04 (see /data/data/Upcycling/consolidation_260904/DESIGN.md).
Old panels -> new page:

  old Fig 8A  ->  A  UreC active-site residues, 6 MICP-complete MAGs x 7 canonical
                     catalytic sites (Table S12); cell colour encodes agreement with the
                     S. pasteurii reference
  old Fig 8B  ->  B  ESMFold UreC backbone agreement with PDB 4CEU chain C (Table S22):
                     TM-score (left sub-axis) and all-residue backbone RMSD (right
                     sub-axis).  The two sub-axes are one lettered panel.
  old S17C    ->  C  PAML codeml M0 omega per urease gene (Table S19b)
  old S17D    ->  D  yn00 pairwise omega by gene and pair class, log10 omega axis,
                     drawn from the per-pair yn00 output rather than the stored medians

Provenance: every number is read from a repository source.  Panel B recomputes nothing
(TM-score and RMSD are the stored measurements) but the MAG labels are stripped of the
"_UreC" suffix carried in the table.  Panel D reapplies the pair filter of the analysis
script (`run_dnds_v3.py`: 0 < omega < 99 and dS > 0.01) and asserts that it reproduces
every count, median and Mann-Whitney P stored in Table S19c.

Sources
  Table_S12_UreC_active_site_residues.csv        7 sites x 6 MAGs + expected residue
  Table_S22_ureC_vs_4CEU_tm.csv                  TM-score and backbone RMSD per MAG
  Table_S19b_codeml_M0_summary.csv               4 genes, codeml M0 omega
  Table_S19c_yn00_hero_vs_rest_summary.csv       stored medians / n / MWU P (asserted)
  research/additional/C3_dnds_codon/yn00_pairwise.csv   the pairwise distribution in D

Colour meanings on this page (one meaning per colour):
  green        observed residue matches the reference (A)
  coral        observed residue differs from the reference (A, legend entry only), and in
               D a pair whose two members are both MICP-complete; in C the coral lollipop
               marks the urease genes of the MICP-complete lineages under study
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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import (HERO, HERO_LT, REST, SPHINGO, PSEUDO, GREEN, GREY, TEXT, AXIS,
                    FS_BODY, FS_STAT, HEROES, hero_col)

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
T1, PH1 = 15.0, 32.0                     # row 1: A and B
T2, PH2 = 72.0, 46.0                     # row 2: C and D
H = T2 + PH2 + 15.0
fig, ax_mm, text_mm, letter = st.page(H)

L_A, X_A, W_A = 4.0, 22.0, 56.0
L_B, X_B, W_B = 92.0, 106.0, 24.0
GAP_B = 34.0
L_C, X_C, W_C = 4.0, 20.0, 28.0
L_D, X_D, W_D = 58.0, 70.0, 99.0

# ---- A: active-site residue match matrix ---------------------------------
letter(L_A, 4, "A")
axA = ax_mm(X_A, T1, W_A, PH1)
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
           loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, handlelength=1.1,
           handleheight=0.9, columnspacing=1.0, fontsize=FS_BODY)

# ---- C: codeml M0 omega per urease gene ----------------------------------
letter(L_C, T2 - 6.0, "C")
axC = ax_mm(X_C, T2, W_C, PH1)
yC = np.arange(len(m0), dtype=float)
axC.hlines(yC, 0, m0.omega_M0, color=HERO, lw=1.0, zorder=2)
axC.scatter(m0.omega_M0, yC, s=20, color=HERO, zorder=3)
for yy, v in zip(yC, m0.omega_M0):
    axC.text(v + m0.omega_M0.max() * 0.06, yy, f"{v:.3f}", ha="left", va="center",
             fontsize=FS_STAT, color=TEXT)
axC.set_yticks(yC)
axC.set_yticklabels(m0.gene, style="italic")
axC.invert_yaxis()
axC.set_xlim(0, m0.omega_M0.max() * 1.75)
axC.set_xlabel("codeml M0 ω")
st.style_axis(axC, left=False)
axC.tick_params(left=False)

# ---- D: yn00 pairwise omega, gene x pair class ---------------------------
letter(L_D, T2 - 6.0, "D")
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
