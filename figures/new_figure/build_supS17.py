"""Suppl Fig S17 - urease gene family selection regime and genome-wide codon usage (C3).

Panels (old files -> new page):
  Figure_S17a_codon_GC3.png       A  genome-wide GC3, MICP-complete vs rest
  Figure_S17b_codon_ENC.png       B  genome-wide effective number of codons, same comparison
  Figure_S17c_codeml_omega.png    C  PAML codeml M0 omega per urease gene
  Figure_S17d_yn00_hero_vs_rest.png D  yn00 pairwise omega by gene and pair class

Panel D is drawn from the per-pair yn00 output rather than from the stored medians, so the
distributions behind the ureG result are visible.  Its omega axis is log10 because the pairwise
estimates span three orders of magnitude: on a linear axis the whiskers of the ureA classes
leave the panel and the ureB / ureC distributions collapse onto the baseline.  The axis limits
are taken from the data, not typed in.  The neutrality line at omega = 1
is drawn for reference.  The pair filter is the one the analysis
script used (`run_dnds_v3.py`: 0 < omega < 99 and dS > 0.01); applying it reproduces every
count, median and P in Table S19c exactly, and the build asserts that.

Sources
  Table_S19a_codon_usage_per_MAG.csv      111 MAGs, GC3 % and ENC
  Table_S19b_codeml_M0_summary.csv        4 genes, codeml M0 omega
  Table_S19c_yn00_hero_vs_rest_summary.csv stored medians / n / MWU P (asserted against)
  research/additional/C3_dnds_codon/yn00_pairwise.csv  595 pairs, the distribution drawn in D

Colour meanings on this page
  A, B: coral = MICP-complete group, grey = the remaining 105 MAGs
  C:    coral, because codeml M0 pools the whole 18-MAG subset and the value describes the
        urease genes of the MICP-complete lineages under study, not a group contrast
  D:    coral = both members of the pair are MICP-complete, grey = neither is,
        light coral = one of each.  Coral therefore means "MICP-complete" on every panel and
        grey means "not MICP-complete" on every panel; no colour carries a second meaning.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import HERO, HERO_LT, REST, TEXT, GREY, FS_STAT

st.setup()
OUT = HERE / "figures"
SUPP = Path(gh.SUPP)
C3 = Path(gh.ADDITIONAL) / "C3_dnds_codon"

# the pair filter of run_dnds_v3.py, reproduced so panel D matches Table S19c
OMEGA_MAX, OMEGA_MIN, DS_MIN = 99.0, 0.0, 0.01
PAIR_CLASSES = [("hero-hero", HERO), ("hero-rest", HERO_LT), ("rest-rest", REST)]

# ------------------------------------------------------------------ A, B: codon usage
cod = pd.read_csv(SUPP / "Table_S19a_codon_usage_per_MAG.csv")
hero = cod.is_hero.astype(bool)
assert sorted(cod.loc[hero, "MAG"]) == sorted(st.HEROES)
stored_cod = pd.read_csv(C3 / "codon_usage_hero_vs_rest.csv").set_index("metric")

codon_panels = []
for col, ylab in [("GC3_pct", "GC at codon position 3 (%)"),
                  ("ENC", "Effective number of codons")]:
    h = cod.loc[hero, col].to_numpy(float)
    r = cod.loc[~hero, col].to_numpy(float)
    p = mannwhitneyu(h, r, alternative="two-sided").pvalue
    s = stored_cod.loc[col]
    assert abs(h.mean() - s.hero_mean) < 1e-9, (col, h.mean(), s.hero_mean)
    assert abs(r.mean() - s.rest_mean) < 1e-9, (col, r.mean(), s.rest_mean)
    assert abs(p - s.MWU_p) < 1e-15, (col, p, s.MWU_p)
    codon_panels.append(dict(ylab=ylab, h=h, r=r, p=p))

# ------------------------------------------------------------------ C: codeml M0
m0 = pd.read_csv(SUPP / "Table_S19b_codeml_M0_summary.csv")
m0 = m0.sort_values("omega_M0").reset_index(drop=True)
assert (m0.omega_M0 < 1).all(), "codeml M0 omega above 1 would not be purifying selection"

# ------------------------------------------------------------------ D: yn00 pairwise
yn = pd.read_csv(C3 / "yn00_pairwise.csv")
n_all = len(yn)
yn = yn[(yn.omega < OMEGA_MAX) & (yn.omega > OMEGA_MIN) & (yn.dS > DS_MIN)]
heroes = set(st.HEROES)
in_h = yn.a.isin(heroes)
in_h_b = yn.b.isin(heroes)
yn = yn.assign(pair=np.where(in_h & in_h_b, "hero-hero",
                             np.where(~in_h & ~in_h_b, "rest-rest", "hero-rest")))
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
TOP1, PH1 = 15.0, 44.0
TOP2 = TOP1 + PH1 + 20.0
PH2 = 52.0
H = TOP2 + PH2 + 14.0
fig, ax_mm, text_mm, letter = st.page(H)

# ---- A, B
for x0, pan, lt in zip([17.0, 72.0], codon_panels, "AB"):
    ax = ax_mm(x0, TOP1, 40.0, PH1)
    groups = [("MICP-\ncomplete", pan["h"]), ("Rest", pan["r"])]
    xs, _ = gh.strip_box(ax, groups, [HERO, REST])
    ax.set_ylabel(pan["ylab"])
    lo = min(pan["h"].min(), pan["r"].min())
    hi = max(pan["h"].max(), pan["r"].max())
    span = hi - lo
    ax.set_ylim(lo - span * 0.08, hi + span * 0.30)
    gh.group_counts(ax, xs, groups, -0.16)
    gh.stat_bracket(ax, xs[0], xs[1], hi + span * 0.12, gh.fmt_p(pan["p"]),
                    drop=span * 0.05)
    letter(x0 - 12.0, TOP1 - 6.0, lt)

# ---- C: codeml M0 omega, horizontal lollipops
axC = ax_mm(140.0, TOP1, 30.0, PH1)
yC = np.arange(len(m0), dtype=float)
axC.hlines(yC, 0, m0.omega_M0, color=HERO, lw=1.0, zorder=2)
axC.scatter(m0.omega_M0, yC, s=20, color=HERO, zorder=3)
for yy, v in zip(yC, m0.omega_M0):
    axC.text(v + 0.004, yy, f"{v:.3f}", ha="left", va="center", fontsize=FS_STAT, color=TEXT)
axC.set_yticks(yC)
axC.set_yticklabels(m0.gene, style="italic")
axC.invert_yaxis()
axC.set_xlim(0, m0.omega_M0.max() * 1.55)
axC.set_xlabel("codeml M0 ω")
st.style_axis(axC, left=False)
axC.tick_params(left=False)
letter(129.0, TOP1 - 6.0, "C")

# ---- D: yn00 pairwise omega, gene x pair class
axD = ax_mm(17.0, TOP2, 149.0, PH2)
rng = np.random.default_rng(gh.JITTER_SEED)
step = 0.26
centres = np.arange(len(GENES), dtype=float)
for gi, g in enumerate(GENES):
    for ci, (pc, col) in enumerate(PAIR_CLASSES):
        v = yn.loc[(yn.gene == g) & (yn.pair == pc), "omega"].to_numpy(float)
        x = centres[gi] + (ci - 1) * step
        bp = axD.boxplot([v], positions=[x], widths=step * 0.72, showfliers=False,
                         patch_artist=True, whis=(0, 100), zorder=1)
        bp["boxes"][0].set(facecolor="white", edgecolor=col, linewidth=0.8)
        for key in ("whiskers", "caps"):
            for a in bp[key]:
                a.set(color=col, linewidth=0.8)
        bp["medians"][0].set(color=col, linewidth=1.4)
        axD.scatter(x + rng.uniform(-step * 0.26, step * 0.26, len(v)), v, s=7,
                    facecolor=col, edgecolor="none", alpha=0.6, zorder=2)
        axD.text(x, -0.05, f"{len(v)}", ha="center", va="top", fontsize=FS_STAT,
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
axD.text(-0.6, -0.05, "pairs", ha="right", va="top", fontsize=FS_STAT, color=TEXT,
         transform=axD.get_xaxis_transform())
for gi, g in enumerate(GENES):
    gh.stat_bracket(axD, centres[gi] - step, centres[gi] + step, yn.omega.max() * 2.2,
                    gh.fmt_p(yn_p[g]), drop=yn.omega.max() * 0.7)

PAIR_LABEL = {"hero-hero": "MICP-complete × MICP-complete",
              "hero-rest": "MICP-complete × rest", "rest-rest": "rest × rest"}
handles = [Patch(facecolor=c, edgecolor="none", label=f"{PAIR_LABEL[pc]} pairs")
           for pc, c in PAIR_CLASSES]
axD.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=3,
           frameon=False, handlelength=1.1, columnspacing=1.8)
letter(5.0, TOP2 - 6.0, "D")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S17")
