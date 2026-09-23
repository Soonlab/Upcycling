"""Suppl Fig S5 (SUBMISSION_v4, 2026-09-23) - growth rate, read-mapping coverage and codon usage.

Derived from build_v2_supS5.py.  The former panels B/C drew a SPAdes contig-coverage proxy
per MAG and split it "by waste source"; the source labels were in fact the binning tools
(Perez-Valera and Elhottova, 2025), and the original authors deposited read-mapping
genome coverage for every MAG (CoverM; Zenodo 10.5281/zenodo.15309541).  The v4 page
therefore shows:

  A  gRodon2 predicted minimum doubling time, MICP-complete vs rest   (unchanged)
  B  read-mapping genome coverage of the MAG in its own sample, MICP-complete vs rest,
     from the deposited metadata (MAG_provenance.tsv)                 (replaces old B, C)
  C  GC3, MICP-complete vs rest                                       (old D)
  D  effective number of codons, MICP-complete vs rest                (old E)

Colour meanings: A blue/orange = the two lineages, grey = other MAG; B-D coral = MICP-complete,
grey = the remaining 105 MAGs.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import HERO, REST, SPHINGO, PSEUDO, TEXT, AXIS, FS_BODY, FS_STAT, hero_col, PROVENANCE

st.setup()
OUT = HERE / "figures_v4"
SUPP = Path(gh.SUPP)
C3 = Path(gh.ADDITIONAL) / "C3_dnds_codon"
RNG = np.random.default_rng(0)


# ================================================================== A: gRodon2 doubling time
gro = pd.read_csv(SUPP / "Table_S16_gRodon_growth_rates_per_MAG.csv")
is_hero = gro.group == "MICP_complete"
gro_h = gro[is_hero]
gro_r = gro[~is_hero]
assert set(gro_h.MAG) <= set(st.HEROES)
p_grow = mannwhitneyu(gro_h.d_hours, gro_r.d_hours, alternative="two-sided").pvalue
assert abs(gro_h.d_hours.median() - 1.06) < 0.005
assert abs(gro_r.d_hours.median() - 1.10) < 0.005
assert abs(p_grow - 0.58) < 0.005
dropped = [m for m in st.HEROES if m not in set(gro_h.MAG)]
assert dropped, "the legend records which MICP-complete MAGs the filter removed"
print(f"  gRodon2: {len(gro_h)} of {len(st.HEROES)} MICP-complete MAGs past the filter "
      f"(removed {', '.join(dropped)}), {len(gro_r)} others")

# ================================================================== B: read-mapping coverage
prov = pd.read_csv(PROVENANCE, sep="\t")
assert len(prov) == 111 and prov.MAG.is_unique
hero = prov.is_candidate.astype(bool)
assert sorted(prov.loc[hero, "MAG"]) == sorted(st.HEROES)
COV = "genome_coverage_x"
cov_h = prov.loc[hero, COV].to_numpy(float)
cov_r = prov.loc[~hero, COV].to_numpy(float)
assert (cov_h > 0).all() and (cov_r > 0).all()
p_cov = mannwhitneyu(cov_h, cov_r, alternative="two-sided").pvalue
print(f"  read-mapping coverage: MICP-complete median {np.median(cov_h):.1f}x  rest median "
      f"{np.median(cov_r):.1f}x  MWU P = {p_cov:.3f}")

# ================================================================== D, E: codon usage
cod = pd.read_csv(SUPP / "Table_S19a_codon_usage_per_MAG.csv")
cod_hero = cod.is_hero.astype(bool)
assert sorted(cod.loc[cod_hero, "MAG"]) == sorted(st.HEROES)
stored_cod = pd.read_csv(C3 / "codon_usage_hero_vs_rest.csv").set_index("metric")

codon_panels = []
for col, ylab in [("GC3_pct", "GC at codon position 3 (%)"),
                  ("ENC", "Effective number of codons")]:
    h = cod.loc[cod_hero, col].to_numpy(float)
    r = cod.loc[~cod_hero, col].to_numpy(float)
    p = mannwhitneyu(h, r, alternative="two-sided").pvalue
    row = stored_cod.loc[col]
    assert abs(h.mean() - row.hero_mean) < 1e-9, (col, h.mean(), row.hero_mean)
    assert abs(r.mean() - row.rest_mean) < 1e-9, (col, r.mean(), row.rest_mean)
    assert abs(p - row.MWU_p) < 1e-15, (col, p, row.MWU_p)
    codon_panels.append(dict(ylab=ylab, h=h, r=r, p=p))

# ================================================================== page
TOP1, PH1 = 15.0, 46.0
TOP2, PH2 = 96.0, 42.0
H = 152.0
fig, ax_mm, text_mm, letter = st.page(H)

# ------------------------------------------------------------------ A: doubling time
axA = ax_mm(18.0, TOP1, 70.0, PH1)
letter(6.0, 8.0, "A")
bp = axA.boxplot([gro_h.d_hours.values, gro_r.d_hours.values], positions=[0, 1], widths=0.5,
                 showfliers=False, patch_artist=True)
for box in bp["boxes"]:
    box.set(facecolor="white", edgecolor=AXIS, linewidth=0.7)
for part in ("whiskers", "caps", "medians"):
    for ln in bp[part]:
        ln.set(color=AXIS, linewidth=0.7)

jit = RNG.uniform(-0.13, 0.13, len(gro_r))
axA.scatter(1 + jit, gro_r.d_hours, s=5, c=REST, alpha=0.85, linewidths=0, zorder=3)
hx = np.linspace(-0.18, 0.18, len(gro_h))
axA.scatter(hx, gro_h.d_hours, s=12, c=[hero_col(m) for m in gro_h.MAG], zorder=4,
            linewidths=0)
axA.set_xticks([0, 1])
axA.set_xticklabels([f"MICP-\ncomplete\nn = {len(gro_h)}", f"rest\nn = {len(gro_r)}"])
axA.set_ylabel("predicted minimum doubling time (h)")
axA.set_xlim(-0.6, 1.6)
axA.set_yscale("log")
axA.set_yticks([0.5, 1, 2, 5, 10, 20])
axA.get_yaxis().set_major_formatter(ScalarFormatter())
lo = float(min(gro_h.d_hours.min(), gro_r.d_hours.min()))
hi = float(max(gro_h.d_hours.max(), gro_r.d_hours.max()))
axA.set_ylim(lo * 0.75, hi * 2.2)
axA.text(0.5, 0.99, f"P = {p_grow:.2f}", ha="center", va="top", fontsize=FS_STAT,
         color=TEXT, transform=axA.transAxes)
st.style_axis(axA)
axA.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, color=SPHINGO,
                           label="Sphingobacterium"),
                    Line2D([], [], marker="o", ls="", ms=3, color=PSEUDO,
                           label="Pseudomonas_E"),
                    Line2D([], [], marker="o", ls="", ms=3, color=REST, label="other MAG")],
           loc="upper left", bbox_to_anchor=(0.0, -0.28), ncol=3, frameon=False,
           fontsize=FS_STAT, handletextpad=0.4, columnspacing=1.4)

# ------------------------------------------------------------------ B: coverage by group
axB = ax_mm(106.0, TOP1, 70.0, PH1)
letter(94.0, 8.0, "B")
grpB = [("MICP-\ncomplete", cov_h), ("Rest", cov_r)]
xsB, _ = gh.strip_box(axB, grpB, [HERO, REST], log=True)
axB.set_ylabel("Read-mapping genome coverage (×)")
axB.set_ylim(1.0, 3000.0)
axB.set_yticks([1, 10, 100, 1000])
gh.group_counts(axB, xsB, grpB, -0.21)
gh.stat_bracket(axB, xsB[0], xsB[1], 1100.0, gh.fmt_p(p_cov), drop=250.0)

# ------------------------------------------------------------------ C, D: codon usage
for x0, xl, pan, lt in zip([18.0, 106.0], [6.0, 94.0], codon_panels, "CD"):
    ax = ax_mm(x0, TOP2, 70.0, PH2)
    letter(xl, TOP2 - 7.0, lt)
    groups = [("MICP-\ncomplete", pan["h"]), ("Rest", pan["r"])]
    xs, _ = gh.strip_box(ax, groups, [HERO, REST])
    ax.set_ylabel(pan["ylab"])
    lo = min(pan["h"].min(), pan["r"].min())
    hi = max(pan["h"].max(), pan["r"].max())
    span = hi - lo
    ax.set_ylim(lo - span * 0.08, hi + span * 0.30)
    gh.group_counts(ax, xs, groups, -0.22)
    gh.stat_bracket(ax, xs[0], xs[1], hi + span * 0.12, gh.fmt_p(pan["p"]),
                    drop=span * 0.05)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S5")
