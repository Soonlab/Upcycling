"""Suppl Fig S5 (consolidated) - growth rate, in-situ abundance proxy and codon usage.

Consolidation of three old supplementary pages onto one page, per
consolidation_260904/DESIGN.md.

Panels (old page -> new panel)
  A  old S12   predicted minimum doubling time by group (gRodon2), log y axis
  B  old S19A  length-weighted mean contig coverage per MAG, MICP-complete vs rest
  C  old S19B  the same coverage by waste source, the six MICP-complete MAGs marked inside
               their own source column
  D  old S17A  genome-wide GC at codon position 3, MICP-complete vs rest
  E  old S17B  effective number of codons, same comparison

The codeml M0 and yn00 panels of old S17 (C and D) are NOT on this page; they move to a
main figure.

Two axis choices are deliberate and are kept from the old pages.  Panel A is logarithmic
because the rest group spans 0.4-17.4 h and on a linear axis its upper outliers compress
every MICP-complete point onto the baseline; the group n stays in the tick label.  Panels B
and C are logarithmic because coverage spans two orders of magnitude.

SOURCE LABELLING (the old S19 defect has since been corrected upstream)
When old S19 was built, the `source` column of Table S21a attached "sheep" to the M-prefixed
MAGs and "swine" to the S-prefixed MAGs - the opposite of the sample-code key used throughout
the manuscript (C cattle, M swine, S sheep, V poultry) and of the independent per-MAG source
column in Table S9a - and Table S21b inherited the swap.  Both tables were rewritten on
2026-09-04 and now agree with the identifier rule, so old S19 can no longer satisfy its own
swap-aware assertions.  Panel C still derives the waste source from the MAG identifier, which
can be checked against the identifier itself, and still asserts it against Table S9a.  The
mapping onto the stored per-source table is resolved from the data rather than assumed: the
stored row whose count, mean, median and standard deviation reproduce a recomputed group is
found, and the build asserts that the row found carries that group's own name.  The drawn
panel is identical to old S19 panel B, which always used the identifier rule.

Sources
  Table_S16_gRodon_growth_rates_per_MAG.csv        A (85 MAGs past the ribosomal anchor filter)
  Table_S21a_abundance_proxy_per_MAG.csv           B, C
  Table_S21b_abundance_proxy_per_source.csv        C (recomputed and asserted against)
  Table_S9a_PCoA_coordinates.csv                   C (independent waste source, asserted against)
  Table_S19a_codon_usage_per_MAG.csv               D, E
  research/additional/C3_dnds_codon/codon_usage_hero_vs_rest.csv  D, E (asserted against)

Colour meanings on this page
  A:    blue = Sphingobacterium MICP-complete MAG, orange = Pseudomonas_E MICP-complete MAG,
        grey = one of the remaining MAGs
  B, D, E: coral = MICP-complete group, grey = the remaining 105 MAGs
  C:    four waste-source colours (cattle brown, swine pink, sheep green, poultry purple) for
        the distributions and a black outline for the six MICP-complete MAGs.  Coral is not
        reused in C, so the MICP-complete MAGs are marked by outline rather than by fill,
        their fill stating which source they came from.
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
from _style import HERO, REST, SPHINGO, PSEUDO, TEXT, AXIS, FS_BODY, FS_STAT, hero_col

st.setup()
OUT = HERE / "figures_v2"
SUPP = Path(gh.SUPP)
C3 = Path(gh.ADDITIONAL) / "C3_dnds_codon"
RNG = np.random.default_rng(0)

COV = "length_weighted_cov"
SOURCE_ORDER = ["cattle", "swine", "sheep", "poultry"]
PREFIX_SOURCE = {"C": "cattle", "M": "swine", "S": "sheep", "V": "poultry"}

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

# ================================================================== B, C: coverage proxy
d = pd.read_csv(SUPP / "Table_S21a_abundance_proxy_per_MAG.csv")
per_source = pd.read_csv(SUPP / "Table_S21b_abundance_proxy_per_source.csv").set_index("source")

hero = d.is_hero.astype(bool)
assert sorted(d.loc[hero, "MAG"]) == sorted(st.HEROES)
assert set(d.source) == set(SOURCE_ORDER), sorted(set(d.source))

# waste source taken from the MAG identifier, checked against the independent Table S9a
d["waste_source"] = d.MAG.str[0].map(PREFIX_SOURCE)
assert d.waste_source.notna().all()
s9 = pd.read_csv(SUPP / "Table_S9a_PCoA_coordinates.csv")
s9_src = dict(zip(s9.MAG, s9.Source.str.lower()))
assert all(s9_src[m] == w for m, w in zip(d.MAG, d.waste_source)), "S9a disagrees with the prefix rule"
n_swapped = int((d.source != d.waste_source).sum())
print(f"  Table S21a source column disagrees with the MAG identifier for {n_swapped} of {len(d)} MAGs "
      f"(the sheep/swine swap old S19 corrected is gone from the table)")

cov_h = d.loc[hero, COV].to_numpy(float)
cov_r = d.loc[~hero, COV].to_numpy(float)
p_cov = mannwhitneyu(cov_h, cov_r, alternative="two-sided").pvalue

# the stored per-source table is recomputed from the per-MAG rows; the stored row that
# belongs to each recomputed group is identified by its own statistics and then asserted to
# carry that group's name, so a relabelling in either direction fails the build
assert set(per_source.index) == set(SOURCE_ORDER), sorted(per_source.index)
for s in SOURCE_ORDER:
    v = d.loc[d.waste_source == s, COV]
    match = [name for name, row in per_source.iterrows()
             if len(v) == row["count"] and abs(v.mean() - row["mean"]) < 1e-6
             and abs(v.median() - row["median"]) < 1e-6
             and abs(v.std(ddof=1) - row["std"]) < 1e-6]
    assert len(match) == 1, (s, match)
    assert match[0] == s, f"Table S21b labels the {s} group as {match[0]}"
print(f"  group means: MICP-complete {cov_h.mean():.1f}x  rest {cov_r.mean():.1f}x  "
      f"MWU P = {p_cov:.4f}")

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
TOP2, PH2 = 90.0, 42.0
H = 146.0
fig, ax_mm, text_mm, letter = st.page(H)

# ------------------------------------------------------------------ A: doubling time
axA = ax_mm(18.0, TOP1, 36.0, PH1)
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
           loc="upper left", bbox_to_anchor=(-0.34, -0.28), ncol=2, frameon=False,
           fontsize=FS_STAT, handletextpad=0.4, columnspacing=1.4)

# ------------------------------------------------------------------ B: coverage by group
axB = ax_mm(70.0, TOP1, 32.0, PH1)
letter(58.0, 8.0, "B")
grpB = [("MICP-\ncomplete", cov_h), ("Rest", cov_r)]
xsB, _ = gh.strip_box(axB, grpB, [HERO, REST], log=True)
axB.set_ylabel("Length-weighted contig coverage (×)")
axB.set_ylim(1.0, 3000.0)
axB.set_yticks([1, 10, 100, 1000])
gh.group_counts(axB, xsB, grpB, -0.21)
gh.stat_bracket(axB, xsB[0], xsB[1], 1100.0, gh.fmt_p(p_cov), drop=250.0)

# ------------------------------------------------------------------ C: coverage by source
axC = ax_mm(118.0, TOP1, 58.0, PH1)
letter(106.0, 8.0, "C")
grpC = [(s, d.loc[d.waste_source == s, COV].to_numpy(float)) for s in SOURCE_ORDER]
xsC, jitC = gh.strip_box(axC, grpC, [st.SOURCE[s] for s in SOURCE_ORDER], log=True)
for si, s in enumerate(SOURCE_ORDER):
    in_src = (d.waste_source == s).to_numpy()
    hero_in_src = hero.to_numpy()[in_src]
    if not hero_in_src.any():
        continue
    axC.scatter(jitC[si][hero_in_src], d.loc[in_src, COV].to_numpy()[hero_in_src],
                s=26, facecolor="none", edgecolor=TEXT, linewidth=0.9, zorder=4)
axC.set_ylabel("Length-weighted contig coverage (×)")
axC.set_ylim(1.0, 3000.0)
axC.set_yticks([1, 10, 100, 1000])
gh.group_counts(axC, xsC, grpC, -0.16)
axC.legend(handles=[Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
                           markeredgecolor=TEXT, markersize=5,
                           label="MICP-complete MAG")],
           loc="upper left", bbox_to_anchor=(0.0, -0.28), frameon=False,
           fontsize=FS_STAT, handletextpad=0.4)

# ------------------------------------------------------------------ D, E: codon usage
for x0, xl, pan, lt in zip([20.0, 108.0], [6.0, 94.0], codon_panels, "DE"):
    ax = ax_mm(x0, TOP2, 56.0, PH2)
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
