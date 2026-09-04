"""Consolidated Fig 4 - trait-module enrichment and the alkali-tolerance signature.

Consolidation of 2026-09-04 (see /data/data/Upcycling/consolidation_260904/DESIGN.md).
Old panels -> new page:

  old Fig 6   ->  A  permutation forest of every trait subcategory whose observed fold
                     change exceeds 1, ranked by fold change, with the 2,000-iteration
                     bootstrap 95 % CI and a stat column carrying the fold change and the
                     BH-FDR q value of the 10,000-iteration one-sided permutation test
  old Fig 4C  ->  B  MICP-critical trait modules, mean hits per 10^3 CDS, MICP-complete
                     vs rest (old Fig 4 panels A and B, the DRAM heat maps, move to the
                     supplementary consolidation and are NOT drawn here)
  old S11     ->  C  alkaliphile signature: Mrp copy count, Nha copy count, proteome pI
                     median and proteome acidic-pI fraction, MICP-complete vs rest, with
                     the two-sided Mann-Whitney P bound to each block

Provenance: no value is typed in.  Panel A re-derives the fold change of every plotted row
as mean(hits per 10^3 CDS | MICP-complete) / mean(hits per 10^3 CDS | rest) from the raw
counts of Table S2a and asserts it against Fold_change in Table S2c, so the drawn point is
a recomputed value.  Panel B recomputes the per-10^3-CDS normalisation from the same raw
counts, asserts it against the stored normalisation of Table S2b, and asserts both group
means against Table S2c.  Panel C recomputes the group means and Mann-Whitney P from
Table S15a and asserts the Mrp fold difference quoted in the manuscript.

Sources
  Table_S2a_trait_module_counts.csv               raw Bakta keyword hits and CDS_total
  Table_S2b_trait_module_per1kCDS.csv             stored normalisation (asserted against)
  Table_S2c_permutation_statistics.csv            fold change, bootstrap CI, permutation q
  Table_S15a_alkaliphile_signature_per_MAG.csv    Mrp / Nha counts and proteome pI

Colour meanings on this page:
  coral        q < 0.05 in A, the MICP-complete group in B
  light coral  q < 0.10 in A
  grey         n.s. in A, the remaining 105 MAGs in B and C, and the neutral
               fold-change = 1 reference line in A
  blue         a Sphingobacterium MICP-complete MAG (points in C)
  orange       a Pseudomonas_E MICP-complete MAG (points in C)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import (HERO, REST, SIG_05, SIG_10, SIG_NS, SPHINGO, PSEUDO, GREY, TEXT,
                    AXIS, FS_BODY, FS_STAT, HEROES, hero_col)

st.setup()
OUT = HERE / "figures_v2"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

Q_STRONG, Q_WEAK = 0.05, 0.10          # the two BH-FDR tiers drawn in A
RNG = np.random.default_rng(0)         # fixed jitter for C, identical on every rebuild

# the MICP-critical trait modules of panel B (category names, not data)
CRIT = [("Alkaline_Osmo::Mrp_complex", "Mrp Na+/H+ antiporter"),
        ("Alkaline_Osmo::Na_H_antiporter", "nhaA-C antiporter"),
        ("Alkaline_Osmo::oxidative", "Oxidative-stress defence"),
        ("Alkaline_Osmo::compatible_solute", "Compatible-solute synthesis"),
        ("CAZyme_proxy::glycoside_hydrolase", "Glycoside hydrolase"),
        ("CAZyme_proxy::carb_binding", "Carbohydrate-binding module"),
        ("Ammonia_N::urea_transport", "Urea transport, urtABC"),
        ("Ammonia_N::GS_GOGAT", "GS-GOGAT, glnA gltBD"),
        ("Biofilm_EPS::quorum", "Quorum sensing"),
        ("Biofilm_EPS::cellulose", "Cellulose synthesis, bcs")]

# the four alkaliphile-signature features of panel C (category names, not data)
FEATURES = [("Mrp_count", "Mrp copies"),
            ("Nha_count", "Nha copies"),
            ("pI_median", "proteome pI (median)"),
            ("pI_acidic_frac", "acidic-pI fraction")]

# ------------------------------------------------------------------ data: A and B
stats = pd.read_csv(SUPP / "Table_S2c_permutation_statistics.csv")
counts = pd.read_csv(SUPP / "Table_S2a_trait_module_counts.csv").set_index("Sample")
per1k_stored = pd.read_csv(SUPP / "Table_S2b_trait_module_per1kCDS.csv").set_index("Sample")

sub_cols = [c for c in counts.columns if "::" in c]
per1k = counts[sub_cols].div(counts["CDS_total"], axis=0) * 1000.0
assert np.allclose(per1k.values, per1k_stored.loc[per1k.index, sub_cols].values)
is_hero = counts["Hero"].astype(bool)
assert set(counts.index[is_hero]) == set(HEROES)
N_HERO, N_REST = int(is_hero.sum()), int((~is_hero).sum())

# A: recompute the fold change of every row of the stored permutation table
recomputed = {}
for _, r in stats.iterrows():
    key = f"{r.Category}::{r.Subcategory}"
    if key not in per1k.columns:
        raise SystemExit(f"{key} missing from Table S2a")
    h = per1k.loc[is_hero, key].mean()
    q = per1k.loc[~is_hero, key].mean()
    recomputed[key] = (h / q) if q > 0 else np.nan
    # the table rounds the means to 3 decimals; compare on that grid
    assert abs(round(h, 3) - r.Hero_mean_per1kCDS) < 1e-9, (key, h, r.Hero_mean_per1kCDS)
    assert abs(round(q, 3) - r.Rest_mean_per1kCDS) < 1e-9, (key, q, r.Rest_mean_per1kCDS)
    if np.isfinite(recomputed[key]):
        assert abs(recomputed[key] - r.Fold_change) < 5e-3, (key, recomputed[key], r.Fold_change)

df = stats[stats.Fold_change > 1].copy()
df["fc"] = [recomputed[f"{c}::{s}"] for c, s in zip(df.Category, df.Subcategory)]
ci = df.Fold_change_CI95.str.strip("[]").str.split(", ", expand=True).astype(float)
df["lo"], df["hi"] = ci[0].values, ci[1].values
df = df.sort_values("fc", ascending=True).reset_index(drop=True)
N_ROW = len(df)


def tier(q):
    return SIG_05 if q < Q_STRONG else (SIG_10 if q < Q_WEAK else SIG_NS)


df["col"] = [tier(q) for q in df.Permutation_q_BH]

# B: MICP-critical modules, group means asserted against the permutation table
keys = [k for k, _ in CRIT]
crit_lbl = [l for _, l in CRIT]
hero_mean = per1k.loc[is_hero, keys].mean()
rest_mean = per1k.loc[~is_hero, keys].mean()
perm = stats.copy()
perm["key"] = perm.Category + "::" + perm.Subcategory
perm = perm.set_index("key")
for k in keys:
    assert abs(round(hero_mean[k], 3) - perm.loc[k, "Hero_mean_per1kCDS"]) < 1e-9, k
    assert abs(round(rest_mean[k], 3) - perm.loc[k, "Rest_mean_per1kCDS"]) < 1e-9, k

# ------------------------------------------------------------------ data: C
alk = pd.read_csv(SUPP / "Table_S15a_alkaliphile_signature_per_MAG.csv")
alk_hero = alk.group == "MICP_complete"
assert sorted(alk.MAG[alk_hero]) == sorted(HEROES)
n_hero_alk, n_rest_alk = int(alk_hero.sum()), int((~alk_hero).sum())
assert (n_hero_alk, n_rest_alk) == (N_HERO, N_REST)
mrp_fold = alk.loc[alk_hero, "Mrp_count"].mean() / alk.loc[~alk_hero, "Mrp_count"].mean()
assert abs(mrp_fold - 11.7) < 0.05, mrp_fold          # the value quoted in the manuscript

# ------------------------------------------------------------------ layout
ROW_MM = 4.6
L_LAB, W_AX = 46.0, 90.0                 # A: label column, plot width
X_AX = L_LAB + 2.0
X_FC, X_Q = X_AX + W_AX + 5.0, X_AX + W_AX + 21.0
TOP = 12.0
H_AX = N_ROW * ROW_MM

T2 = TOP + H_AX + 26.0                   # top of the second row (B and C)
XB, WB = 34.0, 40.0                      # B: axis box
H_B = len(CRIT) * 5.6
XB_STAT = XB + WB + 4.0

XC = [107.0, 147.0]                      # C: 2 x 2 grid of feature blocks
TC = [T2, T2 + 38.0]
WC, HC = 29.0, 25.0

LEG_X = 95.0                             # left edge of the shared point legend
H = TC[1] + HC + 24.0

fig, ax_mm, text_mm, letter = st.page(H)

# ---- A: permutation forest ------------------------------------------------
ax = ax_mm(X_AX, TOP, W_AX, H_AX)
y = np.arange(N_ROW)
ax.errorbar(df.fc, y, xerr=[df.fc - df.lo, df.hi - df.fc], fmt="none",
            ecolor=GREY, elinewidth=0.6, capsize=1.6, capthick=0.6, zorder=2)
ax.scatter(df.fc, y, s=22, c=df.col, edgecolor="white", linewidth=0.4, zorder=3)
ax.axvline(1.0, ls="--", color=GREY, lw=0.7, zorder=1)

ax.set_xscale("log")
ax.set_xlim(0.8, 60)
# the default log locator emits decade ticks outside the drawn range, whose labels then
# land off the canvas; pin the ticks to the range instead
ax.xaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 50]))
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_ylim(-0.8, N_ROW - 0.2)
ax.set_yticks(y)
ax.set_yticklabels([f"{c.replace('_', ' ')}: {s.replace('_', ' ')}"
                    for c, s in zip(df.Category, df.Subcategory)], fontsize=FS_BODY)
ax.set_xlabel("Fold change, MICP-complete / rest (log scale)")
ax.tick_params(axis="y", length=0)
st.style_axis(ax, left=False, bottom=True)

# stat columns, bound to the rows they sit beside
for x, head in ((X_FC, "FC"), (X_Q, "q")):
    text_mm(x, TOP - 1.5, head, fontsize=FS_STAT, ha="left", va="bottom", color=TEXT)
for i, r in df.iterrows():
    yy = TOP + H_AX - (i + 0.5) * ROW_MM
    text_mm(X_FC, yy, f"{r.fc:.2f}", fontsize=FS_STAT, ha="left", va="center")
    q = r.Permutation_q_BH
    text_mm(X_Q, yy, f"{q:.3f}" if q >= 0.001 else f"{q:.0e}",
            fontsize=FS_STAT, ha="left", va="center")

# only the tiers that actually occur among the drawn rows carry a colour on this page
tiers = [(SIG_05, f"q < {Q_STRONG:g}"), (SIG_10, f"q < {Q_WEAK:g}"), (SIG_NS, "n.s.")]
used = set(df.col)
ax.legend(handles=[Line2D([0], [0], marker="o", ls="", ms=4, mec="white", mew=0.4,
                          mfc=c, label=l) for c, l in tiers if c in used],
          loc="lower right", bbox_to_anchor=(1.0, -0.005), fontsize=FS_BODY,
          handletextpad=0.4, borderpad=0.2, labelspacing=0.25)

text_mm(X_AX, TOP + H_AX + 10.5, f"MICP-complete n = {N_HERO}   rest n = {N_REST}   "
                                 f"error bars, bootstrap 95 % CI",
        fontsize=FS_STAT, color=TEXT)
letter(4, 4, "A")

# ---- B: MICP-critical trait modules --------------------------------------
axB = ax_mm(XB, T2, WB, H_B)
yb = np.arange(len(keys))
bh = 0.38
axB.barh(yb - bh / 2, hero_mean[keys].values, bh, color=HERO)
axB.barh(yb + bh / 2, rest_mean[keys].values, bh, color=REST)
axB.set_yticks(yb)
axB.set_yticklabels(crit_lbl, fontsize=FS_BODY)
axB.invert_yaxis()
axB.set_xscale("log")
# a Unicode superscript keeps the label in the page font; mathtext would emit a
# second font-family into the SVG and break the Arial-first rule
axB.set_xlabel("Gene hits per 10³ CDS (log scale)")
axB.tick_params(axis="y", length=0)
st.style_axis(axB, left=False, bottom=True)
axB.xaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10]))
axB.xaxis.set_minor_locator(NullLocator())
axB.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
axB.set_xlim(min(rest_mean.min(), hero_mean.min()) * 0.5, 100)

# the two group means sit in a stat column rather than beside the bars, where the two
# values of a near-tied row would collide
text_mm(XB_STAT, T2 - 1.5, "MICP-c.", fontsize=FS_STAT, ha="left", va="bottom")
text_mm(XB_STAT + 9.5, T2 - 1.5, "rest", fontsize=FS_STAT, ha="left", va="bottom")
for yi, k in enumerate(keys):
    yy = T2 + (yi + 0.5) * (H_B / len(keys))
    text_mm(XB_STAT, yy, f"{hero_mean[k]:.2f}", fontsize=FS_STAT, ha="left", va="center",
            color=HERO)
    text_mm(XB_STAT + 9.5, yy, f"{rest_mean[k]:.2f}", fontsize=FS_STAT, ha="left",
            va="center", color=REST)
# the legend sits above the axis; inside the panel it would cover the shortest bars
axB.legend(handles=[Patch(color=HERO, label=f"MICP-complete n = {N_HERO}"),
                    Patch(color=REST, label=f"Rest n = {N_REST}")],
           loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=1, fontsize=FS_BODY,
           handlelength=1.2, handletextpad=0.4, borderpad=0.2, labelspacing=0.25)
letter(4, T2 - 6.0, "B")

# ---- C: alkaliphile signature --------------------------------------------
for k, (col, ylab) in enumerate(FEATURES):
    axc = ax_mm(XC[k % 2], TC[k // 2], WC, HC)
    hero_v = alk.loc[alk_hero, col].values
    rest_v = alk.loc[~alk_hero, col].values
    bp = axc.boxplot([hero_v, rest_v], positions=[0, 1], widths=0.55,
                     showfliers=False, patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor="white", edgecolor=AXIS, linewidth=0.7)
    for part in ("whiskers", "caps", "medians"):
        for ln in bp[part]:
            ln.set(color=AXIS, linewidth=0.7)
    for xi, vals, cols in ((0, hero_v, [hero_col(m) for m in alk.MAG[alk_hero]]),
                           (1, rest_v, [REST] * n_rest_alk)):
        jit = RNG.uniform(-0.16, 0.16, len(vals))
        axc.scatter(xi + jit, vals, s=5, c=cols, alpha=0.85, linewidths=0, zorder=3)
    p = mannwhitneyu(hero_v, rest_v, alternative="two-sided").pvalue
    lo = min(float(np.min(hero_v)), float(np.min(rest_v)))
    hi = max(float(np.max(hero_v)), float(np.max(rest_v)))
    span = hi - lo
    axc.set_ylim(lo - 0.08 * span, hi + 0.30 * span)
    axc.set_xlim(-0.6, 1.6)
    axc.set_xticks([0, 1])
    axc.set_xticklabels([f"MICP-complete\nn = {n_hero_alk}", f"rest\nn = {n_rest_alk}"])
    txt = f"P = {p:.1e}" if p < 0.01 else f"P = {p:.2f}"
    axc.text(0.5, hi + 0.20 * span, txt, ha="center", va="center", fontsize=FS_STAT,
             color=TEXT)
    axc.set_ylabel(ylab)
    # the default locator also emits ticks outside the drawn range; their labels are still
    # laid out, off the axis, where they can collide with a neighbouring block
    ylo, yhi = axc.get_ylim()
    axc.set_yticks([t for t in axc.get_yticks() if ylo <= t <= yhi])
    st.style_axis(axc)

letter(97.0, T2 - 6.0, "C")
fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, color=SPHINGO,
                           label="Sphingobacterium"),
                    Line2D([], [], marker="o", ls="", ms=3, color=PSEUDO,
                           label="Pseudomonas_E"),
                    Line2D([], [], marker="o", ls="", ms=3, color=REST,
                           label="other MAG")],
           loc="upper left", ncol=3, frameon=False, fontsize=FS_BODY,
           bbox_to_anchor=(LEG_X / st.PAGE_W_MM,
                           1 - (TC[1] + HC + 12.0) * st.MM / (H * st.MM)),
           handletextpad=0.4, columnspacing=1.6)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig4")
