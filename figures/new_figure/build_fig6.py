"""Main Fig 6 - permutation-tested trait-module enrichment, MICP-complete vs rest.

One composed page, 180 mm wide.  Forest of every trait subcategory whose observed fold
change exceeds 1, ranked by fold change, with the 2,000-iteration bootstrap 95 % CI and a
stat column carrying the fold change and the BH-FDR q value of the 10,000-iteration
one-sided permutation test.

Source
  Table_S2c_permutation_statistics.csv   fold change, bootstrap CI, permutation q (BH)
  Table_S2a_trait_module_counts.csv      raw hits and CDS_total, used to RE-DERIVE the
                                         fold change that the table stores
The fold change of every plotted row is recomputed as
mean(hits per 10^3 CDS | MICP-complete) / mean(hits per 10^3 CDS | rest) from Table S2a and
asserted against Fold_change in Table S2c, so the drawn point is a recomputed value.

Colour meanings on this page:
  coral       q < 0.05   (enriched in the MICP-complete group)
  light coral q < 0.10
  grey        n.s.  and the neutral fold-change = 1 reference line
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import SIG_05, SIG_10, SIG_NS, GREY, TEXT, FS_BODY, FS_STAT

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

Q_STRONG, Q_WEAK = 0.05, 0.10          # the two BH-FDR tiers drawn

# ------------------------------------------------------------------ data
stats = pd.read_csv(SUPP / "Table_S2c_permutation_statistics.csv")
counts = pd.read_csv(SUPP / "Table_S2a_trait_module_counts.csv").set_index("Sample")

sub_cols = [c for c in counts.columns if "::" in c]
per1k = counts[sub_cols].div(counts["CDS_total"], axis=0) * 1000.0
is_hero = counts["Hero"].astype(bool)
N_HERO, N_REST = int(is_hero.sum()), int((~is_hero).sum())

# recompute the fold change of every row of the stored table
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

# ------------------------------------------------------------------ layout
ROW_MM = 5.0
L_LAB, W_AX = 46.0, 92.0                 # label column, plot width
X_AX = L_LAB + 2.0
X_FC, X_Q = X_AX + W_AX + 5.0, X_AX + W_AX + 21.0
TOP = 12.0
H_AX = N_ROW * ROW_MM
H = TOP + H_AX + 19.0

fig, ax_mm, text_mm, letter = st.page(H)
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
handles = [Line2D([0], [0], marker="o", ls="", ms=4, mec="white", mew=0.4, mfc=c, label=l)
           for c, l in tiers if c in used]
ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, -0.005),
          fontsize=FS_BODY, handletextpad=0.4, borderpad=0.2, labelspacing=0.25)

text_mm(X_AX, H - 6.0, f"MICP-complete n = {N_HERO}   rest n = {N_REST}   "
                       f"error bars, bootstrap 95 % CI", fontsize=FS_STAT, color=TEXT)
letter(4, 4, "A")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig6")
