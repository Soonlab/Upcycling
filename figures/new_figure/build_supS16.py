"""Suppl Fig S16 - anti-phage defense systems and CRISPR arrays per MAG (C2).

Panels (old files -> new page):
  Figure_S16_defense.png -> A   DefenseFinder systems per MAG, MICP-complete vs rest
  Figure_S16_crispr.png  -> B   minced CRISPR arrays per MAG, MICP-complete vs rest

The finding is a null: every one of the 111 MAGs carries zero detected defense systems and
zero CRISPR arrays.  Drawn honestly - the y-axis starts at 0 and the point cloud sits on
the zero line - rather than rescaled to manufacture spread.  The count of MAGs above zero
is printed per column so a reader sees the null is a real 0 / n, not a plotting failure.

Sources
  Table_S18a_defense_per_MAG.csv   per-MAG counts (111 rows)
  Table_S18b_defense_hero_vs_rest.csv  stored group means and MWU P (asserted against)

Colour meanings on this page
  coral = MICP-complete group (n = 6)      grey = the remaining 105 MAGs
No other colour is used, so neither carries a second meaning.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import HERO, REST, TEXT, FS_STAT

st.setup()
OUT = HERE / "figures"
SUPP = Path(gh.SUPP)

# ------------------------------------------------------------------ data
d = pd.read_csv(SUPP / "Table_S18a_defense_per_MAG.csv")
stored = pd.read_csv(SUPP / "Table_S18b_defense_hero_vs_rest.csv").set_index("metric")

hero = d.is_hero.astype(bool)
assert int(hero.sum()) == len(st.HEROES), hero.sum()
assert sorted(d.loc[hero, "MAG"]) == sorted(st.HEROES)

METRICS = [("n_defense_systems", "Defense systems per MAG"),
           ("n_crispr_arrays", "CRISPR arrays per MAG")]

panels = []
for col, ylab in METRICS:
    h = d.loc[hero, col].to_numpy(float)
    r = d.loc[~hero, col].to_numpy(float)
    p = mannwhitneyu(h, r, alternative="two-sided").pvalue
    s = stored.loc[col]
    assert abs(h.mean() - s.hero_mean) < 1e-9, (col, h.mean(), s.hero_mean)
    assert abs(r.mean() - s.rest_mean) < 1e-9, (col, r.mean(), s.rest_mean)
    assert abs(p - s.MWU_p) < 1e-9, (col, p, s.MWU_p)
    panels.append(dict(col=col, ylab=ylab, h=h, r=r, p=p,
                       n_pos_h=int((h > 0).sum()), n_pos_r=int((r > 0).sum())))

# ------------------------------------------------------------------ page
H = 59.0
fig, ax_mm, text_mm, letter = st.page(H)

LEFT = [16.0, 104.0]
TOP = 14.0
W, PH = 68.0, 34.0

for (x0, pan, lt) in zip(LEFT, panels, "AB"):
    ax = ax_mm(x0, TOP, W, PH)
    groups, colours = [("MICP-complete", pan["h"]), ("Rest", pan["r"])], [HERO, REST]
    xs, _ = gh.strip_box(ax, groups, colours)
    ax.set_ylabel(pan["ylab"])
    ax.set_ylim(-0.08, 1.0)
    ax.set_yticks([0, 0.5, 1.0])
    gh.group_counts(ax, xs, groups, -0.15)
    # count above zero, bound to each column: shows the null is 0 / n, not a missing series
    for x, n_pos, n_tot in zip(xs, [pan["n_pos_h"], pan["n_pos_r"]],
                               [len(pan["h"]), len(pan["r"])]):
        ax.text(x, 0.10, f"{n_pos} / {n_tot} > 0", ha="center", va="bottom",
                fontsize=FS_STAT, color=TEXT)
    gh.stat_bracket(ax, xs[0], xs[1], 0.80, gh.fmt_p(pan["p"]), drop=0.05)
    letter(x0 - 12.5, TOP - 6.0, lt)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S16")
