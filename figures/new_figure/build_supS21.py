"""Suppl Fig S21 - antiSMASH 7 secondary-metabolism profile (C1).

Panels (old file Figure_S21.png -> new page):
  A  total BGC regions per MAG, MICP-complete vs rest
  B  class-level mean BGC regions per MAG, MICP-complete vs rest, with the MWU P per class

Panel B is ranked by MWU P so the two enriched classes lead, and every class is drawn as a
pair of horizontal bars because the value labels of a vertical dot plot collide at this row
count.  Classes are included when either group carries them at an appreciable rate:
hero_mean > 0 or rest_mean > 0.05 regions per MAG.  The cut removes only classes that are
near-absent in both groups and never removes a class with a P below 0.5; the number of
classes dropped is printed at build time and recorded in the journal.

Sources
  Table_S23a_antismash_per_MAG.csv     111 MAGs x per-class region counts
  Table_S23b_antismash_hero_vs_rest.csv stored group means and MWU P (recomputed, asserted)
Every mean and every P drawn here is recomputed from the per-MAG table and asserted against
the stored summary; the `n_regions` row of the summary is the panel A statistic, not a class,
so it is excluded from panel B.

Colour meanings on this page
  coral = MICP-complete group (n = 6)      grey = the remaining 105 MAGs
No other colour is used.
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
from _style import HERO, REST, TEXT, LIGHT, FS_STAT

st.setup()
OUT = HERE / "figures"
SUPP = Path(gh.SUPP)

HERO_MIN = 0.0    # keep a class if the MICP-complete group carries it at all
REST_MIN = 0.05   # or if the rest carry it above this rate (regions per MAG)

# ------------------------------------------------------------------ data
d = pd.read_csv(SUPP / "Table_S23a_antismash_per_MAG.csv")
stored = pd.read_csv(SUPP / "Table_S23b_antismash_hero_vs_rest.csv").set_index("metric")

hero = d.is_hero.astype(bool)
assert sorted(d.loc[hero, "MAG"]) == sorted(st.HEROES)
assert int((~hero).sum()) == len(d) - len(st.HEROES)

# ---- panel A statistic
hA = d.loc[hero, "n_regions"].to_numpy(float)
rA = d.loc[~hero, "n_regions"].to_numpy(float)
pA = mannwhitneyu(hA, rA, alternative="two-sided").pvalue
sA = stored.loc["n_regions"]
assert abs(hA.mean() - sA.hero_mean) < 1e-9 and abs(rA.mean() - sA.rest_mean) < 1e-9
assert abs(pA - sA.MWU_p) < 1e-12, (pA, sA.MWU_p)

# ---- panel B: every class column, recomputed and checked against the stored summary
classes = [c for c in d.columns if c.startswith("BGC_")]
rows = []
for c in classes:
    h = d.loc[hero, c].to_numpy(float)
    r = d.loc[~hero, c].to_numpy(float)
    p = mannwhitneyu(h, r, alternative="two-sided").pvalue
    s = stored.loc[c]
    assert abs(h.mean() - s.hero_mean) < 1e-9, (c, h.mean(), s.hero_mean)
    assert abs(r.mean() - s.rest_mean) < 1e-9, (c, r.mean(), s.rest_mean)
    assert abs(p - s.MWU_p) < 1e-12 * max(1.0, abs(s.MWU_p)) + 1e-15, (c, p, s.MWU_p)
    rows.append(dict(cls=c.replace("BGC_", ""), hero=h.mean(), rest=r.mean(), p=p))
B = pd.DataFrame(rows)

keep = (B.hero > HERO_MIN) | (B.rest > REST_MIN)
dropped = B[~keep]
assert dropped.p.min() > 0.5, dropped.sort_values("p").head()
print(f"  classes: {len(B)} tested | {int(keep.sum())} drawn | {len(dropped)} dropped "
      f"(min P among dropped {dropped.p.min():.3f})")
B = B[keep].sort_values("p").reset_index(drop=True)

# the two enriched classes named in the manuscript must survive and keep their stored P
for cls in ("T3PKS", "RRE-containing"):
    got = B.loc[B.cls == cls]
    assert len(got) == 1, cls
    assert abs(got.p.iloc[0] - stored.loc[f"BGC_{cls}"].MWU_p) < 1e-15, cls

# ------------------------------------------------------------------ page
ROW_MM = 4.2
PH_B = ROW_MM * len(B) + 4.0
TOP_A, PH_A = 15.0, 42.0
TOP_B = TOP_A + PH_A + 16.0
H = TOP_B + PH_B + 12.0

fig, ax_mm, text_mm, letter = st.page(H)

# ---- A
axA = ax_mm(17.0, TOP_A, 46.0, PH_A)
groups, colours = [("MICP-complete", hA), ("Rest", rA)], [HERO, REST]
xs, _ = gh.strip_box(axA, groups, colours)
axA.set_ylabel("BGC regions per MAG")
top = max(hA.max(), rA.max())
axA.set_ylim(-0.5, top * 1.30)
gh.group_counts(axA, xs, groups, -0.12)
gh.stat_bracket(axA, xs[0], xs[1], top * 1.08, gh.fmt_p(pA), drop=top * 0.05)
letter(5.0, TOP_A - 6.0, "A")

# ---- B
axB = ax_mm(40.0, TOP_B, 88.0, PH_B)
y = gh.paired_hbars(axB, B.cls.tolist(), B.hero.to_numpy(), B.rest.to_numpy())
axB.set_xlabel("Mean BGC regions per MAG")
axB.set_xlim(0, max(B.hero.max(), B.rest.max()) * 1.12)
axB.legend(loc="lower left", bbox_to_anchor=(0.0, 1.005), frameon=False, ncol=2,
           handlelength=1.1, borderaxespad=0.0, columnspacing=1.6)

# stat column, one P per row, bound to the row it labels
XP = 1.03
axB.text(XP, -0.9, "MWU P", transform=axB.get_yaxis_transform(), ha="left", va="center",
         fontsize=FS_STAT, color=TEXT)
for yy, p in zip(y, B.p):
    txt = gh.fmt_p(p).replace("P = ", "")
    axB.text(XP, yy, txt, transform=axB.get_yaxis_transform(), ha="left", va="center",
             fontsize=FS_STAT, color=TEXT)

letter(28.0, TOP_B - 6.0, "B")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S21")
