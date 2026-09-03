"""Suppl Fig S4 - CAZyme profile: keyword proxy, dbCAN classes, dbCAN families.

Panels
  A  keyword CAZy proxy, genus-aggregated  (Table_S2b, CAZyme_proxy:: block)
  B  dbCAN class abundance, MICP-complete vs rest, with the permutation statistics
     (Table_S6b per-MAG hits per 10^3 CDS; Table_S6d fold change / q)
  C  dbCAN family abundance, MICP-complete vs rest, top families by MICP-complete mean

Provenance note for panel C.  The shipped Table_S6c (family counts) contains only the
105 non-MICP-complete MAGs, which is why the superseded Figure_S4c rendered as an empty
axis - its hero mean was taken over an empty set.  The six MICP-complete MAGs were
annotated by direct hmmsearch against dbCAN (research/revision/dbcan_direct/<MAG>.tbl,
best HMM per protein), the same route that produced the published hero CLASS counts.
This script re-parses those tables with the recipe of research/revision/04b_dbcan_final.py
and asserts that the per-class sums reproduce Table_S6a exactly for all six MAGs, so the
family panel rests on the same annotation as the published class panel.  The hero/rest
method asymmetry (hmmsearch vs DRAM cazy_best_hit) is inherited from that published
analysis and is unchanged here.

Colour meanings on this page:
  white -> green   mean hits per 10^3 CDS, genus x subcategory (panel A only)
  coral            the MICP-complete group (n = 6): bars, dots and genus labels
  grey             the remaining 105 MAGs: bars, dots, and the lollipop connector
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _supp_traits as tr

st.setup()
OUT = HERE / "figures"
SUPP = tr.SUPP
BASE = Path("/data/data/Upcycling")
TBL_DIR = BASE / "research/revision/dbcan_direct"
CDS_SRC = BASE / "research/extra/gene_category_counts.csv"
CLASSES = ["GH", "GT", "PL", "CE", "AA", "CBM"]
N_FAMILIES = 14

# ------------------------------------------------------------------ panel A
df = tr.load_traits()
cols_a = [c for c in df.columns if c.startswith("CAZyme_proxy::")]
table_a, n_mags = tr.genus_aggregate(df, cols_a)
for g in table_a.index:
    assert np.allclose(table_a.loc[g].values, df[df["Genus"] == g][cols_a].mean().values)

# ------------------------------------------------------------------ panel B
per1k = pd.read_csv(SUPP / "Table_S6b_dbCAN_class_per1kCDS.csv").set_index("fasta")
stats_b = pd.read_csv(SUPP / "Table_S6d_dbCAN_hero_vs_rest.csv").set_index("Class")
is_hero = per1k.index.isin(st.HEROES)
hero_mean = per1k[is_hero][CLASSES].mean()
rest_mean = per1k[~is_hero][CLASSES].mean()
assert is_hero.sum() == len(st.HEROES)
for c in CLASSES:                      # published means and fold changes must reproduce
    assert abs(hero_mean[c] - stats_b.loc[c, "Hero_mean_per1kCDS"]) < 5e-3
    assert abs(rest_mean[c] - stats_b.loc[c, "Rest_mean_per1kCDS"]) < 5e-3
    assert abs(hero_mean[c] / rest_mean[c] - stats_b.loc[c, "Fold_change"]) < 5e-3

# ------------------------------------------------------------------ panel C
def best_hmm_per_protein(path):
    best = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        target, hmm, evalue = f[0], f[2], float(f[4])
        if target not in best or evalue < best[target][1]:
            best[target] = (hmm, evalue)
    return best


FAMILY_RE = re.compile(r"(GH|GT|PL|CE|AA|CBM)(\d+.*)?$")
# dbCAN HMM names carry a subfamily suffix (GT2_Glycos_transf_2) while the DRAM-derived
# family table for the other 105 MAGs carries the bare family (GT2).  Both sides are
# reduced to the bare family so the two groups are counted in the same unit.
FAMILY_KEY = re.compile(r"^(GH|GT|PL|CE|AA|CBM)\d+")


def family_key(name):
    m = FAMILY_KEY.match(name)
    return m.group(0) if m else None


hero_fam = defaultdict(lambda: defaultdict(int))
hero_cls = defaultdict(lambda: defaultdict(int))
for mag in st.HEROES:
    for _, (hmm, _ev) in best_hmm_per_protein(TBL_DIR / f"{mag}.tbl").items():
        name = hmm.replace(".hmm", "")
        m = FAMILY_RE.match(name)
        if m:
            key = family_key(name)
            if key:
                hero_fam[mag][key] += 1
            hero_cls[mag][m.group(1)] += 1

class_counts = pd.read_csv(SUPP / "Table_S6a_dbCAN_class_counts.csv").set_index("fasta")
hero_fam_df = pd.DataFrame(hero_fam).T.fillna(0).astype(int)
for mag in st.HEROES:                  # the re-parse must reproduce the published counts
    for c in CLASSES:
        assert hero_cls[mag].get(c, 0) == class_counts.loc[mag, c], (mag, c)
    assert hero_fam_df.loc[mag].sum() == class_counts.loc[mag, "Total"]
assert all(family_key(f) == f for f in hero_fam_df.columns)

rest_fam_df = pd.read_csv(SUPP / "Table_S6c_dbCAN_family_counts.csv").set_index("fasta")
assert not set(rest_fam_df.index) & set(st.HEROES)
# the DRAM-derived columns are already bare families; keep only those, so the two groups
# are compared family-for-family
rest_fam_df = rest_fam_df[[c for c in rest_fam_df.columns if family_key(c) == c]]
cds = pd.read_csv(CDS_SRC).set_index("Sample")["CDS_total"]
hero_n1k = hero_fam_df.div(cds.loc[hero_fam_df.index], axis=0) * 1000.0
rest_n1k = rest_fam_df.div(cds.loc[rest_fam_df.index], axis=0) * 1000.0
# per-class sums of the normalised family table must reproduce the panel-B class means
for c in CLASSES:
    fams = [f for f in hero_n1k.columns if FAMILY_KEY.match(f).group(1) == c]
    assert abs(hero_n1k[fams].sum(axis=1).mean() - hero_mean[c]) < 5e-3

fam_all = sorted(set(hero_n1k.columns) | set(rest_n1k.columns))
hero_fam_mean = hero_n1k.reindex(columns=fam_all, fill_value=0).mean()
rest_fam_mean = rest_n1k.reindex(columns=fam_all, fill_value=0).mean()
n_hero_with = (hero_n1k.reindex(columns=fam_all, fill_value=0) > 0).sum()
top = hero_fam_mean.sort_values(ascending=False).head(N_FAMILIES).index.tolist()

# ------------------------------------------------------------------ page
ROW_A = 6.0
A_TOP, A_LEFT, A_W = 34.0, 36.0, 108.0
A_H = ROW_A * table_a.shape[0]
B_TOP, B_LEFT, B_W, B_H = A_TOP + A_H + 26.0, 22.0, 96.0, 46.0
C_TOP, C_LEFT, C_W = B_TOP + B_H + 24.0, 44.0, 78.0
C_H = 4.6 * N_FAMILIES
H = C_TOP + C_H + 12.0

fig, ax_mm, text_mm, letter = st.page(H)

axA = ax_mm(A_LEFT, A_TOP, A_W, A_H)
caxA = ax_mm(A_LEFT + A_W + 6.0, A_TOP, 3.0, min(30.0, A_H))
letter(4, 4, "A")
tr.draw_genus_heatmap(fig, axA, caxA, table_a, n_mags)
text_mm(A_LEFT + A_W + 15.5, A_TOP + min(30.0, A_H) / 2,
        "log₁₀(mean hits per 10³ CDS + 1)", rotation=90, ha="center",
        va="center", fontsize=st.FS_BODY)

# --- B: class means, MICP-complete vs rest
axB = ax_mm(B_LEFT, B_TOP, B_W, B_H)
letter(4, B_TOP - 8, "B")
y = np.arange(len(CLASSES))
bh = 0.36
axB.barh(y - bh / 2, hero_mean[CLASSES].values, height=bh, color=st.HERO, zorder=3)
axB.barh(y + bh / 2, rest_mean[CLASSES].values, height=bh, color=st.REST, zorder=3)
for i, c in enumerate(CLASSES):
    axB.text(hero_mean[c] + 0.4, i - bh / 2, f"{hero_mean[c]:.2f}", va="center",
             ha="left", fontsize=st.FS_STAT, color=st.HERO)
    axB.text(rest_mean[c] + 0.4, i + bh / 2, f"{rest_mean[c]:.2f}", va="center",
             ha="left", fontsize=st.FS_STAT, color=st.TEXT)
axB.set_yticks(y)
axB.set_yticklabels(CLASSES, fontsize=st.FS_BODY)
axB.invert_yaxis()
axB.set_xlabel("mean hits per 10³ CDS")
axB.set_xlim(0, float(hero_mean.max()) * 1.18)
st.style_axis(axB, left=False)
axB.tick_params(axis="y", length=0)
axB.legend(handles=[Line2D([], [], color=st.HERO, lw=4,
                           label=f"MICP-complete (n = {int(is_hero.sum())})"),
                    Line2D([], [], color=st.REST, lw=4,
                           label=f"rest (n = {int((~is_hero).sum())})")],
           loc="lower right", fontsize=st.FS_BODY, handlelength=1.2)

# stat columns bound to the class rows
sx = [B_LEFT + B_W + 12.0, B_LEFT + B_W + 30.0, B_LEFT + B_W + 46.0]
head_y = B_TOP - 3.0
for x, head in zip(sx, ["fold\nchange", "permutation\nq (BH)", "MWU\nP"]):
    text_mm(x, head_y, head, ha="center", va="bottom", fontsize=st.FS_STAT,
            color=st.AXIS, linespacing=1.25)
for i, c in enumerate(CLASSES):
    yy = B_TOP + (i + 0.5) * (B_H / len(CLASSES))
    q = stats_b.loc[c, "Permutation_q_BH"]
    text_mm(sx[0], yy, f"{stats_b.loc[c, 'Fold_change']:.2f}", ha="center",
            va="center", fontsize=st.FS_STAT)
    text_mm(sx[1], yy, f"{q:.4f}", ha="center", va="center", fontsize=st.FS_STAT,
            color=st.HERO if q < 0.05 else st.TEXT)
    text_mm(sx[2], yy, f"{stats_b.loc[c, 'MannWhitney_greater_p']:.4f}",
            ha="center", va="center", fontsize=st.FS_STAT)

# --- C: family means, MICP-complete vs rest
axC = ax_mm(C_LEFT, C_TOP, C_W, C_H)
letter(4, C_TOP - 8, "C")
yc = np.arange(len(top))
for i, f in enumerate(top):
    axC.plot([rest_fam_mean[f], hero_fam_mean[f]], [i, i], color=st.LIGHT, lw=2.4,
             solid_capstyle="round", zorder=1)
axC.scatter(rest_fam_mean[top].values, yc, s=16, color=st.REST, zorder=3)
axC.scatter(hero_fam_mean[top].values, yc, s=16, color=st.HERO, zorder=4)
axC.set_yticks(yc)
axC.set_yticklabels(top, fontsize=st.FS_BODY)
axC.invert_yaxis()
axC.set_xlabel("mean hits per 10³ CDS")
axC.set_xlim(-0.05, float(hero_fam_mean[top].max()) * 1.12)
st.style_axis(axC, left=False)
axC.tick_params(axis="y", length=0)
axC.legend(handles=[Line2D([], [], marker="o", ls="", color=st.HERO,
                           label="MICP-complete mean"),
                    Line2D([], [], marker="o", ls="", color=st.REST,
                           label="rest mean")],
           loc="lower right", fontsize=st.FS_BODY, handlelength=1.0)

cx = [C_LEFT + C_W + 12.0, C_LEFT + C_W + 28.0, C_LEFT + C_W + 46.0]
for x, head in zip(cx, ["MICP-complete\nmean", "rest\nmean", "MICP-complete\nMAGs with family"]):
    text_mm(x, C_TOP - 3.0, head, ha="center", va="bottom", fontsize=st.FS_STAT,
            color=st.AXIS, linespacing=1.25)
for i, f in enumerate(top):
    yy = C_TOP + (i + 0.5) * (C_H / len(top))
    text_mm(cx[0], yy, f"{hero_fam_mean[f]:.2f}", ha="center", va="center",
            fontsize=st.FS_STAT, color=st.HERO)
    text_mm(cx[1], yy, f"{rest_fam_mean[f]:.2f}", ha="center", va="center",
            fontsize=st.FS_STAT)
    text_mm(cx[2], yy, f"{int(n_hero_with[f])} / {len(st.HEROES)}", ha="center",
            va="center", fontsize=st.FS_STAT)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S4")
