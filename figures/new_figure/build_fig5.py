"""Main Fig 5 - novel-species delineation of S13 and S16 within Sphingobacterium.

One composed page, 180 mm wide, three panels:
  A  within-panel novelty screen: closest-GTDB-reference ANI per MAG, plus the MAGs for
     which GTDB-Tk returned no species-level ANI at all, listed as their own block
  B  reciprocal-best mmseqs2 AAI of S13 and S16 against the other in-panel Sphingobacterium
  C  skani ANI of the six study Sphingobacterium MAGs against the 63 RefSeq genomes of the
     genus, one block per study MAG

Sources
  Table_S8_novelty_ANI_screen.csv               ANI to closest GTDB reference, novelty flag
  Table_S4b_AAI_S13_S16.csv                     AAI, reciprocal-best hit counts
  Table_S10a_ext_Sphingobacterium_ANI_matrix.csv  6 study MAGs x 63 RefSeq genomes
  Table_S10b_ext_Sphingobacterium_novelty.csv   nearest reference per study MAG (asserted)

Two corrections against the shipped figure (see _job/journal_main_4_7.md)
  1. The shipped panel a is titled "21 candidates < 95 % ANI".  No MAG in Table S8 has an
     ANI below 95: the 21 novelty candidates are exactly the 21 MAGs with NO species-level
     ANI, S13 and S16 among them.  They are drawn here as a named block instead of being
     dropped, and the count of MAGs actually below 95 % is what the axis shows.
  2. The shipped panel c scatters "all 63 RefSeq genomes" per MAG.  skani reports only the
     references it aligns; the matrix stores 0.0 for the rest, and those zeros were plotted
     outside the y-limits, so they were invisible rather than absent.  Only the reported
     hits (ANI >= 80) are drawn here, and each block states how many there are.

Colour meanings on this page:
  coral   a MICP-complete MAG (its point in A, its block header in C) and the < 95 % species
          verdict of a bar in C
  green   >= 95 % species-level match, in C
  blue    the Sphingobacterium lineage: dark for the S13 query block and light for the S16
          query block in B.  Both queries are Sphingobacterium, so the two shades separate
          the query blocks and each block names its query; they do not encode a lineage
          contrast the way blue vs orange does on Fig 8
  grey    a MAG that is not MICP-complete, and every threshold rule
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, REST, GREEN, GREY, TEXT, LIGHT, FS_BODY, FS_STAT, HEROES

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

ANI_SPECIES = 95.0      # species boundary, ANI
AAI_SPECIES = 95.0
AAI_GENUS = 70.0
ANI_REPORTED = 80.0     # skani stores 0.0 for pairs it did not align

# ------------------------------------------------------------------ data
nov = pd.read_csv(SUPP / "Table_S8_novelty_ANI_screen.csv")
has_ani = nov.ANI.notna()
ranked = nov[has_ani].sort_values("ANI").reset_index(drop=True)
no_ani = nov[~has_ani].reset_index(drop=True)
# the novelty flag is exactly "no species-level ANI"; state it by assertion, not in prose
assert set(nov.user_genome[nov.Novel_sp_candidate]) == set(no_ani.user_genome)
assert (ranked.ANI >= ANI_SPECIES).all()
N_BELOW = int((ranked.ANI < ANI_SPECIES).sum())

aai = pd.read_csv(SUPP / "Table_S4b_AAI_S13_S16.csv")
ani_m = pd.read_csv(SUPP / "Table_S10a_ext_Sphingobacterium_ANI_matrix.csv", index_col=0)
ext_nov = pd.read_csv(SUPP / "Table_S10b_ext_Sphingobacterium_novelty.csv").set_index("MAG")

our = [i for i in ani_m.index if i.startswith("OUR_")]
refs = [c for c in ani_m.columns if c.startswith("REF_")]
N_REF = len(refs)
blocks = []
for o in sorted(our, key=lambda s: s.replace("OUR_", "")):
    mag = o.replace("OUR_", "")
    v = ani_m.loc[o, refs]
    hit = v[v >= ANI_REPORTED].sort_values(ascending=False)
    # the nearest-reference row of Table S10b must be reproduced by the matrix
    assert abs(hit.iloc[0] - ext_nov.loc[mag, "Nearest_ANI"]) < 1e-9, mag
    assert (hit.iloc[0] < ANI_SPECIES) == bool(ext_nov.loc[mag, "Novel_species_candidate"]), mag
    blocks.append((mag, hit))

# ------------------------------------------------------------------ layout
TOP = 10.0
H_A = 46.0
W_A = 78.0
X_A = 16.0
X_A2 = X_A + W_A + 14.0          # the no-ANI block sits beside the scatter, same panel
T_B = TOP + H_A + 20.0
H_B1 = 22.0
T_C = T_B + 2 * H_B1 + 14.0
ROW_C = 4.4
H_C_BLOCK = 5 * ROW_C
H = T_C + 2 * (H_C_BLOCK + 20.0) + 12.0

fig, ax_mm, text_mm, letter = st.page(H)

# ---------------------------------------------------------------- panel A
axA = ax_mm(X_A, TOP, W_A, H_A)
is_hero = ranked.user_genome.isin(HEROES).values
axA.scatter(np.arange(len(ranked))[~is_hero], ranked.ANI[~is_hero], s=5, color=REST,
            linewidth=0, zorder=2)
axA.scatter(np.arange(len(ranked))[is_hero], ranked.ANI[is_hero], s=20, color=HERO,
            edgecolor="white", linewidth=0.4, zorder=3)
hero_idx = list(np.where(is_hero)[0])
for k, i in enumerate(hero_idx):
    prev_close = k > 0 and (i - hero_idx[k - 1]) < 6
    dy = -7 if prev_close else 3
    axA.annotate(ranked.user_genome[i], (i, ranked.ANI[i]), xytext=(3, dy),
                 textcoords="offset points", fontsize=FS_STAT, color=HERO)
axA.axhline(ANI_SPECIES, ls="--", lw=0.7, color=GREY)
axA.text(len(ranked) - 1, ANI_SPECIES, f"{ANI_SPECIES:g} % species boundary", fontsize=FS_STAT,
         color=GREY, ha="right", va="bottom")
axA.set_xlabel(f"MAGs with a species-level ANI, ranked (n = {len(ranked)})")
axA.set_ylabel("ANI to closest GTDB reference (%)")
axA.set_xticks([])
st.style_axis(axA, left=True, bottom=True)
axA.spines["bottom"].set_visible(False)

# the 21 MAGs GTDB-Tk left without a species-level ANI, as an identifier block
axA2 = ax_mm(X_A2, TOP, 62.0, H_A)
axA2.axis("off")
text_mm(X_A2, TOP - 1.0, f"no species-level ANI  (n = {len(no_ani)})", fontsize=FS_STAT,
        ha="left", va="bottom", color=TEXT)
names = list(no_ani.user_genome)
NCOL = 4
per = int(np.ceil(len(names) / NCOL))
for k, name in enumerate(names):
    col, row = k // per, k % per
    hero = name in HEROES
    text_mm(X_A2 + 1.5 + col * 15.0, TOP + 3.0 + row * 4.2, name, fontsize=FS_BODY,
            ha="left", va="top", color=HERO if hero else TEXT)

# ---------------------------------------------------------------- panel B
QCOL = {"S13": st.SPHINGO, "S16": st.SPHINGO_LT}
for idx, q in enumerate(["S13", "S16"]):
    sub = aai[aai.Query == q].dropna(subset=["AAI"]).sort_values("AAI")
    ax = ax_mm(X_A + 14.0, T_B + idx * H_B1, 108.0, H_B1 - 6.0)
    y = np.arange(len(sub))
    ax.barh(y, sub.AAI.values, height=0.62, color=QCOL[q], linewidth=0)
    lbl = [f"{t}  {s.replace('Sphingobacterium ', 'S. ')}" if isinstance(s, str) and s
           else f"{t}" for t, s in zip(sub.Target, sub.Target_species)]
    ax.set_yticks(y)
    ax.set_yticklabels(lbl, fontsize=FS_BODY)
    for yi, v in zip(y, sub.AAI.values):
        ax.text(v + 0.4, yi, f"{v:.2f}", va="center", ha="left", fontsize=FS_STAT)
    ax.axvline(AAI_SPECIES, ls="--", lw=0.7, color=GREY)
    ax.axvline(AAI_GENUS, ls=":", lw=0.7, color=GREY)
    ax.set_xlim(55, 100)
    ax.tick_params(axis="y", length=0)
    st.style_axis(ax, left=False, bottom=True)
    text_mm(X_A, T_B + idx * H_B1 - 1.0, f"query {q}", fontsize=FS_STAT, ha="left",
            va="bottom", color=QCOL[q])
    if idx == 1:
        ax.set_xlabel("Reciprocal-best AAI (%)")
    else:
        ax.set_xticklabels([])
        ax.legend(handles=[Line2D([0], [0], ls="--", lw=0.7, color=GREY,
                                  label=f"{AAI_SPECIES:g} % species"),
                           Line2D([0], [0], ls=":", lw=0.7, color=GREY,
                                  label=f"{AAI_GENUS:g} % genus")],
                  loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, fontsize=FS_STAT,
                  handlelength=1.4, handletextpad=0.4, borderpad=0.2, columnspacing=1.2)

# ---------------------------------------------------------------- panel C
CW, CH = 48.0, H_C_BLOCK
for k, (mag, hit) in enumerate(blocks):
    col, row = k % 3, k // 3
    x = X_A + col * (CW + 3.0)
    t = T_C + row * (CH + 20.0)
    ax = ax_mm(x + 19.0, t, CW - 21.0, CH)
    y = np.arange(len(hit))
    cols = [GREEN if v >= ANI_SPECIES else REST for v in hit.values]
    ax.barh(y, hit.values, height=0.62, color=cols, linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels([r.replace("REF_", "") for r in hit.index], fontsize=FS_STAT)
    for yi, v in zip(y, hit.values):
        ax.text(v + 0.3, yi, f"{v:.2f}", va="center", ha="left", fontsize=FS_STAT)
    ax.axvline(ANI_SPECIES, ls="--", lw=0.7, color=GREY)
    ax.set_xlim(84, 102)
    ax.set_ylim(-0.7, max(4, len(hit)) - 0.3)
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    st.style_axis(ax, left=False, bottom=True)
    ax.set_xlabel("ANI (%)" if row == 1 else "")
    if row == 0:
        ax.set_xticklabels([])
    text_mm(x, t - 1.0, f"{mag}   {len(hit)} of {N_REF} RefSeq hits", fontsize=FS_STAT,
            ha="left", va="bottom", color=HERO if mag in HEROES else TEXT)

handles = [Line2D([0], [0], marker="s", ls="", ms=4, mfc=GREEN, mec="none",
                  label=f"ANI >= {ANI_SPECIES:g} %"),
           Line2D([0], [0], marker="s", ls="", ms=4, mfc=REST, mec="none",
                  label=f"ANI < {ANI_SPECIES:g} %"),
           Line2D([0], [0], ls="--", lw=0.7, color=GREY,
                  label=f"{ANI_SPECIES:g} % species boundary")]
fig.legend(handles=handles, loc="lower right",
           bbox_to_anchor=(1 - 6 / 180, 4 / H), fontsize=FS_BODY,
           handletextpad=0.4, borderpad=0.2, labelspacing=0.25)

letter(4, 4, "A")
letter(4, T_B - 5.0, "B")
letter(4, T_C - 6.0, "C")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig5")
