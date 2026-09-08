"""Graphical abstract (revision of 2026-09-09) - one composed page, 190 x 105 mm.

Replaces the four-text-box flow of 2026-04 (scripts/make_graphical_abstract.py in the
public repo), whose animal glyphs did not render and whose content was prose.  Every
element here is drawn from a repository source and no number is typed in:

  left    donut of the 111 MAGs by waste source            Table_S9a_PCoA_coordinates.csv
  centre  radial bar of the MICP module score (0-8) of     _micp_presence.presence()
          every MAG, the six MICP-complete MAGs in coral,   Table_S15a (group flag)
          and a schematic of the ure cluster + cah
  right   finding tiles: lineages, active-site and fold    Table_S12, Table_S22,
          conservation, urease contigs on MGE, Mrp          Table_S17b, Table_S15a,
          prevalence, novel species ANI, MGnify rarity      Table_S10b, Table_S14a
  bottom  the MICP reaction as an icon chain ending in a CaCO3 crystal and a brick

Colours follow _style: coral = MICP-complete, blue / orange = the two lineages,
green = present / product, source palette for the donut, gene palette of Fig 2 for the
operon schematic.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import (FancyBboxPatch, FancyArrow, FancyArrowPatch, Wedge,
                                Polygon, Rectangle, Circle)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _style as st
from _style import (HERO, REST, SPHINGO, PSEUDO, GREEN, GREY, TEXT, AXIS, LIGHT, SOURCE,
                    HEROES, LINEAGE, hero_col)
from _micp_presence import presence

st.setup()
OUT = HERE / "figures_v2"
SUB = Path("/data/data/Upcycling/SUBMISSION_v2")
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
GENE_COL = {"ureA": "#4C72B0", "ureB": "#DD8452", "ureC": "#55A868", "ureD": "#C44E52",
            "ureE": "#8172B3", "ureF": "#937860", "ureG": "#DA8BC3", "cah": "#CCB974"}
GENES = ["ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG", "cah"]
TM_SAME_FOLD = 0.5
ANI_SPECIES = 95.0

# ------------------------------------------------------------------ data
coords = pd.read_csv(SUPP / "Table_S9a_PCoA_coordinates.csv").set_index("MAG")
n_src = coords.Source.value_counts()
SRC_ORDER = ["Cattle", "Swine", "Sheep", "Poultry"]
N_MAG = int(n_src.sum())
assert N_MAG == 111 and set(n_src.index) == set(SRC_ORDER)

grp = pd.read_csv(SUPP / "Table_S15a_alkaliphile_signature_per_MAG.csv", index_col=0)
heroes = sorted(grp.index[grp["group"] == "MICP_complete"])
assert heroes == sorted(HEROES)
mrp_fold = grp.loc[heroes, "Mrp_count"].mean() / grp.loc[grp.index.difference(heroes), "Mrp_count"].mean()
n_lin = pd.Series(LINEAGE).value_counts()

pres = presence(verbose=False)[GENES]
score = pres.sum(axis=1)
n_all8 = int((score == 8).sum())
assert len(score) == N_MAG

sites = pd.read_csv(SUPP / "Table_S12_UreC_active_site_residues.csv")
n_site_match = int((sites[HEROES].values == sites.expected.values[:, None]).sum())
n_site_cells = len(sites) * len(HEROES)

tm = pd.read_csv(SUPP / "Table_S22_ureC_vs_4CEU_tm.csv")
n_fold = int((tm.tm_norm_ref > TM_SAME_FOLD).sum())

ov = pd.read_csv(SUPP / "Table_S17b_ureCah_vs_MGE_overlap.csv").set_index("MAG").loc[HEROES]
n_ure_mge = int(ov.urease_core_MGE_contamination.sum())

nov = pd.read_csv(SUPP / "Table_S10b_ext_Sphingobacterium_novelty.csv").set_index("MAG")
novel = [m for m in ("S13", "S16") if nov.loc[m, "Novel_species_candidate"]]
assert novel == ["S13", "S16"]
ani_txt = " / ".join(f"{nov.loc[m, 'Nearest_ANI']:.1f}" for m in novel)

mg = pd.read_csv(SUPP / "Table_S14a_mgnify_catalog_summary.csv")
n_pool = int(mg.n_species_clusters.sum())
pct_pool = 100 * mg.n_MICP_gene_complete.sum() / n_pool

# ------------------------------------------------------------------ page
W, H = 190.0, 105.0
fig, ax_mm, text_mm, letter = st.page(H, width_mm=W)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)          # y grows downwards, in mm
ax.set_aspect("equal")
ax.axis("off")

BG = "#F7F5F0"
ax.add_patch(Rectangle((0, 0), W, H, facecolor=BG, edgecolor="none", zorder=0))


def panel(x, y, w, h, fc="white", ec=LIGHT, lw=0.6, r=2.5, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, lw=lw, zorder=z))


BOUNDS = []          # (text artist, right edge in mm it must not cross)


def T(x, y, s, size=7, right=None, **kw):
    kw.setdefault("ha", "left")
    kw.setdefault("va", "center")
    kw.setdefault("color", TEXT)
    t = ax.text(x, y, s, fontsize=size, zorder=5, **kw)
    if right is not None:
        BOUNDS.append((t, right))
    return t


# ---- title band
T(W / 2, 6.5, "Livestock-waste metagenomes yield two convergent lineages with a "
              "vertically inherited $\\it{ure}$–$\\it{cah}$ module", 10.5, ha="center",
  fontweight="bold")
T(W / 2, 12.0, "a genomic chassis for alkali-tolerant microbially induced carbonate "
               "precipitation (MICP)", 8.5, ha="center", color=AXIS)

ROW_T, ROW_H = 17.0, 61.0

# ================================================================ left: the input
panel(4, ROW_T, 46, ROW_H)
T(27, ROW_T + 5.0, "111 MAGs from four waste streams", 7.5, ha="center", fontweight="bold")
cx, cy, r_out, r_in = 27.0, ROW_T + 27.0, 15.0, 9.5
a0 = 90.0
for s_ in SRC_ORDER:
    frac = n_src[s_] / N_MAG
    a1 = a0 - 360 * frac
    ax.add_patch(Wedge((cx, cy), r_out, a1, a0, width=r_out - r_in,
                       facecolor=SOURCE[s_.lower()], edgecolor="white", lw=0.8, zorder=3))
    am = np.deg2rad((a0 + a1) / 2)
    rr = (r_out + r_in) / 2
    T(cx + rr * np.cos(am), cy - rr * np.sin(am), f"{n_src[s_]}", 7, ha="center",
      color="white", fontweight="bold")
    a0 = a1
T(cx, cy - 1.5, f"{N_MAG}", 14, ha="center", fontweight="bold")
T(cx, cy + 3.5, "MAGs", 7, ha="center", color=AXIS)
for k, s_ in enumerate(SRC_ORDER):
    yk = ROW_T + 47.5 + (k // 2) * 5.0
    xk = 10.0 + (k % 2) * 21.0
    ax.add_patch(Circle((xk, yk), 1.4, facecolor=SOURCE[s_.lower()], edgecolor="none",
                        zorder=4))
    T(xk + 2.6, yk, f"{s_.lower()}  n = {n_src[s_]}", 6.5)
T(27, ROW_T + 58.0, "cattle · swine · sheep · poultry slurries", 6, ha="center", color=GREY)

# ================================================================ centre: the module
panel(56, ROW_T, 62, ROW_H)
T(87, ROW_T + 5.0, "MICP module score across the panel", 7.5, ha="center",
  fontweight="bold")
# radial bar: one wedge per MAG, ordered by score then identifier, heroes in coral
order = score.sort_values(ascending=False).index.tolist()
order = sorted(order, key=lambda m: (-score[m], m not in heroes, m))
ccx, ccy = 71.0, ROW_T + 26.0
R0, RSC = 6.6, 1.30                  # inner radius and mm per score unit
pitch = 360.0 / N_MAG
R_MAX = R0 + 8 * RSC
for i, m in enumerate(order):
    a_hi = 90.0 - i * pitch
    a_lo = a_hi - pitch * 0.85
    # a faint full-length wedge behind every MAG, so the ring reads as 111 slots and a
    # score of 0 is an empty slot rather than a gap in the data
    ax.add_patch(Wedge((ccx, ccy), R_MAX, a_lo, a_hi, width=R_MAX - R0,
                       facecolor="#EDEDED", edgecolor="none", zorder=2))
    if score[m] == 0:
        continue
    col = HERO if m in heroes else (GREEN if score[m] == 8 else "#B4B4B4")
    ax.add_patch(Wedge((ccx, ccy), R0 + score[m] * RSC, a_lo, a_hi, width=score[m] * RSC,
                       facecolor=col, edgecolor="none", zorder=3))
ax.add_patch(Circle((ccx, ccy), R_MAX, facecolor="none", edgecolor=LIGHT, lw=0.5,
                    zorder=4))
ax.add_patch(Circle((ccx, ccy), R0, facecolor="white", edgecolor=LIGHT, lw=0.5, zorder=4))
T(ccx, ccy - 1.4, f"{n_all8}", 9.5, ha="center", fontweight="bold", color=GREEN)
T(ccx, ccy + 2.2, "of 111", 5.2, ha="center", color=AXIS)
# key beside the ring
kx, ky = 90.0, ROW_T + 11.0
for k, (col, lab) in enumerate(((HERO, f"MICP-complete, n = {len(heroes)}"),
                                (GREEN, f"$\\it{{ureA–G}}$ + $\\it{{cah}}$, n = {n_all8}"),
                                ("#C9C9C9", "score 0–7"))):
    ax.add_patch(Rectangle((kx, ky + k * 4.6 - 1.2), 3.2, 2.4, facecolor=col,
                           edgecolor="none", zorder=4))
    T(kx + 4.2, ky + k * 4.6, lab, 6.0, right=117.0)
T(kx, ky + 15.0, "bar length = genes present", 5.4, color=GREY, right=117.0)
T(kx, ky + 18.3, "centre = MAGs with all eight", 5.4, color=GREY, right=117.0)
# lineages
T(kx, ky + 25.0, "two convergent lineages", 6.5, fontweight="bold", right=117.0)
T(kx, ky + 29.4, f"$\\it{{Sphingobacterium}}$ × {int(n_lin['Sphingobacterium'])}", 6.0,
  color=SPHINGO, right=117.0)
T(kx, ky + 33.2, f"$\\it{{Pseudomonas}}$_E × {int(n_lin['Pseudomonas_E'])}", 6.0,
  color=PSEUDO, right=117.0)

# operon schematic along the bottom of the centre panel
oy = ROW_T + 50.0
ox = 60.5
gw = {"ureA": 3.0, "ureB": 3.0, "ureC": 7.5, "ureD": 4.0, "ureE": 3.0, "ureF": 3.5,
      "ureG": 3.5}
x = ox
for g in GENES[:-1]:
    ax.add_patch(FancyArrow(x, oy, gw[g], 0, width=2.4, head_width=3.4, head_length=1.2,
                            length_includes_head=True, facecolor=GENE_COL[g],
                            edgecolor=AXIS, lw=0.3, zorder=4))
    dy = -3.2 if GENES.index(g) % 2 == 0 else -6.4
    ax.plot([x + gw[g] / 2, x + gw[g] / 2], [oy - 1.5, oy + dy + 1.0], color=GREY,
            lw=0.3, zorder=3)
    T(x + gw[g] / 2, oy + dy, g, 5, ha="center", style="italic")
    x += gw[g] + 0.5
x += 2.5
ax.plot([x - 2.0, x + 0.5], [oy, oy], color=GREY, lw=0.6, ls=":", zorder=3)
ax.add_patch(FancyArrow(x + 1.0, oy, 5.0, 0, width=2.4, head_width=3.4, head_length=1.2,
                        length_includes_head=True, facecolor=GENE_COL["cah"],
                        edgecolor=AXIS, lw=0.3, zorder=4))
T(x + 3.5, oy - 3.2, "cah", 5, ha="center", style="italic")
T(87, oy + 5.2, f"{n_ure_mge}/{len(heroes)} $\\it{{ure}}$ contigs on a mobile element", 6.0,
  ha="center", color=AXIS)
T(87, oy + 8.6, "vertically inherited in both lineages", 6.0, ha="center", color=AXIS)

# ================================================================ right: the findings
panel(124, ROW_T, 62, ROW_H)
T(155, ROW_T + 5.0, "What the six MICP-complete MAGs share", 7.5, ha="center",
  fontweight="bold")
tiles = [
    (f"{n_site_match}/{n_site_cells}", "UreC active-site residues",
     "identical to $\\it{S. pasteurii}$", GREEN),
    (f"{mrp_fold:.1f}×", "Mrp Na⁺/H⁺ antiporter prevalence",
     "the alkali-tolerance signature", HERO),
    (f"{pct_pool:.2f} %", f"of {n_pool:,} livestock species",
     "clusters carry the module", AXIS),
    ("S13 · S16", "novel $\\it{Sphingobacterium}$",
     f"ANI {ani_txt} % to RefSeq", SPHINGO),
]
TILE_H, TILE_P = 10.6, 11.7
ty = ROW_T + 9.0
for big, line1, line2, col in tiles:
    panel(128, ty, 54, TILE_H, fc="white", ec=col, lw=0.9, r=1.8, z=2)
    ax.add_patch(Rectangle((128, ty), 1.6, TILE_H, facecolor=col, edgecolor="none",
                           zorder=3))
    T(132, ty + TILE_H / 2, big, 9.6, fontweight="bold", color=col)
    T(150, ty + 3.6, line1, 5.8, right=181.0)
    T(150, ty + 7.2, line2, 5.4, color=AXIS, right=181.0)
    ty += TILE_P
T(155, ty + 2.2, f"urease fold kept {n_fold}/{len(tm)} · ESMFold TM > {TM_SAME_FOLD:g}", 5.8,
  ha="center", color=AXIS)

# arrows between the three panels
for xa in (50.5, 118.5):
    ax.add_patch(FancyArrowPatch((xa, ROW_T + ROW_H / 2), (xa + 5.0, ROW_T + ROW_H / 2),
                                 arrowstyle="-|>", mutation_scale=9, lw=1.4, color=AXIS,
                                 zorder=6))

# ================================================================ bottom: the reaction
BY = 87.5
panel(4, BY - 6.0, 182, 16.0, fc="#EEF3EE", ec=GREEN, lw=0.7)
T(8, BY - 2.8, "waste-coupled MICP", 6.5, fontweight="bold", color=GREEN)


def molecule(x, label, col, r=4.0):
    ax.add_patch(Circle((x, BY + 1.5), r, facecolor="white", edgecolor=col, lw=1.0,
                        zorder=4))
    T(x, BY + 1.5, label, 6.2, ha="center", color=col, fontweight="bold")


def step(x0, x1, label):
    ax.add_patch(FancyArrowPatch((x0, BY + 1.5), (x1, BY + 1.5), arrowstyle="-|>",
                                 mutation_scale=7, lw=1.0, color=AXIS, zorder=4))
    T((x0 + x1) / 2, BY - 2.6, label, 5.3, ha="center", color=AXIS, style="italic")


molecule(44, "urea", AXIS)
step(49, 60, "urease")
molecule(67, "NH₄⁺", PSEUDO, r=3.6)
T(72.5, BY + 1.5, "+", 7, ha="center")
molecule(78, "CO₂", AXIS, r=3.6)
step(83, 95, "carbonic\nanhydrase")
molecule(102, "CO₃²⁻", SPHINGO, r=4.2)
T(110.5, BY + 1.5, "+ Ca²⁺", 6.5, ha="center", fontweight="bold")
step(117, 129, "pH ↑ · nucleation")
# CaCO3 rhombohedron
crystal = np.array([[133, BY + 5.0], [139, BY + 5.0], [142, BY - 1.0], [136, BY - 1.0]])
ax.add_patch(Polygon(crystal, closed=True, facecolor=GREEN, edgecolor="white", lw=0.8,
                     zorder=4))
ax.add_patch(Polygon(crystal + np.array([2.0, -2.5]), closed=True, facecolor="#5FA27A",
                     edgecolor="white", lw=0.8, zorder=3))
T(145, BY + 1.5, "CaCO₃", 7, fontweight="bold", color=GREEN)
step(157, 165, "")
# brick wall icon
BX = 167.5
for rr_ in range(3):
    for c_ in range(3):
        off = -1.8 if rr_ % 2 else 0.0
        ax.add_patch(Rectangle((BX + c_ * 3.9 + off, BY - 3.2 + rr_ * 3.0), 3.6, 2.6,
                               facecolor="#C9A27A", edgecolor="white", lw=0.5, zorder=4))
T((BX + 5.6), BY + 8.0, "biocement", 6.2, ha="center", fontweight="bold",
  color="#8E6C3A")

# a label that runs past its panel is a failure st.audit cannot see (there is no other
# text to collide with), so every registered label is measured against its right edge
fig.canvas.draw()
r = fig.canvas.get_renderer()
scale = W / fig.canvas.get_width_height()[0]
for t, right in BOUNDS:
    x1 = t.get_window_extent(renderer=r).x1 * scale
    assert x1 <= right, (t.get_text()[:40], round(x1, 1), right)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Graphical_abstract", dpi=300)
for ext in ("png", "pdf", "svg"):
    (SUB / f"03_Graphical_abstract.{ext}").write_bytes((OUT / f"Graphical_abstract.{ext}").read_bytes())
print(f"  {N_MAG} MAGs | all-8 {n_all8} | sites {n_site_match}/{n_site_cells} | fold {n_fold} | "
      f"Mrp {mrp_fold:.1f}x | MGnify {pct_pool:.2f} % of {n_pool} | ANI {ani_txt}")
