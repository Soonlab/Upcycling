"""Main Fig 2 - MICP module completeness by genus and per-gene prevalence.

Panels (reading order left to right):
  A  MICP module score (ureA-G + cah, 0-8) per GTDB-Tk genus, box + jittered points,
     the six MICP-complete MAGs overplotted as coral rings
  B  per-gene prevalence, MICP-complete group (n = 6) versus the remaining 105 MAGs

Sources
  MAGs_FASTA_files/bakta_results/*/*.tsv  gene presence, via _micp_presence.presence():
                                    a CDS-only keyword scan of all 111 annotations
  Table_S1d_GTDB_Tk_classification.tsv   genus per MAG (all 111)
  Table_S15a_alkaliphile_signature_per_MAG.csv   the MICP-complete / rest group flag

Provenance note (see _job/JOURNAL.md).  Neither earlier source for this panel is usable.
pangenome_work/MICP_Pangenome_Final_Summary.csv, used by the shipped figure, covers only
100 of the 111 MAGs.  Table_S1a_ace_samples_list.csv, used by the first version of this
rebuild, lists only 45 MAGs - 13 of the 66 it omits carry a ure CDS, so the non-MICP-
complete prevalences came out far too low - and it counts the Bakta 5_ureB_sRNA
non-coding feature as a copy of the beta subunit.  Both panels are therefore built from
the annotations directly; _micp_presence documents the defects and asserts the
relationship to S1a.

Colour meanings on this page (one meaning per colour):
  coral  = the MICP-complete group (box outline of the two lineage genera in A, bars in B)
  grey   = every other MAG / the rest group
  black  = medians and individual MAG points
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import HERO, REST, TEXT, GREY, LIGHT, FS_BODY, FS_STAT
from _micp_presence import presence

st.setup()
OUT = HERE / "figures"
SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")

GENES = ["ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG", "cah"]

# ------------------------------------------------------------------ data
tax = pd.read_csv(SUPP / "Table_S1d_GTDB_Tk_classification.tsv", sep="\t",
                  index_col="user_genome")
genus = (tax["classification"].str.extract(r"g__([^;]*)")[0]
         .replace("", np.nan).fillna("Unclassified"))
PANEL = sorted(genus.index)

grp = pd.read_csv(SUPP / "Table_S15a_alkaliphile_signature_per_MAG.csv", index_col=0)
heroes = sorted(grp.index[grp["group"] == "MICP_complete"])
assert heroes == sorted(st.HEROES), heroes
assert len(genus) == len(grp) == 111, (len(genus), len(grp))
n_hero, n_rest = len(heroes), len(PANEL) - len(heroes)

pres = presence().reindex(PANEL)[GENES]
assert pres.notna().all().all()
score = pres.sum(axis=1)
# five of the six MICP-complete MAGs carry the full module; S26 has no protein-coding
# urease beta subunit in its Bakta annotation (see _micp_presence)
assert (score[heroes] == len(GENES)).sum() == 5, score[heroes].to_dict()
assert pres.loc["S26", "ureB"] == 0 and score["S26"] == len(GENES) - 1

is_hero = pd.Series(pres.index.isin(heroes), index=pres.index)
prev_hero = pres[is_hero.values].mean() * 100
prev_rest = pres[~is_hero.values].mean() * 100
assert len(pres[is_hero.values]) == n_hero
assert len(pres[~is_hero.values]) == n_rest

df = pd.DataFrame({"genus": genus.reindex(PANEL), "score": score})
gcount = df["genus"].value_counts()
TOP_N = 9  # style constant: genera drawn individually; the remainder are pooled
top = list(gcount.head(TOP_N).index)
df["grp"] = np.where(df["genus"].isin(top), df["genus"], "Other genera")
order = (df.groupby("grp")["score"].mean().sort_values(ascending=False).index.tolist())

# ------------------------------------------------------------------ page
H = 78.0
fig, ax_mm, text_mm, letter = st.page(H)

LEFT, TOP = 34.0, 11.0
axA = ax_mm(LEFT, TOP, 62.0, 50.0)
axB = ax_mm(LEFT + 80.0, TOP, 58.0, 50.0)
letter(4.0, 6.0, "A")
letter(LEFT + 68.0, 6.0, "B")

# ---- A: module score by genus (horizontal; genus names would collide when rotated)
order = order[::-1]  # highest mean at the top of a horizontal axis
data = [df.loc[df["grp"] == g, "score"].values for g in order]
bp = axA.boxplot(data, positions=range(len(order)), widths=0.62, patch_artist=True,
                 vert=False,
                 medianprops=dict(color=TEXT, lw=0.9),
                 whiskerprops=dict(color=GREY, lw=0.7),
                 capprops=dict(color=GREY, lw=0.7),
                 flierprops=dict(marker="", markersize=0))
for patch in bp["boxes"]:
    patch.set_facecolor(LIGHT)
    patch.set_edgecolor(GREY)
    patch.set_linewidth(0.7)

rng = np.random.default_rng(0)  # jitter only; carries no data
for i, g in enumerate(order):
    sub = df[df["grp"] == g]
    y = i + rng.normal(0, 0.13, size=len(sub))
    hero_mask = sub.index.isin(heroes)
    axA.scatter(sub["score"].values[~hero_mask], y[~hero_mask], s=5,
                color=TEXT, alpha=0.55, lw=0, zorder=3)
    if hero_mask.any():
        axA.scatter(sub["score"].values[hero_mask], y[hero_mask], s=26,
                    facecolor="none", edgecolor=HERO, lw=1.1, zorder=4)

axA.set_yticks(range(len(order)))
axA.set_yticklabels([f"{g} ({int((df['grp'] == g).sum())})" for g in order],
                    fontsize=FS_BODY)
axA.set_xlabel("MICP module score (max 8)")
axA.set_xlim(-0.4, 8.6)
axA.set_xticks([0, 2, 4, 6, 8])
axA.set_ylim(-0.7, len(order) - 0.3)
st.style_axis(axA)
ring = Patch(facecolor="none", edgecolor=HERO, lw=1.1,
             label=f"MICP-complete (n = {n_hero})")
axA.legend(handles=[ring], loc="upper left", bbox_to_anchor=(0.0, -0.16),
           fontsize=FS_STAT, handlelength=1.0, borderpad=0.2)

# ---- B: per-gene prevalence
x = np.arange(len(GENES))
w = 0.38
axB.bar(x - w / 2, prev_hero[GENES].values, w, color=HERO,
        label=f"MICP-complete (n = {n_hero})")
axB.bar(x + w / 2, prev_rest[GENES].values, w, color=REST,
        label=f"Rest (n = {n_rest})")
for xi, v in zip(x - w / 2, prev_hero[GENES].values):
    axB.text(xi, v + 1.5, f"{v:.0f}", ha="center", va="bottom", fontsize=FS_STAT,
             color=TEXT)
for xi, v in zip(x + w / 2, prev_rest[GENES].values):
    axB.text(xi, v + 1.5, f"{v:.0f}", ha="center", va="bottom", fontsize=FS_STAT,
             color=TEXT)
axB.set_xticks(x)
axB.set_xticklabels(GENES, rotation=45, ha="right", fontstyle="italic")
axB.set_ylabel("MAGs with the gene (%)")
axB.set_ylim(0, 112)
axB.set_yticks([0, 25, 50, 75, 100])
st.style_axis(axB)
axB.legend(loc="upper right", bbox_to_anchor=(1.02, 1.12), fontsize=FS_STAT,
           handlelength=1.0, borderpad=0.2)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig2")
