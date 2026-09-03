"""Main Fig 3 - gene order of the ure cluster on the main ure contig of each
MICP-complete MAG.

One synteny track per MAG (six tracks, top to bottom, ordered by the number of ure genes
recovered on the contig).  Arrows are drawn to scale from the Bakta CDS coordinates;
arrow direction is the coding strand; colour is gene identity.

Sources
  MAGs_FASTA_files/bakta_results/<MAG>/<MAG>.gff3   CDS coordinates, strand, product
  Table_S1c_hero_cluster_audit.csv                  the contig of record for each MAG,
                                                    the number of ure genes on it and the
                                                    cluster span; every track is asserted
                                                    against this table
  Table_S3a_HGT_ureCah_cluster.csv                  mobile-element count and regional GC
                                                    of the same window (stat column)

The gene classifier is the recipe of scripts/01_main_figures.py (Bakta product string ->
ure subunit / carbonic anhydrase), re-applied here; the contig drawn is the one named in
Table S1c rather than the densest-cluster search of the old script, because for S23 the
search ties at four distinct ure genes and picks contig_151 while the table of record
names contig_220 (recorded in the journal).

Colour meanings on this page (one meaning per colour): each of the eight MICP genes has
its own colour, and grey is any other CDS in the window.  No other colour is used.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrow, Patch
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import TEXT, GREY, AXIS, FS_BODY, FS_STAT

st.setup()
OUT = HERE / "figures"
BASE = Path("/data/data/Upcycling")
SUPP = BASE / "SUBMISSION/Supplementary_tables"
BAKTA = BASE / "MAGs_FASTA_files/bakta_results"

GENES = ["ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG", "cah"]
# gene identity palette: eight distinct hues plus grey for everything else
GENE_COL = {"ureA": "#4C72B0", "ureB": "#DD8452", "ureC": "#55A868", "ureD": "#C44E52",
            "ureE": "#8172B3", "ureF": "#937860", "ureG": "#DA8BC3", "cah": "#CCB974",
            "other": "#D0D0D0"}
FLANK_BP = 2000  # style constant: window padding drawn either side of the cluster


def classify(product, gene):
    p = (product or "").lower()
    g = (gene or "").lower()
    if "urease subunit alpha" in p:
        return "ureC"
    if "urease subunit beta" in p:
        return "ureB"
    if "urease subunit gamma" in p:
        return "ureA"
    if g in ("ured", "uree", "uref", "ureg"):
        return g[:3] + g[3].upper()
    for k in ("ured", "uree", "uref", "ureg"):
        if k in p:
            return k[:3] + k[3].upper()
    if "carbonic anhyd" in p:
        return "cah"
    return "other"


def read_gff(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            rows.append(dict(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                             cls=classify(attrs.get("product", ""), attrs.get("gene", ""))))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ data
audit = pd.read_csv(SUPP / "Table_S1c_hero_cluster_audit.csv", index_col=0)
hgt = pd.read_csv(SUPP / "Table_S3a_HGT_ureCah_cluster.csv", index_col=0)
assert sorted(audit.index) == sorted(st.HEROES), sorted(audit.index)

tracks = {}
for mag in audit.index:
    ann = read_gff(BAKTA / mag / f"{mag}.gff3")
    ctg = audit.loc[mag, "main_contig"]
    ure = ann[(ann.contig == ctg) & (ann.cls.str.startswith("ure"))]
    n_distinct = ure.cls.nunique()  # the audit counts distinct ure genes, not CDS copies
    assert n_distinct == int(audit.loc[mag, "ure_genes_on_main_contig"]), \
        (mag, n_distinct, audit.loc[mag, "ure_genes_on_main_contig"])
    span_kb = (ure.end.max() - ure.start.min()) / 1000
    assert abs(span_kb - float(audit.loc[mag, "cluster_span_kb_main"])) < 0.02, \
        (mag, span_kb, audit.loc[mag, "cluster_span_kb_main"])
    lo, hi = ure.start.min() - FLANK_BP, ure.end.max() + FLANK_BP
    win = ann[(ann.contig == ctg) & (ann.end >= lo) & (ann.start <= hi)].sort_values("start")
    tracks[mag] = dict(ctg=ctg, span=span_kb, lo=lo, win=win,
                       n_ure=n_distinct, length_kb=(hi - lo) / 1000,
                       mge=int(hgt.loc[mag, "MobileElements"]),
                       dgc=float(hgt.loc[mag, "DeltaGC"]))

order = sorted(tracks, key=lambda m: (-tracks[m]["n_ure"], -tracks[m]["span"]))

# ------------------------------------------------------------------ page
ROW_H, GAP = 13.0, 6.0
TOP, LEFT, PLOT_W = 12.0, 26.0, 118.0
STAT_X = LEFT + PLOT_W + 6.0
H = TOP + len(order) * (ROW_H + GAP) + 15.0
fig, ax_mm, text_mm, letter = st.page(H)
letter(4.0, 6.0, "A")

# stat column headers (table column headers, not a title)
text_mm(STAT_X, TOP - 4.6, "ure genes", fontsize=FS_STAT, ha="left", color=TEXT)
text_mm(STAT_X + 13.0, TOP - 4.6, "MGE", fontsize=FS_STAT, ha="left", color=TEXT)
text_mm(STAT_X + 21.5, TOP - 4.6, "Δ GC", fontsize=FS_STAT, ha="left", color=TEXT)

LEVELS = (0.42, 1.02, 1.62)  # label rows above the track, used only when labels collide


def place_levels(items, kb_per_mm):
    """Assign each gene label the lowest row that clears the previous label on that row.

    Label width is estimated from the character count at the 6.5 pt stat size, converted
    to kb through the track's own scale, so a tight three-gene cluster inside a 32 kb
    window staggers instead of printing on top of itself.
    """
    last = [-1e9] * len(LEVELS)
    out = []
    for xc, text in items:
        half = len(text) * FS_STAT * 0.353 * 0.52 * kb_per_mm / 2
        for lv in range(len(LEVELS)):
            if xc - half > last[lv] + 0.15 * kb_per_mm:
                last[lv] = xc + half
                out.append(lv)
                break
        else:
            last[0] = xc + half
            out.append(0)
    return out


for i, mag in enumerate(order):
    t = tracks[mag]
    top = TOP + i * (ROW_H + GAP)
    ax = ax_mm(LEFT, top, PLOT_W, ROW_H)
    kb_per_mm = t["length_kb"] / PLOT_W
    labelled = [((r.start + r.end) / 2000 - t["lo"] / 1000, r.cls)
                for _, r in t["win"].iterrows() if r.cls != "other"]
    levels = place_levels(labelled, kb_per_mm)
    for _, r in t["win"].iterrows():
        x0 = (r.start - t["lo"]) / 1000
        x1 = (r.end - t["lo"]) / 1000
        col = GENE_COL.get(r.cls, GENE_COL["other"])
        if r.strand == "-":
            xs, dx = x1, -(x1 - x0)
        else:
            xs, dx = x0, (x1 - x0)
        head = min(0.30, abs(dx) * 0.45)
        ax.add_patch(FancyArrow(xs, 0, dx, 0, width=0.40, head_width=0.58,
                                head_length=head, length_includes_head=True,
                                facecolor=col, edgecolor=AXIS, lw=0.35))
    for (xc, text), lv in zip(labelled, levels):
        y = LEVELS[lv]
        ax.plot([xc, xc], [0.24, y - 0.04], color=GREY, lw=0.4, zorder=1)
        ax.text(xc, y, text, ha="center", va="bottom", fontsize=FS_STAT,
                style="italic", color=TEXT)
    # each track is scaled to its own window: a shared scale would squeeze the five
    # short clusters into a quarter of the page beside the 32 kb M1 window
    ax.set_xlim(-0.3, t["length_kb"] + 0.3)
    ax.set_ylim(-0.7, LEVELS[-1] + 0.62)
    ax.set_yticks([])
    step = max(1.0, round(t["length_kb"] / 6))
    ax.set_xticks(np.arange(0, t["length_kb"] + 0.01, step))
    st.style_axis(ax, left=False)
    if mag == order[-1]:
        ax.set_xlabel("Position in the drawn window (kb)")
    # bare identifiers to the left of each track
    text_mm(LEFT - 2.0, top + ROW_H / 2 - 1.8, mag, ha="right", va="center",
            fontsize=FS_BODY, fontweight="bold", color=TEXT)
    text_mm(LEFT - 2.0, top + ROW_H / 2 + 1.4, t["ctg"], ha="right", va="center",
            fontsize=FS_STAT, color=GREY)
    # stat columns bound to the row
    text_mm(STAT_X, top + ROW_H / 2 - 0.2, f"{t['n_ure']} / 7", ha="left", va="center",
            fontsize=FS_STAT, color=TEXT)
    text_mm(STAT_X + 13.0, top + ROW_H / 2 - 0.2, f"{t['mge']}", ha="left", va="center",
            fontsize=FS_STAT, color=TEXT)
    text_mm(STAT_X + 21.5, top + ROW_H / 2 - 0.2, f"{t['dgc']:.1f}", ha="left",
            va="center", fontsize=FS_STAT, color=TEXT)

# cah is deliberately absent from the key: Table S1c records cah_on_main_contig = False
# for all six MAGs, so no cah arrow is drawn and a key entry would have no mark
drawn = [g for g in GENES if any((t["win"].cls == g).any() for t in tracks.values())]
handles = [Patch(facecolor=GENE_COL[g], edgecolor=AXIS, lw=0.35, label=g) for g in drawn] + \
          [Patch(facecolor=GENE_COL["other"], edgecolor=AXIS, lw=0.35, label="other CDS")]
fig.legend(handles=handles, loc="lower center", ncol=9, fontsize=FS_STAT, frameon=False,
           bbox_to_anchor=(0.5, 0.008), handlelength=1.1, columnspacing=1.1,
           handletextpad=0.4)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig3")
