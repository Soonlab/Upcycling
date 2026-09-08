"""Consolidated main Fig 2 - architecture and dosage of the ure-cah cluster.

Consolidation of 2026-09-04 (consolidation_260904/DESIGN.md): old main Fig 3, old Suppl
Fig S14 and panel B of old Suppl Fig S15 become one page.  Panel A of old S15 (geNomad
plasmid/virus contigs per MAG, whole panel) is not carried here; it moves to the
consolidated supplementary page.

Panels (reading order, top to bottom):
  A  gene order of the ure cluster on the main ure contig of each MICP-complete MAG,
     one synteny track per MAG                                            [old Fig 3]
  B  MICP pathway gene dosage, the six MICP-complete MAGs by urease subunit, carbonic
     anhydrase family and calcium / carbonate / cation transport gene     [old Fig S14]
  C  per-MAG cross-check of the urease and carbonic-anhydrase contigs against the
     geNomad mobile-element calls: one square per urease-core contig and per carbonic-
     anhydrase contig, filled coral where geNomad flagged that contig as a plasmid, and
     beside it the MAG's total plasmid- and virus-flagged contig counts as bars
                                                                          [old Fig S15B]

Revision of 2026-09-09: all type on this page is set 1.3x larger than the page body
(FS_A_* constants; second round extended it from A to B and C) with the stat block tightened; panel B spans the same width as the
synteny tracks of A, its stat columns and colour bar sit in the same right-hand block as
A's stat columns; panel C, previously a table of numbers, is drawn as contig squares plus
count bars.

Sources
  MAGs_FASTA_files/bakta_results/<MAG>/<MAG>.gff3   CDS coordinates, strand, product (A)
  Table_S1c_hero_cluster_audit.csv                  the contig of record for each MAG,
                                                    the number of ure genes on it and the
                                                    cluster span; every track in A is
                                                    asserted against this table
  Table_S3a_HGT_ureCah_cluster.csv                  mobile-element count and regional GC
                                                    of the same window (stat columns of A)
  Table_S15b_stoichiometry_per_MAG.csv              gene copy numbers and the two stored
                                                    completeness calls (B)
  Table_S17a_genomad_summary_per_MAG.csv            per-MAG plasmid/virus contig counts,
                                                    used to cross-check S17b (C)
  Table_S17b_ureCah_vs_MGE_overlap.csv              the urease / CA contig lists and the
                                                    two contamination counts (C)

Method notes carried over from the source scripts.
  A  the gene classifier is the recipe of scripts/01_main_figures.py (Bakta product
     string -> ure subunit / carbonic anhydrase), re-applied here; the contig drawn is
     the one named in Table S1c rather than the densest-cluster search of the old script,
     because for S23 the search ties at four distinct ure genes and picks contig_151
     while the table of record names contig_220 (recorded in the journal).
  B  the Ca_pathway call is recomputed here (Ca_transporter or Ca_ATPase present) and
     asserted against the stored flag, which is what makes the "at least one Ca-handling
     gene in 5 of 6" statement checkable: S26 carries neither, and its five
     CO3_transporter copies do not enter the stored call.
  C  the contig lists of S17b are split on the comma and every contig becomes one
     square; a square is coral when its contig is named in the matching *_on_plasmid /
     *_on_virus column.  The number of coral squares per MAG is asserted against the two
     stored contamination counts, and the per-MAG plasmid contig counts of S17a and S17b
     are asserted equal.  Only CA_on_plasmid is non-empty in the source (S26, contig_66).

Colour meanings on this page (one meaning per colour):
  gene identity palette (A)  one hue per MICP gene, grey for any other CDS in the window
  green intensity (B)        gene copy number; the stat columns of B are text, not
                             colour, so nothing competes with the copy-number scale
  blue / orange (B, C)       Sphingobacterium and Pseudomonas_E MAG identifiers
  coral (C)                  a urease or CA contig that geNomad flagged as a plasmid
  dark grey / hatched (C)    plasmid-flagged and virus-flagged contig totals per MAG
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrow, Patch, Rectangle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
from _style import (HERO, GREEN, TEXT, GREY, AXIS, LIGHT, FS_BODY, FS_STAT,
                    HEROES, hero_col)

st.setup()
OUT = HERE / "figures_v2"
BASE = Path("/data/data/Upcycling")
SUPP = BASE / "SUBMISSION/Supplementary_tables"
BAKTA = BASE / "MAGs_FASTA_files/bakta_results"

GENES = ["ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG", "cah"]
# gene identity palette: eight distinct hues plus grey for everything else
GENE_COL = {"ureA": "#4C72B0", "ureB": "#DD8452", "ureC": "#55A868", "ureD": "#C44E52",
            "ureE": "#8172B3", "ureF": "#937860", "ureG": "#DA8BC3", "cah": "#CCB974",
            "other": "#D0D0D0"}
FLANK_BP = 2000  # style constant: window padding drawn either side of the cluster

# panel B column groups (category names, not data)
URE = ["ureA", "ureB", "ureC", "ureD_H", "ureE", "ureF", "ureG"]
CA = ["cah_alphaCA", "canA_gammaCA", "cynT_betaCA", "CA_generic"]
ION = ["Ca_transporter", "Ca_ATPase", "CO3_transporter", "Na_H_antiporter_Mrp",
       "K_transport"]
COLS = URE + CA + ION
NICE = {"ureD_H": "ureD/H", "cah_alphaCA": "cah (α-CA)", "canA_gammaCA": "canA (γ-CA)",
        "cynT_betaCA": "cynT (β-CA)", "CA_generic": "CA generic",
        "Ca_transporter": "Ca transporter", "Ca_ATPase": "Ca ATPase",
        "CO3_transporter": "CO₃ transporter", "Na_H_antiporter_Mrp": "Mrp Na⁺/H⁺",
        "K_transport": "K⁺ transport"}


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


# ------------------------------------------------------------------ data: panel A
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

# ------------------------------------------------------------------ data: panel B
dos = pd.read_csv(SUPP / "Table_S15b_stoichiometry_per_MAG.csv")
dos = dos[dos.group == "MICP_complete"].set_index("MAG").loc[HEROES]
recomputed_ca = ((dos.Ca_transporter > 0) | (dos.Ca_ATPase > 0)).astype(int)
assert (recomputed_ca == dos.Ca_pathway).all()   # the stored call is Ca transporter/ATPase
mat = dos[COLS].values.astype(int)

# ------------------------------------------------------------------ data: panel C
gsum = pd.read_csv(SUPP / "Table_S17a_genomad_summary_per_MAG.csv")
assert sorted(gsum.MAG[gsum.group == "MICP_complete"]) == sorted(HEROES)
ov = pd.read_csv(SUPP / "Table_S17b_ureCah_vs_MGE_overlap.csv").set_index("MAG").loc[HEROES]


def n_contigs(cell):
    return 0 if pd.isna(cell) else len(str(cell).split(","))


def split(cell):
    return [] if pd.isna(cell) else [c.strip() for c in str(cell).split(",")]


squares = {}       # MAG -> list of (group, contig, flagged_as)
for mag in HEROES:
    r = ov.loc[mag]
    flagged = {c: "plasmid" for c in split(r.urease_on_plasmid) + split(r.CA_on_plasmid)}
    flagged.update({c: "virus" for c in split(r.urease_on_virus) + split(r.CA_on_virus)})
    ure = [("urease", c, flagged.get(c)) for c in split(r.urease_core_contigs)]
    ca = [("CA", c, flagged.get(c)) for c in split(r.CA_contigs)]
    assert sum(f is not None for _, _, f in ure) == int(r.urease_core_MGE_contamination)
    assert sum(f is not None for _, _, f in ca) == int(r.CA_MGE_contamination)
    squares[mag] = ure + ca
n_ure_max = max(sum(g == "urease" for g, _, _ in v) for v in squares.values())
n_ca_max = max(sum(g == "CA" for g, _, _ in v) for v in squares.values())
assert not any(f == "virus" for v in squares.values() for _, _, f in v)
# the per-MAG plasmid/virus counts of the two tables must agree
assert (ov.n_plasmid_contigs.values ==
        gsum.set_index("MAG").loc[HEROES].n_plasmid_contigs.values).all()

# ------------------------------------------------------------------ page geometry
FS_A_BODY, FS_A_STAT = FS_BODY * 1.3, FS_STAT * 1.3      # panel A type, 1.3x the page
FS_A_TITLE = st.FS_TITLE * 1.3
FS_A_TAB = 9.5                                            # A's stat block, a step larger

ROW_H, GAP = 11.5, 5.0           # GAP holds a row of 9 pt tick labels
TOP, LEFT, PLOT_W = 12.0, 26.0, 112.0
STAT_X = LEFT + PLOT_W + 4.0                              # left edge of the stat block
STAT_C = [STAT_X + 8.0, STAT_X + 22.0, STAT_X + 32.0]     # column centres (A)
A_END = TOP + len(order) * (ROW_H + GAP)
A_LEG_Y = A_END + 5.5            # gene key of panel A
B_LET_Y = A_LEG_Y + 9.0
B_TOP = B_LET_Y + 7.0
B_H = 24.0
C_LET_Y = B_TOP + B_H + 28.0     # clears the rotated column labels of B
C_TOP = C_LET_Y + 7.0
C_H = 29.0
H = C_TOP + C_H + 13.0
assert H <= 235.0, H             # single-page ceiling

fig, ax_mm, text_mm, letter = st.page(H)


def fy(y_mm):
    return 1.0 - y_mm / H


# ================================================================== A: synteny
letter(4.0, 6.0, "A")

# stat column headers (table column headers, not a title)
for xc, head in zip(STAT_C, ("ure genes", "MGE", "Δ GC")):
    text_mm(xc, TOP - 5.2, head, fontsize=FS_A_TAB, ha="center", color=TEXT)

LEVELS = (0.42, 1.60, 2.78)  # label rows above the track (one 8.5 pt line apart at
                             # ROW_H), used only when labels collide


def place_levels(items, kb_per_mm):
    """Assign each gene label the lowest row that clears the previous label on that row.

    Label width is estimated from the character count at the panel's label size,
    converted to kb through the track's own scale, so a tight three-gene cluster inside
    a 32 kb window staggers instead of printing on top of itself.
    """
    last = [-1e9] * len(LEVELS)
    out = []
    for xc, text in items:
        half = len(text) * FS_A_STAT * 0.353 * 0.52 * kb_per_mm / 2
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
        ax.add_patch(FancyArrow(xs, 0, dx, 0, width=0.50, head_width=0.72,
                                head_length=head, length_includes_head=True,
                                facecolor=col, edgecolor=AXIS, lw=0.35))
    for (xc, text), lv in zip(labelled, levels):
        y = LEVELS[lv]
        ax.plot([xc, xc], [0.24, y - 0.04], color=GREY, lw=0.4, zorder=1)
        ax.text(xc, y, text, ha="center", va="bottom", fontsize=FS_A_STAT,
                style="italic", color=TEXT)
    # each track is scaled to its own window: a shared scale would squeeze the five
    # short clusters into a quarter of the page beside the 32 kb M1 window
    ax.set_xlim(-0.3, t["length_kb"] + 0.3)
    ax.set_ylim(-0.75, LEVELS[-1] + 0.75)
    ax.set_yticks([])
    step = max(1.0, round(t["length_kb"] / 6))
    ax.set_xticks(np.arange(0, t["length_kb"] + 0.01, step))
    ax.tick_params(axis="x", labelsize=FS_A_BODY)
    st.style_axis(ax, left=False)
    if mag == order[-1]:
        ax.set_xlabel("Position in the drawn window (kb)", fontsize=FS_A_TITLE)
    # bare identifiers to the left of each track
    text_mm(LEFT - 2.0, top + ROW_H / 2 - 2.2, mag, ha="right", va="center",
            fontsize=FS_A_BODY, fontweight="bold", color=TEXT)
    text_mm(LEFT - 2.0, top + ROW_H / 2 + 1.8, t["ctg"], ha="right", va="center",
            fontsize=FS_A_STAT, color=GREY)
    # stat columns bound to the row
    for xc, val in zip(STAT_C, (f"{t['n_ure']} / 7", f"{t['mge']}", f"{t['dgc']:.1f}")):
        text_mm(xc, top + ROW_H / 2 - 0.2, val, ha="center", va="center",
                fontsize=FS_A_TAB, color=TEXT)

# cah is deliberately absent from the key: Table S1c records cah_on_main_contig = False
# for all six MAGs, so no cah arrow is drawn and a key entry would have no mark
drawn = [g for g in GENES if any((t["win"].cls == g).any() for t in tracks.values())]
handles = [Patch(facecolor=GENE_COL[g], edgecolor=AXIS, lw=0.35, label=g) for g in drawn] + \
          [Patch(facecolor=GENE_COL["other"], edgecolor=AXIS, lw=0.35, label="other CDS")]
fig.legend(handles=handles, loc="upper center", ncol=9, fontsize=FS_A_STAT, frameon=False,
           bbox_to_anchor=((LEFT + PLOT_W / 2) / st.PAGE_W_MM, fy(A_LEG_Y)),
           handlelength=1.1, columnspacing=1.1, handletextpad=0.4)

# ================================================================== B: gene dosage
letter(4.0, B_LET_Y, "B")
axB = ax_mm(LEFT, B_TOP, PLOT_W, B_H)          # the same width as the tracks of A
axB.imshow(mat, cmap=st.seq_cmap("copies", hi=GREEN), vmin=0, vmax=mat.max(),
           aspect="auto")
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        axB.text(j, i, mat[i, j], ha="center", va="center", fontsize=FS_A_BODY,
                 color="white" if mat[i, j] > mat.max() / 2 else TEXT)
axB.set_xticks(range(len(COLS)))
axB.set_xticklabels([NICE.get(c, c) for c in COLS], fontsize=FS_A_BODY, rotation=90)
axB.set_yticks(range(len(HEROES)))
axB.set_yticklabels(HEROES, fontsize=FS_A_BODY)
for tick, mag in zip(axB.get_yticklabels(), HEROES):
    tick.set_color(hero_col(mag))
axB.tick_params(length=0)
for s_ in axB.spines.values():
    s_.set_visible(False)

# stat columns bound to the rows, in the same right-hand block as A's stat columns
B_STAT_C = [STAT_X + 5.0, STAT_X + 16.0]
for xc, head, vals in zip(B_STAT_C, ("urease\ncore", "Ca\npathway"),
                          (dos.urease_core_complete, dos.Ca_pathway)):
    text_mm(xc, B_TOP - 7.5, head, fontsize=FS_A_STAT, ha="center", va="center",
            color=TEXT)
    for i, v in enumerate(vals):
        text_mm(xc, B_TOP + (i + 0.5) * B_H / len(HEROES), str(int(v)),
                fontsize=FS_A_STAT, ha="center", va="center", color=TEXT)

cax = ax_mm(STAT_X + 23.5, B_TOP + 1.0, 2.4, B_H - 2.0)
cb = fig.colorbar(axB.images[0], cax=cax)
cb.set_label("gene copies", fontsize=FS_A_BODY)
cb.ax.tick_params(labelsize=FS_A_BODY, length=2)
cb.outline.set_visible(False)

# ================================================================== C: MGE cross-check
letter(4.0, C_LET_Y, "C")
SQ = 0.78                        # square side in slot units
GAP_SLOT = 1.6                   # empty slots between the urease and the CA group
n_slots = n_ure_max + GAP_SLOT + n_ca_max
W_SQ = 68.0                      # width of the square block (mm)
axC = ax_mm(LEFT, C_TOP, W_SQ, C_H)
axC.set_xlim(-0.6, n_slots + 0.2)
axC.set_ylim(len(HEROES) - 0.4, -1.5)
axC.axis("off")
FLAG_COL = {None: GENE_COL["other"], "plasmid": HERO, "virus": AXIS}
for i, mag in enumerate(HEROES):
    k_ure = k_ca = 0
    for grp, ctg, flag in squares[mag]:
        if grp == "urease":
            x = k_ure
            k_ure += 1
        else:
            x = n_ure_max + GAP_SLOT + k_ca
            k_ca += 1
        axC.add_patch(Rectangle((x - SQ / 2, i - SQ / 2), SQ, SQ,
                                facecolor=FLAG_COL[flag], edgecolor=AXIS, lw=0.35))
    axC.text(-0.55, i, mag, ha="right", va="center", fontsize=FS_A_BODY,
             color=hero_col(mag))
# short block headers over the two contig groups
axC.text((n_ure_max - 1) / 2, -0.95, "urease-core contigs", ha="center", va="center",
         fontsize=FS_A_STAT, color=TEXT)
axC.text(n_ure_max + GAP_SLOT + (n_ca_max - 1) / 2, -0.95, "CA contigs", ha="center",
         va="center", fontsize=FS_A_STAT, color=TEXT)
for t in axC.texts:
    t.set_clip_on(False)

# geNomad-flagged contig totals per MAG, as bars on the same rows
axC2 = ax_mm(LEFT + W_SQ + 8.0, C_TOP, PLOT_W - W_SQ - 8.0, C_H)
yb = np.arange(len(HEROES))
hb = 0.36
axC2.barh(yb - hb / 2, ov.n_plasmid_contigs.values, hb, color=AXIS, label="plasmid")
axC2.barh(yb + hb / 2, ov.n_virus_contigs.values, hb, facecolor="white", edgecolor=AXIS,
          lw=0.5, hatch="////", label="virus")
xmax_c = float(max(ov.n_plasmid_contigs.max(), ov.n_virus_contigs.max()))
# a zero count has no bar and needs no label; labelling it would print on top of the
# label of the other bar of the same row
for yy, v in zip(yb - hb / 2, ov.n_plasmid_contigs.values):
    if v > 0:
        axC2.text(v + xmax_c * 0.02, yy, f"{int(v)}", ha="left", va="center",
                  fontsize=FS_A_STAT)
for yy, v in zip(yb + hb / 2, ov.n_virus_contigs.values):
    if v > 0:
        axC2.text(v + xmax_c * 0.02, yy, f"{int(v)}", ha="left", va="center",
                  fontsize=FS_A_STAT)
axC2.set_ylim(len(HEROES) - 0.4, -1.5)
axC2.set_yticks([])
axC2.set_xlim(0, xmax_c * 1.18)
axC2.set_xlabel("geNomad-flagged contigs per MAG", fontsize=FS_A_TITLE)
axC2.tick_params(axis="x", labelsize=FS_A_BODY)
st.style_axis(axC2, left=False)
axC2.legend(handles=[Patch(facecolor=AXIS, label="plasmid-flagged"),
                     Patch(facecolor="white", edgecolor=AXIS, hatch="////",
                           label="virus-flagged")],
            loc="upper right", bbox_to_anchor=(1.0, 1.06), ncol=1, fontsize=FS_A_STAT,
            handlelength=1.1, handleheight=0.9, borderpad=0.2, labelspacing=0.3)
# key for the squares, on the panel-letter row so it sits above the block it explains
fig.legend(handles=[Patch(facecolor=GENE_COL["other"], edgecolor=AXIS, lw=0.35,
                          label="contig, not MGE-flagged"),
                    Patch(facecolor=HERO, edgecolor=AXIS, lw=0.35,
                          label="contig flagged as plasmid (" + ", ".join(
                              f"{m} {c}" for m in HEROES for _, c, f in squares[m]
                              if f is not None) + ")")],
           loc="upper left", bbox_to_anchor=((LEFT + 2.0) / st.PAGE_W_MM, fy(C_LET_Y - 0.5)),
           ncol=2, fontsize=FS_A_STAT, frameon=False, handlelength=1.1, columnspacing=1.2,
           handletextpad=0.4)

print(f"  page height {H:.1f} mm")
st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig2")
