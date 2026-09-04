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
     geNomad mobile-element calls                                         [old Fig S15B]

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
  C  the contig counts are recomputed by splitting the comma-joined contig lists of S17b;
     the two contamination columns are read as stored.  The flag columns
     urease_on_plasmid / urease_on_virus / CA_on_virus are empty for every MAG in the
     source and carry no information, so they are not drawn; CA_on_plasmid is non-empty
     for S26 alone and is represented by that MAG's CA contamination count.  The per-MAG
     plasmid contig counts of S17a and S17b are asserted equal.

Colour meanings on this page (one meaning per colour):
  gene identity palette (A)  one hue per MICP gene, grey for any other CDS in the window
  green intensity (B)        gene copy number; the stat columns of B are text, not
                             colour, so nothing competes with the copy-number scale
  blue / orange (C)          Sphingobacterium and Pseudomonas_E MAG identifiers
  coral (C)                  a non-zero mobile-element contamination count
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrow, Patch

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


tab_cols = [("urease contigs", [n_contigs(v) for v in ov.urease_core_contigs]),
            ("CA contigs", [n_contigs(v) for v in ov.CA_contigs]),
            ("plasmid contigs", list(ov.n_plasmid_contigs)),
            ("virus contigs", list(ov.n_virus_contigs)),
            ("urease on MGE", list(ov.urease_core_MGE_contamination)),
            ("CA on MGE", list(ov.CA_MGE_contamination))]
# the per-MAG plasmid/virus counts of the two tables must agree
assert (ov.n_plasmid_contigs.values ==
        gsum.set_index("MAG").loc[HEROES].n_plasmid_contigs.values).all()

# ------------------------------------------------------------------ page geometry
ROW_H, GAP = 11.5, 4.5
TOP, LEFT, PLOT_W = 12.0, 26.0, 118.0
STAT_X = LEFT + PLOT_W + 6.0
A_END = TOP + len(order) * (ROW_H + GAP)
A_LEG_Y = A_END + 5.0            # gene key of panel A
B_LET_Y = A_LEG_Y + 10.0
B_TOP = B_LET_Y + 8.0
B_H = 22.0
C_LET_Y = B_TOP + B_H + 29.0     # clears the rotated column labels of B
C_TOP = C_LET_Y + 7.0
C_H = 34.0
H = C_TOP + C_H + 3.0
assert H <= 235.0, H             # single-page ceiling

fig, ax_mm, text_mm, letter = st.page(H)


def fy(y_mm):
    return 1.0 - y_mm / H


# ================================================================== A: synteny
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
fig.legend(handles=handles, loc="upper center", ncol=9, fontsize=FS_STAT, frameon=False,
           bbox_to_anchor=(0.5, fy(A_LEG_Y)), handlelength=1.1, columnspacing=1.1,
           handletextpad=0.4)

# ================================================================== B: gene dosage
letter(4.0, B_LET_Y, "B")
axB = ax_mm(24.0, B_TOP, 104.0, B_H)
axB.imshow(mat, cmap=st.seq_cmap("copies", hi=GREEN), vmin=0, vmax=mat.max(),
           aspect="auto")
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        axB.text(j, i, mat[i, j], ha="center", va="center", fontsize=FS_BODY,
                 color="white" if mat[i, j] > mat.max() / 2 else TEXT)
axB.set_xticks(range(len(COLS)))
axB.set_xticklabels([NICE.get(c, c) for c in COLS], fontsize=FS_BODY, rotation=90)
axB.set_yticks(range(len(HEROES)))
axB.set_yticklabels(HEROES)
for tick, mag in zip(axB.get_yticklabels(), HEROES):
    tick.set_color(hero_col(mag))
axB.tick_params(length=0)
for s in axB.spines.values():
    s.set_visible(False)

# stat columns bound to the rows
for dx, head, vals in ((len(COLS) + 0.4, "urease\ncore", dos.urease_core_complete),
                       (len(COLS) + 2.6, "Ca\npathway", dos.Ca_pathway)):
    axB.text(dx, -1.15, head, fontsize=FS_STAT, ha="left", va="center", color=TEXT)
    for i, v in enumerate(vals):
        axB.text(dx, i, str(int(v)), fontsize=FS_STAT, ha="left", va="center", color=TEXT)
for t in axB.texts:
    t.set_clip_on(False)

cb = fig.colorbar(axB.images[0], ax=axB, fraction=0.03, pad=0.22, aspect=12)
cb.set_label("gene copies", fontsize=FS_BODY)
cb.ax.tick_params(labelsize=FS_BODY, length=2)
cb.outline.set_visible(False)

# ================================================================== C: MGE cross-check
letter(4.0, C_LET_Y, "C")
axC = ax_mm(50.0, C_TOP, 84.0, C_H)
axC.set_xlim(0, len(tab_cols))
axC.set_ylim(len(HEROES), -1.4)
axC.axis("off")
for j, (head, vals) in enumerate(tab_cols):
    axC.text(j + 0.5, -1.0, head.replace(" ", "\n", 1), ha="center", va="center",
             fontsize=FS_STAT, color=TEXT)
    for i, v in enumerate(vals):
        axC.text(j + 0.5, i, str(int(v)), ha="center", va="center", fontsize=FS_BODY,
                 color=HERO if (int(v) > 0 and "on MGE" in head) else TEXT)
for i, mag in enumerate(HEROES):
    axC.text(-0.15, i, mag, ha="right", va="center", fontsize=FS_BODY,
             color=hero_col(mag))
    axC.axhline(i + 0.5, color=LIGHT, lw=0.5, xmin=0.0, xmax=1.0)
axC.axhline(-0.5, color=AXIS, lw=0.7)
for t in axC.texts:
    t.set_clip_on(False)

print(f"  page height {H:.1f} mm")
st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig2")
