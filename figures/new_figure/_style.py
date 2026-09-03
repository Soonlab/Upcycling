"""Shared style for the Upcycling SVG figure rebuild (main Fig 1-8, Suppl Fig S1-S21).

Method carried over from the Rectal_Organoid svg_supp rebuild (2026-08-20/21) and the
SNUH_KMJ new_figure rebuild (2026-09-03/04):

  * matplotlib, one composed page per figure, laid out in millimetres
  * SVG keeps text as text (svg.fonttype='none'); after saving, every font-family
    declaration is rewritten to "Arial, 'Liberation Sans', 'DejaVu Sans', sans-serif"
  * PDF embeds TrueType (fonttype 42); PNG at 200 dpi
  * uniform type scale: body / tick labels 7 pt, axis titles 8 pt,
    panel letters 10 pt bold, stat annotations 6.5 pt
  * spines left+bottom only
  * NO panel titles and NO explanatory sentences on the figure itself - those belong
    to the manuscript legend.  Only element-bound labels survive: axis and tick labels,
    legends, per-datapoint value labels and stat columns, short block headers that
    separate data blocks, table column headers, and direction cues.
  * panel letters left-aligned on a small number of shared x positions
  * no data value is hard-coded in a build script; everything is loaded or recomputed
    from a repository source file (style constants and category names excepted)

Palette (one colour, one meaning per page):
  HERO / REST      MICP-complete group (coral) vs the remaining 105 MAGs (grey)
  SPHINGO / PSEUDO the two convergent lineages, Sphingobacterium (blue) and
                   Pseudomonas_E (orange) - the convention already used by main Fig 8
  SOURCE           waste source, four categorical colours that do not reuse any of the
                   above (cattle brown, swine pink, sheep green, poultry purple)
  GREEN            "present / match / complete" in presence-absence and match heat maps
  SIG_*            significance tiers on forest plots (q<0.05 coral, q<0.10 light coral,
                   n.s. grey) - the same coral as HERO because both mean "enriched in
                   the MICP-complete group"
Sequential heat maps use seq_cmap() (white -> GREEN) so that blue and orange stay
reserved for the two lineages on every page.
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------------------------------------------ palette
HERO     = "#C53E1F"   # MICP-complete group (n = 6)
REST     = "#8A8A8A"   # remaining 105 MAGs
SPHINGO  = "#1F4E78"   # Sphingobacterium lineage (S13, S16, S23, C22)
PSEUDO   = "#D9772B"   # Pseudomonas_E lineage (M1, S26)
GREEN    = "#2F7D4F"   # present / match / complete
GREY     = "#6E6E6E"   # neutral reference lines, thresholds
LIGHT    = "#E6E6E6"   # very light rules and grids
TEXT     = "#222222"   # all text and panel letters
AXIS     = "#444444"
HERO_LT  = "#EBC5BC"   # light coral (q < 0.10 tier, secondary hero series)
SPHINGO_LT = "#A6C8E0"
PSEUDO_LT  = "#F2CBA5"
SOURCE = {"cattle": "#8E6C3A", "swine": "#D46A9E", "sheep": "#2E8B57", "poultry": "#7B4F9D"}
SIG_05, SIG_10, SIG_NS = HERO, HERO_LT, REST

# the six MICP-complete MAGs and their lineage (category names, not data)
HEROES = ["S13", "S16", "S23", "C22", "M1", "S26"]
LINEAGE = {"S13": "Sphingobacterium", "S16": "Sphingobacterium", "S23": "Sphingobacterium",
           "C22": "Sphingobacterium", "M1": "Pseudomonas_E", "S26": "Pseudomonas_E"}
LINEAGE_COL = {"Sphingobacterium": SPHINGO, "Pseudomonas_E": PSEUDO}


def hero_col(mag):
    return LINEAGE_COL[LINEAGE[mag]]


FS_BODY  = 7.0
FS_TITLE = 8.0
FS_PANEL = 10.0
FS_STAT  = 6.5

MM = 1 / 25.4  # mm -> inch
PAGE_W_MM = 180.0  # two-column journal width


def seq_cmap(name="white_green", hi=GREEN, lo="#FFFFFF"):
    return LinearSegmentedColormap.from_list(name, [lo, hi])


def diverging_cmap(name="blue_orange", lo=SPHINGO, hi=PSEUDO, mid="#F7F7F5"):
    return LinearSegmentedColormap.from_list(name, [lo, mid, hi])


def setup():
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # Liberation Sans is metric-compatible with Arial; DejaVu covers stray glyphs
        "font.family": ["Liberation Sans", "DejaVu Sans"],
        "font.size": FS_BODY,
        "axes.labelsize": FS_TITLE,
        "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_BODY,
        "ytick.labelsize": FS_BODY,
        "legend.fontsize": FS_BODY,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": TEXT,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "text.color": TEXT,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "legend.frameon": False,
    })


def style_axis(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    if not left:
        ax.tick_params(left=False)
    if not bottom:
        ax.tick_params(bottom=False)


def page(height_mm, width_mm=PAGE_W_MM):
    """Create the figure plus mm-based placement helpers.

    Returns (fig, ax_mm, text_mm, letter).  ax_mm(l, t, w, h) places an Axes by its
    top-left corner in millimetres from the top-left of the page; text_mm(x, y, s)
    places figure text at a mm coordinate; letter(x, y, 'A') places a panel letter.
    """
    fw, fh = width_mm * MM, height_mm * MM
    fig = plt.figure(figsize=(fw, fh))

    def ax_mm(l, t, w, h):
        return fig.add_axes([l * MM / fw, 1 - (t + h) * MM / fh,
                             w * MM / fw, h * MM / fh])

    def text_mm(x, y, s, **kw):
        kw.setdefault("fontsize", FS_BODY)
        kw.setdefault("color", TEXT)
        kw.setdefault("va", "top")
        return fig.text(x * MM / fw, 1 - y * MM / fh, s, **kw)

    def letter(x, y, s):
        return fig.text(x * MM / fw, 1 - y * MM / fh, s, fontsize=FS_PANEL,
                        fontweight="bold", color=TEXT, ha="left", va="top")

    return fig, ax_mm, text_mm, letter


# any family matplotlib may emit (Liberation Sans, DejaVu Sans from mathtext, or a
# fallback list) is rewritten to the Arial-first stack
_FONT_RE = re.compile(
    r"font-family:\s*(?:'?Liberation Sans'?|'?DejaVu Sans'?|sans-serif)"
    r"(?:,\s*(?:'?Liberation Sans'?|'?DejaVu Sans'?|sans-serif))*")


def save(fig, outdir, stem, dpi=200):
    """Save SVG + PDF + PNG; rewrite the SVG font-family to Arial-first."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    svg = outdir / f"{stem}.svg"
    fig.savefig(svg)
    fig.savefig(outdir / f"{stem}.pdf")
    fig.savefig(outdir / f"{stem}.png", dpi=dpi)
    txt = svg.read_text()
    txt = _FONT_RE.sub(
        "font-family: Arial, 'Liberation Sans', 'DejaVu Sans', sans-serif", txt)
    svg.write_text(txt)
    plt.close(fig)
    print(f"saved {stem}.svg/.pdf/.png in {outdir}")


# ------------------------------------------------------------------ self-audit
def _bbox(artist, renderer):
    """Window extent of an artist.

    For a Text carrying a bbox patch (the white boxes used in node diagrams), the drawn
    footprint is the PATCH, which is larger than the glyph extent that
    Text.get_window_extent reports.  Measuring only the glyphs lets two boxed labels sit
    visibly on top of each other while the checker reports no overlap.
    """
    try:
        bb = artist.get_window_extent(renderer=renderer)
    except Exception:
        return None
    patch = getattr(artist, "get_bbox_patch", lambda: None)()
    if patch is not None:
        try:
            artist.update_bbox_position_size(renderer)
            pb = patch.get_window_extent(renderer=renderer)
            if pb.width > 0 and pb.height > 0:
                bb = bb.union([bb, pb])
        except Exception:
            pass
    return bb if bb.width > 0 and bb.height > 0 else None


def audit(fig, pad=1.0, verbose=True):
    """Geometric self-audit of a composed page.

    Uses the real matplotlib renderer, so a label that wraps or overflows its slot is
    measured at its true size - the failure mode that a shape-box checker misses.

    Reports: text-vs-text overlaps, text running outside the canvas, and Axes that
    overlap each other.  Returns a dict; empty lists mean the page passes.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()

    texts = []
    for t in fig.findobj(match=lambda o: o.__class__.__name__ == "Text"):
        s = (t.get_text() or "").strip()
        if not s or not t.get_visible():
            continue
        bb = _bbox(t, r)
        if bb is None:
            continue
        texts.append((s, bb))

    overlaps = []
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = texts[i][1], texts[j][1]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > pad and oy > pad:
                overlaps.append((texts[i][0][:42], texts[j][0][:42],
                                 round(ox, 1), round(oy, 1)))

    outside = [(s[:42], round(bb.x0, 1), round(bb.x1, 1), round(bb.y0, 1), round(bb.y1, 1))
               for s, bb in texts
               if bb.x0 < -pad or bb.y0 < -pad or bb.x1 > W + pad or bb.y1 > H + pad]

    axes = [(ax, _bbox(ax, r)) for ax in fig.axes]
    ax_overlaps = []
    for i in range(len(axes)):
        for j in range(i + 1, len(axes)):
            a, b = axes[i][1], axes[j][1]
            if a is None or b is None:
                continue
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > pad and oy > pad:
                ax_overlaps.append((i, j, round(ox, 1), round(oy, 1)))

    res = dict(text_overlaps=overlaps, text_outside=outside, axes_overlaps=ax_overlaps,
               n_text=len(texts))
    if verbose:
        n = len(overlaps) + len(outside) + len(ax_overlaps)
        print(f"  audit: {len(texts)} text items | overlaps {len(overlaps)} | "
              f"outside {len(outside)} | axes overlaps {len(ax_overlaps)} "
              f"-> {'PASS' if n == 0 else 'FAIL'}")
        for o in overlaps[:12]:
            print(f"    overlap: {o[0]!r} x {o[1]!r}  ({o[2]}x{o[3]} px)")
        for o in outside[:12]:
            print(f"    outside: {o[0]!r} x[{o[1]},{o[2]}] y[{o[3]},{o[4]}]")
        for o in ax_overlaps[:12]:
            print(f"    axes {o[0]} x axes {o[1]} overlap ({o[2]}x{o[3]} px)")
    return res


_SENTENCE = re.compile(r"[A-Za-z]{3,}[^.;]*[.;]")


def prose_scan(fig, max_words=9):
    """Flag figure text that reads as a sentence rather than an element-bound label.

    The 2026-08-20/21 review rounds removed every explanatory sentence, method
    paragraph and stat footnote from the figures.  This catches a regression.
    """
    fig.canvas.draw()
    hits = []
    for t in fig.findobj(match=lambda o: o.__class__.__name__ == "Text"):
        s = (t.get_text() or "").strip()
        if not s or not t.get_visible():
            continue
        words = s.split()
        if len(words) > max_words and _SENTENCE.search(s):
            hits.append(s[:80])
    if hits:
        print(f"  prose_scan: {len(hits)} sentence-like item(s) -> FAIL")
        for h in hits:
            print(f"    prose: {h!r}")
    else:
        print("  prose_scan: no sentence-like text -> PASS")
    return hits
