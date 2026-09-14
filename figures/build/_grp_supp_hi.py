"""Plotting helpers shared by build_supS16.py .. build_supS21.py.

Layout and drawing only - this module holds no data value.  Kept out of `build_*.py`
so that `_job/check_no_hardcoded.py` scans the scripts that carry the science.
"""

import numpy as np

import _style as st
from _style import HERO, REST, GREY, TEXT, LIGHT, FS_BODY, FS_STAT

JITTER_SEED = 20260904  # fixed so the point cloud is identical on every rebuild
SUPP = "/data/data/Upcycling/SUBMISSION/Supplementary_tables"
ADDITIONAL = "/data/data/Upcycling/research/additional"


SUPERSCRIPT = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")


def fmt_p(p):
    """P value as it should read on a figure: 2 significant digits, sci notation below 1e-3.

    The exponent uses Unicode superscript characters rather than matplotlib mathtext.
    Mathtext is typeset in the mathtext font and is emitted into the SVG with its own
    font-family, which breaks the one-font Arial-first rule that _job/verify_outputs.py
    enforces; Liberation Sans and Arial both carry these superscript glyphs.
    """
    if p >= 1e-3:
        return f"P = {p:.3g}"
    m, e = f"{p:.1e}".split("e")
    return f"P = {m} × 10{str(int(e)).translate(SUPERSCRIPT)}"


def strip_box(ax, groups, colours, jitter_w=0.28, box_w=0.42, point_s=9, log=False,
              value_fmt=None):
    """Box + jittered points, one column per group.

    `groups` is a list of (label, values).  Returns (column x positions, per-group arrays of
    the jittered x coordinate of every input value) so that a later overlay can be drawn at
    the same x as the point it marks.
    The box shows median and inter-quartile range; whiskers span the full range.
    """
    rng = np.random.default_rng(JITTER_SEED)
    xs = np.arange(len(groups), dtype=float)
    jitters = []
    for x, (label, vals), col in zip(xs, groups, colours):
        v = np.asarray(vals, dtype=float)
        finite = np.isfinite(v)
        v = v[finite]
        if len(v) == 0:
            jitters.append(np.full(len(finite), np.nan))
            continue
        bp = ax.boxplot([v], positions=[x], widths=box_w, showfliers=False,
                        patch_artist=True, whis=(0, 100), zorder=1)
        bp["boxes"][0].set(facecolor="white", edgecolor=col, linewidth=0.8)
        for key in ("whiskers", "caps"):
            for a in bp[key]:
                a.set(color=col, linewidth=0.8)
        bp["medians"][0].set(color=col, linewidth=1.4)
        jx = x + rng.uniform(-jitter_w, jitter_w, len(v))
        full = np.full(len(finite), np.nan)
        full[finite] = jx
        jitters.append(full)
        ax.scatter(jx, v, s=point_s, facecolor=col, edgecolor="none", alpha=0.55,
                   zorder=2, clip_on=True)
    ax.set_xticks(xs)
    ax.set_xticklabels([g[0] for g in groups])
    if log:
        ax.set_yscale("log")
    st.style_axis(ax)
    ax.set_xlim(xs[0] - 0.65, xs[-1] + 0.65)
    return xs, jitters


def stat_bracket(ax, x0, x1, y, text, colour=TEXT, drop=0.0):
    """A comparison bracket with a short stat string above it."""
    ax.plot([x0, x0, x1, x1], [y - drop, y, y, y - drop], lw=0.7, color=colour,
            clip_on=False, solid_capstyle="butt")
    ax.text((x0 + x1) / 2, y, text, ha="center", va="bottom", fontsize=FS_STAT,
            color=colour, clip_on=False)


def group_counts(ax, xs, groups, y, fmt="n = {n}"):
    """Sample size under each column, bound to the column."""
    for x, (_, vals) in zip(xs, groups):
        n = int(np.sum(np.isfinite(np.asarray(vals, dtype=float))))
        ax.text(x, y, fmt.format(n=n), ha="center", va="top", fontsize=FS_STAT,
                color=TEXT, transform=ax.get_xaxis_transform(), clip_on=False)


def hero_rest_groups(df, col, hero_mask, hero_label="MICP-complete", rest_label="Rest"):
    return ([(hero_label, df.loc[hero_mask, col].to_numpy()),
             (rest_label, df.loc[~hero_mask, col].to_numpy())],
            [HERO, REST])


def paired_hbars(ax, rows, hero_vals, rest_vals, bar_h=0.36):
    """Horizontal paired bars, MICP-complete above rest, one pair per row."""
    y = np.arange(len(rows), dtype=float)
    ax.barh(y - bar_h / 2, hero_vals, height=bar_h, color=HERO, edgecolor="none",
            zorder=2, label="MICP-complete (n = 6)")
    ax.barh(y + bar_h / 2, rest_vals, height=bar_h, color=REST, edgecolor="none",
            zorder=2, label="Rest (n = 105)")
    ax.set_yticks(y)
    ax.set_yticklabels(rows)
    ax.invert_yaxis()
    ax.grid(axis="x", color=LIGHT, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    st.style_axis(ax, left=False)
    ax.tick_params(left=False)
    return y


def value_labels(ax, xs, vals, fmt="{v:.2f}", dy=0.0, fontsize=FS_STAT, colour=TEXT,
                 ha="center", va="bottom"):
    for x, v in zip(xs, vals):
        ax.text(x, v + dy, fmt.format(v=v), ha=ha, va=va, fontsize=fontsize,
                color=colour, clip_on=False)


def threshold_line(ax, y, label, colour=GREY, x=0.99):
    ax.axhline(y, color=colour, lw=0.8, ls="--", zorder=1)
    ax.text(x, y, label, transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=FS_STAT, color=colour)
