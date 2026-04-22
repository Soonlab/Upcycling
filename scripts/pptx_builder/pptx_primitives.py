"""
Higher-level chart primitives built on pptx_style.py.

Each primitive produces a complete native pptx chart (axes, data marks,
labels, legend) on a slide.  All text is editable textboxes; all marks are
shapes/connectors so the reviewer / co-author can double-click and edit
anything in place without ungrouping.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

from pptx.dml.color import RGBColor
from pptx.util import Inches, Pt

from pptx_style import (AxisSpec, OKABE_ITO, SPHINGO, PSEUDO, REST, HERO_MIX,
                        HERO_SET, SPHINGO_HEROES, PSEUDO_HEROES, SOURCE_COLOR,
                        add_text, add_line, add_rect, add_circle,
                        add_arrow_shape, data_to_plot_x, data_to_plot_y,
                        draw_axis, hero_color_for, nice_log_ticks,
                        panel_label, sig_star)


# =============================================================================
# Helper utilities
# =============================================================================

def median(seq: Sequence[float]) -> float:
    s = sorted(seq)
    n = len(s)
    if n == 0:
        return 0.0
    if n % 2 == 1:
        return s[n // 2]
    return 0.5 * (s[n // 2 - 1] + s[n // 2])


def quartiles(seq: Sequence[float]) -> tuple[float, float, float]:
    s = sorted(seq)
    n = len(s)
    if n == 0:
        return 0.0, 0.0, 0.0
    mid = n // 2
    q1 = median(s[:mid])
    q3 = median(s[mid + n % 2:])
    return q1, median(s), q3


def tukey_whiskers(seq: Sequence[float]) -> tuple[float, float]:
    """Return (lower_whisker, upper_whisker) at 1.5 IQR."""
    if not seq:
        return 0.0, 0.0
    q1, _, q3 = quartiles(seq)
    iqr = q3 - q1
    lo = q1 - 1.5 * iqr
    hi = q3 + 1.5 * iqr
    lower = min((v for v in seq if v >= lo), default=q1)
    upper = max((v for v in seq if v <= hi), default=q3)
    return lower, upper


# =============================================================================
# Bar chart (vertical)
# =============================================================================

def draw_vertical_bars(shapes, ax: AxisSpec, categories: Sequence[str],
                       values: Sequence[float], *,
                       colors: Optional[Sequence[RGBColor]] = None,
                       bar_outline: RGBColor = OKABE_ITO["black"],
                       value_labels: Optional[Sequence[str]] = None,
                       value_label_size: int = 8) -> None:
    n = len(categories)
    slot = ax.width / max(n, 1)
    bar_w = slot * 0.6
    for i, (cat, val) in enumerate(zip(categories, values)):
        col = (colors[i] if colors else OKABE_ITO["blue"])
        cx = ax.x0 + slot * (i + 0.5)
        x0 = cx - bar_w / 2
        y1 = data_to_plot_y(val, ax)
        y0 = data_to_plot_y(ax.ylim[0], ax)
        if y1 > y0:
            y1, y0 = y0, y1
        add_rect(shapes, x0, y1, bar_w, y0 - y1, fill=col,
                 outline=bar_outline, outline_weight_pt=0.4)
        if value_labels:
            add_text(shapes, x0 - 0.2, y1 - 0.25, bar_w + 0.4, 0.2,
                     value_labels[i], size=value_label_size, align="center")

    # X tick labels
    if ax.x_tick_labels is None:
        ax.x_tick_labels = list(categories)
        ax.x_ticks = [ax.xlim[0] + (i + 0.5) * (ax.xlim[1] - ax.xlim[0]) / max(n, 1)
                      for i in range(n)]
    draw_axis(shapes, ax)


# =============================================================================
# Horizontal bar chart
# =============================================================================

def draw_horizontal_bars(shapes, ax: AxisSpec, categories: Sequence[str],
                         values_a: Sequence[float],
                         values_b: Optional[Sequence[float]] = None,
                         *, color_a: RGBColor = OKABE_ITO["purple"],
                         color_b: RGBColor = REST,
                         label_a: str = "Series A",
                         label_b: str = "Series B",
                         annotations: Optional[Sequence[str]] = None,
                         bold_rows: Optional[set[str]] = None,
                         category_colors: Optional[dict[str, RGBColor]] = None
                         ) -> None:
    n = len(categories)
    slot = ax.height / max(n, 1)
    bar_h = slot * 0.35 if values_b is not None else slot * 0.6
    for i, cat in enumerate(categories):
        yc = ax.y0 + slot * (i + 0.5)
        x0 = ax.x0
        if values_b is not None:
            # series A on top (above centre)
            ya_top = yc - bar_h - 0.015
            width_a = data_to_plot_x(values_a[i], ax) - x0
            add_rect(shapes, x0, ya_top, max(width_a, 0.001), bar_h,
                     fill=color_a, outline=OKABE_ITO["black"],
                     outline_weight_pt=0.3)
            yb_top = yc + 0.015
            width_b = data_to_plot_x(values_b[i], ax) - x0
            add_rect(shapes, x0, yb_top, max(width_b, 0.001), bar_h,
                     fill=color_b, outline=OKABE_ITO["black"],
                     outline_weight_pt=0.3)
        else:
            y_top = yc - bar_h / 2
            width = data_to_plot_x(values_a[i], ax) - x0
            add_rect(shapes, x0, y_top, max(width, 0.001), bar_h,
                     fill=color_a, outline=OKABE_ITO["black"],
                     outline_weight_pt=0.3)

        # annotation
        if annotations and annotations[i]:
            ann_x = max(
                data_to_plot_x(values_a[i], ax),
                data_to_plot_x(values_b[i] if values_b else 0, ax),
            ) + 0.05
            add_text(shapes, ann_x, yc - 0.12, 2.0, 0.2, annotations[i],
                     size=8, color=color_a, bold=True)

        # category label
        lbl_color = (category_colors.get(cat) if category_colors
                     else OKABE_ITO["black"])
        lbl_bold = bool(bold_rows and cat in bold_rows)
        add_text(shapes, ax.x0 - 1.6, yc - 0.12, 1.55, 0.22,
                 cat, size=9, align="right", color=lbl_color,
                 bold=lbl_bold)

    draw_axis(shapes, ax, draw_x=True, draw_y=True)

    # Legend
    if values_b is not None:
        lx = ax.x0 + ax.width - 2.4
        ly = ax.y0 - 0.4
        add_rect(shapes, lx, ly, 0.18, 0.14, fill=color_a)
        add_text(shapes, lx + 0.22, ly - 0.02, 1.2, 0.2, label_a, size=9)
        add_rect(shapes, lx + 1.2, ly, 0.18, 0.14, fill=color_b)
        add_text(shapes, lx + 1.42, ly - 0.02, 1.0, 0.2, label_b, size=9)


# =============================================================================
# Grouped (paired) vertical bar
# =============================================================================

def draw_grouped_bars(shapes, ax: AxisSpec, categories: Sequence[str],
                      groups: Sequence[Sequence[float]], *,
                      group_labels: Sequence[str],
                      group_colors: Sequence[RGBColor],
                      value_label_digits: int = 2,
                      reference_line: Optional[tuple[float, str, RGBColor]] = None
                      ) -> None:
    n = len(categories)
    k = len(groups)
    slot = ax.width / max(n, 1)
    total_bar = slot * 0.7
    bar_w = total_bar / k
    for i in range(n):
        for j in range(k):
            col = group_colors[j]
            val = groups[j][i]
            cx = ax.x0 + slot * (i + 0.5) - total_bar / 2 + bar_w * (j + 0.5)
            x0 = cx - bar_w / 2
            y1 = data_to_plot_y(val, ax)
            y_base = data_to_plot_y(ax.ylim[0], ax)
            top, bot = sorted([y1, y_base])
            add_rect(shapes, x0, top, bar_w, bot - top, fill=col,
                     outline=OKABE_ITO["black"], outline_weight_pt=0.3)
            # value label on top
            add_text(shapes, x0 - 0.1, top - 0.22, bar_w + 0.2, 0.18,
                     f"{val:.{value_label_digits}f}", size=7,
                     align="center")

    # Category labels under bars
    if ax.x_tick_labels is None:
        ax.x_tick_labels = list(categories)
        ax.x_ticks = [ax.xlim[0] + (i + 0.5) * (ax.xlim[1] - ax.xlim[0]) / max(n, 1)
                      for i in range(n)]
    draw_axis(shapes, ax)

    # Reference line
    if reference_line is not None:
        ref_val, ref_label, ref_col = reference_line
        yref = data_to_plot_y(ref_val, ax)
        add_line(shapes, ax.x0, yref, ax.x0 + ax.width, yref,
                 color=ref_col, weight_pt=0.9, dash="dash")
        add_text(shapes, ax.x0 + ax.width - 2.2, yref - 0.22, 2.2, 0.2,
                 ref_label, size=9, color=ref_col, align="right", bold=True)

    # Legend at top-right (inside panel to avoid being clipped)
    k = len(group_labels)
    lx = ax.x0 + ax.width - 2.0
    ly = ax.y0 + 0.05
    for j, (lbl, col) in enumerate(zip(group_labels, group_colors)):
        add_rect(shapes, lx, ly + j * 0.34, 0.18, 0.14, fill=col)
        add_text(shapes, lx + 0.22, ly + j * 0.34 - 0.04, 1.75, 0.32, lbl,
                 size=9)


# =============================================================================
# Dual-axis bar (e.g., TM-score + RMSD)
# =============================================================================

def draw_dual_axis_bars(shapes, ax_left: AxisSpec, ax_right: AxisSpec,
                        categories: Sequence[str],
                        values_left: Sequence[float],
                        values_right: Sequence[float], *,
                        left_label: str, right_label: str,
                        color_left_fn,  # fn(cat) -> RGBColor
                        color_right: RGBColor = OKABE_ITO["black"],
                        reference_left: Optional[tuple[float, str]] = None
                        ) -> None:
    n = len(categories)
    slot = ax_left.width / max(n, 1)
    bar_w = slot * 0.3

    for i, cat in enumerate(categories):
        cx = ax_left.x0 + slot * (i + 0.5)
        # left bar
        x0l = cx - bar_w - 0.02
        y1l = data_to_plot_y(values_left[i], ax_left)
        y0l = data_to_plot_y(ax_left.ylim[0], ax_left)
        top, bot = sorted([y1l, y0l])
        add_rect(shapes, x0l, top, bar_w, bot - top, fill=color_left_fn(cat),
                 outline=OKABE_ITO["black"], outline_weight_pt=0.3)
        add_text(shapes, x0l - 0.1, top - 0.22, bar_w + 0.2, 0.18,
                 f"{values_left[i]:.3f}", size=7, align="center")
        # right bar (hatched via light fill + outline)
        x0r = cx + 0.02
        y1r = data_to_plot_y(values_right[i], ax_right)
        y0r = data_to_plot_y(ax_right.ylim[0], ax_right)
        top, bot = sorted([y1r, y0r])
        add_rect(shapes, x0r, top, bar_w, bot - top,
                 fill=OKABE_ITO["lightgrey"], outline=color_right,
                 outline_weight_pt=0.4)
        add_text(shapes, x0r - 0.1, top - 0.22, bar_w + 0.2, 0.18,
                 f"{values_right[i]:.2f}", size=7, align="center")

    # Left axis (draw_axis takes care of ticks)
    if ax_left.x_tick_labels is None:
        ax_left.x_tick_labels = list(categories)
        ax_left.x_ticks = [ax_left.xlim[0] + (i + 0.5) *
                           (ax_left.xlim[1] - ax_left.xlim[0]) / max(n, 1)
                           for i in range(n)]
    draw_axis(shapes, ax_left, draw_y=True, draw_x=True)

    # Right spine + right-axis ticks (manual because draw_axis only does left)
    add_line(shapes, ax_left.x0 + ax_left.width, ax_left.y0,
             ax_left.x0 + ax_left.width, ax_left.y0 + ax_left.height,
             weight_pt=0.75)
    if ax_right.y_ticks is not None:
        for tval, tlbl in zip(ax_right.y_ticks,
                              ax_right.y_tick_labels or
                              [f"{v:g}" for v in ax_right.y_ticks]):
            yp = data_to_plot_y(tval, ax_right)
            add_line(shapes, ax_left.x0 + ax_left.width, yp,
                     ax_left.x0 + ax_left.width + 0.06, yp, weight_pt=0.6)
            add_text(shapes, ax_left.x0 + ax_left.width + 0.08, yp - 0.1,
                     0.5, 0.2, tlbl, size=8, align="left")
    if ax_right.y_title:
        add_text(shapes, ax_left.x0 + ax_left.width + 0.6,
                 ax_left.y0 + ax_left.height / 2 - 0.15,
                 0.7, 0.3, ax_right.y_title, size=10, align="left",
                 bold=True, rotation=-90)

    # reference line
    if reference_left is not None:
        val, txt = reference_left
        yref = data_to_plot_y(val, ax_left)
        add_line(shapes, ax_left.x0, yref, ax_left.x0 + ax_left.width, yref,
                 color=OKABE_ITO["black"], weight_pt=0.9, dash="dash")
        add_text(shapes, ax_left.x0 + 0.05, yref - 0.22, 3.5, 0.2,
                 txt, size=9, color=OKABE_ITO["black"], bold=True)

    # Legend
    lx = ax_left.x0 + 0.1
    ly = ax_left.y0 - 0.4
    add_rect(shapes, lx, ly, 0.18, 0.14, fill=OKABE_ITO["blue"])
    add_text(shapes, lx + 0.22, ly - 0.02, 2.6, 0.2, left_label, size=9)
    add_rect(shapes, lx + 2.7, ly, 0.18, 0.14,
             fill=OKABE_ITO["lightgrey"], outline=color_right)
    add_text(shapes, lx + 2.92, ly - 0.02, 2.0, 0.2, right_label, size=9)


# =============================================================================
# Boxplot with jitter
# =============================================================================

def draw_box_with_jitter(shapes, ax: AxisSpec, labels: Sequence[str],
                         data: Sequence[Sequence[float]], *,
                         box_colors: Sequence[RGBColor],
                         point_colors: Optional[Sequence[Sequence[RGBColor]]] = None,
                         n_below: bool = True,
                         box_width_frac: float = 0.4,
                         p_annotation: Optional[str] = None,
                         bracket: bool = True) -> None:
    n = len(labels)
    slot = ax.width / max(n, 1)
    for i, (lbl, values) in enumerate(zip(labels, data)):
        cx = ax.x0 + slot * (i + 0.5)
        if not values:
            continue
        q1, med, q3 = quartiles(values)
        wlo, whi = tukey_whiskers(values)
        y_q1 = data_to_plot_y(q1, ax)
        y_q3 = data_to_plot_y(q3, ax)
        y_med = data_to_plot_y(med, ax)
        y_wlo = data_to_plot_y(wlo, ax)
        y_whi = data_to_plot_y(whi, ax)

        box_w = slot * box_width_frac
        x0 = cx - box_w / 2

        # Box
        add_rect(shapes, x0, y_q3, box_w, y_q1 - y_q3,
                 fill=box_colors[i], outline=OKABE_ITO["black"],
                 outline_weight_pt=0.6)
        # Median
        add_line(shapes, x0, y_med, x0 + box_w, y_med,
                 color=OKABE_ITO["black"], weight_pt=1.1)
        # Whiskers
        add_line(shapes, cx, y_q1, cx, y_wlo, weight_pt=0.5)
        add_line(shapes, cx, y_q3, cx, y_whi, weight_pt=0.5)
        # Whisker caps
        add_line(shapes, cx - 0.04, y_wlo, cx + 0.04, y_wlo,
                 weight_pt=0.5)
        add_line(shapes, cx - 0.04, y_whi, cx + 0.04, y_whi,
                 weight_pt=0.5)

        # Jittered points
        import random
        rng = random.Random(17 + i)
        for jidx, v in enumerate(values):
            jx = cx + (rng.random() - 0.5) * box_w * 0.75
            yv = data_to_plot_y(v, ax)
            pt_col = (point_colors[i][jidx] if point_colors
                      else OKABE_ITO["black"])
            add_circle(shapes, jx, yv, 0.028, fill=pt_col,
                       outline=OKABE_ITO["black"], outline_weight_pt=0.2)

        # Category label
        add_text(shapes, cx - 0.8, ax.y0 + ax.height + 0.08, 1.6, 0.22,
                 lbl, size=9, align="center")
        if n_below:
            add_text(shapes, cx - 0.8, ax.y0 + ax.height + 0.3, 1.6, 0.20,
                     f"n = {len(values)}", size=8, align="center",
                     color=OKABE_ITO["grey"])

    draw_axis(shapes, ax, draw_x=True, draw_y=True)

    # P-value bracket (n == 2 case)
    if p_annotation and n == 2 and bracket:
        x_left = ax.x0 + slot * 0.5
        x_right = ax.x0 + slot * 1.5
        y_bracket = ax.y0 + 0.25
        add_line(shapes, x_left, y_bracket, x_right, y_bracket,
                 weight_pt=0.6)
        add_line(shapes, x_left, y_bracket, x_left, y_bracket + 0.1,
                 weight_pt=0.6)
        add_line(shapes, x_right, y_bracket, x_right, y_bracket + 0.1,
                 weight_pt=0.6)
        add_text(shapes, (x_left + x_right) / 2 - 1.0, y_bracket - 0.26,
                 2.0, 0.22, p_annotation, size=9, align="center",
                 bold=True)


# =============================================================================
# Scatter plot
# =============================================================================

def draw_scatter(shapes, ax: AxisSpec,
                 xs: Sequence[float], ys: Sequence[float],
                 *, colors: Sequence[RGBColor],
                 radii_in: Sequence[float],
                 outline: RGBColor = OKABE_ITO["black"],
                 outline_weight_pt: float = 0.25,
                 reference_h: Optional[tuple[float, str]] = None,
                 draw_axes: bool = True) -> None:
    for x, y, col, r in zip(xs, ys, colors, radii_in):
        xp = data_to_plot_x(x, ax)
        yp = data_to_plot_y(y, ax)
        add_circle(shapes, xp, yp, r, fill=col, outline=outline,
                   outline_weight_pt=outline_weight_pt)
    if draw_axes:
        draw_axis(shapes, ax)

    if reference_h is not None:
        val, txt = reference_h
        yref = data_to_plot_y(val, ax)
        add_line(shapes, ax.x0, yref, ax.x0 + ax.width, yref,
                 color=OKABE_ITO["vermilion"], weight_pt=0.8, dash="dash")
        add_text(shapes, ax.x0 + 0.05, yref - 0.22, 3.0, 0.2,
                 txt, size=8, color=OKABE_ITO["vermilion"], bold=True)


# =============================================================================
# Heatmap as native table of filled rects + textbox in each cell
# =============================================================================

def draw_categorical_heatmap(shapes, x0: float, y0: float,
                             cell_w: float, cell_h: float,
                             row_labels: Sequence[str],
                             col_labels: Sequence[str],
                             values: Sequence[Sequence],
                             cell_color_fn,  # fn(row, col, value) -> RGBColor
                             cell_text_fn=None,  # fn(row, col, value) -> str
                             cell_text_color_fn=None,
                             row_label_colors: Optional[dict[str, RGBColor]] = None,
                             row_label_bold: Optional[set[str]] = None,
                             col_label_rotation: int = 0,
                             cell_font_size: int = 9,
                             row_label_size: int = 10,
                             col_label_size: int = 10,
                             col_label_bold: Optional[set[str]] = None) -> None:
    for ri, row in enumerate(row_labels):
        for ci, col in enumerate(col_labels):
            val = values[ri][ci]
            col_fill = cell_color_fn(row, col, val)
            xc = x0 + ci * cell_w
            yc = y0 + ri * cell_h
            add_rect(shapes, xc, yc, cell_w, cell_h, fill=col_fill,
                     outline=OKABE_ITO["lightgrey"],
                     outline_weight_pt=0.25)
            if cell_text_fn is not None:
                txt = cell_text_fn(row, col, val)
                if txt:
                    col_txt = (cell_text_color_fn(row, col, val)
                               if cell_text_color_fn
                               else OKABE_ITO["black"])
                    add_text(shapes, xc + 0.02, yc + cell_h / 2 - 0.12,
                             cell_w - 0.04, 0.22, txt,
                             size=cell_font_size, align="center",
                             bold=True, color=col_txt)

    # Row labels (left)
    for ri, row in enumerate(row_labels):
        lc = (row_label_colors.get(row) if row_label_colors
              else OKABE_ITO["black"])
        lb = bool(row_label_bold and row in row_label_bold)
        add_text(shapes, x0 - 1.0, y0 + ri * cell_h + cell_h / 2 - 0.12,
                 0.95, 0.22, row, size=row_label_size,
                 align="right", color=lc, bold=lb)
    # Col labels (above)
    ncols = len(col_labels)
    for ci, col in enumerate(col_labels):
        lb = bool(col_label_bold and col in col_label_bold)
        add_text(shapes, x0 + ci * cell_w - 0.2,
                 y0 - 0.45, cell_w + 0.4, 0.4, col,
                 size=col_label_size, align="center", bold=lb,
                 rotation=col_label_rotation)


# =============================================================================
# Forest plot (horizontal rows, point + error bar)
# =============================================================================

def draw_forest_plot(shapes, ax: AxisSpec, rows: Sequence[str],
                     point_vals: Sequence[float],
                     err_lows: Sequence[float], err_highs: Sequence[float],
                     *, point_colors: Sequence[RGBColor],
                     q_values: Sequence[float],
                     reference_x: float = 1.0,
                     ref_label: str = "FC = 1 (no enrichment)") -> None:
    n = len(rows)
    row_h = ax.height / max(n, 1)
    # Reference vertical line
    xref = data_to_plot_x(reference_x, ax)
    add_line(shapes, xref, ax.y0, xref, ax.y0 + ax.height,
             color=OKABE_ITO["grey"], weight_pt=0.8, dash="dash")
    add_text(shapes, xref + 0.05, ax.y0 - 0.32, 3.0, 0.24, ref_label,
             size=8, color=OKABE_ITO["grey"])

    for i, row in enumerate(rows):
        yc = ax.y0 + row_h * (i + 0.5)
        x_pt = data_to_plot_x(point_vals[i], ax)
        x_lo = data_to_plot_x(err_lows[i], ax)
        x_hi = data_to_plot_x(err_highs[i], ax)

        # Error bar
        add_line(shapes, x_lo, yc, x_hi, yc, weight_pt=0.6)
        add_line(shapes, x_lo, yc - 0.04, x_lo, yc + 0.04, weight_pt=0.6)
        add_line(shapes, x_hi, yc - 0.04, x_hi, yc + 0.04, weight_pt=0.6)
        # Point
        add_circle(shapes, x_pt, yc, 0.05, fill=point_colors[i],
                   outline=OKABE_ITO["black"], outline_weight_pt=0.3)

        # Row label (left)
        add_text(shapes, ax.x0 - 2.4, yc - 0.12, 2.3, 0.22, row,
                 size=9, align="right")

        # Right annotation
        ann = (f"FC = {point_vals[i]:.2f}   "
               f"[{err_lows[i]:.2f}, {err_highs[i]:.2f}]   "
               f"q = {q_values[i]:.1e}")
        add_text(shapes, ax.x0 + ax.width + 0.05, yc - 0.12, 4.2, 0.22,
                 ann, size=8)

    draw_axis(shapes, ax, draw_x=True, draw_y=False)


# =============================================================================
# Synteny track (multiple rows of ure arrows)
# =============================================================================

def draw_synteny_track(shapes, x_in: float, y_in: float, track_w: float,
                       track_h: float, features: Sequence[dict],
                       *, track_label: str,
                       track_label_color: RGBColor = OKABE_ITO["black"],
                       span_kb: float = 40.0) -> None:
    """features: list of {"gene":str, "start_kb":float, "end_kb":float,
                          "strand":"+/-", "color":RGBColor}."""
    # Track baseline
    add_line(shapes, x_in, y_in + track_h / 2, x_in + track_w,
             y_in + track_h / 2, weight_pt=0.5, color=OKABE_ITO["grey"])
    # Track label (left)
    add_text(shapes, x_in - 1.2, y_in + track_h / 2 - 0.1, 1.15, 0.2,
             track_label, size=10, align="right", bold=True,
             color=track_label_color)

    for ft in features:
        s = ft["start_kb"] / span_kb
        e = ft["end_kb"] / span_kb
        xa = x_in + s * track_w
        xe = x_in + e * track_w
        arrow_w = max(0.02, xe - xa)
        add_arrow_shape(shapes, xa, y_in, arrow_w, track_h,
                        fill=ft["color"], reverse=ft.get("strand") == "-")
        # Gene label inside arrow (centred) if wide enough
        if arrow_w >= 0.3:
            add_text(shapes, xa, y_in + track_h / 2 - 0.1, arrow_w, 0.2,
                     ft["gene"], size=7, align="center", bold=True,
                     color=OKABE_ITO["black"])
