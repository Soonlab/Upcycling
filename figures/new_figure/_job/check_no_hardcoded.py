#!/usr/bin/env python
"""Rule 4 check: no data value may be hard-coded in a build script.

Style constants (positions in mm, font sizes, axis limits, tick lists, alpha, linewidths)
are legitimate.  What must not appear is a measurement: a correlation, a P value, a count
that should have been read from a file.  This flags the literals that look like data and
leaves the judgement to a reader, rather than pretending to decide automatically.

Heuristics for a suspicious literal:
  * scientific notation with a negative exponent (a P value written out)
  * a float with 3+ decimal places (a coefficient or an effect size)
Literals inside an obvious layout call (ax_mm, set_xlim, set_ylim, set_xticks, set_yticks,
figsize, linspace, arange) or on a line that assigns a style constant are exempt.

Run: /home/soon/miniconda3/bin/python _job/check_no_hardcoded.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LAYOUT = re.compile(r"\b(ax_mm|set_x?lim|set_y?lim|set_xticks|set_yticks|set_xbound|"
                    r"set_ybound|figsize|linspace|arange|MaxNLocator|set_bounds|margins|"
                    r"bbox_to_anchor|transAxes|alpha|lw|linewidth|width|height|pad|"
                    r"fontsize|labelpad|columnspacing|handletextpad|labelspacing|borderpad|"
                    r"solid_capstyle|zorder|ms|s=)\b")
SCI = re.compile(r"\b\d+(?:\.\d+)?[eE]-\d+\b")
DEC3 = re.compile(r"(?<![\w.])\d+\.\d{3,}(?![\w])")


def scan(path):
    hits = []
    for i, line in enumerate(path.read_text().splitlines(), 1):
        code = line.split("#", 1)[0]
        if not code.strip() or code.lstrip().startswith(('"""', "'''")):
            continue
        if LAYOUT.search(code):
            continue
        found = SCI.findall(code) + DEC3.findall(code)
        # a tolerance in an assertion is a check, not data
        if found and ("assert" in code or "abs(" in code):
            continue
        if found:
            hits.append((i, line.strip()[:100], found))
    return hits


def main():
    scripts = sorted(ROOT.glob("build_*.py"))
    total = 0
    for s in scripts:
        hits = scan(s)
        total += len(hits)
        mark = "clean" if not hits else f"{len(hits)} literal(s) to review"
        print(f"{s.name:24s} {mark}")
        for ln, txt, found in hits:
            print(f"    L{ln}: {txt}")
            print(f"         -> {found}")
    print(f"\n{len(scripts)} scripts scanned, {total} literal(s) flagged for review")
    return 0


if __name__ == "__main__":
    sys.exit(main())
