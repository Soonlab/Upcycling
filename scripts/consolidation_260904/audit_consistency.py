"""Cross-artefact consistency audit for the consolidated submission.

Checks the manuscript, the legends, the figure pages and the three table workbooks
against each other. Prints one line per check and exits non-zero if any check fails.
"""
import re
import sys
from pathlib import Path

import openpyxl

BASE = Path("/data/data/Upcycling")
MAN = BASE / "consolidation_260904/manuscript/01_Manuscript.md"
LEG = BASE / "consolidation_260904/manuscript/02_Figure_legends.md"
FIGDIR = BASE / "new_figure/figures_v2"
TABDIR = BASE / "SUBMISSION/Supplementary_tables_v2"
MAX_H_MM = 235.0
PAGE_W_MM = 180.0

man = MAN.read_text()
leg = LEG.read_text()
fails = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  |  {detail}" if detail else ""))
    if not ok:
        fails.append(name)


# body text only: strip the legend-style trailing tables and the reference list
body = man.split("## References")[0]

# ---------------------------------------------------------------- 1 figure callouts
fig_calls = set()
for m in re.finditer(r"Fig(?:\.|ure)\s+(S?)(\d+)([a-z](?:\s*,\s*[a-z])*)?", body):
    supp, num, panels = m.group(1), m.group(2), m.group(3)
    stem = f"Fig_S{num}" if supp else f"Fig{num}"
    if panels:
        for pl in re.findall(r"[a-z]", panels):
            fig_calls.add((stem, pl))
    else:
        fig_calls.add((stem, None))

expected_stems = [f"Fig{i}" for i in range(1, 6)] + [f"Fig_S{i}" for i in range(1, 6)]
called_stems = {s for s, _ in fig_calls}
check("every main/supp figure is cited in the body",
      set(expected_stems) <= called_stems,
      f"uncited: {sorted(set(expected_stems) - called_stems) or 'none'}")
check("no callout to a figure outside Fig1-5 / Fig_S1-S5",
      called_stems <= set(expected_stems),
      f"stray: {sorted(called_stems - set(expected_stems)) or 'none'}")

# panels actually lettered in the legends
leg_panels = {}
for block in re.split(r"\n### (?=Fig)", leg)[1:]:
    m = re.match(r"Fig(?:ure)?\s+(S?)(\d+)", block)
    if not m:
        continue
    stem = f"Fig_S{m.group(2)}" if m.group(1) else f"Fig{m.group(2)}"
    leg_panels[stem] = {x.lower() for x in re.findall(r"\*\*\(([A-E])\)\*\*", block)}
bad = [(s, p) for s, p in fig_calls if p and p not in leg_panels.get(s, set())]
check("every cited figure panel exists in that figure's legend",
      not bad, f"missing: {bad or 'none'}")

# ---------------------------------------------------------------- 2 legends <-> files
check("one legend per figure, ten in total",
      set(leg_panels) == set(expected_stems),
      f"legends for {sorted(leg_panels)}")
missing_files = [s for s in expected_stems
                 if not all((FIGDIR / f"{s}.{e}").exists() for e in ("png", "svg", "pdf"))]
check("every figure exists as png + svg + pdf", not missing_files,
      f"missing: {missing_files or 'none'}")

# ---------------------------------------------------------------- 3 page geometry
def svg_mm(p):
    head = p.read_text()[:600]
    w = re.search(r'width="([\d.]+)pt"', head)
    h = re.search(r'height="([\d.]+)pt"', head)
    return (float(w.group(1)) * 25.4 / 72, float(h.group(1)) * 25.4 / 72) if w and h else (None, None)

geo = {s: svg_mm(FIGDIR / f"{s}.svg") for s in expected_stems if (FIGDIR / f"{s}.svg").exists()}
badw = {s: g for s, g in geo.items() if g[0] is None or abs(g[0] - PAGE_W_MM) > 0.5}
badh = {s: round(g[1], 1) for s, g in geo.items() if g[1] and g[1] > MAX_H_MM}
check(f"every page is {PAGE_W_MM:.0f} mm wide", not badw, f"off: {badw or 'none'}")
check(f"every page fits within {MAX_H_MM:.0f} mm of height", not badh,
      f"too tall: {badh or 'none'}")
print("      heights (mm): " + ", ".join(f"{s} {geo[s][1]:.1f}" for s in expected_stems if s in geo))

# ---------------------------------------------------------------- 4 table callouts
wb = {}
for tag, fname in (("S1", "Table_S1_per_MAG_measurements.xlsx"),
                   ("S2", "Table_S2_comparative_statistics.xlsx"),
                   ("S3", "Table_S3_reference_panels_and_methods.xlsx")):
    wb[tag] = set(openpyxl.load_workbook(TABDIR / fname, read_only=True).sheetnames)

sheet_calls = set(re.findall(r"sheets?\s+((?:S\d+\.\d+(?:\s*[–-]\s*S\d+\.\d+)?)(?:\s+and\s+S\d+\.\d+)?)", man + leg))
named = set()
for s in sheet_calls:
    named |= set(re.findall(r"S\d+\.\d+", s))
    rng = re.match(r"S(\d+)\.(\d+)\s*[–-]\s*S(\d+)\.(\d+)", s)
    if rng and rng.group(1) == rng.group(3):
        named |= {f"S{rng.group(1)}.{i}" for i in range(int(rng.group(2)), int(rng.group(4)) + 1)}

all_sheets = {sh.split("_")[0]: tag for tag, shs in wb.items() for sh in shs if sh != "README"}
unknown = sorted(n for n in named if n not in all_sheets)
check("every cited workbook sheet exists", not unknown, f"unknown: {unknown or 'none'}")

tab_calls = set(re.findall(r"Table\s+(S?\d+)\b", body))
check("no callout to a table outside Table 1, 2, S1, S2, S3",
      tab_calls <= {"1", "2", "S1", "S2", "S3"},
      f"stray: {sorted(tab_calls - {'1','2','S1','S2','S3'}) or 'none'}")
check("all five tables are cited in the body",
      {"1", "2", "S1", "S2", "S3"} <= tab_calls,
      f"uncited: {sorted({'1','2','S1','S2','S3'} - tab_calls) or 'none'}")

# ---------------------------------------------------------------- 5 section refs
sec_defined = set(re.findall(r"^#{2,3} (\d+(?:\.\d+)?)", man, re.M))
sec_cited = set(re.findall(r"§\s*(\d+\.\d+)", man))
dangling = sorted(s for s in sec_cited if s not in sec_defined)
check("every § cross-reference resolves to a section that exists",
      not dangling, f"dangling: {dangling or 'none'}")

# ---------------------------------------------------------------- 6 references
ref_block = man.split("## References")[1].split("## Table 1")[0]
entries = re.findall(r"^\d+\.\s+([A-Z][A-Za-z\u00C0-\u024F\-']+)", ref_block, re.M)
surnames = set(entries)
NAME = r"[A-Z][A-Za-z\u00C0-\u024F\-']+"
cited = {m.group(1) for m in re.finditer(
    rf"({NAME})(?:\s+and\s+{NAME}|,\s+{NAME}(?:,\s+{NAME})*)?(?:\s+et al\.?)?,\s*\d{{4}}", body)}
orphan = sorted(s for s in surnames if s not in cited)
missing_ref = sorted(c for c in cited if c not in surnames)
check("every citation in the body has a reference entry",
      not missing_ref, f"no entry: {missing_ref or 'none'}")
check("no reference entry is uncited", not orphan, f"orphan: {orphan or 'none'}")

# ---------------------------------------------------------------- 7 word count
lines = man.split("\n")
s = next(i for i, l in enumerate(lines) if l.startswith("## 1. Introduction"))
e = next(i for i, l in enumerate(lines) if l.startswith("## CRediT"))
w = len(" ".join(lines[s:e]).split())
banner = re.search(r"\*\*Word count \(body Intro→Conclusions\):\*\* ([\d,]+)", man)
check("body is within the 7,000-word budget", w <= 7000, f"{w} words")
check("the stated word count matches the text",
      banner and int(banner.group(1).replace(",", "")) == w,
      f"banner {banner.group(1) if banner else '?'} vs actual {w}")

# ---------------------------------------------------------------- 8 display-item count
check("main figures <= 5", len([s for s in expected_stems if not s.startswith("Fig_S")]) <= 5)
check("supplementary figures <= 5", len([s for s in expected_stems if s.startswith("Fig_S")]) <= 5)
n_tab = len(re.findall(r"^## Table \d", man, re.M)) + len(wb)
check("tables <= 5 in total", n_tab <= 5, f"{n_tab} tables")

print()
print(f"{'ALL CHECKS PASS' if not fails else str(len(fails)) + ' CHECK(S) FAILED: ' + ', '.join(fails)}")
sys.exit(1 if fails else 0)
