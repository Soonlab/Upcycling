"""Cross-artefact consistency audit - Microbial Genomics copy (2026-09-23): panel callouts are capital letters
(matching the figures and legends), and the body ends at "## Conflicts of interest" (Society back-matter order).
Run with UPCYCLING_MAN_DIR=../_build/audit_view.

Reads the figures shipped in the package (Figures/) rather than new_figure/figures_v2.

Checks the manuscript, the legends, the figure pages and the three table workbooks
against each other. Prints one line per check and exits non-zero if any check fails.
"""
import re
import sys
from pathlib import Path

import openpyxl

BASE = Path("/data/data/Upcycling")
# UPCYCLING_MAN_DIR=/data/data/Upcycling/SUBMISSION_v2 audits the submission master instead of the 09-04 copy
import os
MAN_DIR = Path(os.environ.get("UPCYCLING_MAN_DIR", Path(__file__).resolve().parent.parent))
MAN = MAN_DIR / "01_Manuscript.md"
LEG = MAN_DIR / "02_Figure_legends.md"
FIGDIR = MAN_DIR / "Figures"
TABDIR = BASE / "SUBMISSION/Supplementary_tables_v2"
# since 2026-09-23 the master ships its own workbooks (S1 = reference panels, S2 = per-MAG, S3 = statistics; sheets S1A ...)
if (MAN_DIR / "Supplementary_tables").is_dir() and list((MAN_DIR / "Supplementary_tables").glob("Table_S1_*.xlsx")):
    TABDIR = MAN_DIR / "Supplementary_tables"
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
for m in re.finditer(r"Fig(?:\.|ure)\s+(S?)(\d+)([A-Ea-e](?:\s*,\s*[A-Ea-e])*)?", body):
    supp, num, panels = m.group(1), m.group(2), m.group(3)
    stem = f"Fig_S{num}" if supp else f"Fig{num}"
    if panels:
        for pl in re.findall(r"[A-Ea-e]", panels.lower()):
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
for tag in ("S1", "S2", "S3"):
    fname = sorted(TABDIR.glob(f"Table_{tag}_*.xlsx"))
    assert len(fname) == 1, f"expected one Table_{tag}_*.xlsx in {TABDIR}, found {fname}"
    wb[tag] = set(openpyxl.load_workbook(fname[0], read_only=True).sheetnames)

named = set()
# legacy scheme "sheet S1.1", "sheets S3.7–S3.8"
for s_ in set(re.findall(r"sheets?\s+((?:S\d+\.\d+(?:\s*[–-]\s*S\d+\.\d+)?)(?:\s+and\s+S\d+\.\d+)?)", man + leg)):
    named |= set(re.findall(r"S\d+\.\d+", s_))
    rng = re.match(r"S(\d+)\.(\d+)\s*[–-]\s*S(\d+)\.(\d+)", s_)
    if rng and rng.group(1) == rng.group(3):
        named |= {f"S{rng.group(1)}.{i}" for i in range(int(rng.group(2)), int(rng.group(4)) + 1)}
# lettered scheme since 2026-09-23: "Table S1A", "Table S3G,H", "(S1A–S1B)"
for m in re.finditer(r"(?<!Fig\. )\bS(\d)([A-Z])((?:\s*,\s*[A-Z]\b)*)(?:\s*[–-]{1,2}\s*S\1([A-Z]))?", man + leg):  # MGen copy: skip "Fig. S4B" figure-panel callouts
    d, a, more, b = m.group(1), m.group(2), m.group(3), m.group(4)
    named |= {f"S{d}{chr(c)}" for c in range(ord(a), ord(b) + 1)} if b else {f"S{d}{a}"}
    named |= {f"S{d}{x}" for x in re.findall(r"[A-Z]", more or "")}

all_sheets = {sh.split("_")[0]: tag for tag, shs in wb.items() for sh in shs if sh != "README"}
unknown = sorted(n for n in named if n not in all_sheets)
check("every cited workbook sheet exists", not unknown, f"unknown: {unknown or 'none'}")

tab_calls = set(re.findall(r"Table\s+(S?\d+)[A-Z]?\b", body))  # "Table S2A" counts as S2
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
# author-date list since 2026-09-19 ("Family, I., ... Year. Title"); the numbered form is still accepted
entries = re.findall(r"^(?:\d+\.\s+)?([A-Z][A-Za-z\u00C0-\u024F\-']+),? [A-Z].*?\b(\d{4})\b", ref_block, re.M)
check("the reference list was parsed", len(entries) > 0, f"{len(entries)} entries")
surnames = {e[0] for e in entries}
NAME = r"[A-Z][A-Za-z\u00C0-\u024F\-']+"
cite_body = body.split("## 1. Introduction", 1)[1] if "## 1. Introduction" in body else body  # skip title page / affiliations
cited = {m.group(1) for m in re.finditer(
    rf"({NAME})(?:\s+and\s+{NAME}|,\s+{NAME}(?:,\s+{NAME})*)?(?:\s+et al\.?)?,\s*\d{{4}}", cite_body)}
orphan = sorted(s for s in surnames if s not in cited)
missing_ref = sorted(c for c in cited if c not in surnames)
check("every citation in the body has a reference entry",
      not missing_ref, f"no entry: {missing_ref or 'none'}")
check("no reference entry is uncited", not orphan, f"orphan: {orphan or 'none'}")

# ---------------------------------------------------------------- 7 word count
lines = man.split("\n")
s = next(i for i, l in enumerate(lines) if l.startswith("## 1. Introduction"))
e = next(i for i, l in enumerate(lines) if l.startswith("## Conflicts of interest"))
w = len(" ".join(lines[s:e]).split())
banner = re.search(r"\*\*Word count \(body Intro→Conclusions\):\*\* ([\d,]+)", man)
# 7,000 was a self-imposed target; on 2026-09-19 the author chose complete tool citations over it (journal has no body limit)
check("body is within the word budget (7,000 target + citation allowance)", w <= 7050, f"{w} words")
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
