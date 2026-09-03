#!/usr/bin/env python
"""Rule check over every rendered figure in new_figure/figures/.

This runs on the SAVED files, not on a live matplotlib figure, so it catches a script
whose in-process audit passed but whose output was never regenerated.  Checks:

  1. the SVG, PDF and PNG of each stem all exist
  2. the page is 180 mm wide
  3. every font-family declaration is Arial-first
  4. no sentence-like text survived onto the page
  5. panel letters are present
  6. the expected set of stems (Fig1-Fig8, Fig_S1-Fig_S21) is complete

Run: /home/soon/miniconda3/envs/dram_env/bin/python _job/verify_outputs.py
"""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

FIG = Path(__file__).resolve().parent.parent / "figures"
SENTENCE = re.compile(r"[A-Za-z]{3,}[^.;]*[.;]")
LETTER = re.compile(r"^[A-Z]$")
EXPECTED = [f"Fig{i}" for i in range(1, 9)] + [f"Fig_S{i}" for i in range(1, 22)]


def mm(v):
    v = v.strip()
    for suf, f in (("mm", 1.0), ("pt", 25.4 / 72), ("px", 25.4 / 96), ("in", 25.4)):
        if v.endswith(suf):
            return float(v[:-len(suf)]) * f
    return float(v) * 25.4 / 96


def check(stem):
    problems = []
    svg = FIG / f"{stem}.svg"
    for ext in ("svg", "pdf", "png"):
        if not (FIG / f"{stem}.{ext}").exists():
            problems.append(f"missing .{ext}")
    if not svg.exists():
        return problems, {}
    raw = svg.read_text()
    fams = set(re.findall(r"font-family:\s*([^;\"'}]+)", raw))
    bad = [f for f in fams if not f.strip().startswith("Arial")]
    if bad:
        problems.append(f"font-family not Arial-first: {sorted(bad)[:3]}")
    root = ET.fromstring(raw)
    w = mm(root.get("width", "0"))
    h = mm(root.get("height", "0"))
    if abs(w - 180.0) > 0.6:
        problems.append(f"page width {w:.1f} mm, expected 180")
    texts = []
    for t in root.iter("{http://www.w3.org/2000/svg}text"):
        s = "".join(t.itertext()).strip()
        if s:
            texts.append(s)
    prose = [s for s in texts if len(s.split()) > 9 and SENTENCE.search(s)]
    if prose:
        problems.append(f"{len(prose)} sentence-like item(s): {prose[0][:60]!r}")
    letters = [s for s in texts if LETTER.match(s)]
    if not letters:
        problems.append("no panel letters found")
    return problems, dict(w=w, h=h, n_text=len(texts), n_letters=len(letters))


def main():
    present = {p.stem for p in FIG.glob("*.svg")}
    fails = 0
    print(f"{'figure':12s} {'w x h (mm)':>12s} {'text':>5s} {'letters':>8s}  status")
    for stem in EXPECTED:
        if stem not in present:
            print(f"{stem:12s} {'-':>12s} {'-':>5s} {'-':>8s}  MISSING")
            fails += 1
            continue
        problems, info = check(stem)
        status = "PASS" if not problems else "FAIL"
        fails += bool(problems)
        dims = f"{info.get('w', 0):.0f} x {info.get('h', 0):.0f}" if info else "-"
        print(f"{stem:12s} {dims:>12s} {info.get('n_text', 0):5d} "
              f"{info.get('n_letters', 0):8d}  {status}")
        for p in problems:
            print(f"    - {p}")
    extra = sorted(present - set(EXPECTED))
    if extra:
        print(f"unexpected stems (not in the 29-figure set): {extra}")
    print(f"\n{len(EXPECTED)} figures expected, {fails} missing or with problems")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
