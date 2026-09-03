# new_figure — SVG/PNG rebuild of the Upcycling figure set

Every manuscript figure (main Fig 1–8, Suppl Fig S1–S21) rebuilt as an editable vector page,
applying the figure method settled during the Rectal_Organoid `svg_supp` work (2026-08-20/21)
and used for SNUH_KMJ `new_figure` (2026-09-03/04). The final deliverable is the **PNG**
(200 dpi, 180 mm wide); SVG (text as text) and PDF (embedded TrueType) sit beside it.

Numbering is unchanged from `SUBMISSION/` (no consolidation), so `02_Figure_legends.md`
callouts still apply. Supplementary figures that were several files (S4 a/b/c, S16, S17 a–d,
S19) are one composed page each.

```
_style.py                 shared palette, type scale, mm layout helpers, self-audit
build_fig{1..8}.py        main figures            -> figures/Fig{N}.{svg,pdf,png}
build_supS{1..21}.py      supplementary figures   -> figures/Fig_S{N}.{svg,pdf,png}
_job/TASK.md              the build contract (format, allowed text, provenance, colour, checks)
_job/JOURNAL.md           merged build journal (sources, recomputations, drops, discrepancies)
_job/journal_*.md         per-group journals as written by the build groups
_job/verify_outputs.py    rule check over the rendered files
_job/check_no_hardcoded.py numeric-literal scan of the builders
_job/run_all.sh           rebuild everything, then verify
```

## Rebuild

```bash
cd /data/data/Upcycling/new_figure
PY=/home/soon/miniconda3/envs/dram_env/bin/python     # has Biopython for the two trees
for f in build_fig*.py build_supS*.py; do $PY "$f"; done
$PY _job/verify_outputs.py
# or:  bash _job/run_all.sh
```

Scripts are independent; order does not matter. Each script runs `st.audit` (text overlaps,
text outside the canvas, overlapping Axes) and `st.prose_scan` (sentence-like text) before saving.

## What changed against the 2026-04 PNG/PPTX set

- No panel titles and no explanatory text; only element-bound labels remain.
- One page per figure at 180 mm, uniform 7/8/10 pt type scale, spines left+bottom.
- Every number is loaded or recomputed from a repository source and asserted where a stored
  value exists; unsourced numbers were dropped and journaled.
- One colour carries one meaning per page (coral = MICP-complete, grey = rest, blue/orange =
  Sphingobacterium/Pseudomonas_E, green = present/match, four source colours).
- Fig 3 tracks are parsed from the Bakta gff3 (the 06-22 pptx re-typeset had used inline
  constants); Fig 7B is a recomputed trait-module PCoA (the pptx version reused panel A's
  coordinates).
