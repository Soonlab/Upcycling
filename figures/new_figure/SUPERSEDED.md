# Superseded by `../new_figure_v2/` (2026-09-04)

These 29 builders produced the 29-page figure set that the 2026-09-04 consolidation replaced
with 10 pages. They are kept for provenance. The current set is built by
`../new_figure_v2/build_v2_*.py`, which import the same `_style.py` and the same helpers from
this directory.

**Three of these scripts no longer run**, because the supplementary tables they assert against
were re-exported later the same day:

| script | why it aborts |
|---|---|
| `build_fig4.py` | asserts the shipped DRAM product table still differs from the intact distillate on the six MICP-complete rows; the table has since been re-exported and the two are now identical |
| `build_supS4.py` | asserts the dbCAN family table contains no MICP-complete rows; it now covers all 111 MAGs |
| `build_supS19.py` | carries a swine/sheep label-swap map for the coverage table; the table has since been re-emitted with the labels already corrected |

Each aborts before drawing, so no stale output was produced. The v2 builders keep the same
data routes and replace each stale assertion with one that matches the corrected tables.
