# Consolidated submission package — 2026-09-04 (text: authors' revision of 2026-09-17, ported 2026-09-23)

## 2026-09-23 — authors' revised text ported

`01_Manuscript.md` and `02_Figure_legends.md` now carry the authors' revised text
(`Manuscript_revised.docx`, saved 2026-09-17; `submission_v3/Manuscript.docx` is the same
text without the English author block). The port is `port_revised_260923.py`; the
pre-port files are `*.pre_revised_260923.*`. Body 4,320 words, abstract 198 words,
57 cited references (all DOI-verified; the 09-19 corrections re-applied to the new
sentences). The supplementary workbooks follow the revised text's numbering — **S1 =
reference panels and methods (sheets S1A–S1O), S2 = per-MAG measurements (S2A–S2T), S3 =
comparative statistics (S3A–S3M)**; the 09-04 workbooks (old S1/S2/S3 order, `S1.1` sheet
names) are in `_superseded_260923/`. Both audits pass on the ported master (20 structural
checks, 109 numeric checks). Decisions, claim-strength comparison and the open author items
are in `REVISED_PORT_REPORT_260923.md`. The sections below describe the 09-04 package and
remain valid for the figures and the data.

This package replaces the 2026-09-04 morning package in `SUBMISSION/`. The science is
unchanged; the display items and the manuscript length were consolidated to the limits
requested, and two data errors found during the consolidation were corrected.

## What changed

| | before | now |
|---|---|---|
| main figures | 8 | **5** |
| supplementary figures | 21 | **5** |
| tables | 3 main + 23 supplementary (52 files) | **2 main + 3 supplementary workbooks** |
| body length | 9,731 words | **6,997 words** |

No analysis was re-run and no result was dropped. Every panel of the 29 old pages is either
carried into one of the 10 new pages or was an explicitly labelled duplicate of a main panel
(old Fig S9, S13, S20, S21). All 52 shipped supplementary files are carried into the three
workbooks, verified programmatically. The old-to-new map is
`../consolidation_260904/MAPPING.md`.

## Contents

* `01_Manuscript.md` / `.docx` — body 7,005 words, 60 references (author-date list, every entry DOI-verified 2026-09-19; see `references/README_EndNote.md`), 2 main tables.
* `references/` — EndNote library (`.ris`), EndNote-ready manuscript with temporary citations, and the scripts that verify references against CrossRef and regenerate the list.
* `02_Figure_legends.md` / `.docx` — 10 figure legends and 5 table legends.
* `Figures/` — `Fig1`–`Fig5` and `Fig_S1`–`Fig_S5`, each as `.png` (200 dpi, the print
  deliverable), `.svg` and `.pdf`. Every page is 180 mm wide and fits within 235 mm of
  height, so each is one printed page.
* `Main_tables/` — `Table_1` and `Table_2` as editable Excel files (one sheet each: title, table, legend and footnote), exported from the markdown tables in `01_Manuscript.md` by `build_main_tables_xlsx.py`; re-run it after any edit to either table in the manuscript.
* `Supplementary_tables/` — since 2026-09-23: `Table_S1` (reference panels and method support,
  15 sheets S1A–S1O), `Table_S2` (per-MAG measurements, 20 sheets S2A–S2T), `Table_S3`
  (comparative statistics, 13 sheets S3A–S3M), each with a README sheet naming the shipped file behind every sheet; plus the
  two raw artefacts too large or too unstructured to be a sheet, the full DRAM metabolism
  workbook and the IQ-TREE SH/AU report.
* `03_Graphical_abstract` — 190 x 105 mm, `.png` (300 dpi), `.svg` and `.pdf`, rebuilt
  2026-09-09 from `../new_figure/build_v2_graphical_abstract.py`; every number in it is
  loaded or recomputed from the same supplementary sources as the figures.
* `04_`/`05_` cover letters — carried over unchanged.
* `rebuild_docx.sh` — regenerates both `.docx` from the markdown. Run it after any `.md` edit.

## Two data corrections made during this consolidation

1. **The Mrp antiporter result is a prevalence, not a dosage.** The `Mrp_count` column of the
   per-MAG table is binary presence. Mrp is detected in 2 of the 6 MICP-complete MAGs
   (S13, S16) and 3 of the remaining 105, so the 11.7-fold figure is a prevalence ratio over
   small counts. The manuscript had described it as "0.33 versus 0.029 copies". Corrected in
   §3.10, §4.3, §4.5, the abstract, the Highlights, the Conclusions and the Fig. 4 legend.
2. **The minimum species-level ANI is 95.08 %, not 95.06 %.** The true minimum is M6 against
   *Stenotrophomonas indicatrix*. Corrected in §3.1 and the Fig. 5 legend.

## Figure revision of 2026-09-09

Three review rounds on the 2026-09-04 pages; the science and every number are unchanged.
Main Fig 1: the tree and the gene-presence matrix are one radial panel (A), so the old
panels C and D are now B and C — the manuscript callouts and the legend were renumbered.
The graphical abstract was redrawn from the same sources as the figures (the 2026-04
version was four text boxes whose animal glyphs did not render).  Fig 2C (was a table of
counts) and Fig 5A right (was a list of identifiers) are drawn as
graphics; Fig 3A shows the catalytic residues with their flanking alignment columns;
Fig 2 and Fig 4/5 use larger type; panel heights and widths were aligned on every page.
The pre-revision pages are kept in `../new_figure/figures_v2/_pre260909/`; the changes
are logged in `../new_figure/_job/JOURNAL.md` (entries dated 2026-09-09).

## Verification

Two audits are shipped alongside, both re-runnable:

* `../consolidation_260904/audit_consistency.py` — 19 structural checks: figure and table
  callouts resolve, every cited panel exists in its legend, one legend per figure, page
  geometry, workbook sheets, section cross-references, citations against the reference list,
  and the display-item and word-count budgets. **All pass.**
* `../consolidation_260904/audit_numbers.py` — 90 numeric checks recomputing every headline
  value in the manuscript from the shipped supplementary source. **All agree.**

Run both with `/home/soon/miniconda3/envs/dram_env/bin/python`.

## Still blocked on the authors

Unchanged from the previous package: author names, affiliations and ORCIDs; funder and grant
number; the NCBI BioProject accession replacing `PRJNA-XXXXXXX` in two places; a Zenodo DOI;
three suggested reviewers; the journal decision, after which one cover letter is deleted; and
the decision on whether S26 stays in the analysis set.
