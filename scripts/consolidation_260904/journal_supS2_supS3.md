# Build note — consolidated Suppl Fig S2 and Fig S3 (2026-09-04)

Scripts (new, nothing existing was modified):

* `/data/data/Upcycling/new_figure/build_v2_supS2.py` → `new_figure/figures_v2/Fig_S2.{svg,pdf,png}`
* `/data/data/Upcycling/new_figure/build_v2_supS3.py` → `new_figure/figures_v2/Fig_S3.{svg,pdf,png}`

Interpreter: `/home/soon/miniconda3/envs/dram_env/bin/python`. Both scripts were run to
completion; both end in `st.audit` + `st.prose_scan` + `st.save`.

---

## Fig S2 — DRAM metabolic context and dbCAN CAZyme validation

| new panel | old source | content |
|---|---|---|
| A | `build_fig4.py` panel A | genus-aggregated mean completeness of the 24 fractional DRAM modules, 9 genera with n ≥ 2 |
| B | `build_fig4.py` panel B | genus-aggregated prevalence of the 14 CAZy / nitrogen-metabolism present-absent DRAM modules |
| C | `build_supS4.py` panel B | dbCAN class abundance (GH, GT, PL, CE, AA, CBM), mean hits per 10³ CDS, MICP-complete vs rest, with the fold-change / permutation-q / MWU-P stat columns |
| D | `build_supS4.py` panel C | dbCAN family detail, top 14 families by MICP-complete mean, with the mean / mean / MAGs-with-family stat columns |

**Not on this page.** `build_fig4.py` panel C (MICP-critical trait modules) moves to a main
figure, per DESIGN.md. `build_supS4.py` panel A (keyword CAZy proxy) belongs to consolidated
Fig S1.

**DESIGN.md says "old S4C/D" for panel D — there is no panel D in `build_supS4.py`.** That
script draws exactly three panels (A keyword proxy, B dbCAN classes, C dbCAN families). New
panel D is old S4C alone; nothing was dropped and no second half exists.

Layout: A and B are the top row (heat maps sharing the white→green ramp, one colour bar
each); C and D are the bottom row, side by side, sharing one coral/grey key placed in the
free band under C (an in-axes legend sat on top of the shortest bars in C and the tied dots
in D).

**Final page height: 234.2 mm × 180 mm**, inside the 235 mm one-page ceiling; the script
carries `assert H <= 235.0, H` immediately after `H` is computed, as the Fig 1 and Fig 2
builders do, so the ceiling cannot be lost in a later edit.

The first build came out at 237.9 mm. Two changes brought it to 234.2 mm, with no panel,
row or column dropped:

* heat-map row pitch `ROW_H` 5.0 → **4.8 mm**, the pitch consolidated Fig S1 settled on
  (9 genus rows in each of A and B; −1.8 mm on the page, and the measured label band moves
  up with the heat-map bottom rather than being squeezed). A 7 pt genus label is ~3.0 mm
  tall, so the rows stay clear of each other and the audit confirms it.
* bottom margin under the C / D row 14.0 → **12.0 mm** (−2.0 mm); panel D's x-axis title
  still clears the page edge by ~4 mm.

The rotated module names under A and B are not budgeted with a constant, and were **not**
shrunk. The script draws the two heat maps once on a throwaway page, measures the deepest
tick-label bottom with the real renderer, closes it, and builds the real page from that
measurement (`label_depth_mm`, re-asserted on the final figure). The old `build_fig4.py`
reserved a hard 96 mm (`LAB_H = 96.0`), which is about 16 mm more than the labels actually
need at this width.

Audit: `audit: 173 text items | overlaps 0 | outside 0 | axes overlaps 0 -> PASS`,
`prose_scan: no sentence-like text -> PASS`.

---

## Fig S3 — biosafety, mobile elements and defence systems

| new panel | old source | content |
|---|---|---|
| A | `build_supS8.py` (whole page) | abricate hits per MAG, four database blocks (CARD, VFDB, ResFinder, PlasmidFinder), MICP-complete vs rest, two-sided Mann-Whitney P per block |
| B | `build_supS15.py` panel A | geNomad plasmid-flagged and virus-flagged contigs per MAG, MICP-complete vs rest |
| C | `build_supS16.py` panel A | DefenseFinder anti-phage systems per MAG |
| D | `build_supS16.py` panel B | minced CRISPR arrays per MAG |

**Not on this page.** `build_supS15.py` panel B (the per-MICP-complete-MAG urease / CA vs
geNomad mini table) moves to a main figure. Its source Table S17b is still read here for the
one cross-table check that touches a quantity panel B draws (per-MAG plasmid contig counts
of S17a and S17b must agree).

Layout: two rows of a single four-column grid (`LEFT = 22`, `W = 30`, `GAP = 10`), so all
eight box-and-point columns share a width and a baseline. Row 1 is panel A's four database
blocks; row 2 is B (two columns), C (one) and D (one). Panel letters therefore run
A (row 1) → B, C, D (row 2, left to right).

C and D keep the honest null rendering of the old script: y held at 0–1, ticks 0 / 0.5 / 1,
the whole point cloud on the zero line, `0 / 6 > 0` and `0 / 105 > 0` printed under each
column, and the recomputed `P = 1` bracket. Nothing was rescaled.

**Final page height: 148.0 mm × 180 mm.**

Audit: `audit: 93 text items | overlaps 0 | outside 0 | axes overlaps 0 -> PASS`,
`prose_scan: no sentence-like text -> PASS`.

---

## Discrepancies between the source tables and what the old scripts drew

### 1. Table S5a has been re-exported — `build_fig4.py` no longer runs (blocking for the old script)

`build_fig4.py` carried the source-fidelity note that
`SUBMISSION/Supplementary_tables/Table_S5a_DRAM_product.tsv` held all-zero rows for the six
MICP-complete MAGs, and asserted that fact:

```python
assert set(diff[diff > 0].index) == set(HEROES)
assert (shipped.loc[HEROES, num_cols].sum(axis=1) == 0).all()
```

**That re-export has since landed.** The shipped table is now cell-for-cell identical to the
intact distillate `/data/pangenome_work/dram_distilled/product.tsv` (111 × 98, zero differing
cells, identical file size). Consequence: **`build_fig4.py` now aborts** at that assertion —
verified by running it (`AssertionError: []`, line 69; it fails before drawing, so
`figures/Fig4.*` from 02:54 are untouched).

Handling in `build_v2_supS2.py`: the **data-loading route is unchanged** — panels A and B are
still built from the intact distillate, exactly as `build_fig4.py` was. The assertion was
inverted to the check that is now correct: the shipped Table S5a must reproduce the intact
distillate cell for cell, including the six hero rows, and those hero rows must be non-empty.
No drawn value changes.

### 2. Table S6c has been re-exported and now contains all 111 MAGs — `build_supS4.py` no longer runs

`build_supS4.py` asserted `not set(rest_fam_df.index) & set(HEROES)`, because the shipped
`Table_S6c_dbCAN_family_counts.csv` held only the 105 non-MICP-complete MAGs (which is why
the superseded `Figure_S4c` rendered as an empty axis). The file was re-exported on
2026-09-04 03:17 and now has **111 rows including the six MICP-complete MAGs**, with hero
row totals matching Table S6a (`S13 158, S16 249, S23 264, C22 135, M1 60, S26 75`).
Consequence: **`build_supS4.py` now aborts** at that assertion (verified, line 115;
`figures/Fig_S4.*` untouched).

Handling in `build_v2_supS2.py`: the hero families are still recomputed by re-parsing the
direct hmmsearch tables `research/revision/dbcan_direct/<MAG>.tbl` (the documented
provenance route), and every old assert is kept — per-class sums reproduce Table S6a, family
totals reproduce Table S6a `Total`, per-class sums of the normalised family table reproduce
the panel-C class means. The stale "heroes absent from S6c" assert was replaced by the
stronger check that the re-parse reproduces the **newly published hero rows of Table S6c
family for family** (verified: 0 absolute difference across 192 families × 6 MAGs). The rest
group is now `Table S6c` minus the six heroes (105 rows, asserted). **No drawn value
changed**: the new panel D reproduces the old `Figure_S4c` numbers exactly (GT2 2.33/1.66,
GH29 1.84/0.08, … GH177 0.73/0.06, and the same MAGs-with-family counts).

### 3. Old Fig S8 was drawn 4 mm off the page

`build_supS8.py` used `LEFT = 22, W = 33, GAP = 10`, so its fourth block (PlasmidFinder)
spanned 151–184 mm on a 180 mm page: the right edge of that axes was clipped. Little data
was lost (the clipped strip is empty axis beyond the "rest" column), but the page was out of
spec. The new grid is `LEFT = 22, W = 30, GAP = 10` → right edge 172 mm, with
`assert COL_X[-1] + W <= st.PAGE_W_MM` so the defect cannot come back.

### 4. Cosmetic changes made when merging (no data affected)

* Old S15A put the quantity name on each y-axis label; on the merged page the two geNomad
  columns share one y-axis label (`geNomad contigs per MAG`) and the quantity name moves to
  the block header above each column, matching how panel A labels its four database blocks.
* Old S4B/S4C each carried their own in-axes legend; the merged page carries one shared
  coral/grey key for C and D. Old S8 and S15A's lineage key (blue / orange / grey) is
  carried on Fig S3 as one figure-level legend with a fourth entry for the coral
  MICP-complete group used by C and D.
* Old S4B's stat header `permutation\nq (BH)` is shortened to `permut.\nq (BH)` and old
  S4C's `MICP-complete\nMAGs with family` to a three-line `MICP-compl.\nMAGs with\nfamily`,
  so the stat columns fit beside the narrower side-by-side panels.

### 5. Checks run on the outputs

`_job/check_no_hardcoded.py` reports both new scripts **clean**. Both SVGs are 180.0 mm
wide, every `font-family` declaration is Arial-first, and the panel letters A, B, C, D are
present on each page.
