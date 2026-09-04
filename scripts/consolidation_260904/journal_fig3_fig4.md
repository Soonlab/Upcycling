# Consolidated Fig 3 and Fig 4 — build note (2026-09-04)

Interpreter: `/home/soon/miniconda3/envs/dram_env/bin/python`
New builders: `/data/data/Upcycling/new_figure/build_v2_fig3.py`, `build_v2_fig4.py`
Output directory: `/data/data/Upcycling/new_figure/figures_v2/` (SVG + PDF + PNG per stem)
No existing `build_*.py` and no `_style.py` was modified.

## Panel mapping

### Fig 3 — catalytic and structural conservation of urease (stem `Fig3`)

| new | old | content | source |
|---|---|---|---|
| A | Fig 8A (`build_fig8.py`) | UreC active-site residue match matrix, 6 MAGs × 7 catalytic sites | `Table_S12_UreC_active_site_residues.csv` |
| B | Fig 8B (`build_fig8.py`) | ESMFold TM-score (left sub-axis) + backbone RMSD (right sub-axis) vs PDB 4CEU chain C — kept as the two sub-axes it already was, lettered once | `Table_S22_ureC_vs_4CEU_tm.csv` |
| C | S17C (`build_supS17.py`) | PAML codeml M0 ω per urease gene | `Table_S19b_codeml_M0_summary.csv` |
| D | S17D (`build_supS17.py`) | yn00 pairwise ω by gene and pair class, log10 ω axis retained | `research/additional/C3_dnds_codon/yn00_pairwise.csv`, asserted against `Table_S19c` |

Layout: row 1 = A (left) + B (right); row 2 = C (left) + D (right). Reading order A → B → C → D.
S17 panels A and B (GC3, ENC) are **not** on this page; per DESIGN.md they go to the new Fig S5 D/E.

### Fig 4 — trait-module enrichment and the alkali-tolerance signature (stem `Fig4`)

| new | old | content | source |
|---|---|---|---|
| A | Fig 6 (`build_fig6.py`) | permutation forest, all 22 subcategories with FC > 1, bootstrap 95 % CI, FC and q stat columns | `Table_S2a` (raw counts, FC re-derived), `Table_S2c` |
| B | Fig 4 panel **C only** (`build_fig4.py`) | MICP-critical trait modules, mean hits per 10³ CDS, MICP-complete vs rest | `Table_S2a` + `Table_S2b` (normalisation asserted) + `Table_S2c` (group means asserted) |
| C | S11 (`build_supS11.py`) | alkaliphile signature — Mrp copies, Nha copies, proteome median pI, acidic-pI fraction, with Mann-Whitney P per block | `Table_S15a_alkaliphile_signature_per_MAG.csv` |

Layout: A spans the full page width on top (22 rows, none dropped); B is bottom-left; C is
bottom-right, re-laid out from the old 1 × 4 strip into a 2 × 2 grid of blocks so that it
fits beside B. Reading order A → B → C.
Old Fig 4 panels A and B (the two DRAM genus heat maps) are **not** drawn here — they belong
to the supplementary consolidation (new Fig S2 A/B) built separately.

## Provenance / asserts

Every assert of the source scripts is preserved, and two were added.

* Fig 3: `ref_len` single-reference check (Fig 8B); `omega_M0 < 1` (S17C); the full
  `run_dnds_v3.py` pair filter (0 < ω < 99, dS > 0.01) reproduced and asserted against every
  n, median and Mann-Whitney P in Table S19c (S17D). 595 pairs computed → 574 pass the filter.
* Fig 4: the fold change of **every** row of Table S2c is re-derived from the raw Table S2a
  counts as mean(per-10³-CDS | hero) / mean(per-10³-CDS | rest) and asserted against the
  stored `Fold_change` (tolerance 5e-3) and both stored group means (3-dp grid) — the
  `build_fig6.py` re-derivation, kept intact; the per-10³-CDS normalisation is asserted
  against `Table_S2b`; the hero set is asserted against `_style.HEROES` from both
  `Table_S2a` and `Table_S15a`; the Mrp fold difference quoted in the manuscript (11.7×) is
  recomputed and asserted.
  * **added**: `assert (n_hero_alk, n_rest_alk) == (N_HERO, N_REST)` — the S15a group split
    and the S2a `Hero` flag now have to agree (6 / 105). They do.

`_job/check_no_hardcoded.py` heuristics on the two new scripts: Fig 3 clean; Fig 4 flags one
literal, `0.001` in `f"{q:.3f}" if q >= 0.001 else f"{q:.0e}"`, which is the display-format
cutoff carried over verbatim from `build_fig6.py`, not a datum.

## Page geometry and acceptance gate

| page | width | height | audit | prose_scan |
|---|---|---|---|---|
| Fig3 | 180.0 mm | **133.0 mm** | 144 text items, 0 overlaps, 0 outside, 0 axes overlaps → PASS | no sentence-like text → PASS |
| Fig4 | 180.0 mm | **226.2 mm** | 153 text items, 0 overlaps, 0 outside, 0 axes overlaps → PASS | no sentence-like text → PASS |

Layout iterations needed to reach the gate:

1. Fig 3 — the `ω = 1` reference label sat 3 mm past the right edge of the canvas; panel D
   narrowed from 104 mm to 99 mm.
2. Fig 4 — panel C's first column collided with panel B's stat column; B was narrowed
   (34–74 mm) and C moved to 107/147 mm, and the panel letter C moved to x = 97 mm.
3. Fig 4 — **generic trap worth recording**: the default y locator of the four panel-C
   boxplots also emits ticks *outside* the drawn `ylim` (e.g. −0.5 and 9 for `Mrp_count`,
   whose range is 0–1). Matplotlib still lays those labels out, off the axis, where they
   collided with a neighbouring block. Fixed by pruning the ticks to the limits after
   `set_ylim`. The old `build_supS11.py` carries the same phantom ticks; they only escape
   its audit because its 1 × 4 strip leaves nothing beside them to collide with.
4. Fig 4 — panel B's group legend, drawn `loc="lower right"` inside the axes by
   `build_fig4.py`, covers the two shortest bars (`st.audit` compares text against text
   only, so it never flagged this in the old page). Moved above the axis.

## Discrepancies found between a source table and what the old script drew

1. **`Mrp_count` is binary, but the axis says "copies" (old S11 → new Fig 4C).**
   `Table_S15a.Mrp_count` takes only the values 0 and 1 across all 111 MAGs. The old S11
   axis label "Mrp copies" and the manuscript's "11.7×" therefore describe a *prevalence*
   ratio (2/6 = 0.333 vs 3/105 = 0.0286), not a copy-number ratio. Only 2 of the 6
   MICP-complete MAGs (S13, S16 — both Sphingobacterium) carry Mrp at all. Drawn as in the
   source, label unchanged, but the manuscript wording should say prevalence.
2. **`Table_S15a` and `Table_S15b` disagree on Mrp presence for 2 MAGs.**
   `S15a.Mrp_count` = 0 while `S15b.Na_H_antiporter_Mrp` = 1 for C10 and M13 (both "rest").
   The panel is drawn from S15a, as the old S11 was; if S15b is the correct call the rest
   mean rises to 5/105 and the fold drops from 11.7× to ~7.0×. Needs a source decision.
3. **Panel 3A's legend carries a category with no cell.**
   The "differs" (coral) legend entry of old Fig 8A is drawn although every one of the 42
   cells matches the reference — inherited deliberately from `build_fig8.py` so the colour
   key is complete; noted here so it is not read as a plotting bug.
4. **Panel 4A point vs interval provenance.** The plotted fold change is the value recomputed
   from Table S2a; the 95 % CI is the stored bootstrap interval from Table S2c, which was
   generated around the stored `Fold_change`. The two agree to < 5e-3 (asserted), so the
   point can sit off the stored interval centre by that much. Same behaviour as old Fig 6.
5. **Not a discrepancy of this page, but carried over from the split:** `build_fig4.py`
   documents that the shipped `Table_S5a_DRAM_product.tsv` has all-zero rows for the six
   MICP-complete MAGs and that the intact distillate is
   `/data/pangenome_work/dram_distilled/product.tsv`. That defect belongs to the two DRAM
   panels, which are **not** on this page; the supplementary page that inherits them must
   keep the same source-fidelity check and Table S5a must still be re-exported.
