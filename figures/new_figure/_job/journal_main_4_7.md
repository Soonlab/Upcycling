# Build journal — main Fig 4, 5, 6, 7 (group `main_4_7`, 2026-09-04)

Interpreter `/home/soon/miniconda3/envs/dram_env/bin/python` (matplotlib 3.10.8, pandas 1.5.3,
scipy 1.15.2, scikit-bio 0.5.8). Scripts `build_fig4.py` … `build_fig7.py`; stems `Fig4`…`Fig7`
in `new_figure/figures/` as SVG + PDF + PNG (200 dpi).

| page | canvas | audit | prose_scan | verify_outputs |
|---|---|---|---|---|
| Fig4 | 180 × 241 mm | PASS (101 texts, 0 overlaps) | PASS | PASS |
| Fig5 | 180 × 230 mm | PASS (144 texts, 0 overlaps) | PASS | PASS |
| Fig6 | 180 × 141 mm | PASS (79 texts, 0 overlaps) | PASS | PASS |
| Fig7 | 180 × 96 mm | PASS (33 texts, 0 overlaps) | PASS | PASS |

---

## 🔴 Discrepancy 1 (blocking, needs a data fix before submission) — Table S5a is corrupt

`SUBMISSION/Supplementary_tables/Table_S5a_DRAM_product.tsv` carries **all-zero rows for the six
MICP-complete MAGs** (C22, M1, S13, S16, S23, S26): every one of its 32 numeric module columns and
66 boolean module columns is 0 / False for those six genomes.

Cause: `scripts/07_dram_hero_annotate.sh` (and its twin `research/extra/run_dram_hero.sh`) re-ran
DRAM for the six MAGs into `/data/pangenome_work/dram_output/`, then re-distilled into
`dram_output/distillate/`, overwriting the full-panel distillate with a hero-only one. The run was
interrupted (`dram_output/S13.log` ends inside "Getting hits from pfam"; the per-MAG output
directories no longer exist), so the distillate it left behind has the 105 other MAGs from the
earlier run and empty rows for the six it was rebuilding. That file is what was copied into
SUBMISSION as Table S5a.

The intact full-panel distillate survives at **`/data/pangenome_work/dram_distilled/product.tsv`**
(2026-02-09). `build_fig4.py` asserts the relationship rather than asserting it in prose: the two
files are identical for all 105 non-hero MAGs and differ **only** on the six hero rows, which are
non-empty in the intact file (numeric row sums 3.21 – 14.53; 5 – 10 boolean modules each).

Consequences for what was published:
- The shipped Figure 4a is genus-aggregated with a `row_sum > 0` filter, so the four zeroed
  *Sphingobacterium* MAGs were silently dropped and the panel is labelled **"Sphingobacterium
  (n = 2)"**. The genus has **six** MAGs. The rebuilt panel shows n = 6.
- Manuscript §"DRAM distillation …" (line 154 of `01_Manuscript.md`) states that urease and
  nitrogen-metabolism modules "reach full completeness (1.0) in MICP-complete MAGs while averaging
  0.35–0.60 in the remainder". That claim cannot come from Table S5a at all — see Discrepancy 2.

**Action for the user:** re-export Table S5a from `/data/pangenome_work/dram_distilled/product.tsv`
(and Table S5b/S5c likewise — `genome_stats.tsv` differs between the two directories as well).

## 🔴 Discrepancy 2 — DRAM has no urease / CA / antiport / cobalamin module

`product.tsv` has 98 columns and **none** of them matches urease, carbonic anhydrase, Na⁺/H⁺
antiport or cobalamin/B12 (`grep -Ei 'urease|anhyd|antiport|cobalamin|B12'` returns nothing). The
Fig 4 legend's panel-a list ("… nitrogen metabolism, stress response, vitamin (cobalamin / B12)
biosynthesis …") and panel-b list ("urease, carbonic anhydrase, nitrogen metabolism, Na⁺/H⁺
antiport, cobalamin biosynthesis, CAZymes") are not obtainable from DRAM. The shipped panel b in
fact plots **Bakta keyword-scan trait modules**, not DRAM modules, on a log axis.

Rebuild: panel C keeps the Bakta trait modules (the honest source of that panel) and its axis says
so; the legend needs rewording. The 11-module "MICP-critical" DRAM heat map in
`scripts/pptx_builder/build_editable_figures.py::fig4a_dram_heatmap` is **hard-coded invented
data** ("Representative per-genus module completeness … distilled from") — not used here, and it
should not be used anywhere.

## Discrepancy 3 — Fig 4 panel structure changed 2 → 3 panels

DRAM scores 32 modules as a completeness fraction and 66 as present/absent. The shipped panel a
mixed both under one colour bar labelled "Mean module completeness", so a genus-mean of a boolean
column was being read as a completeness fraction. Rebuilt as two panels with two colour bars:

- **A** 24 fractional modules (std > 0 and ≥ 3 MAGs carrying) → "Mean module completeness"
- **B** the 14 present/absent CAZy and nitrogen-metabolism modules → "Fraction of genus with module"
- **C** the Bakta MICP-critical comparison (old panel b)

Dropped by scope, recorded here: the other 16 present/absent modules with variance
(methanogenesis ×4, SCFA ×8, sulfur ×2, other reductases ×2) are in Table S5 but are not part of
what this figure is about. **Legend must be renumbered a/b → A/B/C.**

## Discrepancy 4 — Fig 5a: "21 candidates < 95 % ANI" is wrong

In `Table_S8_novelty_ANI_screen.csv` **no MAG has ANI < 95** (minimum 95.06). The 21
`Novel_sp_candidate = True` rows are exactly the 21 rows with **no species-level ANI at all**
(`ANI` empty), S13 and S16 among them — the script asserts that identity. The shipped panel drops
those 21 rows with `dropna` and then titles the panel with their count, so the panel shows 90
points all above the cutoff under a caption saying 21 are below it. The legend sentence "Twenty-one
of 111 MAGs fall below 95 %, including S13 and S16 (no species-level ANI available)" carries the
same error.

Rebuild: panel A plots the 90 MAGs that have an ANI, and lists the 21 without one as a named block
("no species-level ANI (n = 21)"), heroes in coral. **Legend needs the wording fixed.**

## Discrepancy 5 — Fig 5c: "all 63 RefSeq genomes" are not all there

`Table_S10a` is a 6 × 63 matrix in which skani stores **0.0** for any pair it did not align. Each
study MAG has only 3–5 non-zero entries (C13 4, C22 3, S13 5, S16 4, S23 4, V3 5). The shipped
panel scatters every column and relies on `ylim(72, 101)` to hide the zeros, so "per-reference ANI
(all 63 RefSeq genomes)" in its key describes points that are not drawn. (This is the same class of
defect as the "fake bar" already fixed once in the 2026-06-22 Figure_Master pass.)

Rebuild: only the reported hits (ANI ≥ 80) are drawn, as ranked bars per MAG, and each block header
states "n of 63 RefSeq hits". Every block's top hit is asserted against `Table_S10b`
(`Nearest_ANI`, `Novel_species_candidate`).

## Fig 7b — recomputed, and it reproduces the published numbers exactly

The trait-module ordination has **no stored coordinates** anywhere in the repository (Table S9a is
the pan-genome ordination only). Recomputed from `Table_S2b` with the recipe of
`research/revision/05_PCoA_by_source.py` as modified in
`scripts/revision/regenerate_panelized_figures.py::figure7`: per-module z-score (ddof = 0),
Euclidean distance, classical PCoA, PERMANOVA on waste source, 999 permutations, seed 0.

| quantity | recomputed | published |
|---|---|---|
| pseudo-F, source | 2.7056 | 2.71 |
| p | 0.001 | 0.001 |
| PC1 variance | 16.47 % | 16.5 % |
| PC2 variance | 10.38 % | 10.4 % |

All four are asserted in the script; the p value is stable across seeds 0/1/42. Note the recipe
that reproduces the figure is the **z-scored Euclidean** one; the docstring of the original
`05_PCoA_by_source.py` says "Bray-Curtis on the original counts" in a comment while the code
z-scores — the code is what reproduces the published value.

---

## Recomputation and assertions per figure

**Fig4** — genus from `Table_S1d`; genera with ≥ 2 MAGs, hero-bearing genera first. Panel C means
recomputed from `Table_S2a` (hits ÷ CDS_total × 1000), asserted elementwise against the stored
`Table_S2b`, and each group mean asserted against `Hero_mean_per1kCDS` / `Rest_mean_per1kCDS` in
`Table_S2c` at the 3-decimal grid the table stores.

**Fig5** — novelty-flag ↔ missing-ANI identity asserted; every panel-C block's maximum asserted
against `Table_S10b`.

**Fig6** — every one of the 38 rows of `Table_S2c` re-derived from `Table_S2a`: both group means
asserted at 3 decimals and the fold change asserted to < 5 × 10⁻³ of the stored value. Drawn points
are the recomputed fold changes. CI parsed from the stored `[lo, hi]` strings (bootstrap CI cannot
be re-derived without the 2,000 resamples). 22 rows have fold change > 1 and are drawn.

**Fig7** — panel A uses the stored coordinates and stored pseudo-F / p / variance; panel B is fully
recomputed as above.

## Colour → meaning (checked one meaning per page)

| page | colour | meaning |
|---|---|---|
| Fig4 | white→green ramp | A mean module completeness; B fraction of genus carrying the module (separate bars, each labelled) |
| | coral / grey | MICP-complete genus label and bar / the rest |
| Fig5 | coral | a MICP-complete MAG (point in A, block header in C) |
| | green / grey | ANI ≥ 95 % / < 95 % in C |
| | dark, light blue | S13 and S16 query blocks in B (both *Sphingobacterium*; the two shades separate query blocks, they are not the Fig 8 lineage contrast) |
| | grey dashed / dotted | 95 % species and 70 % genus rules |
| Fig6 | coral / light coral / grey | q < 0.05 / q < 0.10 / n.s.; light coral is omitted from the key because no drawn row falls in that tier |
| Fig7 | brown, pink, green, purple | cattle, swine, sheep, poultry |
| | coral ring + label | a MICP-complete MAG |

## `check_no_hardcoded.py` — flagged literals, all legitimate

- `build_fig4.py` L250 `FixedLocator([0.001, 0.01, 0.1, 1, 10, 100])` — log-axis tick positions.
- `build_fig6.py` L120 `q >= 0.001` — display-format switch between decimal and exponent notation.
- `build_fig7.py` L19–20 — the docstring quoting the published values; the code asserts them.
- `build_fig5.py` — clean.

Thresholds that appear as named constants (`ANI_SPECIES = 95.0`, `AAI_GENUS = 70.0`,
`ANI_REPORTED = 80.0`, `Q_STRONG/Q_WEAK`, `MIN_PER_GENUS`, `MIN_MAGS_WITH_MODULE`) are decision
rules, not measurements, and every one of them is printed on the page it governs.

## Layout notes worth reusing

- `st.audit` compares text against text only. On Fig 4 a colour bar was initially drawn **on top of**
  a long rotated tick label and the audit still passed; the fix measures the deepest tick-label
  bottom with the real renderer and places the colour bars below it.
- A log axis emits decade ticks outside the drawn range whose labels land off-canvas
  (`FixedLocator` + `NullLocator` on the minors fixes it), and on a scatter panel the first x tick
  and first y tick print on top of each other in the spine corner (`MaxNLocator(prune="both")`).
- Mathtext (`10$^3$`) emits a second `font-family` into the SVG and fails the Arial-first check;
  use the Unicode superscript instead.
- MAGs that project onto the same point need their labels fanned; Fig 7 clusters hero markers
  within 28 display px and steps the labels by 7 pt.
