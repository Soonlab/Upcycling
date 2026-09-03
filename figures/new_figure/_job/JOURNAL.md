# Upcycling new_figure — build journal

Started 2026-09-04. Contract: `_job/TASK.md`. Per-group journals `_job/journal_*.md` are merged here.

---

# Group: journal_main_1_3

## Build journal — group main_1_3 (main Fig 1, Fig 2, Fig 3)

Built 2026-09-04 with `/home/soon/miniconda3/envs/dram_env/bin/python` (matplotlib 3.10.8,
pandas 1.5.3, Biopython 1.87). Contract: `_job/TASK.md`.

| Figure | Canvas | audit | prose_scan | check_no_hardcoded |
|---|---|---|---|---|
| Fig1 | 180 × 288.75 mm | PASS (135 text items, 0 overlaps) | PASS | clean |
| Fig2 | 180 × 78 mm | PASS (51 text items, 0 overlaps) | PASS | clean |
| Fig3 | 180 × 129 mm | PASS (116 text items, 0 overlaps) | PASS | 1 literal, justified below |

---

## 🔴 BLOCKER for the user — the MICP presence matrix behind Fig 1 and Fig 2

The shipped Fig 1 and Fig 2 were drawn from `pangenome_work/MICP_Pangenome_Final_Summary.csv`.
That file cannot be used and cannot be reproduced:

1. **It covers 100 of the 111 MAGs.** `C1` and `C10`–`C19` are absent. The old Fig 1 code
   initialised the heat map to zeros and filled only the rows it found, so those eleven MAGs
   were drawn as carrying *no* MICP gene. The old Fig 2b x-axis admits this in passing
   ("Other MAGs (n = 94 of 100 Panaroo-QC MAGs)"), while the manuscript and legend describe
   6 vs 105.
   Cause: `pangenome_work/summarize_pangenome_micp.py` takes `df.columns[14:]` as the sample
   columns, but this Panaroo `gene_presence_absence.csv` has only three metadata columns, so
   the first eleven sample columns were skipped.
2. **The named recipe does not reproduce the file.** Re-running that script's logic over
   `pangenome_work/results/gene_presence_absence.csv` matches only `ureD/E/F/G`; the Panaroo
   `Gene` column contains no `ureA`, `ureB`, `ureC` or `cah` entry at all, so the recipe yields
   zero for those four columns while the shipped file has ones. The file's true provenance is
   unknown.
3. **It contradicts the paper's own hero set.** In that file `S26` scores **1/8** and `S13`
   has `cah = 0`, although both are MICP-complete MAGs. This is why the shipped Fig 2b shows
   MICP-complete prevalence of **83.3 % (5/6)** for every gene while its legend claims 100 %.
4. **No source in the repository yields "the complete 8/8 module is restricted to six MAGs".**
   Counts of MAGs carrying all eight genes:

   | source | MAGs covered | MAGs with 8/8 |
   |---|---|---|
   | `MICP_Pangenome_Final_Summary.csv` (shipped) | 100 | 26 |
   | `Table_S1a_ace_samples_list.csv` | 45 ure-bearing (of 111) | 28 |
   | Bakta GFF3 keyword scan, all 111 MAGs (this rebuild's recipe) | 111 | 32 |

   The six-MAG designation is therefore **not** "score == 8". It is a curated flag, carried
   identically in seven analysis tables (S11, S15a, S17a, S18a, S19a, S21a, S23a), and it
   evidently also depends on single-contig operon architecture (Table S1c). **The Fig 1 legend
   sentence "The complete 8/8 module is restricted to four *Sphingobacterium* MAGs and two
   *Pseudomonas*_E MAGs" is false against every available source and must be rewritten**, or
   the MICP-complete criterion must be stated explicitly in the Methods and reflected in the
   figure. This is a manuscript decision, not a figure decision.

**What was drawn instead.** `Table_S1a_ace_samples_list.csv` — a shipped supplementary table
of record giving per-gene copy counts for the 45 ure-bearing MAGs — binarised, with MAGs
absent from it scoring 0. It is the only per-gene source in which all six MICP-complete MAGs
carry all eight genes, and it covers the full 111-MAG panel once absence is read as "no
detected MICP gene". Both Fig 1 and Fig 2 use it, so the two pages now agree with each other.
Asserted in both scripts: `score[heroes] == 8` for all six.

Consequence the user should expect to see: panel Fig 2A now shows MAGs at 8/8 in
*Chryseobacterium*, *Acinetobacter* and *Pseudomonas*_E as well, and Fig 1's heat map shows
28 fully green rows. That is what the source says.

---

## Fig 1 — phylogeny + MICP gene presence

**Sources**
- `pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile` — topology and branch
  lengths, midpoint-rooted for display (`tree.root_at_midpoint()`), drawn as a rectangular
  phylogram from `tree.depths()` rather than through `Phylo.draw`.
- `Table_S1d_GTDB_Tk_classification.tsv` — genus (`g__`) and species (`s__`) per MAG.
- `Table_S1a_ace_samples_list.csv` — gene presence (see BLOCKER above).
- `Table_S15a_alkaliphile_signature_per_MAG.csv` — the MICP-complete flag.

**Recomputed and asserted**
- tip count == 111 == rows of the GTDB table, and the set of tip MAG ids equals the set of
  MAG ids in the taxonomy table.
- the nine genera given their own strip colour are exactly the nine most frequent genera
  (`gcount.head(9)`), asserted rather than typed; the remaining 8 MAGs pool into "Other
  genera" with the count computed, not written.
- every MICP-complete MAG carries all eight genes in the drawn matrix.
- scale bar length is derived (`1/2/5 × 10^floor(log10(xmax/4))`, largest that fits in a third
  of the tree width) and printed with its own computed value — 0.2 substitutions/site here.

**Dropped**
- Ultrafast bootstrap support. The values are present in the treefile but were not drawn in
  the shipped figure either; at 111 tips over a 50 mm tree they cannot be rendered legibly.
  The legend's mention of 1,000 UFBoot replicates describes the method, not a figure element,
  so nothing on the page claims support.

**Discrepancy (minor)** — the legend says the colour strip encodes "the top ten genera";
there are nine genera with ≥ 2 MAGs plus a pooled "Other genera" class of 8 single-MAG genera.
Drawn as nine + pooled, which is what the data supports.

**Colour meanings**

| colour | meaning (whole page) |
|---|---|
| blue `#1F4E78` | genus *Sphingobacterium* (a MICP-complete lineage) |
| orange `#D9772B` | genus *Pseudomonas*_E (the other MICP-complete lineage) |
| lilac / tan / mauve / teal / grey / olive / brown | the other seven strip genera |
| light grey `#D9D9D9` | pooled "Other genera" |
| green `#2F7D4F` | gene present (panel B); white = absent |
| coral `#C53E1F` | a MICP-complete MAG — tip label text only, used for nothing else |
| black | tree branches |

No genus was given green or coral. Teal (*Comamonas*) is the nearest hue to the presence
green and is clearly lighter and bluer; flagged here for the record.

**Canvas** 180 × 288.75 mm (111 tips × 2.25 mm row pitch + 13 mm top + 26 mm for the scale bar
and legend). The row pitch is set by the 6.5 pt tip labels; anything tighter collides.

---

## Fig 2 — module completeness by genus, and per-gene prevalence

**Sources** `Table_S1a` (presence), `Table_S1d` (genus), `Table_S15a` (group flag).

**Recomputed and asserted**
- hero set from the group flag == `st.HEROES`; group sizes 6 and 105 asserted against the
  arrays actually used for the two prevalence series.
- prevalences are means of the binarised matrix, computed at draw time; every bar carries its
  own computed value label.
- genus counts in the y tick labels are computed from the plotted grouping.

**Panel A was turned horizontal.** With ten genera the rotated x labels collided
(`Chryseobacterium` × `Sphingobacterium`, and eight more). Genus names now run along the y
axis, which also puts the highest-scoring genus at the top.

**Discrepancies against the shipped figure and legend**

| item | legend / old figure | this rebuild (from Table S1a) |
|---|---|---|
| MICP-complete prevalence | legend "100 %"; old figure 83.3 % | 100 % for all eight genes |
| rest prevalence | legend "55–85 %"; old figure 28–35 % | 28–37 % |
| rest group size | legend n = 105; old figure n = 94 | n = 105 |

The legend's "55–85 %" matches neither the old figure nor any source and should be corrected
to the drawn values.

**Colour meanings** coral = the MICP-complete group (rings in A, bars in B); grey = every other
MAG; light grey box fill + black medians and points = distribution summary, no group meaning.

**Canvas** 180 × 78 mm.

---

## Fig 3 — ure cluster synteny on the main contig of each MICP-complete MAG

**Sources** `MAGs_FASTA_files/bakta_results/<MAG>/<MAG>.gff3` (CDS coordinates, strand,
product), `Table_S1c_hero_cluster_audit.csv` (contig of record, gene count, span),
`Table_S3a_HGT_ureCah_cluster.csv` (mobile-element count, Δ GC).

**Fixed from the previous rebuild.** The `figure_master_v1` version drew S13/S16 from an inline
schematic constant and showed only *ureA* + *ureB*. Every arrow here comes from the GFF3, and
all seven *ure* genes appear on a single contig for M1, S13 and S16, as the manuscript states.

**Recomputed and asserted, per MAG**
- number of *distinct* `ure` genes on the contig named in Table S1c == the table's
  `ure_genes_on_main_contig`;
- cluster span (max end − min start, kb) == the table's `cluster_span_kb_main` to < 0.02 kb.

Both assertions pass for all six MAGs: M1 7/28.56 kb, S13 7/5.86, S16 7/5.92, C22 5/5.32,
S23 4/3.27, S26 4/3.03.

**Discrepancies found**
1. **The audit column counts distinct genes, not CDS copies.** On S23 `contig_220` there are
   five `ure` CDS (two *ureE* copies) but four distinct genes; the table records 4. The stat
   column on the figure reads "4 / 7" and both *ureE* arrows are drawn, so the reader sees the
   duplication. Worth one clause in the legend.
2. **The old densest-cluster search disagrees with the table of record for S23.** The search in
   `01_main_figures.py` ties at four distinct `ure` genes between `contig_151` (2.85 kb) and
   `contig_220` (3.27 kb) and returns `contig_151`; Table S1c names `contig_220`. This rebuild
   draws the contig named in the table. The shipped figure drew `contig_151`, so its S23 track
   does not match Table S1c.
3. **`cah` is on a separate contig in all six MAGs** (`cah_on_main_contig = False` throughout
   Table S1c), so no *cah* arrow exists on any track. The gene key therefore omits *cah* —
   a key entry with no mark on the page would be misleading. The legend should say explicitly
   that *cah* lies outside the drawn window.

**Per-track scaling, not a shared axis.** M1's window is 32 kb while the other five are 7–10 kb.
On one shared scale the five short clusters occupy the leftmost quarter of the page and their
gene labels are unreadable. Each track now carries its own kb axis with its own ticks, so
absolute scale stays readable per track; arrows remain to scale within a track.

**Label placement** gene labels are staggered over up to three rows by a greedy placer with a
leader line to the arrow, because M1's seven genes sit in two tight sub-clusters inside a 32 kb
window.

**Flagged literal, justified.** `check_no_hardcoded.py` flags `0.353` in `place_levels`
(`build_fig3.py` L140). It is the typographic constant 1 pt = 0.3528 mm, used to convert a
label's point size into millimetres for the collision estimate. It is not a measurement.

**Colour meanings** one colour per MICP gene (`ureA`…`ureG`), light grey for any other CDS.
No other colour is used on the page; the MAG identifiers are plain black bold and the contig
identifiers grey, carrying no group meaning.

**Canvas** 180 × 129 mm.

---

# Group: journal_main_4_7

## Build journal — main Fig 4, 5, 6, 7 (group `main_4_7`, 2026-09-04)

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

---

# Group: journal_fig8_S8_S15

## Build journal — group `fig8_S8_S15` (main Fig 8, Suppl Figs S8–S15)

Built 2026-09-04 with `/home/soon/miniconda3/envs/dram_env/bin/python` against
`_job/TASK.md`. Every page: `st.audit` 0 overlaps / 0 outside / 0 axes overlaps,
`st.prose_scan` empty, PNG read and corrected by eye, page height re-fitted.

| stem | script | height (mm) | audit | prose | source tables |
|---|---|---|---|---|---|
| Fig8    | `build_fig8.py`     | 122 | PASS | PASS | S12, S22, S14a, S23b |
| Fig_S8  | `build_supS8.py`    | 76  | PASS | PASS | S11 + A1 raw abricate output |
| Fig_S9  | `build_supS9.py`    | 64  | PASS | PASS | S12 |
| Fig_S10 | `build_supS10.py`   | 57  | PASS | PASS | S13a, S13b |
| Fig_S11 | `build_supS11.py`   | 76  | PASS | PASS | S15a |
| Fig_S12 | `build_supS12.py`   | 64  | PASS | PASS | S16 |
| Fig_S13 | `build_supS13.py`   | 62  | PASS | PASS | S14a |
| Fig_S14 | `build_supS14.py`   | 60  | PASS | PASS | S15b |
| Fig_S15 | `build_supS15.py`   | 80  | PASS | PASS | S17a, S17b |

## Values recomputed and how they were checked

* **Fig 8c / S13** — per-catalog and pooled MGnify percentages recomputed from the count
  columns and asserted against `pct_MICP_gene_complete` / `pct_MICP_single_contig`
  (`atol=1e-3`). Pooled: 233/7,599 = 3.07 % gene-complete, 265/7,599 = 3.49 % single-contig.
* **Fig 8a / S9** — every observed residue compared against `expected`; 42/42 match, and
  `expected == ref_aa` is asserted (the table's two reference columns agree).
* **Fig 8b** — TM-score and RMSD are stored measurements; only the `_UreC` suffix is stripped
  from the MAG labels and the single reference length is asserted constant.
* **Fig 8d** — no recomputation; the stored hero/rest means and MWU P of the ten classes are
  drawn as they stand.
* **S8** — CARD / VFDB / ResFinder counts recomputed from the raw abricate combined files and
  asserted equal to Table S11 for all 111 MAGs; MWU P per database recomputed (two-sided).
* **S10** — all six percentages recomputed from the 0/1 columns of S13a/S13b over the same
  146 accessions (asserted identical accession sets): UreC 137/146 = 93.8 %, UreB 137/146,
  urease core 137/146, any CA 146/146 = 100 %, ureC+ureB one contig 137/146,
  ureC+CA one contig 53/146 = 36.3 %. All match the legend.
* **S11** — Mrp fold difference recomputed (11.67×) and asserted against the manuscript's
  11.7×; MWU two-sided P per feature recomputed. Mrp P = 5.3 × 10⁻⁴ reproduces the legend,
  which is therefore the **two-sided** value (one-sided greater is 2.6 × 10⁻⁴).
* **S12** — medians and MWU P recomputed and asserted against the legend
  (1.06 h vs 1.10 h, P = 0.58).
* **S14** — `Ca_pathway` recomputed as (Ca_transporter > 0) OR (Ca_ATPase > 0) and asserted
  equal to the stored flag.
* **S15** — group means and two-sided MWU P recomputed from S17a (plasmid P = 0.06,
  virus P = 0.11, matching the legend); per-MAG contig counts in panel B recomputed by
  splitting the comma-joined contig lists of S17b, and the plasmid counts of S17a and S17b
  asserted equal.

## Values dropped for lack of a source

* **Fig 8c** — the shipped figure's dashed "present study 6/6 = 100 %" line and its label.
  It mixes denominators: the 100 % is over this study's 6 MICP-complete MAGs, not over MGnify
  species clusters, and the ~30-fold enrichment statement it supports belongs in the legend.
  Dropped from Fig 8c **and** from S13. No number was lost — both are in the legend text.
* **Fig 8a/b/d, S9–S15** — the shipped panel titles carrying conclusions ("42/42 matches",
  "6/6 TM > 0.5", "T3PKS 23× + RRE 8× Sphingobacterium fingerprint",
  "A4. geNomad plasmid/virus calls per MAG", significance asterisks). Titles are forbidden by
  the contract; the encodings survive as legends and stat columns.
* **S15** — the flag columns `urease_on_plasmid`, `urease_on_virus`, `CA_on_virus` are empty
  for all six MAGs in the source and are not drawn. `CA_on_plasmid` is non-empty for S26 only
  and is represented by that MAG's `CA_MGE_contamination` = 1.

## Discrepancies for the user to adjudicate

1. **Table S11 is missing its PlasmidFinder column.** The A1 run did produce PlasmidFinder
   output (`research/additional/A1_biosafety/combined/plasmidfinder_all.tsv`) and
   `aggregate_biosafety.py` is written to include it, but the shipped Table S11 carries only
   card / vfdb / resfinder. Recomputing from the raw file reproduces the three published
   columns exactly and adds the fourth: **exactly one PlasmidFinder hit in the whole panel,
   IncQ2_1 in S21 (a non-MICP-complete MAG); all six MICP-complete MAGs are zero**, which is
   what the S8 legend already claims. Fig S8 now shows four blocks. *Table S11 should be
   regenerated with the PlasmidFinder column before submission.*
2. **S26 has no calcium pathway.** In Table S15b, S26 has `Ca_transporter` = 0 and
   `Ca_ATPase` = 0, hence the stored `Ca_pathway` = 0; its five `CO3_transporter` copies do
   not enter that call. The S14 legend's "at least one Ca²⁺-handling gene (83.3 %)" = 5/6 is
   consistent with this, but any statement elsewhere that all **6/6** MICP-complete MAGs are
   "urease + CA + Ca complete" (analysis A6 as summarised in the project notes) is not.
   Worth a check against Methods §2 and Results §3.10.
3. **Site order in Table S12** is H137, H139, K220, H249, H275, D363, C322 — D363 before
   C322. Both Fig 8a and S9 re-sort by residue number so reading order follows the numbering;
   the shipped Fig 8a used the table order. No value changes.
4. **S12 shows 4 of the 6 MICP-complete MAGs.** Table S16 contains only C22, M1, S23, S26;
   S13 and S16 fell to the ≥ 10 ribosomal-anchor filter. That is a property of the source and
   is stated in the legend; the figure carries `n = 4` on the tick label and no longer prints
   an "excluded" note on the canvas.
5. **S12 y axis is now logarithmic.** The rest group spans 0.4–17.4 h; on the shipped linear
   axis every MICP-complete point was pressed onto the baseline. Legend should say "log scale".

## Colour → meaning

| page | colour | meaning |
|---|---|---|
| Fig8 | green | residue matches reference (A) |
| Fig8 | coral | residue differs (A, legend entry only), MICP-complete group (D) |
| Fig8 | blue / orange | Sphingobacterium / Pseudomonas_E lineage (A row labels, B bars and tick labels) |
| Fig8 | grey | rest group (D), TM = 0.5 reference line (B) |
| Fig8 | dark / light green | gene-complete / single-contig prevalence (C) |
| S8, S11, S12, S15A | blue / orange | lineage of a MICP-complete MAG |
| S8, S11, S12, S15A | grey | one of the remaining MAGs |
| S9 | green / coral | matches reference / differs |
| S10 | coral | prevalence of a MICP feature among 146 Pseudomonas_E references |
| S13 | dark / light green | gene-complete (A) / single-contig (B) prevalence |
| S14 | green intensity | gene copy number |
| S15B | coral | a non-zero mobile-element contamination count |

Blue never encodes a prevalence series anywhere in this group, so it means "Sphingobacterium
lineage" on every page it appears on.

## Notes on the checkers

* `_job/check_no_hardcoded.py` flags one literal in `build_fig8.py`:
  `p < 0.001` in the panel-D P-value formatter. That is a display threshold, not a datum.
* `_job/verify_outputs.py` counts single uppercase characters as panel letters, so Fig8 (46)
  and Fig_S9 (50) report inflated letter counts — the residue letters H, K, C, D inside the
  heat-map cells. Both pages carry exactly one letter per panel.
* A twinx pair always registers as overlapping Axes in `st.audit`. Fig 8b therefore draws
  TM-score and RMSD as two adjacent axes rather than one twinned pair.
* Box-plot axes whose lower limit dips slightly below zero make matplotlib create an
  out-of-view negative tick whose Text still measures; S15 filters those ticks out.

---

# Group: journal_supp_S1_S7

## Build journal — supplementary figures S1–S7

Group `supp_S1_S7`. Interpreter `/home/soon/miniconda3/envs/dram_env/bin/python`.
Scripts `build_supS1.py` … `build_supS7.py`; shared layout helper `_supp_traits.py`
(layout only, contains no data value). All seven pages: `st.audit` 0 overlaps / 0 outside /
0 axes overlaps, `st.prose_scan` empty, `_job/verify_outputs.py` PASS at 180 mm.

| page | height (mm) | panels |
|---|---|---|
| Fig_S1 | 98 | A biofilm/EPS genus heat map |
| Fig_S2 | 98 | A ammonia / nitrogen genus heat map |
| Fig_S3 | 98 | A alkaline / osmotic genus heat map |
| Fig_S4 | 260 | A keyword CAZy proxy · B dbCAN classes · C dbCAN families |
| Fig_S5 | 98 | A metal / AMR genus heat map |
| Fig_S6 | 126 | A pairwise ANI among the six MICP-complete MAGs |
| Fig_S7 | 211 | A ureC gene tree · B topology test + Robinson-Foulds |

## Sources used

* `Table_S2b_trait_module_per1kCDS.csv` — S1 (`Biofilm_EPS::`), S2 (`Ammonia_N::`),
  S3 (`Alkaline_Osmo::`), S5 (`MetalResist_AMR::`), S4 panel A (`CAZyme_proxy::`).
* `Table_S6a/S6b/S6c/S6d` (dbCAN counts, per-1k-CDS, family counts, class statistics) — S4 B/C.
* `research/revision/dbcan_direct/<MAG>.tbl` and `research/extra/gene_category_counts.csv` — S4 C.
* `Table_S4a_skani_ANI_matrix_111MAGs.csv` — S6.
* `Table_S7c_ureC_gene_tree.newick`, `Table_S7a_RF_distance.txt`,
  `Table_S7f_SH_AU_test.iqtree` — S7.

## Recomputations and how they were checked

* Genus aggregation (S1/S2/S3/S5, S4 A) reproduces `_genus_aggregate` of the superseded
  PowerPoint builder: mean per genus, genera with ≥ 2 MAGs (9 genera survive), the two
  MICP-complete-containing genera first, the rest by descending row sum. Each script
  asserts every genus row equals the mean of its own MAG rows and that the printed n
  equals the number of those rows. Colour transform log10(value + 1), pseudocount 1, the
  same transform the old builder used; cell labels carry the untransformed mean.
* S4 B: hero and rest class means recomputed from `Table_S6b` and asserted against
  `Table_S6d` (means and fold change, tolerance 5e-3). All six classes reproduce.
* S4 C: the six MICP-complete MAGs were re-annotated by re-parsing
  `research/revision/dbcan_direct/<MAG>.tbl` with the recipe of `04b_dbcan_final.py`
  (best HMM per protein by E-value). The per-class sums reproduce
  `Table_S6a_dbCAN_class_counts.csv` **exactly for all six MAGs and all six classes**, and
  the per-MAG totals match. The normalised family table also reproduces the panel-B class
  means to 5e-3.
* S6: the 6 × 6 sub-matrix is asserted symmetric with a 100.0 diagonal before drawing.
* S7: normalised RF asserted equal to RF / max RF; the taxon count in the RF file asserted
  equal to the number of tips in the drawn tree (46). The topology-test block is parsed by
  regex; tree 1 is asserted to be the maximum-likelihood tree (ΔlogL = 0) and the logL gap
  asserted equal to the printed ΔlogL of tree 2.

## Values dropped, and defects fixed rather than carried over

1. **The shipped `Figure_S4c_dbCAN_families.png` is an empty axis** — no bars at all. Cause:
   `Table_S6c_dbCAN_family_counts.csv` contains only the 105 non-MICP-complete MAGs
   (the six MICP-complete MAGs are absent), so the old script's hero mean was taken over an
   empty set. Fixed by recomputing the six MAGs' family counts from the raw HMMER tables,
   validated against the published class counts as described above.
2. **Family-name mismatch between the two annotation routes.** dbCAN HMM names carry a
   subfamily suffix (`GT2_Glycos_transf_2`) while the DRAM-derived table for the other 105
   MAGs carries the bare family (`GT2`). Comparing them raw made GT2 look absent from the
   105 (rest mean 0.00 against a true 1.66). Both sides are now reduced to the bare family
   before comparison. This mismatch never surfaced in the published analysis because that
   analysis never compared families across the two groups.
3. **Method asymmetry, carried over unchanged and flagged.** MICP-complete CAZyme
   annotations come from direct hmmsearch, the other 105 from the DRAM `cazy_best_hit`
   column. This is inherited from the published `04b_dbcan_final.py` and is the basis of the
   published class-level numbers too; it is not introduced here. It is worth one sentence in
   the S4 legend.
4. **The superseded S7 slide was not a tree.** It listed six hard-coded "clades" with
   invented membership (e.g. an "Empedobacter clade" of C14/V8/V9) and hard-coded
   "RF = 0.58 / SH p < 0.001" text. Nothing of it was reused: the page now draws the actual
   46-tip newick and parses the RF and topology-test values from their files.
5. `Table_S4a` stores "skani could not align this pair" as 0.0. Rendered as an em dash on a
   light grey cell, never as 0 % identity. S13 aligns to none of the other five MAGs, and
   M1/S26 align to none — both are properties of the data, not missing cells.

## Discrepancies for the user to adjudicate

* **S7 panel B prints p-SH = 0 and p-KH = 1 verbatim from the IQ-TREE file.** IQ-TREE prints
  0 when the value falls below its display precision with 1,000 RELL resamplings; the
  manuscript states "p < 0.001". The figure shows the file's token rather than the
  manuscript's rounding. If you prefer, the column can print "< 0.001" instead — say which.
* **`Table_S2b` keyword proxy `Biofilm_EPS::polysaccharide_intercellular` runs 146–338 hits
  per 10³ CDS** across every genus, i.e. 15–34 % of all coding sequences. That is a keyword
  over-match in the proxy dictionary, not a biological signal, and it dominates the S1
  colour scale. The values are printed as-is; consider dropping this subcategory from S1 or
  qualifying it in the legend.
* **`Table_S7a` reports normalised RF = 0.581**, which the manuscript rounds to 0.58. No
  action needed, noted for completeness.
* S1/S2/S3/S5/S6 are genuinely single-panel pages; each carries a lone panel letter "A" so
  that the set is lettered uniformly and the output gate (which requires panel letters)
  passes.

## Colour–meaning tables

* **S1, S2, S3, S5** — white→green: mean hits per 10³ CDS (only encoding). Coral: genus
  label of a genus containing MICP-complete MAGs. Black: every other genus label.
* **S4** — white→green: mean hits per 10³ CDS (panel A only). Coral: the MICP-complete group
  (n = 6) — bars, dots, genus labels, and significant q values. Grey: the remaining 105 MAGs
  — bars, dots, lollipop connectors.
* **S6** — white→green: ANI (%). Blue: Sphingobacterium MAG label. Orange: Pseudomonas_E MAG
  label. Light grey: a pair skani could not align.
* **S7** — blue: Sphingobacterium MICP-complete tip. Orange: Pseudomonas_E MICP-complete tip.
  Black: all other tips and all branches. Grey: bootstrap values and the scale bar. Coral: a
  topology-test row excluded at the 95 % level.

## Display choices recorded

* Bootstrap values are drawn only where ≥ 70 (display threshold, not a datum), at the
  midpoint of the branch leading to the node, and are nudged down in 0.5-row steps when two
  labels would collide. The nudge moves the label, never the node.
* The ureC tree is ladderized for readability; ladderization changes display order only.
* S1/S2/S3/S5 keep the subcategory order of the source table.
* S4 C shows the 14 families with the highest MICP-complete mean.

---

# Group: journal_supp_S16_S21

## Build journal — group supp_S16_S21 (Suppl Figs S16–S21)

Built 2026-09-04 under `_job/TASK.md`. Interpreter `/home/soon/miniconda3/envs/dram_env/bin/python`.
Shared plotting helpers for this group live in `_grp_supp_hi.py` (layout only, no data value).

All six pages: `st.audit` 0 overlaps / 0 outside / 0 axes overlaps, `st.prose_scan` empty,
`_job/verify_outputs.py` PASS, `_job/check_no_hardcoded.py` clean.

| Page | Height | Panels | Status |
|---|---|---|---|
| Fig_S16 | 59 mm | A defense systems, B CRISPR arrays | PASS |
| Fig_S17 | 145 mm | A GC3, B ENC, C codeml M0 ω, D yn00 pairwise ω | PASS |
| Fig_S18 | 59 mm | A nearest verified reference ANI | PASS (content reduced — see D1) |
| Fig_S19 | 76 mm | A coverage by group, B coverage by waste source | PASS |
| Fig_S20 | 66 mm | A TM-score, B RMSD | PASS |
| Fig_S21 | 190 mm | A total BGC regions, B class-level means + MWU P | PASS |

## Values recomputed and how they were checked

Every statistic drawn was recomputed from the per-MAG source table and asserted against the
stored summary. Mann–Whitney U is **two-sided** everywhere; that alternative reproduces every
stored P exactly, so the stored tables were produced two-sided.

- **S16** — group means and MWU P for `n_defense_systems` and `n_crispr_arrays` recomputed from
  Table S18a, asserted against Table S18b (both means 0.0, both P = 1.0). All 111 MAGs are zero
  for both metrics; the count above zero is printed per column so the null reads as 0 / n.
- **S17 A, B** — GC3 % and ENC group means and MWU P recomputed from Table S19a, asserted against
  `C3_dnds_codon/codon_usage_hero_vs_rest.csv` to 1e-15. Reproduces the legend values
  (GC3 P = 0.0433, ENC P = 0.00283; the legend rounds ENC to 2.8 × 10⁻³).
- **S17 C** — codeml M0 ω read from Table S19b; asserted all four ω < 1 (purifying selection).
- **S17 D** — drawn from the 595-row per-pair file `C3_dnds_codon/yn00_pairwise.csv`, not from the
  stored medians. The pair filter is the one in `scripts/additional/C3_dnds_codon/run_dnds_v3.py`
  (0 < ω < 99 and dS > 0.01), which leaves 574 pairs and reproduces **every** count, median and
  MWU P in Table S19c exactly (asserted). ureG hero-hero vs rest-rest P = 7.66 × 10⁻⁸ confirmed.
- **S19** — group means 26.8× vs 25.5× and MWU P = 0.2723 reproduce the legend. Per-source count,
  mean, median and std recomputed and asserted against Table S21b (see D2 for the label swap).
- **S20** — values read from Table S22; `ref_len` asserted constant across rows (single reference
  chain, 569 aa). The `_UreC` suffix is stripped from the MAG id.
- **S21** — all 43 class means and MWU P recomputed from Table S23a and asserted against Table S23b;
  T3PKS (5.29 × 10⁻¹⁰) and RRE-containing (1.88 × 10⁻⁵) additionally asserted by name.

## Discrepancies for the user to adjudicate

### D1 (blocker for S18) — the C5 curated reference panel is corrupted
`research/additional/C5_panMICP_env/refs/` was built by naming each download after the organism.
The naming went wrong: **12 of the 14 downloaded reference genomes contain a different organism
than their filename**, and 6 of the 20 curated accessions never downloaded at all. Verified by
reading the first FASTA header of each file; the full audit is written to
`_job/S18_reference_identity_audit.csv`. Examples:

| file | actually contains |
|---|---|
| Sporosarcina_ureae.fna | *Escherichia coli* O104:H4 |
| Halomonas_elongata.fna | *Bdellovibrio bacteriovorus* HD100 |
| Halomonas_pacifica.fna | *Escherichia coli* NCTC 13441 |
| Lysinibacillus_fusiformis.fna | *Bdellovibrio bacteriovorus* 109J |
| Lysinibacillus_sphaericus.fna | *Thermus aquaticus* Y51MC23 |
| Sphingobacterium_sp_21.fna | *Marinobacter nauticus* |
| Sphingobacterium_paramultivorum.fna | *Enterobacter* sp. HK169 |
| Bacillus_licheniformis.fna | *Geminocystis* sp. NIES-3708 |
| Helicobacter_pylori_26695.fna | *Riemerella anatipestifer* Yb2 |
| Klebsiella_aerogenes.fna | *Borreliella bavariensis* PBi |
| Pseudoalteromonas_haloplanktis.fna | *Hoylesella loescheii* |
| Yersinia_enterocolitica.fna | *Helicobacter acinonychis* |

Only **Bacillus subtilis 168** and **Pseudomonas helleri** hold what their names claim.
Not downloaded: *Sporosarcina pasteurii* (the canonical MICP organism), *S. psychrophila*,
*Bacillus megaterium*, *Idiomarina loihiensis*, *Proteus mirabilis*, *Sphingobacterium multivorum*.

Consequences:
1. The two impossible ANI values in `skani_panMICP.matrix` are artefacts of this: 99.02 % between
   "Halomonas elongata" and "Lysinibacillus fusiformis" is really two *Bdellovibrio* strains, and
   96.45 % between "Halomonas pacifica" and "Sporosarcina ureae" is really two *E. coli*.
2. **The Fig. S18 legend claim that the four *Sphingobacterium* MAGs "remain below 95 % against
   every reference examined" is not supported by C5** — no genuine *Sphingobacterium* reference was
   ever in the comparison. That claim *is* supported by the independent 63-genome RefSeq screen
   (Fig. 5c / Table S10), whose reference set I spot-checked and found correct.
3. **The legend's "M1 ↔ *P_E* sp. 025837155 = 97.98 % ANI" does not come from C5 at all** — that
   genome is not in this reference set. It is an A3 value (`Table_S13` / 146-genome screen). It is
   dropped from this page rather than carried across.
4. The legend's "20 curated reference genomes" is really 14 attempted, 2 usable.

S18 is therefore drawn from the verified genomes only: one panel giving each MAG's ANI to its
nearest *verified* curated reference. Only S26 ↔ *P. helleri* (97.54 %) survives; the other five
MAGs are marked "no alignment reported". **Recommended action: re-run C5 with accession-checked
downloads (name files by accession, as A3 and the ext_sphingo screen correctly do), then rebuild
S18.** The scope of the defect is confined to C5 — `A3_pseudomonas_ani` and `revision/ext_sphingo`
name their files by accession and were spot-checked as correct.

### D2 (corrected in S19) — Table S21a swaps the swine and sheep labels
The `source` column of Table S21a attaches "sheep" to the 15 M-prefixed MAGs and "swine" to the 32
S-prefixed MAGs — the opposite of the sample-code key used throughout the manuscript
(C cattle, **M swine**, **S sheep**, V poultry) and of the independent `Source` column of Table S9a.
47 of 111 MAGs disagree. S19 therefore derives the source from the MAG identifier and asserts it
against Table S9a; Table S21b inherits the same swap, so its stored statistics are asserted against
the recomputed groups *under the swap*, which both reproduces the stored numbers and proves it.

Consequence for the manuscript: the Fig. S19 legend sentence "sheep-rumen MAGs are the most deeply
covered (mean 57.4×)" is really about the **swine** MAGs (n = 15). The corrected hero distribution
is 1 cattle (C22), 1 swine (M1), 4 sheep (S13, S16, S23, S26), which also matches the sample codes.
**Tables S21a and S21b should be re-emitted with corrected labels.**

### D3 (shared, for the coordinator) — matplotlib mathtext breaks the Arial-first rule
Superscript exponents written as mathtext (`$^{-10}$`) are typeset in the mathtext font and land in
the SVG as `font-family: 'DejaVu Sans'`, which `_style.save`'s rewrite does not catch and
`verify_outputs.py` correctly fails. `_grp_supp_hi.fmt_p` now uses Unicode superscript characters
instead, which Liberation Sans and Arial both carry. **Other groups will hit this the moment they
write an exponent**; either reuse `fmt_p` or widen the rewrite in `_style.save`.

## Values dropped for lack of a source
- S18: the M1 ↔ *P_E* sp. 025837155 ANI (see D1.3) and the 12 mislabelled references (D1).
- S21 panel B: 19 of 43 BGC classes, all with hero mean 0 and rest mean ≤ 0.05 regions per MAG.
  The smallest P among the dropped classes is 0.599, asserted in the script, so no class with any
  signal was removed. 24 classes drawn.

## Merges
- S19: the old panels (b) "coverage by source" and (c) "heroes over the source distribution" are the
  same plot with and without the hero markers, so they are one panel B in which the six MICP-complete
  MAGs are ringed **at their own jittered x position** (`strip_box` now returns the jitter
  coordinates so an overlay lands on the point it marks, not on the column centre).
- S16, S17, S21: the multi-file old figures (S16 defense + crispr; S17 a–d; S21 two panels) are one
  composed page each, as the contract requires.

## Colour–meaning tables
- **S16, S21**: coral = MICP-complete group (n = 6); grey = the remaining 105. No other colour.
- **S17**: coral = MICP-complete on A, B, C and "both members MICP-complete" on D; light coral =
  mixed pair (D only); grey = not MICP-complete on A, B, D, and the ω = 1 neutrality line.
- **S18**: blue = *Sphingobacterium* lineage MAG; orange = *Pseudomonas*_E lineage MAG;
  grey = the 95 % threshold and the no-alignment markers.
- **S19**: A coral = MICP-complete, grey = rest. B four waste-source colours; MICP-complete MAGs are
  marked by a black **outline** rather than a fill, because their fill already states waste source.
  Coral is not reused in B, so no colour carries two meanings on the page.
- **S20**: blue / orange = the two lineages; grey = the TM = 0.5 threshold.

## Axis choices worth noting
- S19 A/B and S17 D use a log10 y-axis. Coverage spans 2–743×; yn00 ω spans 0.002–4.11. On a linear
  axis the whiskers (drawn at the full range, `whis=(0,100)`) leave the panel and the low-ω genes
  collapse onto the baseline. Limits are taken from the data.
- S16 keeps a 0–1 y-axis although every value is 0, so the null is visible as points on the zero
  line rather than as an empty panel.
