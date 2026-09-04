# Build note - consolidated Fig 5 and Fig S1 (2026-09-04)

Scripts (new files only; no existing builder, helper or `_style.py` was touched):

* `/data/data/Upcycling/new_figure/build_v2_fig5.py`   -> `new_figure/figures_v2/Fig5.svg/.pdf/.png`
* `/data/data/Upcycling/new_figure/build_v2_supS1.py`  -> `new_figure/figures_v2/Fig_S1.svg/.pdf/.png`

Interpreter for both: `/home/soon/miniconda3/envs/dram_env/bin/python`.

## Panel provenance - main Fig 5

| new panel | old source | content |
|---|---|---|
| A | `build_fig5.py` panel A | closest-GTDB-reference ANI per MAG, ranked, with the 95 % species boundary, plus the block of 21 MAGs with no species-level ANI |
| B | `build_fig8.py` panel C | MGnify livestock species-cluster prevalence of the MICP gene-complete profile and of the single-contig ureC + CA architecture, per biome and pooled |
| C | `build_fig8.py` panel D | antiSMASH BGC class means, MICP-complete vs rest, with the Mann-Whitney P stat column |
| D | `build_fig7.py` panel A | Jaccard PCoA of the Panaroo pan-genome coloured by waste source, MICP-complete ringed, PERMANOVA stat block |

Dropped as instructed: old `build_fig5.py` panels B (S13/S16 reciprocal-best AAI) and C
(skani vs the 63 RefSeq *Sphingobacterium* genomes), and old `build_fig7.py` panel B
(z-scored trait-module ordination). Their numbers survive in Table S4b, Table S10a/b and
Table S2b respectively. Dropping fig7 panel B also removes the only skbio/PERMANOVA
recomputation from this page, so the new script has no permutation step.

Sources read: `Table_S8_novelty_ANI_screen.csv`, `Table_S14a_mgnify_catalog_summary.csv`,
`Table_S23b_antismash_hero_vs_rest.csv`, `Table_S9a_PCoA_coordinates.csv`,
`Table_S9b_PERMANOVA_global.csv`.

Asserts carried over verbatim:
* A - the novelty flag equals exactly "no species-level ANI"; every ranked ANI >= 95.
* B - per-catalog `pct_MICP_gene_complete` and `pct_MICP_single_contig` are **recomputed**
  from `n_MICP_gene_complete` / `n_MICP_single_contig_ureC_CA` over `n_species_clusters`
  and asserted (`atol = 1e-3`) against the stored columns; the pooled bar is then derived
  from the same counts, not from a stored value.
* C - all ten BGC classes must be present in Table S23b.
* D - the `Hero` flag of Table S9a must equal the six MICP-complete MAGs, and the four
  waste sources must be exactly the four expected categories.

`_job/check_no_hardcoded.py` flags one literal in the script, `0.001` on the line
`txt = f"{p:.0e}" if p < 0.001 else f"{p:.2f}"`. It is a P-value *formatting* threshold
carried over unchanged from `build_fig8.py`, not a measurement.

## Panel provenance - Suppl Fig S1

| new panel | old source | Table S2b block |
|---|---|---|
| A | `build_supS1.py` | `Biofilm_EPS::` (10 subcategories) |
| B | `build_supS2.py` | `Ammonia_N::` (6) |
| C | `build_supS3.py` | `Alkaline_Osmo::` (6) |
| D | `build_supS4.py` panel A | `CAZyme_proxy::` (5) |
| E | `build_supS5.py` | `MetalResist_AMR::` (6) |

The dbCAN panels of old `build_supS4.py` (panels B and C) are NOT on this page; per
DESIGN.md they belong to the consolidated Fig S2, built separately.

Construction is unchanged and still goes through `_supp_traits.py`: genus mean over genera
with >= 2 MAGs, the two MICP-complete-containing genera first and the rest by descending
row sum, colour on `log10(value + 1)`, cell labels carrying the untransformed mean. Each
panel re-derives its table from the per-MAG rows of Table S2b and asserts the aggregate
against those rows, exactly as the five old scripts did; the new script additionally
asserts that all five blocks span the same nine genera.

### Colour scale: NOT shared - five colour bars

The old scripts each set `vmax` to the maximum of their own log10 matrix, and those maxima
differ by a factor of five:

| panel | block | log10 vmax | raw max (mean hits per 10^3 CDS) |
|---|---|---|---|
| A | Biofilm_EPS | 2.530 | 338.2 |
| B | Ammonia_N | 0.483 | 2.04 |
| C | Alkaline_Osmo | 1.498 | 30.5 |
| D | CAZyme_proxy | 1.003 | 9.06 |
| E | MetalResist_AMR | 1.403 | 24.3 |

Forcing one scale would flatten panels B and D to near-white, so the per-panel scaling of
the old scripts is preserved and **each panel carries its own colour bar** with its own
`log₁₀(mean hits per 10³ CDS + 1)` caption. The row order also differs between panels,
because the descending-row-sum rank is computed within each block (only the first two
rows, Sphingobacterium and Pseudomonas_E, are fixed).

## Layout and gate results

| page | width | height | audit | prose_scan |
|---|---|---|---|---|
| Fig5 | 180.0 mm | 218.0 mm | 109 text items, 0 overlaps, 0 outside, 0 axes overlaps -> PASS | no sentence-like text -> PASS |
| Fig_S1 | 180.0 mm | 227.6 mm | 410 text items, 0 overlaps, 0 outside, 0 axes overlaps -> PASS | no sentence-like text -> PASS |

Both passed on the rendered output as well: SVG page width 180 mm for both, and the only
font-family declaration in either SVG is the Arial-first stack.

### Fig S1 height revision - 278.0 mm -> 227.6 mm (one-page ceiling 235 mm)

The first lay-out at the old 6.0 mm row pitch came out 278.0 mm, 43 mm over the printed
page. Two levers were used, in the order the coordinator listed; the third was not needed
and would not have helped:

1. **Row pitch 6.0 -> 4.8 mm** (`H_ROW`). Nine genus rows per panel, three panel bands on
   the page, so this alone returns 3 x 10.8 = 32.4 mm. A 7 pt genus label needs about
   3.0 mm of line height and a 6.5 pt cell value about 2.3 mm, so 4.8 mm still clears
   both; the audit confirms no text pair collides.
2. **Column-label wrapping in band A** (`WRAP_CHARS = 15`, applied only where the cell
   pitch is at least `WRAP_MIN_PITCH = 8.0` mm, which is panel A alone at 11.2 mm). The
   34 mm label band of panel A was set by the single label
   "polysaccharide intercellular"; wrapped at its space it becomes two lines of about
   15 mm, so `HDR_A` drops 34 -> 18 mm. `HDR_BC` and `HDR_DE` go 26 -> 25 mm, matching the
   measured label lengths. The wrap inserts a newline at an existing space - the label
   text itself is unchanged, and the half-panel labels are left on one line because their
   6.0 mm cell pitch cannot hold a two-line rotated block without collision.
3. **Re-packing was not applicable.** Panel A carries 10 columns, the most on the page,
   so it cannot join a pair; the existing A / B|C / D|E packing is already the tightest
   three-band arrangement.

Nothing was dropped: all five panels, all nine genus rows in each and every module column
(10 / 6 / 6 / 5 / 6) are still on the page, with the per-panel colour scales and per-panel
colour bars unchanged, the palette and the 7 / 8 / 10 / 6.5 pt type scale unchanged, and
no titles or prose. New height 227.6 mm = 4 + 18 + 43.2 + 8 + 25 + 43.2 + 8 + 25 + 43.2
+ 10. Fig 5 was not touched and stays at 218.0 mm.

Panel letters run in reading order: Fig 5 A (full-width row 1), B (row 2 left), C (row 2
right), D (row 3 left with its legend to the right); Fig S1 A (full-width row 1), B / C
(row 2), D / E (row 3).

## Discrepancies between a source table and what an old script drew

Carried over from `build_fig5.py` (already documented in `_job/journal_main_4_7.md`, still
true of the consolidated panel A):

1. **The shipped panel a was titled "21 candidates < 95 % ANI", which the table does not
   support.** No MAG in Table S8 has an ANI below 95 (the ranked block spans 90 MAGs, all
   >= 95). The 21 novelty candidates are exactly the 21 MAGs for which GTDB-Tk returned no
   species-level ANI at all, S13 and S16 among them. The new panel A draws them as a named
   identifier block and the assert states the equivalence, so the page carries no
   "< 95 %" claim.

Noted for completeness, affecting panels that were dropped rather than redrawn:

2. Old `build_fig5.py` panel C: the Table S10a matrix stores 0.0 for reference pairs skani
   never aligned; the shipped figure plotted those zeros outside the y-limits so they read
   as absent rather than unreported. Panel dropped here; the caveat now only matters for
   Table S10a.
3. Old `build_fig7.py` panel B: the trait-module ordination has **no stored coordinates
   anywhere in the repository** - the old script had to recompute it. With that panel
   dropped, Fig 5 D now uses only the stored pan-genome coordinates of Table S9a and
   recomputes nothing.
4. Old `build_supS4.py` panels B/C: Table S6c (dbCAN family counts) contains only the 105
   non-MICP-complete MAGs, which is why the superseded Figure_S4c rendered empty. Those
   panels are not on this page; whoever builds the consolidated Fig S2 inherits that
   caveat.

No new discrepancy was found in the four Fig 5 panels or the five Fig S1 blocks: every
assert of the old scripts was kept and all of them pass.
