# Build note — consolidated Fig S4 and Fig S5 (2026-09-04)

Scripts (new, nothing existing was modified):

* `/data/data/Upcycling/new_figure/build_v2_supS4.py` → `new_figure/figures_v2/Fig_S4.{svg,pdf,png}`
* `/data/data/Upcycling/new_figure/build_v2_supS5.py` → `new_figure/figures_v2/Fig_S5.{svg,pdf,png}`

Interpreter for both: `/home/soon/miniconda3/envs/dram_env/bin/python`.
Both import `_style as st` (and `_grp_supp_hi as gh`), keep the palette, the 7 / 8 / 10 / 6.5 pt
type scale, spines left+bottom, no panel titles and no explanatory sentences on the page.

---

## Fig S4 — panel provenance

| new panel | old page/panel | content |
|---|---|---|
| A | S6 | pairwise skani ANI among the six MICP-complete MAGs, values printed, em dash where skani reported no alignment |
| B | S7 panel A | ureC maximum-likelihood gene tree, 46 tips, UFBoot ≥ 70 shown, MICP-complete tips coloured by lineage |
| C | S18 panel A | ANI of the six against the organism-verified MICP reference panel (C5 v2) |
| D | S18 panel B | the 8 reported MAG–reference pairs, ranked, against the 95 % species boundary |
| E | S10 | MICP feature prevalence within the 146 Pseudomonas_E reference genomes |

Page: **180 × 230 mm**. Layout: A top-left, B (the tree) in the right column running the full
height of the left column's two stacked panels, C below A, then D and E side by side across the
bottom. Reading order left→right then top→bottom therefore gives A, B, C, D, E.

### How the S18 two-panel question was resolved
**Not merged.** The old S18 panels A and B are kept as two separate lettered panels (C and D)
and the Pseudomonas_E prevalence panel shifted to E — the alternative the task allowed.
Reason: the reference matrix needs a column of full organism binomials (longest
"Sphingobacterium paramultivorum", "Pseudoalteromonas haloplanktis" ≈ 36 mm at 6.5 pt) and the
ranked-pair strip needs its own label column of "MAG vs organism" strings (≈ 39 mm) plus an ANI
axis. Side by side they need ~175 mm of the 180 mm width in a band that also has to sit beside
a 46-tip tree; nothing legible remained. Abbreviating the genus was rejected because the panel
carries four genera beginning with P (Pseudomonas, Pseudoalteromonas, Proteus, Priestia) and
two with S (Sporosarcina, Sphingobacterium).

Only change to the old S18 panel A drawing: it is **transposed** (20 reference rows × 6 MAG
columns instead of 6 × 20) so that the organism names read horizontally in the 84 mm left
column instead of rotated 90°. Same cells, same values, same dashes, same colours.

### Data route — unchanged
* A: `Table_S4a_skani_ANI_matrix_111MAGs.csv`, symmetry and diagonal asserted; a stored 0.0
  is drawn as an em dash (no alignment), never as "0 % identity". **This honesty rendering is
  kept**, in A (em dash) and C (dash), exactly as the old pages did it.
* B: `Table_S7c_ureC_gene_tree.newick`, ladderized, drawn from clade branch lengths.
* C, D: the **rebuilt, organism-verified** C5 panel is used exactly as build_supS18.py did —
  `research/additional/C5_panMICP_env_v2/reference_manifest.tsv` (all 20 rows asserted
  `status == verified`), `skani_panMICP.matrix` (declared size and symmetry asserted), and
  `Table_S20_skani_hero_vs_refs.tsv` asserted pair-by-pair against the matrix, including the
  check that every reported pair is within the MAG's own genus.
* E: `Table_S13a` + `Table_S13b`, same accession set asserted, every percentage recomputed
  from the 0/1 columns.

### Deviation to flag
Old S7 **panel B** (SH/AU topology test + Robinson–Foulds distance vs the species tree) is not
drawn: the task specifies the tree alone for panel B. DESIGN.md describes S4 B as "ureC gene
tree with species-tree congruence", so the congruence numbers are currently unplaced. The
script still reads `Table_S7a_RF_distance.txt` and `Table_S7f_SH_AU_test.iqtree`, still runs
every assertion the old script had, and prints the values, so they remain verified and can be
quoted in the legend:
`RF 50/86 (normalised 0.581), 46 taxa | p-AU ureC tree 1, species tree 2.51e-37 | tree 2 excluded at 95 %`.

---

## Fig S5 — panel provenance

| new panel | old page/panel | content |
|---|---|---|
| A | S12 | predicted minimum doubling time by group (gRodon2), **log y axis**, group n in the tick label |
| B | S19 panel A | length-weighted mean contig coverage per MAG, MICP-complete vs rest |
| C | S19 panel B | the same coverage by waste source, the six MICP-complete MAGs ringed inside their own source column |
| D | S17 panel A | genome-wide GC3, MICP-complete vs rest |
| E | S17 panel B | effective number of codons, same comparison |

Page: **180 × 146 mm**. Row 1 = A, B, C; row 2 = D, E.

Old S17 panels C (codeml M0 ω) and D (yn00 pairwise ω) are **not** on this page — they move to
a main figure built elsewhere, as instructed. Nothing else from S17 was dropped.

Deliberate choices kept: A's logarithmic axis and its "MICP-complete / n = 4", "rest / n = 81"
tick labels; B and C logarithmic; C marks the six by outline rather than a fifth fill colour.

---

## Audit results (acceptance gate)

Both scripts end with `st.audit(fig)`, `st.prose_scan(fig)`, `st.save(fig, OUT, stem)`.

```
Fig_S4  audit: 331 text items | overlaps 0 | outside 0 | axes overlaps 0 -> PASS
        prose_scan: no sentence-like text -> PASS
Fig_S5  audit:  70 text items | overlaps 0 | outside 0 | axes overlaps 0 -> PASS
        prose_scan: no sentence-like text -> PASS
```

Independent check on the saved files (the checks `_job/verify_outputs.py` applies): SVG/PDF/PNG
present for both stems, page width 180.0 mm, every `font-family` Arial-first, panel letters
A–E present, no sentence-like text in the SVG.

`_job/check_no_hardcoded.py` flags two literals in `build_v2_supS4.py`
(`max_depth * 0.035`, `max_depth * 0.012`) — the same two it flags in the old
`build_supS7.py`; both are tree-x layout offsets, not measurements. `build_v2_supS5.py` is clean.

Two layout iterations were needed: the tree's "bootstrap ≥ 70 shown" cue collided with a tip
label (moved up beside the scale bar), and the `n = 6` group count collided with the two-line
"MICP-\ncomplete" tick label on S5 B/D/E (offset moved from −0.16 to −0.21/−0.22 axes units).

---

## Discrepancies between a source table and what an old script drew

1. **🔴 `build_supS19.py` can no longer run.** Its swap-aware assertion block fails:
   `AssertionError: ('swine', 15, 32.0)`. When S19 was rendered (`figures/Fig_S19.png`,
   2026-09-04 02:31) `Table_S21a/b` carried the sheep↔swine label swap the script corrected
   for. Both tables were re-emitted at 03:17 today
   (`SUBMISSION/_revision_260904/tables/TABLE_REEXPORT_LOG.md`) and now agree with the
   sample-code key, so the script's `STORED_LABEL` swap map is stale and its asserts fire.
   *Effect on the drawn panel: none* — S19 always derived the source from the MAG identifier,
   so old S19 panel B and new S5 panel C are the same plot. Only the cross-check route changed.
   In `build_v2_supS5.py` the mapping is no longer assumed in either direction: each recomputed
   group is matched to the stored row that reproduces its count, mean, median and SD, and the
   build asserts that the matched row carries that group's own name — so a relabelling in
   either direction now fails the build. Recomputed groups (all asserted against S21b):
   cattle 22 / 5.85×, swine 15 / 57.43×, sheep 32 / 20.32×, poultry 42 / 28.63×.
   **Action for someone else:** `new_figure/build_supS19.py` still needs fixing or retiring, and
   `SUBMISSION/02_Figure_legends.md` (Fig S19 legend) plus `01_Manuscript.md` §3 still describe
   the swap in the present tense ("the `source` column of the shipped Table S21a assigns
   'sheep' to the M-prefixed MAGs") — that sentence is now false of the shipped table, though
   the "swine, mean 57.4×, n = 15" statement it protects is correct.

2. **Not a defect, but it looks like one.** In S4 panel E the first three bars are all
   137/146 (93.8 %). `Table_S13a` really does have `UreC_alpha == UreB_beta_gamma` on every one
   of the 146 rows, so `urease_core` and `ureC + ureB on one contig` land on the same 137
   genomes. Verified in the table, not a plotting error; a reviewer may still query it.

3. Old S12's filter behaviour is unchanged and still worth stating in the legend: the ≥ 10
   ribosomal-protein anchor filter removed **S13 and S16**, so S5 panel A shows 4 of the 6
   MICP-complete MAGs (4 + 81 = 85 rows in `Table_S16`). The script asserts that some MAGs were
   dropped, so a silent change in the filter would be caught.

4. No discrepancy found for S6, S7, S10, S17A/B or S18: every stored value those scripts
   asserted against still reproduces (ANI matrix symmetry, RF/AU parse, Table S20 pair table,
   `codon_usage_hero_vs_rest.csv` means and P values, the S13a/S13b accession sets).
