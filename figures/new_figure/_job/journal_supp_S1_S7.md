# Build journal — supplementary figures S1–S7

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
