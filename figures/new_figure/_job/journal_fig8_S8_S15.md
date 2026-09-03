# Build journal — group `fig8_S8_S15` (main Fig 8, Suppl Figs S8–S15)

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
