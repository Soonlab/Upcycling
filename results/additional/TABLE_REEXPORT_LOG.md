# Supplementary table re-export — 2026-09-04

Repairs the four defective supplementary tables found during the figure rebuild
(`/data/data/Upcycling/new_figure/_job/JOURNAL.md`).

Script: `reexport_tables.py` (idempotent; `--check` runs every assertion and writes nothing).
Interpreter: `/home/soon/miniconda3/envs/dram_env/bin/python`. Run log: `reexport_run.log`.
Every replaced file was copied to `originals/` first.

## Summary

| Table | Before | After | Defect |
|---|---|---|---|
| S5a DRAM product | 111 × 98, six rows all zero | 111 × 98, six rows populated | crashed hero-only re-run overwrote the distillate |
| S5b DRAM metabolism summary | 1,319,283 B (crashed run) | 1,322,115 B (intact run) | same |
| S5c DRAM genome stats | 2,547 B (crashed run) | 2,547 B (intact run) | same; content was already identical |
| S11 biosafety | 111 × 5 | 111 × 6 | PlasmidFinder column never written |
| S21a abundance per MAG | 111 × 10, 47 rows mislabelled | 111 × 10, labels corrected | swine/sheep swapped in the prefix map |
| S21b abundance per source | 4 × 5, two rows swapped | 4 × 5, recomputed | consequence of S21a |
| S6c dbCAN families | 105 × 181, no MICP-complete MAGs | 111 × 192 | table came from the older DRAM-only run |
| Supplementary_tables.zip | 52 entries | 52 entries, rebuilt | — |

---

## 1. Table S5a / S5b / S5c — DRAM distillate

**Source used:** `/data/pangenome_work/dram_distilled/{product.tsv, metabolism_summary.xlsx,
genome_stats.tsv}` — the intact 2026-02-09 full-panel distillation
(`dram_distilled/distill.log` ends "Completed distillation").

**What was wrong.** All three shipped files came from
`/data/pangenome_work/dram_output/distillate/`. That directory is the product of the
2026-04-12 hero-only re-run: `dram_hero_run.log` records the six MICP-complete MAGs being
re-annotated and then `[distill fail]`, and `dram_output/distill.log` holds the traceback —
`FileNotFoundError: '/data/pangenome_work/dram_output/*/annotations.tsv'`. The annotate step
had already replaced the output directory, so the distillate that survived carries the six
MICP-complete MAGs as all-zero rows.

**Assertions that passed**

* the two distillates share row order, column list and 111 genomes;
* the shipped file's six MICP-complete rows sum to exactly 0 across all 32 numeric columns;
* the intact file's six rows are all non-zero;
* the two files differ in **exactly** those six rows and nowhere else.

**Restored module-completeness sums** (sum over the 32 numeric module columns):
S13 7.73, S16 10.65, S23 11.42, C22 3.21, M1 14.53, S26 14.48.

**Note.** S5c's content was already identical to the intact run when parsed; it was replaced
anyway so that all three S5 files carry one provenance. S5b is a genuinely different workbook.

**Consequence for the manuscript.** The shipped Figure 4a labelled *Sphingobacterium* "n = 2"
because four of its six MAGs were zero rows. Any statement drawn from the shipped S5a about
the MICP-complete MAGs is unreliable and should be re-read off the repaired table.

---

## 2. Table S11 — biosafety counts

**Source used:** `research/additional/A1_biosafety/combined/{card,vfdb,resfinder,plasmidfinder}_all.tsv`,
re-aggregated with the recipe of `aggregate_biosafety.py` (group by the MAG parsed out of
`#FILE`, reindex over the 111 Bakta directories, fill absent MAGs with 0).

**What was wrong.** `aggregate_biosafety.py` iterates over all four databases but skips any
whose combined file is absent. The aggregate was written at 04:48 and
`combined/plasmidfinder_all.tsv` only at 04:58, so the fourth column was silently dropped
from both the analysis file and the shipped table.

**Assertions that passed**

* 111 Bakta directories found;
* re-aggregated CARD, VFDB and ResFinder counts equal the shipped values for all 111 MAGs;
* the `group` column is unchanged.

**Raw totals:** CARD 289 hits / 74 MAGs, VFDB 1,303 / 75, ResFinder 35 / 24,
**PlasmidFinder 1 hit in the whole panel** — `IncQ2_1` (FJ696404, 98.0 % coverage,
82.8 % identity) in **S21**, a non-MICP-complete MAG. **Zero PlasmidFinder hits in all six
MICP-complete MAGs**, which is what the §3.9 text and the Fig. S8 legend already assert;
the table now carries the evidence.

Column order after the fix: `MAG, card, vfdb, resfinder, plasmidfinder, group`.

The upstream `A1_biosafety/biosafety_counts_per_MAG.csv` was refreshed identically so a
future re-export cannot reproduce the omission.

---

## 3. Table S21a / S21b — waste-source labels

**Source used:** the MAG identifier prefix, with the sample-code key stated in
`02_Figure_legends.md` (**C** cattle, **M** swine, **S** sheep, **V** poultry), cross-checked
against the independent `Source` column of Table S9a.

**What was wrong.** `scripts/additional/C6_abundance_proxy/run_abundance_proxy.py` line 26:

```python
return {"C":"cattle","S":"swine","M":"sheep","V":"poultry"}.get(p, "unknown")
```

`S` and `M` are transposed. 47 of 111 MAGs carried the wrong source label. **This line must
be corrected before the C6 analysis is re-run** — the fix here is applied to the outputs.

**Assertions that passed**

* every MAG in S21a is present in Table S9a;
* the corrected labels agree with Table S9a for all 111 MAGs;
* no column other than `source` changed — the coverage numbers are untouched.

**Group sizes:** shipped sheep 15 / swine 32 → corrected **swine 15 / sheep 32**
(cattle 22 and poultry 42 unchanged).

**Per-source length-weighted contig coverage, shipped → corrected**

| source | n | mean (×) | median (×) |
|---|---|---|---|
| cattle | 22 → 22 | 5.85 → 5.85 | 4.67 → 4.67 |
| poultry | 42 → 42 | 28.63 → 28.63 | 9.75 → 9.75 |
| sheep | 15 → 32 | 57.43 → **20.32** | 33.27 → 11.99 |
| swine | 32 → 15 | 20.32 → **57.43** | 11.99 → 33.27 |

**Consequence for the manuscript.** The Fig. S19 legend sentence "sheep-rumen MAGs are the
most deeply covered (mean 57.4×)" is about the **swine** group (n = 15). The corrected sheep
mean is 20.3× (n = 32). MICP-complete MAG sources after the fix: C22 cattle, M1 swine,
S13/S16/S23/S26 sheep — which is what the sample codes say.

The upstream `C6_abundance_proxy/abundance_proxy_{per_MAG,per_source}.csv` were refreshed
identically.

---

## 4. Table S6c — dbCAN family counts

**Sources used:** `research/revision/dbcan_direct/<MAG>.tbl` (HMMER-3 hmmsearch, best HMM per
protein by E-value) for the six MICP-complete MAGs, and the `cazy_best_hit` / `cazy_ids`
columns of `/data/pangenome_work/dram_output/all_annotations.tsv` for the other 105 — the two
routes used by `04b_dbcan_final.py`.

**What was wrong.** S6a, S6b and S6d are outputs of the hero-inclusive `04b_dbcan_final.py`,
but S6c is byte-identical to `research/revision/dbCAN_family_counts.csv`, the **older**
DRAM-only family table covering 105 MAGs. `04b_dbcan_final.py` computes `rest_fam` but never
writes a family table, so the family level was never updated. With no MICP-complete rows in
the table, the shipped Figure S4c averaged over an empty group and rendered blank.

Family names also disagreed between the two routes: the hmmsearch side kept the full HMM name
(`GT2_Glycos_transf_2`) while the DRAM side yields the bare family (`GT2`), so the same family
appeared as two different columns. Both sides are now reduced to the bare family
(`^(GH|GT|PL|CE|AA|CBM)\d+`).

**Assertion that passed:** for each of the six CAZy classes, the per-MAG sum of that class's
family columns equals the class count in **Table S6a**, for all 111 MAGs. This pins the new
family table to the shipped class table on both annotation routes.

**Shape:** 105 × 181 → **111 × 192**. Example of what the normalisation recovers — GT2:
MICP-complete mean 9.83 copies, rest mean 6.73; the shipped table had no MICP-complete rows
at all.

A copy is kept beside the other hero-inclusive artefacts as
`research/revision/dbCAN_final_family_counts.csv`.

**Method asymmetry (inherited, not a defect):** the six MICP-complete MAGs are annotated by
direct hmmsearch and the other 105 by DRAM's `cazy_best_hit`. This is the published design of
`04b_dbcan_final.py` and applies to S6a/S6b/S6d as well; it belongs in the legend.

---

## 5. Audit of the remaining tables

Every table was checked for the three defect classes seen above: rows missing for the six
MICP-complete MAGs, a mis-mapped categorical column, and a column present in the raw output
but absent from the table. Each per-MAG table's shape and column list were compared with the
raw analysis file named in `02_Figure_legends.md` §Supplementary tables.

**Column/shape comparison against the named source: all 17 comparable tables match exactly**
(S11 matched its source too — the defect was upstream, in the source itself, and both were
repaired).

**No further defects found.** Tables whose MICP-complete count is below 6/6 were each checked
and are legitimate:

| Table | MICP-complete rows | Why this is correct |
|---|---|---|
| S16 gRodon | 4 / 6 | S13 and S16 fall below the ≥ 10 ribosomal-anchor filter; stated in the legend |
| S22 ESMFold TM | 0 / 6 by name | keys are `S13_UreC` etc., all six present |
| S10a ext. ANI matrix | 0 / 6 by name | keys are `OUR_S13` etc.; all six study *Sphingobacterium* present |
| S10b, S1b | 4 / 6 | *Sphingobacterium*-only tables by design; M1 and S26 are *Pseudomonas*_E |
| S13a, S13b | 0 / 6 | rows are 146 NCBI *Pseudomonas*_E reference accessions, not study MAGs |
| S14b | 0 / 6 | rows are 7,599 MGnify species-cluster representatives |
| S1a | 45 rows | the *ure*-bearing subset of the panel; all six present |

**Left unfixed / for the user**

1. `scripts/additional/C6_abundance_proxy/run_abundance_proxy.py` line 26 still has the
   swapped prefix map. Outputs are corrected; the code is not (outside this group's scope).
2. `pangenome_work/summarize_pangenome_micp.py` and its output
   `MICP_Pangenome_Final_Summary.csv` remain unrepaired. That file covers 100 of 111 MAGs and
   cannot be reproduced from the Panaroo matrix it names; it is **not** a supplementary table,
   so it was out of scope here, but it is the blocker recorded at the top of the figure
   journal and it drives the "8/8 module restricted to six MAGs" claim.
3. Table S1c reports C22 with 2,193 contigs while the repaired S5c reports 463 scaffolds for
   the same MAG. These are different metrics from different tools (assembly contigs vs DRAM
   scaffolds) and are not necessarily in conflict, but the two numbers appear in the same
   submission and may draw a reviewer's eye.

## 6. Effect on the rebuilt figures

None expected. The 2026-09-04 figure rebuild had already bypassed all four defects by reading
the true sources: Fig 4 from the intact distillate, Fig S8 from the raw abricate output,
Fig S19 by deriving the source from the MAG identifier, and Fig S4c by re-parsing
`dbcan_direct/`. The repaired tables now agree with what the figures already show.
The figures were not re-run by this group.

---

## Second pass, 2026-09-04 (group "consolidate")

Two further tables were re-exported after the *ureB* sRNA/CDS defect was traced during the
figure rebuild. Both fixes live in the same `reexport_tables.py` (functions `fix_s1a`,
`fix_s15b`), which remains idempotent and supports `--check`.

### Table S1a — per-gene copy counts

| | shipped | re-exported |
|---|---|---|
| MAGs covered | 45 | **111** |
| feature types counted | every Bakta feature | **CDS only** |
| MAGs carrying all eight genes | 28 | **27** |

Two defects, both confirmed against the Bakta annotations:

1. The keyword scan matched the gene-name column over **every feature type**, so the Bakta
   `5_ureB_sRNA` non-coding RNA (RFAM RF02514, product "5' ureB small RNA") was counted as a
   copy of the urease β-subunit. The feature occurs 58 times panel-wide. For **S26** (0 CDS +
   2 sRNA) and **S11** (0 CDS + 1 sRNA) it was the *only* *ureB* feature, so the table
   recorded the gene as present where no protein-coding β-subunit is annotated. S26 is one of
   the six MICP-complete MAGs, and this is the origin of its long-standing disagreement with
   Table S1c (`ureABC_all_present = False`) and the pangenome summary (1/8).
2. The table listed only 45 of the 111 MAGs. Thirteen of the omitted MAGs carry a real *ure*
   CDS (V9 carries seven of the eight genes), so reading absence as "no MICP gene" understated
   every prevalence in the non-MICP-complete group.

Assertions: the all-feature rule reproduces the shipped presence/absence for all 45 covered
MAGs; excluding non-CDS features changes only *ureB*; exactly S11 and S26 lose a gene.

Known residual: *ureA*/*ureC* copy counts differ from the shipped table for four MAGs (C13,
S23, S26, V3), where a feature carries a gene name and a product naming a different subunit,
so α and γ swap between the two routes. Presence is unaffected.

### Table S15b — stoichiometry

The `ureB` column was byte-identical to the all-feature scan, so the CDS-only value is well
defined. It was re-exported on that basis, and the two derived flags that are a pure function
of it were recomputed (both definitions verified against the shipped table first:
`urease_core_complete = ureA & ureB & ureC`, `MICP_stoich_complete = urease_core & any_CA`).

- 52 MAGs have a lower *ureB* copy count once the sRNA is excluded.
- Derived flags change for **S11** and **S26** only: `urease_core_complete` 1 → 0 and
  `MICP_stoich_complete` 1 → 0 for both.
- `ureA` and `ureG` use a broader keyword rule than the shared scan (75 and 17 MAGs differ
  respectively) and were **left as shipped**; only the documented *ureB* defect was repaired.

### Knock-on effects handled

- `new_figure/_micp_presence.py` asserted that Table S1a had 45 rows and matched the
  all-feature scan. That assertion now fails against the re-export, which would have broken
  every rebuild of Fig 1 and Fig 2. The check was made version-aware: it pins the CDS-only
  scan against the re-exported table when the table has 111 rows, and retains the original
  45-MAG check otherwise.
- Fig 1, Fig 2 and Fig S14 rebuilt against the corrected tables; `verify_outputs.py` 29/29 PASS.
- `Supplementary_tables.zip` rebuilt: 52 entries, 0.92 MB.

Originals of both tables are in `originals/`.
