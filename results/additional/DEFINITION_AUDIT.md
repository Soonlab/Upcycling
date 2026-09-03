# Audit — what defines the six "MICP-complete" MAGs?

2026-09-04. Prompted by the 29-figure rebuild journal
(`/data/data/Upcycling/new_figure/_job/JOURNAL.md`), which found that the manuscript's
statement that mapping *ureA–G* + *cah* onto the phylogeny "identified a MICP-complete set of
six MAGs" is not supported by any source in the repository.

## Question

§3.2 and the Fig. 1 legend present the six MAGs (S13, S16, S23, C22, M1, S26) as the outcome
of a gene-presence criterion. Is there a rule over the available data that returns exactly
those six?

## Method

Twelve candidate rules were evaluated over every per-MAG source in the repository: a direct
keyword scan of all 111 Bakta GFF3 annotations (`ure` subunit and carbonic-anhydrase product
strings, recording the contig of each hit), Table S1a (per-gene copy counts), Table S15b
(stoichiometry flags), the shipped `pangenome_work/MICP_Pangenome_Final_Summary.csv`, and each
of those restricted to the two lineages that the six MAGs occupy. Scan script:
`_revision_260904/scan_micp_architecture.py`.

## Result — no rule reproduces the six

| # | Candidate rule | MAGs selected | Reproduces the six? | Extra | Missing |
|---|---|---|---|---|---|
| A | Bakta GFF3: all 8 genes present (score 8/8) | 32 | no | 27 MAGs | S26 |
| B | Bakta GFF3: *ureA–G* present + *cah* present | 32 | no | 27 MAGs | S26 |
| C | Bakta GFF3: *ureA–G* on one shared contig | 26 | no | 23 MAGs | C22, S23, S26 |
| D | Bakta GFF3: *ureA–G* on one contig + *cah* | 26 | no | 23 MAGs | C22, S23, S26 |
| E | Bakta GFF3: *ureC* and *cah* on the same contig | 1 | no | M13 | all six |
| F | Table S1a: all 8 genes with copy number > 0 | 28 | no | 22 MAGs | — |
| G | Table S15b: urease core + accessory + any CA | 33 | no | 27 MAGs | — |
| H | Rule G + a Ca-handling gene | 7 | no | M9, V3 | S26 |
| I | Shipped pangenome summary: `Total_Score == 8` | 26 | no | 22 MAGs | S13, S26 |
| J | Rule A restricted to the two lineages | 21 | no | 16 MAGs | S26 |
| K | Rule F restricted to the two lineages | 17 | no | 11 MAGs | — |
| L | Rule C restricted to the two lineages | 16 | no | 13 MAGs | C22, S23, S26 |

The closest any rule comes is H (7 MAGs: the six minus S26, plus M9 and V3).

## The S26 problem

S26 fails every gene-content rule because **its Bakta annotation contains no urease β-subunit**.
The complete set of urease and carbonic-anhydrase features in `S26.gff3`:

| contig | product | n |
|---|---|---|
| contig_306, contig_85 | Urease subunit alpha (*ureC*) | 2 |
| contig_306 | Urease subunit gamma (*ureA*) | 1 |
| — | **Urease subunit beta (*ureB*)** | **0** |
| contig_137, contig_140 | Urease accessory protein UreD | 2 |
| contig_140, contig_25 | Urease accessory protein UreE | 2 |
| contig_140, contig_25 | Urease accessory protein UreF | 2 |
| contig_140, contig_25 | Urease accessory protein UreG | 2 |
| contig_189 | Urease accessory protein UreH-like | 1 |
| contig_19, contig_66, contig_27 | Carbonic anhydrase (α ×2, γ ×1) | 3 |

Three independent records agree that S26 is not urease-gene-complete:

* the direct Bakta scan above (no *ureB*);
* `Table_S1c_hero_cluster_audit.csv`, which records `ureABC_all_present = False` and
  `ureABCDEFG_all_present = False` for S26 alone among the six;
* the shipped `MICP_Pangenome_Final_Summary.csv`, which scores S26 **1/8**.

Two tables nevertheless report `ureB = 2` for S26 — Table S1a and Table S15b. Neither is
traceable to the Bakta annotation, and Table S15b additionally reports `ureA = 7` for S26 where
the annotation carries one. These two tables are aggregated by a looser keyword route than the
product-string match used here; the discrepancy is recorded rather than resolved, because both
routes are shipped and they disagree.

S26 also has `Ca_transporter = 0` and `Ca_ATPase = 0` in Table S15b, so its stored `Ca_pathway`
flag is 0. It carries five CO₃-transporter copies, which do not enter that call.

## Provenance of the six

The list is a hard-coded literal in every analysis script that uses it — at least fourteen
occurrences, all identical:

```
research/revision/01_cluster_and_quality.py:18       HERO=["C22","M1","S13","S16","S23","S26"]
research/revision/05_PCoA_by_source.py:24            HERO = {"C22","M1","S13","S16","S23","S26"}
research/additional/A1_biosafety/aggregate_biosafety.py:7   HERO = {"S13","S16","S23","C22","M1","S26"}
research/additional/A5_alkaliphile/…:22              HERO = {"S13","S16","S23","C22","M1","S26"}
scripts/02_trait_module_scans.py:38                  HERO = {"S13","S16","S23","S26","C22","M1"}
…
```

No script derives it. The seven per-MAG supplementary tables that carry a group flag
(S11, S15a, S17a, S18a, S19a, S21a, S23a) all encode the same six, because all seven were
written from the same literal.

## Conclusion adopted in the manuscript

The six-MAG set is a **curated shortlist**, not the output of a threshold. It was assembled
from the *Sphingobacterium* and *Pseudomonas*_E MAGs that combine a near-complete urease
operon, a carbonic anhydrase, and cluster contiguity good enough to inspect (Table S1c) — but
no single one of those conditions, nor any conjunction of them tested here, returns exactly
these six, and S26 is retained despite lacking a detectable *ureB*.

Accordingly the manuscript now:

1. states the criterion as applied, and states plainly that it is a curated shortlist rather
   than a reproducible threshold (new Methods §2.2a);
2. drops the claim that gene mapping "identified" the six or that the 8/8 module is
   "restricted" to them — 28–32 MAGs carry all eight genes depending on the source (§3.2,
   Fig. 1 sentence, Abstract, Highlights, Conclusions);
3. records the S26 *ureB* gap and the source disagreement in Results §3.2 and in
   Limitations §4.4;
4. keeps the six as the analysis set, because every downstream analysis was run on them, and
   frames them as chassis candidates rather than as an exhaustive detection result.
