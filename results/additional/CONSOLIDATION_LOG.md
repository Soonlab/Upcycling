# Consolidation pass — 2026-09-04 (group "consolidate")

Runs after five separate editing passes (figure rebuild, table re-export, manuscript
corrections, legend rewrite, *ureB* fix). Backups of both prose files as they stood at the
start of this pass are in `consolidate_backup/`.

---

## Job 1 — corrected C5 results propagated into the prose

Source: `research/additional/C5_panMICP_env_v2/C5_RERUN_LOG.md` and `reference_manifest.tsv`.

The earlier passes had described the C5 panel as *withdrawn* — correct at the time, but the
panel has since been rebuilt with all 20 references organism-verified, so every statement
about it was out of date.

| location | before | after |
|---|---|---|
| Methods §2.23 | panel "failed quality control", results "reported only in reduced form", two verified genomes, re-acquisition pending | panel composition listed; the accession-list fault stated (2 of 20 accessions named the intended taxon, 5 do not exist, 13 unrelated); rebuild procedure and organism verification described; the two renamed taxa and the one unfillable slot recorded |
| Results §3.11 | "contributes only one usable comparison"; a sentence deferring the below-95 % claim to §3.8 | 8 of 120 pairs align, all eight values given; S23/C22 above the boundary explained as the expected positive control; **new headline** — no MICP-complete MAG shares measurable identity with any canonical ureolytic reference, *S. pasteurii* DSM 33 included |
| Discussion §4.4 | "one supporting analysis was withdrawn at quality control" | "had to be rebuilt at quality control", with the provenance lesson stated as a recommendation |
| Fig S18 legend | one-panel reduced figure; reference-integrity paragraph | two panels described (6 × 20 matrix, ranked reported pairs); provenance paragraph rewritten; the unsupported "below 95 % against every reference" claim removed; title changed to "organism-verified" |

The M1 ↔ *P_E* sp. 025837155 = 97.98 % value was already correctly attributed to the A3
screen in both files; left as is.

---

## Job 2 — Table S1a re-exported, Table S15b corrected

Added `fix_s1a()` and `fix_s15b()` to `tables/reexport_tables.py` (same idempotent,
`--check`-capable structure; originals in `tables/originals/`). Detail and assertions are in
`tables/TABLE_REEXPORT_LOG.md`, second-pass section. Summary:

- **Table S1a** 45 MAGs → **111**, all-feature counts → **CDS only**, all-eight count 28 → **27**.
  S26 and S11 are the only MAGs that lose a gene; both had the non-coding `5_ureB_sRNA`
  (Rfam RF02514) as their sole *ureB* feature.
- **Table S15b** *ureB* re-exported CDS-only. 52 MAGs have a lower copy count; derived flags
  change for **S11 and S26 only** (`urease_core_complete` and `MICP_stoich_complete` 1 → 0).
  `ureA`/`ureG` use a broader keyword rule than the shared scan and were left as shipped.
- `Supplementary_tables.zip` rebuilt (52 entries, 0.92 MB).

**Knock-on defect caught and fixed.** `new_figure/_micp_presence.py` asserted that Table S1a
had 45 rows and matched the all-feature scan. The re-export made that assertion fail, which
would have broken every future rebuild of Fig 1 and Fig 2. The check is now version-aware:
against a 111-row table it pins the CDS-only scan exactly, and it retains the original 45-MAG
check otherwise. Fig 1, Fig 2 and Fig S14 rebuilt; `verify_outputs.py` 29/29 PASS.

---

## Job 3 — consistency sweep

Fixed:

1. **Five supplementary tables were never cited.** S3, S4, S6, S7 and S9 existed and had
   legends but no callout in the body; two of them were instead referenced by raw filename
   (`HGT_ureCah_cluster.csv`, `dbCAN_final_MICP-complete_vs_rest_class.csv`). Callouts added at
   §3.3 (S3), §3.5 (S6), §3.7 (S7), §3.7a (S9) and §3.8 (S4), replacing the bare filenames.
   All of S1–S23 are now cited.
2. **S26 stoichiometry.** With the corrected Table S15b, S26 fails the urease-core flag as
   well as the calcium-pathway flag. §3.10 said the sole exception was calcium; it now records
   both, and the Fig S14 legend no longer claims "all six MAGs carry a complete urease core".
3. **Table state described as unfixed.** Nine notes across both files still said tables were
   corrupt, omitted a column, swapped labels, or "should be re-exported". All were rewritten
   to the corrected state: Methods §2.2a, Results §3.2 (twice), and the Fig 4, Fig S8, Fig S14,
   Fig S19 legends plus the Table S1, S5, S6, S11 and S21 entries.
4. **Abstract word count.** The editorial note claimed 386 words; a whitespace count gives 391
   (and 249, not 245, for the mSystems version). Both corrected, with a caution that the
   400-word margin is thin.

Checks run, all passing:

| check | result |
|---|---|
| one legend per figure, Fig 1–8 and S1–S21 | PASS, in numerical order |
| every Table S1–S23 cited in the manuscript | PASS |
| every cited table present in `Supplementary_tables/` | PASS |
| every figure filename in the legends exists on disk | PASS, 29 stems |
| no stale phrase from the superseded state remains (13 probes) | PASS |
| Fig 1 legend defers to Methods §2.2a, not §2.7 | PASS |
| Fig 1 legend and Methods agree that 27 MAGs carry all eight genes | PASS |
| `new_figure/_job/verify_outputs.py` | 29/29 PASS |

Citation cross-check: 63 reference entries, 43 distinct surnames cited author–year. Five
apparent orphans (Mukherjee, Skolnick, Standley, Tiedje, Yu) are all second authors of
"X and Y, year" citations whose reference entries are first-author-only — false positives, no
action. Twenty reference entries are never cited by author–year; the manuscript mixes
author–year prose citations with a numbered Vancouver list, so these are reachable only by
number. Pre-existing, unchanged by this pass, and flagged below.

Word counts: **abstract 392** (MB limit 400), **alternative abstract 249** (mSystems limit
250), **body Introduction→Conclusions 9,848**.

---

## Left unresolved

1. **Body length.** 9,848 words against the ~8,000-word Microbial Biotechnology guidance.
   mSystems does not cap body text. This is a journal-choice question, not a factual one.
2. **Citation style is mixed.** Prose uses author–year; the reference list is numbered
   Vancouver, and 20 entries are never reached by an author–year callout. Reconciling this
   needs the target journal's style, so it was not touched.
3. **Table S1c versus Table S5c contig counts.** S1c reports 2,193 contigs for C22, the
   repaired S5c reports 463 scaffolds. Different tools and different metrics, both in the
   package; a note in whichever table the authors prefer would settle it.
4. **`MICP_Pangenome_Final_Summary.csv` is still unrepaired** (covers 100 of 111 MAGs). No
   figure or table now depends on it and Methods §2.2a says so explicitly, but the file remains
   in the working tree.
5. **`Sphingobacterium` sp. 21** has no genome under that strain designation, so the C5 panel
   carries 20 verified references with that slot filled by other *Sphingobacterium* type
   strains rather than by the originally intended entry.
