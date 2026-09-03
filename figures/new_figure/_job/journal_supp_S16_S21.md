# Build journal — group supp_S16_S21 (Suppl Figs S16–S21)

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
