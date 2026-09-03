# C5 pan-MICP environment ANI — re-run with accession-verified references (2026-09-04)

Supersedes `research/additional/C5_panMICP_env/`. Outputs: this directory, plus the re-exported
`SUBMISSION/Supplementary_tables/Table_S20_skani_hero_vs_refs.tsv` and `new_figure/figures/Fig_S18.*`.

## 1. What was wrong

The figure rebuild flagged the C5 reference panel as containing organisms other than the ones
its filenames claimed. The re-run establishes that **the fault was in the accession list, not in
the download**. Every accession in `curated_refs.tsv` was queried against the NCBI datasets API:

| intended organism | accession used in v1 | what that accession actually is |
|---|---|---|
| Sporosarcina pasteurii | GCF_000189475.1 | does not exist |
| Sporosarcina psychrophila | GCF_000293575.1 | does not exist |
| Sporosarcina ureae | GCF_000299475.1 | *Escherichia coli* O104:H4 str. 2009EL-2071 |
| Bacillus subtilis 168 | GCF_000009045.1 | **correct** |
| Bacillus megaterium | GCF_002209305.2 | does not exist |
| Bacillus licheniformis | GCF_001548095.1 | *Geminocystis* sp. NIES-3708 |
| Lysinibacillus sphaericus | GCF_001399775.1 | *Thermus aquaticus* Y51MC23 |
| Lysinibacillus fusiformis | GCF_000691605.1 | *Bdellovibrio bacteriovorus* |
| Halomonas pacifica | GCF_900119685.1 | *Escherichia coli* |
| Halomonas elongata | GCF_000196175.1 | *Bdellovibrio bacteriovorus* HD100 |
| Idiomarina loihiensis | GCF_002288455.1 | does not exist |
| Pseudoalteromonas haloplanktis | GCF_000613485.1 | *Hoylesella loescheii* (suppressed) |
| Helicobacter pylori 26695 | GCF_001077795.1 | *Riemerella anatipestifer* Yb2 |
| Klebsiella aerogenes | GCF_000196215.1 | *Borreliella bavariensis* PBi |
| Proteus mirabilis | GCF_000648555.1 | does not exist |
| Yersinia enterocolitica | GCF_000009305.1 | *Helicobacter acinonychis* str. Sheeba |
| Sphingobacterium sp. 21 | GCF_000284615.1 | *Marinobacter nauticus* ATCC 49840 |
| Sphingobacterium multivorum | GCF_000620325.1 | *Salinispora tropica* CNY012 |
| Sphingobacterium paramultivorum | GCF_001719105.1 | *Enterobacter* sp. HK169 |
| Pseudomonas helleri | GCF_001043025.1 | **correct** (DSM 29165) |

Two of twenty accessions named the intended taxon. The C5 reference panel as published was
therefore fictitious: no *Sporosarcina*, no *Sphingobacterium* and no ureolytic clinical reference
was ever in the comparison, and the ANI values between mislabelled entries (99.02 % "Halomonas
elongata" vs "Lysinibacillus fusiformis", 96.45 % "Halomonas pacifica" vs "Sporosarcina ureae")
are two strains of one species compared under two wrong names.

The six MICP-complete MAG files were checked and are sound: `HERO_*.fna` are byte-identical
(md5) to `MAGs_FASTA_files/{MAG}.fasta.gz`, with matching contig counts.

## 2. How the panel was rebuilt

`resolve_accessions.py` → `download_refs.py` → `run_skani_v2.sh` → `export_table_S20.py`.

1. Each intended **taxon name** resolved against the NCBI datasets API, preferring the RefSeq
   reference genome, else the best-assembled current RefSeq genome (`resolved_accessions.tsv`).
2. Downloaded **by accession**, file named by accession. The organism NCBI reports for the
   accession is compared against the intended taxon before the file is accepted
   (`reference_manifest.tsv`, with strain designation).
3. skani with the original parameters: `triangle --full-matrix`, and
   `dist -q HERO_* -r <all> --min-af 30`.

**20 of 20 references are organism-verified.** Notes on three entries:

- *Sporosarcina pasteurii* — recovered as **GCF_041295575.1** (DSM 33, complete genome). The
  canonical MICP organism was absent from the published panel entirely.
- *Halomonas pacifica* → **Bisbaumannia pacifica** (GCF_007989625.1) and *Bacillus megaterium* →
  **Priestia megaterium** (GCF_006094495.1): current NCBI names for the same tax_id, confirmed
  against the taxonomy endpoint, not assumed from the name.
- *Sphingobacterium* sp. 21 — no genome exists under that strain designation. The slot was not
  filled with an unrelated genome; the panel instead carries the two *Sphingobacterium* type
  strains the manuscript actually compares against (*S. multivorum* GCF_016894225.1,
  *S. paramultivorum* GCF_014109745.1), plus *P. helleri* NF1a alongside DSM 29165.

Accessions replaced are recorded per row in Table S20 (`Superseded_accession_in_v1`).

## 3. Results, old versus new

skani reports a pair only where the genomes align (roughly ANI > 80 % at `--min-af 30`).
Eight of 120 MAG × reference pairs are reported.

| MAG | reference (verified) | new ANI | v1 |
|---|---|---|---|
| S23 | *Sphingobacterium paramultivorum* BIGb0170 | **98.96** | absent from panel |
| S26 | *Pseudomonas helleri* NF1a | **97.88** | not in panel |
| S26 | *Pseudomonas helleri* DSM 29165 | **97.54** | 97.54 (reproduced exactly) |
| C22 | *Sphingobacterium paramultivorum* | **96.71** | absent from panel |
| C22 | *Sphingobacterium multivorum* | **94.18** | absent from panel |
| S16 | *Sphingobacterium multivorum* | **93.85** | absent from panel |
| S23 | *Sphingobacterium multivorum* | **91.99** | absent from panel |
| S16 | *Sphingobacterium paramultivorum* | **91.70** | absent from panel |

S13 and M1 align with no reference in the panel.

Two independent cross-checks of the rebuilt panel against the ext_sphingo screen behind Fig. 5c,
which used a separate, correctly-named 63-genome reference set: **S16 vs *S. multivorum* = 93.85 %**
and **S23 vs *S. paramultivorum* = 98.96 %** are identical to the values reported there. The
rebuilt C5 panel and the Fig. 5c panel now agree to two decimals.

## 4. Effect on the manuscript's C5 claims

**Claim: S26 ↔ *P. helleri* DSM 29165 = 97.54 % — SURVIVES, unchanged.** The value is reproduced
exactly against the verified DSM 29165 genome, and the current *P. helleri* reference NF1a gives
97.88 %, the same species-level assignment.

**Claim: the four *Sphingobacterium* MICP-complete MAGs remain below 95 % against every reference
examined — DOES NOT SURVIVE as written, and should be removed from the S18 legend.** Against
genuine *Sphingobacterium* references, S23 reaches 98.96 % and C22 reaches 96.71 %, both above the
species boundary. This is not a new problem for the paper: S23 and C22 have GTDB species
assignments (*S. paramultivorum*, *S. siyangense*) and Fig. 5c already reports them at 98.96 % and
99.16 % as the pipeline's positive control. The novelty claim rests on **S13 and S16**, and it
holds — S16 reaches only 93.85 %, and S13 aligns with nothing in the panel. The v1 legend
generalised a statement about S13/S16 to all four MAGs and was never supported by C5 data.

**Note on a value quoted in the S18 legend:** M1 ↔ *P_E* sp. 025837155 = 97.98 % is an **A3**
result (146-genome *Pseudomonas*_E screen), not a C5 one. M1 aligns with no genome in the C5
panel. The legend should attribute it to A3 or drop it.

**New, and worth stating:** none of the six MICP-complete MAGs shares measurable whole-genome
identity with any canonical ureolytic reference — *Sporosarcina pasteurii* included, now that it
is actually in the panel. Only same-genus references align. This is direct support for the
paper's framing of these MAGs as non-canonical MICP chassis, and it is a statement C5 could not
make before because the canonical references were never present.

## 5. Files

```
resolve_accessions.py        taxon -> current RefSeq accession
resolved_accessions.tsv      resolution result, 19 of 20 taxa (Sphingobacterium sp. 21 unresolved)
download_refs.py             download by accession + organism verification
reference_manifest.tsv       accession, intended taxon, reported organism, strain, size, status
run_skani_v2.sh              skani triangle + dist, original parameters
skani_panMICP.matrix(.af)    26 x 26 all-vs-all (20 references + 6 MAGs)
skani_hero_vs_refs.tsv       reported pairs
export_table_S20.py          writes SUBMISSION/.../Table_S20_skani_hero_vs_refs.tsv
refs/                        GCF_*.fna (accession-named) + HERO_*.fna
logs/                        accession_status.tsv, per-accession API reports, download zips, skani log
Table_S20_OLD_backup.tsv     the superseded table
build_supS18_v1_backup.py    the superseded figure script
```
