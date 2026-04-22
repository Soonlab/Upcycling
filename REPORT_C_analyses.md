# C1-C6 reviewer-defense additions (extended panel)

**Date**: 2026-04-20 (initial) / 2026-04-21 (C1 re-run + C3 v3 + C4 HF ESMFold) / **2026-04-22 (integrated into main manuscript)**

All in-silico, no wet-lab.

**Integration status (2026-04-22)**: all C1–C6 analyses — alongside the earlier A1–A7+B reviewer-defense layer — have been **integrated into the initial submission body** as Methods §2.11–2.24, Results §3.9–3.12, Discussion §4.6–4.10, supplementary Figures S8–S21, and Tables S11–S23. The revision-card document (`proposed_supplementary_additions.md`) is retained as a historical planning artefact and is superseded by the submission package at `/data/data/Upcycling/SUBMISSION/` (manuscript itself is gitignored until publication).


## C1 antiSMASH BGC scan

`results/additional/C1_antismash/antismash_per_MAG.csv` + `_hero_vs_rest.csv`

First run (2026-04-20) failed because the MITE database was missing from
`antismash_env` (ValueError: No matching database in location .../databases/mite).
On 2026-04-21 we ran `download-antismash-databases` (9.4 GB installed) and
re-launched the scan over all 111 MAGs in `antismash_env` (8 jobs × 4 threads).

### Total BGC load
Hero mean 6.2 regions/MAG vs rest 6.9 — heroes are neither BGC-rich nor
BGC-poor (MWU not significant).

### Top class-level enrichments (Mann-Whitney U)

| metric | hero_mean | rest_mean | MWU p |
|--------|-----------|-----------|-------|
| BGC_T3PKS          | 0.67 | 0.029 | **5.3e-10** |
| BGC_RRE-containing | 0.83 | 0.105 | **1.9e-05** |
| BGC_NRP-metallophore | 0.00 | 0.46  | 0.082 (depleted in heroes) |
| BGC_resorcinol     | 0.00 | 0.24  | 0.18 |
| BGC_NAGGN          | 0.33 | 0.14  | 0.21 |

### Per-hero BGC profile

```
MAG  n_regions  classes
C22  3          RRE-containing;T3PKS;terpene
M1   7          NAGGN;RiPP-like;betalactone;redox-cofactor;terpene-precursor
S13  5          RRE-containing;T3PKS;arylpolyene;terpene
S16  5          RRE-containing;T3PKS;arylpolyene;terpene;terpene-precursor
S23  5          RRE-containing;T3PKS;arylpolyene;betalactone;terpene-precursor
S26  12         NAGGN;NRPS;RiPP-like;arylpolyene;betalactone;hydrogen-cyanide;...
```

**Interpretation**: The four Sphingobacterium heroes (S13/S16/S23/C22) carry a
characteristic **T3PKS + RRE-containing + arylpolyene + terpene** secondary-
metabolism fingerprint — an apparently lineage-specific signature consistent
with published Sphingobacterium metabolomes. The two Pseudomonas_E heroes
(M1/S26) show a distinct, richer profile (NAGGN / NRPS / RiPP-like /
betalactone / hydrogen-cyanide), typical of environmental Pseudomonas.
The MICP-complete trait is therefore superimposed on a metabolically active
background, not a stripped-down chassis — relevant for downstream engineering
(metabolic load, precursor competition).

## C2 DefenseFinder + minced CRISPR

```
              metric  hero_mean  rest_mean  MWU_p
0  n_defense_systems        0.0        0.0    1.0
1    n_crispr_arrays        0.0        0.0    1.0
```

Interpretation: both groups are anti-phage-naive. This is a **null result** in the
between-group sense but a positive one for chassis engineering - no programmable
systems to displace.

## C3 dN/dS + codon usage

### Codon usage (111 MAGs, whole-proteome)

```
    metric  hero_mean  rest_mean     MWU_p
0  GC3_pct  51.413327  70.065052  0.043288
1      ENC  49.988579  39.004205  0.002830
```

Hero MAGs carry a low-GC codon background (Sphingobacterium dominates the hero
set) with reduced codon bias (higher ENC).

### codeml M0 single-ω (v3: 18-MAG subset + FastTree topology)

```
gene  n_seqs  aln_len  omega_M0  tree_dN  tree_dS  lnL       kappa
ureA  18      4215     0.08708   6.2701   72.0004  -25402.09 1.41
ureB  17      1749     0.05942   2.3414   39.4026  -6372.31  1.37
ureC  18      1875     0.02605   1.3080   50.2195  -16463.71 1.50
ureG  18      828      0.04082   2.3150   56.7136  -7609.11  1.36
```

v1/v2 runs failed because PAML's star-tree MAXNSONS limit is exceeded with
>16 leaves. v3 uses 6 heroes + 12 stratified non-hero reps and a FastTree
bifurcating topology. All four urease genes are under **strong purifying
selection** (ω << 1); ureC (catalytic α) is the most constrained.

### yn00 pairwise ω partitioned by group

```
gene  hero_hero_n  hh_median  rest_rest_n  rr_median  hero_rest_n  hr_median  MWU_hh_vs_rr_p
ureA  14           0.767      65           0.717      71           0.697      0.832
ureB  8            0.132      65           0.086      60           0.146      0.382
ureC  12           0.137      62           0.121      71           0.146      0.572
ureG  12           0.312      62           0.074      72           0.196      7.66e-08
```

**ureG** shows an ~4× elevation in hero-hero ω vs rest-rest ω (MWU p = 7.7e-8),
while ureA/B/C are indistinguishable across buckets. Interpretation: **relaxed
purifying selection** on the GTPase accessory within the two convergent
MICP-complete lineages, the catalytic subunits remain conserved.

## C4 ESMFold UreC vs PDB 4CEU

First run (2026-04-20) failed with `ModuleNotFoundError: openfold` because
fair-esm's ESMFold v1 imports openfold, which isn't pip-installable without
CUDA build steps.

2026-04-21 workaround: load ESMFold via HuggingFace Transformers
(`EsmForProteinFolding` from `facebook/esmfold_v1`) which bundles the required
openfold code. The RTX 5090 (sm_120) is not supported by the installed
torch 2.5.1+cu121 (supports up to sm_90), so inference runs on CPU (fp32).
Predictions are then aligned to 4CEU chain C with `tmtools.tm_align` and
TM-score (normalised to prediction and to reference) + RMSD are reported for
each of the 6 MICP-complete UreCs.

### Results (all 6 heroes folded successfully)

| MAG  | pred_len | TM-norm(pred) | TM-norm(4CEU) | RMSD (Å) |
|------|----------|---------------|----------------|----------|
| S13  | 616      | 0.576         | 0.620          | 4.17     |
| S16  | 616      | 0.570         | 0.613          | 4.31     |
| S23  | 616      | 0.569         | 0.612          | 4.34     |
| C22  | 573      | 0.673         | **0.678**      | **3.52** |
| M1   | 566      | 0.600         | 0.597          | 4.18     |
| S26  | 572      | 0.596         | 0.599          | 4.04     |

**All 6 MICP-complete hero UreCs produce predictions with TM-score > 0.5 and
backbone RMSD < 4.4 Å against the *S. pasteurii* urease α-subunit**, indicating
that the (β/α)₈ TIM-barrel fold is conserved across both the novel
Sphingobacterium lineage (S13/S16/S23/C22) and the divergent Pseudomonas_E
lineage (M1/S26). This corroborates the MSA-based active-site conservation
result (A2: 42/42 residue match) at the whole-chain structural level.

Output: `results/additional/C4_esmfold/ureC_vs_4CEU_tm.csv` +
`pdb/{HERO}.pdb`. Figure: `figures/additional/C4_esmfold_superposition.png`.

## C5 Pan-MICP environment ANI

```
Ref_file	Query_file	ANI	Align_fraction_ref	Align_fraction_query	Ref_name	Query_name
HERO_C22.fna	HERO_C22.fna	100.00	99.52	99.52	...
HERO_S23.fna	HERO_C22.fna	96.78	46.22	59.76	...
HERO_S16.fna	HERO_C22.fna	93.54	29.23	37.89	...
HERO_S26.fna	HERO_S26.fna	100.00	99.77	99.77	...
Pseudomonas_helleri.fna	HERO_S26.fna	97.54	83.59	89.36	NZ_JYLD01000010.1 ...
HERO_S23.fna	HERO_S23.fna	100.00	...
```

Key external hit: **S26 ↔ *P. helleri* DSM 29165 = 97.54% ANI** - same
assignment as the A3 screen with 146 NCBI genomes. Sphingobacterium heroes
remain <95% against all environmental refs, reinforcing the novel-species
designation.

## C6 In-situ abundance proxy

(SPAdes contig `cov_X` length-weighted; raw FASTQ unavailable.)

```
                metric  hero_mean  rest_mean  hero_median  rest_median     MWU_p
0  length_weighted_cov  26.756106  25.545793     5.751760     9.046683  0.272327
1           median_cov  23.963773  24.645675     5.040609     8.572691  0.333487
2             mean_cov  25.546721  24.735012     5.306006     9.178483  0.314320
```

```
    source  count       mean     median        std
0   cattle     22   5.852739   4.666129   3.624958
1  poultry     42  28.631291   9.748761  47.308926
2    sheep     15  57.427073  33.271047  62.777102
3    swine     32  20.317635  11.993405  24.865754
```

Heroes carry **no abundance advantage** over rest (MWU p ≈ 0.3), so the MICP
trait is not an artefact of outlier-dominant populations. Sheep-rumen MAGs are
the most deeply covered.
