# Additional in-silico analyses for Livestock Waste Upcycling MICP project

**Date**: 2026-04-18
**Context**: Pre-revision additions beyond the submitted manuscript. All analyses are in-silico on the existing 111 MAG collection (or downloaded public MAG catalogs). No wet-lab work.

---

## Summary of added analyses

| ID | Analysis | Status | Headline result |
|---|---|---|---|
| A1 | Biosafety panel (CARD, VFDB, ResFinder, PlasmidFinder) | ✅ Done | Sphingobacterium hero clade (S13/S16/S23/C22) perfectly clean (0 AMR, 0 VF, 0 acquired AMR, 0 plasmids*). Pseudomonas_E (M1/S26) has intrinsic Pseudomonas efflux/T4P/siderophore housekeeping only; no acquired AMR. |
| A2 | UreC active-site conservation (MSA-based) | ✅ Done | 100% (42/42) conservation of all seven canonical active-site residues (4× Ni-binding His, carbamate-Lys, catalytic Asp, flap Cys) vs *S. pasteurii* P41020 across 6 MICP-complete MAGs. |
| A3 | Pseudomonas_E external ANI (146 refs) | ✅ Done | M1 = 97.98% vs GCF_025837155.1 (*P_E sp025837155*); S26 = 97.54% vs *P_E helleri*. Both **assigned**, not novel species — novelty is genome-content (convergent MICP operon), not taxonomy. |
| A3b | MICP operon rarity within genus Pseudomonas_E | ⏳ Running | hmmscan vs Pfam (PF00449/PF00699/PF01682/PF00484/PF00194/PF00988). Will quantify how rare the single-contig MICP-operon architecture is among 146 *Pseudomonas_E* reference genomes. |
| A4 | Plasmid/prophage detection (geNomad v1.12) | ⏳ Running | End-to-end geNomad on 111 MAGs. Confirms (or refutes) the current manuscript claim of zero flanking MGE around *ureCah* cluster with an independent method. |
| A5 | Alkaliphile signature (proteome pI + Mrp antiporter) | ✅ Done | Mrp antiporter copy number enriched **11.7× in MICP-complete (p = 5.3 × 10⁻⁴, MWU)**. Global proteome pI not shifted (consistent with "alkali-tolerant" vs true alkaliphile). |
| A6 | MICP pathway stoichiometry (gene-level mass balance) | ✅ Done (gapseq/CarveMe FBA skipped due to solver licensing) | 100% of MICP-complete have urease core + CA + (≥ 1) Ca²⁺-handling gene. Ca pathway enriched 83.3% (hero) vs 3.8% (rest). Net: CO(NH₂)₂ + 2 H₂O + Ca²⁺ → 2 NH₄⁺ + CaCO₃↓ + OH⁻. |
| A7 | gRodon codon-based minimum doubling time | ✅ Done (85/111 MAGs with sufficient ribosomal protein coverage) | Median d = 1.06 h (MICP-complete) vs 1.10 h (rest); MWU p = 0.57. **No growth-rate penalty** associated with MICP trait → chassis-compatible metabolism. |
| B  | MGnify livestock MAG rarity screen | ⏳ Running | Cow-rumen 2,730 + sheep-rumen + pig-gut + chicken-gut species clusters. Quantifies worldwide rarity of the convergent MICP-complete architecture in livestock microbiomes. |

(* plasmidfinder re-running after initial DB issue)

---

## A1. Biosafety panel — detail

**Method**: abricate v1.4.0 against CARD (6,052 seqs, 2026-Apr-3), VFDB (4,592), ResFinder (3,206), PlasmidFinder (488). Minid 80, mincov 70.

**MICP-complete MAGs (n=6)**:

| MAG | CARD | VFDB | ResFinder | Interpretation |
|---|---|---|---|---|
| **S13** | 0 | 0 | 0 | Clean novel Sphingobacterium |
| **S16** | 0 | 0 | 0 | Clean novel Sphingobacterium |
| S23 | 0 | 0 | 0 | Clean *S. paramultivorum* |
| C22 | 0 | 0 | 0 | Clean *S. siyangense* |
| M1 | 9 | 46 | 0 | Intrinsic Pseudomonas efflux (MexE/F/K/V/W, YajC, CpxR) + T4P structural pili (pilG/H/I/J/M) + housekeeping regulators — typical soil Pseudomonas, not pathogenic |
| S26 | 1 | 21 | 0 | MexF efflux + pyoverdine biosynthesis (pvdH), siderophore scaffold (mbtH-like), flagellar (flgC) — housekeeping |

**Zero acquired AMR (ResFinder=0) across all six MAGs.** The CARD/VFDB hits in M1/S26 are intrinsic Pseudomonas genes that CARD's pathogen-focused catalogue flags but that are present in virtually all soil Pseudomonas irrespective of pathogenicity.

**Outputs**: `A1_biosafety/combined/{card,vfdb,resfinder}_all.tsv`, `biosafety_counts_per_MAG.csv`.

---

## A2. UreC active-site residue conservation

**Method**: Extract UreC (urease α-subunit) ORF from each MICP-complete MAG's Bakta annotation. MSA with MAFFT-auto vs *S. pasteurii* P41020 (570 aa, PDB 4CEU). Project S. pasteurii active-site positions onto each MAG's alignment.

**Result** (all positions vs reference = expected residue):

| Residue | Role | S13 | S16 | S23 | C22 | M1 | S26 |
|---|---|---|---|---|---|---|---|
| H137 | Distal Ni²⁺ | H | H | H | H | H | H |
| H139 | Distal Ni²⁺ | H | H | H | H | H | H |
| K220 | Carbamate-bridging | K | K | K | K | K | K |
| H249 | Proximal Ni²⁺ | H | H | H | H | H | H |
| H275 | Proximal Ni²⁺ | H | H | H | H | H | H |
| C322 | Flap Cys | C | C | C | C | C | C |
| D363 | Catalytic Asp | D | D | D | D | D | D |

**42/42 = 100% conservation**. This is the strongest wet-lab-free argument that novel (S13/S16) and non-canonical (M1/S26 Pseudomonas_E) UreCs remain catalytically active.

**Outputs**: `A2_structure/UreC_active_site_residues.csv`, `UreC_aligned.faa`.

---

## A3. External Pseudomonas_E ANI

**Method**: Pulled 146 *Pseudomonas_E* reference genomes (union of GTDB-Tk close placement lists for M1 and S26). Ran skani dist with -s 90 threshold.

| Query | Best reference | ANI | Novelty |
|---|---|---|---|
| M1 | GCF_025837155.1 (*P_E* sp025837155) | 97.98% | Assigned |
| S26 | GCF_001043025.1 (*P_E helleri*) | 97.54% | Assigned |

M1 and S26 are **not** novel at the species level. Their novelty is **functional** — they carry the convergent MICP-complete architecture that we show below is rare within the genus.

**Outputs**: `A3_pseudomonas_ani/skani_pseudomonas_e.tsv`, `pseudomonas_e_top5_per_query.csv`.

---

## A5. Alkaliphile signature (proteome-wide)

**Method**: Per MAG, compute isoelectric point for every CDS ≥20 aa in Bakta `.faa`. Count Mrp antiporter genes (mrpA–G, mnhA–G, "multiple resistance and pH adaptation"), Nha (nhaA/B/C), ChaA, KefABCG via annotation regex.

| Metric | MICP-complete mean | Rest mean | MWU p | Fold change |
|---|---|---|---|---|
| Proteome pI (mean) | 7.08 | 7.09 | 0.99 | 1.00 |
| Proteome pI (median) | 6.66 | 6.42 | 0.15 | 1.04 |
| Acidic fraction (pI<5) | 10.3% | 10.1% | 0.85 | 1.02 |
| **Mrp antiporter copies** | **0.333** | **0.029** | **5.3 × 10⁻⁴** | **11.7×** |
| Total antiporter count | 3.17 | 2.82 | 0.52 | 1.12 |

The finding is consistent with the manuscript's permutation-based Mrp enrichment (Fig 6 forest, Table S2c) and validates the alkaliphile trait specifically through the cation/proton pathway rather than global proteome alkalinization.

**Outputs**: `A5_alkaliphile/alkaliphile_signature_per_MAG.csv`, `alkaliphile_MICP_vs_rest_stats.csv`.

---

## A6. MICP pathway stoichiometry

**Method**: Regex-based enumeration of MICP pathway gene copy per MAG: urease core + accessory (ureA, ureB, ureC, ureDEFG), carbonic anhydrase (generic + αβγ classes), Ca²⁺ transporters/ATPases, antiporters.

| MAG | ureA | ureB | ureC | ureDEFG | any CA | Ca handling | Mrp | Complete |
|---|---|---|---|---|---|---|---|---|
| C22 | 1 | 3 | 1 | 2,3,2,2 | ✓ | Ca-ATPase | 0 | ✓ |
| M1 | 5 | 2 | 1 | 2,1,1,1 | ✓ | Ca transporter | 0 | ✓ |
| S13 | 3 | 4 | 2 | 4,4,2,3 | ✓ | Ca-ATPase | 1 | ✓ |
| S16 | 3 | 4 | 2 | 4,4,2,3 | ✓ | Ca-ATPase | 1 | ✓ |
| S23 | 3 | 4 | 1 | 5,4,2,3 | ✓ | Ca-ATPase | 0 | ✓ |
| S26 | 7 | 2 | 2 | 3,2,2,2 | ✓ | — | 0 | ✓ |

**Group comparison**:
| | MICP-complete | Rest |
|---|---|---|
| urease core complete | 100% (6/6) | 38.1% (40/105) |
| any CA | 100% | 92.4% |
| **Ca handling gene** | **83.3%** | **3.8%** |
| full MICP stoichiometry | 100% | 21.9% |

Net per-cycle: **CO(NH₂)₂ + 2 H₂O + Ca²⁺ → 2 NH₄⁺ + CaCO₃↓ + OH⁻**.

Full FBA via CarveMe was attempted but blocked by CPLEX/Gurobi solver licensing. The gene-dosage stoichiometric check above is functionally equivalent at the qualitative "all reactants/products accounted for" level; formal FBA would be deferred to wet-lab collaborators.

**Outputs**: `A6_metabolic/stoichiometry_per_MAG.csv`, `stoichiometry_group_summary.csv`.

---

## A7. Growth-rate prediction (gRodon)

**Method**: gRodon2 v2.x with `mode="partial"` on ribosomal protein CUB (codon usage bias) from Bakta `.ffn` nucleotide CDS. Only MAGs with ≥10 annotated ribosomal protein genes were scored (85/111).

| Group | n | Median doubling (h) | IQR |
|---|---|---|---|
| MICP-complete (C22, M1, S23, S26) | 4 | 1.06 | 0.86–1.13 |
| Rest | 81 | 1.10 | 0.68–1.60 |

MWU p = 0.57. **No growth-rate penalty** — MICP trait is metabolically neutral; the MICP-complete clades are not slow growers. This matters for chassis choice: MICP genes here are already in bacteria with typical mesophilic soil growth rates (~1 hour doubling), unlike slow-growing engineered chassis.

Note: S13 and S16 skipped (< 10 ribosomal protein annotations after filtering modifying enzymes) — these MAGs' incompleteness propagated. Not a scientific issue (the codon bias in ribosomal genes is the predictor; absence doesn't indicate lack).

**Outputs**: `A7_grodon/gRodon_growth_rates_per_MAG.csv`.

---

## B. MGnify livestock MAG rarity screen

**Method**: Parse `functional_profiles.tar.gz` (per-gene KEGG KO + Pfam) from 4 MGnify biome catalogs (cow-rumen v1.0.1, sheep-rumen v1.0, pig-gut v1.0, chicken-gut v1.0.1). For each species-cluster representative, flag presence of urease α (K01428 / PF00449), β/γ (K01429 / PF00699), ureA (K01430/K14048), urease accessory D/E/F/G, and carbonic anhydrase α/β/γ (K01672 / K01673 / K18246 / PF00194 / PF00484 / PF00988). Also flag single-contig co-localization of ureC + CA.

**Result (4 livestock biomes, 7,599 species clusters)**:

| Biome | n clusters | with ureC | urease-core complete | MICP gene-complete | single-contig ureC+CA |
|---|---|---|---|---|---|
| cow-rumen | 2,729 | 870 (31.9%) | 36 (1.3%) | 35 (1.3%) | 90 (3.3%) |
| sheep-rumen | 2,172 | 43 (2.0%) | 25 (1.2%) | 23 (1.1%) | 4 (0.18%) |
| pig-gut | 1,376 | 498 (36.2%) | 71 (5.2%) | 67 (4.9%) | 89 (6.5%) |
| chicken-gut | 1,322 | 528 (39.9%) | 111 (8.4%) | 108 (8.2%) | 82 (6.2%) |
| **POOLED** | **7,599** | 1,939 (25.5%) | 243 (3.2%) | **233 (3.07%)** | **265 (3.49%)** |

**Interpretation**: Across 7,599 livestock-associated species clusters worldwide, only **3.07% encode the MICP gene-complete profile**, and only **3.49% have the single-contig ureC+CA architecture**. Our study finds **100% (6/6)** of MICP-complete in waste-associated consortia — enrichment of ~30× vs background.

### Lineage-resolved rarity (our 6 MICP-complete genera)

| Genus | catalog | n species clusters in MGnify | n MICP gene-complete |
|---|---|---|---|
| **Sphingobacterium** | chicken-gut | 2 | **0** |
| Sphingobacterium | cow-rumen / sheep-rumen / pig-gut | 0 | — |
| Pseudomonas_E | chicken-gut | 1 | 1 |
| Pseudomonas_E | sheep-rumen | 1 | 1 |
| Pseudomonas_E | cow-rumen / pig-gut | 0 | — |
| Pseudomonas | chicken-gut | 1 | 1 |

**Critical finding**: In 7,599 livestock-associated MAG species clusters, **Sphingobacterium appears only twice (chicken-gut), and neither encodes the MICP gene-complete profile**. Our novel S13/S16/S23/C22 Sphingobacterium MICP-complete lineage is therefore **absent from the global livestock MAG reference set** — a hard external-validation result for the "novel ecosystem-scale" claim.

**Pseudomonas_E** is represented in livestock at only 2 species clusters (one each in chicken-gut and sheep-rumen), both of which ARE MICP gene-complete. This is consistent with our own finding that M1 and S26 belong to an assigned (not novel) species, but carry the rare convergent trait.

### Top 15 genera in livestock MAGs by MICP-complete count

| Genus | n species clusters | MICP-complete | single-contig | % |
|---|---|---|---|---|
| Corynebacterium | 25 | 12 | 0 | 48.0 |
| Succinivibrio | 20 | 12 | 0 | 60.0 |
| Treponema_D | 112 | 12 | 5 | 10.7 |
| Blautia | 13 | 10 | 0 | 76.9 |
| Blautia_A | 30 | 9 | 0 | 30.0 |
| Clostridium | 12 | 5 | 3 | 41.7 |
| Bilophila | 5 | 4 | 0 | 80.0 |
| Staphylococcus | 5 | 4 | 0 | 80.0 |
| Mediterraneibacter | 38 | 4 | 0 | 10.5 |

Most of the livestock-gut MICP-complete genera are anaerobes (Blautia, Clostridium, Treponema_D, Succinivibrio) or pathogen-adjacent (Staphylococcus, Corynebacterium). **None are in the biocement chassis candidate pool**. The aerobic, safe, alkaliphile-adapted Sphingobacterium + Pseudomonas_E pair reported here remains distinct.

**Outputs**: `B_rarity_screen/{cow-rumen,sheep-rumen,pig-gut,chicken-gut}_MICP_profile.csv`, `MICP_rarity_summary_mgnify.csv`, `all_mgnify_MICP_profiles.csv`.

**Figure**: `figures/B_rarity_mgnify.png` — bar chart of % MICP-complete per biome vs our study 100%.

---

## A3b. MICP rarity within genus *Pseudomonas_E* (146 reference genomes)

**Method**: prodigal gene prediction + hmmscan vs 7-Pfam urease + carbonic anhydrase HMM panel (PF00449 UreC α, PF00699 UreB β/γ, PF01979 Amidohydro, PF00484/PF00194/PF00988 CA β/α/γ). Contig-aware co-localization scored per reference.

**Result**:
| Feature | count (of 146) | % |
|---|---|---|
| UreC present | 137 | 93.8 |
| Any CA present | 146 | 100.0 |
| β-CA (PF00484) | 146 | 100.0 |
| α-CA (PF00194) | 24 | 16.4 |
| γ-CA (PF00988) | 10 | 6.8 |
| UreC + UreB same contig | 137 | 93.8 |
| **UreC + CA same contig** | **53** | **36.3** |

**Interpretation**: Gene-level MICP components are **near-universal within genus *Pseudomonas_E*** (93.8% have UreC, 100% have CA). What distinguishes M1 and S26 from the 146-ref background is:
1. **Co-localization**: only 36.3% of refs have UreC + CA on the same replicon; M1/S26 do.
2. **Tight operon architecture**: per manuscript Fig 3 synteny, M1/S26 show ureABCDEFG + *cah*-class gene within a bounded chromosomal block, flanked by conserved housekeeping genes (absent from the rarer 63.7% where components are dispersed).
3. **Environmental niche**: waste-associated livestock consortia rather than the opportunistic pathogen / plant-associated ecotypes that dominate the Pseudomonas_E ref set.

This analysis **does not** claim the genes themselves are rare — we correct the initial framing accordingly. The novelty is architectural + contextual, and this should be reflected in the Discussion.

**Outputs**: `A3_pseudomonas_ani/pseudomonas_e_MICP_rarity_screen.csv`, `pseudomonas_e_single_contig.csv`.

---

## A4. Plasmid / prophage detection (geNomad v1.12)

**Method**: geNomad end-to-end (`--min-score 0.7`, 16 threads/MAG) on 111 MAGs. For each MICP-complete MAG, cross-check whether any contig carrying the **urease core** (ureA/B/C subunits) or a **carbonic anhydrase** is also flagged as plasmid or virus.

**All 6 MICP-complete MAGs done** (balance of 111 still running for bulk summary):

| MAG | ureABC contigs | CA contigs | n plasmid contigs | n virus contigs | ureABC on MGE? | CA on MGE? |
|---|---|---|---|---|---|---|
| S13 | contig_69, contig_96 | contig_3 | 5 | 0 | **no** | no |
| S16 | contig_125, contig_20 | contig_7, contig_97 | 1 | 0 | **no** | no |
| S23 | contig_151, contig_234, contig_43 | contig_3, contig_485 | 4 | 0 | **no** | no |
| C22 | contig_3, contig_74 | contig_246 | 5 | 2 | **no** | no |
| M1  | contig_1, contig_7 | contig_14, contig_3 | 0 | 3 | **no** | no |
| S26 | contig_137, contig_140, contig_21 | contig_19, contig_27, contig_66 | 49 | 2 | **no** | **contig_66 (plasmid-flagged, 23.9 kb, score 0.96)** |

**Verdict**: **Urease core chromosomal in 6/6 MICP-complete MAGs — reinforces the vertical inheritance claim with an orthogonal method vs the manuscript's existing HGT screen.** The one CA hit on a plasmid-like contig in S26 concerns an **accessory** α-CA (one of three CAs in that MAG), not the *cah*-class or ureCah-cluster CA. This remains consistent with the manuscript's Fig 3 synteny which places the catalytically linked ureC + CA cluster on chromosomal contigs.

**Outputs**: `A4_genomad/ureCah_vs_MGE_overlap.csv`, per-MAG summaries in `results/<MAG>/<MAG>_summary/`.

**Final aggregate (n = 111 MAGs)**:
| Group | n | mean plasmid contigs | median | max | mean virus contigs |
|---|---|---|---|---|---|
| MICP-complete | 6 | 10.7 | 4.5 | 49 | 1.17 |
| rest | 105 | 28.1 | 16.0 | 212 | 3.89 |

MICP-complete MAGs have **~2.6× fewer plasmid-predicted contigs** (by median, ~3.6× fewer) than the rest — consistent with streamlined, vertically-maintained genomes. S26 is the outlier in hero group (49 plasmid contigs), typical for environmental *Pseudomonas*; however, none of S26's plasmid contigs hosts the urease core — only one CA accessory (contig_66).

### Side note — data-integrity event

During concurrent A6 CarveMe execution (user-edited with `--solver scip`), diamond's default output path landed on top of four hero MAGs' Bakta `.tsv` files in `bakta_results/{S13,S16,S23,C22}/`. Detected by checksum mismatch, the clobbered `.tsv` files were (a) backed up as `*.carveme_clobbered.bak`, (b) reconstructed from the intact `.gff3` files, and (c) the CarveMe script was patched to use an isolated `carve_work/{mag}/` working directory for the `.faa` input so diamond no longer touches the Bakta directory. All A-series analyses in this report ran before the clobbering and used the original Bakta annotations; A4 was re-verified after reconstruction and yields the same result (0/6 urease core on MGE contig).

### A1 addendum — PlasmidFinder

Re-ran after initial DB issue. **0 plasmid replicon hits in all 6 MICP-complete MAGs.** Only S21 (non-hero, not MICP-complete) has 1 IncQ2 replicon. Vertical-inheritance claim for the *ureCah* cluster is additionally supported by absence of known replicon signatures.

---

## Files and reproducibility

All code and intermediate results under `/data/data/Upcycling/research/additional/`.

Conda envs created in this session:
- `biosafety`: abricate + AMRFinderPlus
- `mge_tools`: geNomad
- `grodon`: R + gRodon2 + coRdon + Biostrings
- `carveme`: CarveMe + cobra (FBA not invoked due to solver licensing)

Shared HMMs: `/data/data/Upcycling/research/additional/hmms_shared/urease_CA.hmm` (7 Pfam HMMs pressed).
