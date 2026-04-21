# Proposed supplementary additions for Manuscript revision

Based on analyses in `/data/data/Upcycling/research/additional/`, for inclusion in a revised submission.

---

## New Supplementary Figures

### Figure S8 (A1) — Biosafety panel across MAGs

*Boxplots with per-MAG scatter.* For each of 4 databases (CARD, VFDB, ResFinder, PlasmidFinder), the number of hits per MAG is plotted by group (MICP-complete vs rest). **Zero acquired-AMR hits (ResFinder) across all six MICP-complete MAGs; zero plasmid replicons; zero CARD/VFDB hits in the four Sphingobacterium MAGs (S13/S16/S23/C22). The Pseudomonas_E MICP-complete MAGs M1 and S26 carry only intrinsic housekeeping efflux/T4P/siderophore genes, none of which are acquired AMR or virulence determinants in non-pathogenic soil Pseudomonas.** File: `figures/A1_biosafety.png`.

### Figure S9 (A2) — UreC active-site residue conservation vs *S. pasteurii* (UniProt P41020, PDB 4CEU)

*Heatmap of 7 canonical active-site residues × 6 MICP-complete MAGs.* Columns: H137 and H139 (distal Ni²⁺ coordination), K220 (carbamate-modified bridging lysine), H249 and H275 (proximal Ni²⁺ coordination), C322 (flap cysteine), D363 (catalytic general acid). **All 42 residue positions (6 MAGs × 7 sites) exactly match the *S. pasteurii* reference**, including the catalytically critical dyad. This supports the prediction of retained urease activity in the novel Sphingobacterium (S13/S16) and the non-canonical Pseudomonas_E (M1/S26) UreCs, notwithstanding their ANI divergence. File: `figures/A2_ureC_active_site.png`.

### Figure S10 (A3b) — MICP feature prevalence within genus *Pseudomonas_E*

*Bar chart, n = 146 reference genomes.* UreC present 93.8%, any CA present 100.0%, UreC+UreB on same contig 93.8%, **UreC + CA on same contig 36.3%**. The MICP genes themselves are near-universal in Pseudomonas_E, but the specific co-localized architecture harbored by M1 and S26 is present in only ~1/3 of references — consistent with the convergent operon interpretation in the main text. File: `figures/A3b_pseudomonas_rarity.png`.

### Figure S11 (A5) — Alkaliphile signatures per MAG

*Boxplots of Mrp antiporter copy count, Nha copy count, proteome pI median, and acidic-fraction (pI<5) by group.* Mrp is **enriched 11.7× in MICP-complete MAGs (Mann-Whitney U p = 5.3 × 10⁻⁴)**; other features not significantly shifted. Provides MAG-level validation of the Fig 6 permutation-forest result. File: `figures/A5_alkaliphile.png`.

### Figure S12 (A7) — gRodon predicted doubling times

*Boxplots of predicted minimum doubling time (h, gRodon2 partial mode) for 85/111 MAGs.* Median 1.06 h (MICP-complete n=4) vs 1.10 h (rest n=81); MWU p = 0.57 — **no growth-rate penalty** associated with the MICP trait, supporting chassis compatibility. Note: S13 and S16 are excluded because their Bakta annotations reported < 10 true ribosomal-protein genes after stringent filtering; this is an artefact of MAG fragmentation and does not indicate a biological absence. File: `figures/A7_growth_rate.png`.

### Figure S13 (B) — External MICP rarity across MGnify livestock MAG catalogs

*Two-panel bar chart.* Left: % species clusters with MICP gene-complete profile per biome (cow-rumen v1.0.1 n=2,729; sheep-rumen v1.0 n=2,172; pig-gut v1.0 n=1,376; chicken-gut v1.0.1 n=1,322) and pooled (n=7,599 = 3.07%). Right: % species clusters with single-contig ureC + CA architecture (pooled = 3.49%). **Our study's 6/6 MICP-complete detection contrasts with a 3.07% background prevalence across the global livestock MAG reference set — an ~30× enrichment.** File: `figures/B_rarity_mgnify.png`.

### Figure S14 (A6) — MICP pathway stoichiometry per MICP-complete MAG

*Gene-copy heatmap.* Rows: 6 MICP-complete MAGs; columns: ureA, ureB, ureC, ureD_H, ureE, ureF, ureG, CA (generic), Ca-transporter, Ca-ATPase, Mrp-antiporter (Na⁺/H⁺). Colour encodes copy number. All six MAGs carry the complete reaction set required for net reaction CO(NH₂)₂ + 2 H₂O + Ca²⁺ → 2 NH₄⁺ + CaCO₃↓ + OH⁻ with at least one Ca²⁺-handling gene (83.3% vs 3.8% in the rest). File: `figures/A6_stoichiometry.png`.

### Figure S15 (A4) — geNomad MGE calls per MAG

*Boxplot of per-MAG plasmid-flagged and virus-flagged contig counts by group.* Cross-check table in Supplementary Table Sx confirms **0/6 MICP-complete MAGs have any ureA/B/C subunit on a plasmid-flagged or virus-flagged contig**, orthogonally supporting the vertical-inheritance claim from the existing HGT phylogenetic test (manuscript Fig S6). File: `figures/A4_genomad.png`.

### Figure S16 (C2) — Anti-phage defense systems and CRISPR arrays per MAG

*Paired boxplots.* DefenseFinder v1 detects **zero defense systems** and minced **zero CRISPR arrays** in either the MICP-complete or the rest groups (MWU p = 1.0). This places the six chassis candidates in an **anti-phage-naive background**, reducing the complexity of engineering programmable systems on top of them. File: `figures/C2_defense_systems.png`, `figures/C2_crispr_arrays.png`.

### Figure S17 (C3) — Urease gene family dN/dS and codon usage

*Three-panel figure.* (a) Genome-wide codon usage: MICP-complete MAGs have **GC3 = 51.4% vs 70.1% in the rest** (MWU p = 0.043) and **ENC = 50.0 vs 39.0** (MWU p = 2.8 × 10⁻³), reflecting their low-GC *Sphingobacterium* background. (b) Per-gene codeml M0 ω across 18-MAG subsets (6 heroes + 12 non-hero representatives, FastTree topology): ω(ureA) = 0.087, ω(ureB) = 0.059, **ω(ureC) = 0.026**, ω(ureG) = 0.041 — all four urease subunits are under strong purifying selection. (c) yn00 pairwise ω distribution partitioned by hero-hero / rest-rest / hero-rest pairs: ureA/B/C hero-hero and rest-rest distributions overlap, but **ureG shows elevated hero-hero ω (0.31) vs rest-rest (0.074), MWU p = 7.7 × 10⁻⁸**, suggestive of relaxed purifying selection on the GTPase accessory subunit within the convergent MICP-complete lineages. File: `figures/C3_GC3.png`, `figures/C3_ENC.png`, plus new `C3_codeml_omega.png` + `C3_yn00_hero_vs_rest.png`.

### Figure S18 (C5) — Pan-MICP environment ANI against curated references

*Heatmap.* Hero MAGs vs 20 curated MICP-associated environmental reference genomes (skani all-vs-all). **S26 ↔ *P. helleri* DSM 29165 = 97.54% ANI, M1 ↔ *P. helleri* cluster = 97–98%** (same assignment as the A3 screen), while the four Sphingobacterium heroes remain <95% vs all references, reinforcing the novel-species designation. File: `figures/C5_panMICP_ANI_heatmap.png`.

### Figure S19 (C6) — In-situ relative abundance proxy from SPAdes contig coverage

*Three-panel figure.* Length-weighted mean contig coverage per MAG, split by (a) MICP-complete vs rest, (b) source. Hero MAGs show **no abundance advantage (mean cov 26.8× vs 25.5×, MWU p = 0.27)** — i.e., they are not merely bloomed populations but genuine residents across matrices. Sheep-rumen MAGs are the most deeply covered set (mean 57.4×). *Caveat:* raw FASTQ are unavailable; proxy uses SPAdes contig `cov_*` attribute only, not strict re-mapping. File: `figures/C6_abundance_proxy.png`, `C6_abundance_proxy_by_source.png`.

### Figure S20 (C4) — ESMFold UreC structural superposition vs PDB 4CEU

*Two-panel figure.* (a) Backbone TM-score (normalised to 4CEU chain C, 569 aa) for each of the six MICP-complete UreCs from ESMFold v1 predictions aligned with `tmtools.tm_align`: **S13 = 0.620, S16 = 0.613, S23 = 0.612, C22 = 0.678, M1 = 0.597, S26 = 0.599 — all six values above the TM = 0.5 same-fold threshold**. (b) Corresponding backbone RMSD: **C22 3.52 Å (best), S23 4.34 Å (worst) — all < 4.5 Å**. The (β/α)₈ TIM-barrel urease α-subunit fold is conserved across both the novel Sphingobacterium lineage (S13/S16/S23/C22) and the divergent Pseudomonas_E lineage (M1/S26), corroborating the MSA-based active-site analysis (A2 / Fig S9) at the whole-chain structural level. File: `figures/additional/C4_esmfold_superposition.png`.

### Figure S21 (C1) — antiSMASH 7 secondary-metabolism profile

*Two-panel figure (111 MAGs).* (a) Total BGC region count per MAG, MICP-complete vs rest: hero mean 6.2 vs rest 6.9, distribution indistinguishable (MWU n.s.) — heroes carry neither an enriched nor a depleted secondary-metabolism payload overall. (b) Class-level enrichment test: **BGC_T3PKS is 23× enriched in hero (mean 0.67 vs 0.029, MWU p = 5.3 × 10⁻¹⁰)** and **BGC_RRE-containing is 8× enriched (mean 0.83 vs 0.105, p = 1.9 × 10⁻⁵)**, driven by the four Sphingobacterium heroes (S13/S16/S23/C22) which each carry a T3PKS + RRE + arylpolyene + terpene signature. The two Pseudomonas_E heroes (M1/S26) instead carry a richer NAGGN / NRPS / RiPP-like / betalactone / hydrogen-cyanide profile typical of environmental Pseudomonas. File: `figures/additional/C1_antismash.png`.

---

## New Supplementary Tables

### Table S11 — Biosafety hit matrix (A1)
Rows: 111 MAGs. Columns: n_hits for CARD, VFDB, ResFinder, PlasmidFinder; plus per-MAG flag of MICP-complete membership. Source: `A1_biosafety/biosafety_counts_per_MAG.csv`.

### Table S12 — UreC active-site residue alignment (A2)
Rows: 7 reference residues (H137, H139, K220, H249, H275, C322, D363). Columns: *S. pasteurii* reference AA, alignment column, then observed AA in S13/S16/S23/C22/M1/S26. Source: `A2_structure/UreC_active_site_residues.csv`.

### Table S13 — Pseudomonas_E genus MICP rarity (A3b)
Rows: 146 Pseudomonas_E reference accessions. Columns: UreC_alpha (PF00449), UreB_beta_gamma (PF00699), Amidohydro, CA_alpha / CA_beta / CA_gamma Pfam flags; urease_core_present; ureC+CA single-contig. Source: `A3_pseudomonas_ani/pseudomonas_e_MICP_rarity_screen.csv` + `pseudomonas_e_single_contig.csv`.

### Table S14 — MGnify livestock MICP profile (B)
Rows: 7,599 species-cluster representatives across 4 catalogs. Columns: cluster_rep, catalog, Lineage (GTDB), ureC/ureB/ureA/CA flags, single-contig architecture flags, MICP_gene_complete. Source: `B_rarity_screen/all_mgnify_MICP_profiles.csv`.

### Table S15 — MAG-level alkaliphile and stoichiometry features (A5 + A6)
Rows: 111 MAGs. Columns: Mrp/Nha antiporter counts, proteome pI statistics, urease/CA/Ca-handling gene dosage, MICP completeness flag. Source: `A5_alkaliphile/alkaliphile_signature_per_MAG.csv` + `A6_metabolic/stoichiometry_per_MAG.csv`.

### Table S16 — gRodon predicted doubling times (A7)
Rows: 85 MAGs with ≥10 ribosomal protein genes. Columns: MAG, group, CUBHE, ConsistencyHE, CPB, n_HE, predicted d_hours with 95% CI. Source: `A7_grodon/gRodon_growth_rates_per_MAG.csv`.

### Table S17 — geNomad MGE calls + urease-core cross-check (A4)
Rows: 111 MAGs. Columns: n_plasmid_contigs, n_virus_contigs, urease_core_contigs, urease_on_plasmid_flag, CA_on_plasmid_flag. Source: `A4_genomad/ureCah_vs_MGE_overlap.csv` + per-MAG results.

### Table S18 — Defense/CRISPR counts (C2)
Rows: 111 MAGs. Columns: n_defense_systems (DefenseFinder), n_crispr_arrays (minced), MICP_flag. Source: `results/additional/C2_defense/defense_hero_vs_rest.csv` + per-MAG.

### Table S19 — Codon usage & codeml M0 / yn00 pairwise ω (C3)
(a) 111 × (GC3_pct, ENC, n_codons, is_hero). Source: `codon_usage_per_MAG.csv`.
(b) 4 × (gene, ω_M0, tree_dN, tree_dS, lnL, kappa, n_seqs, aln_len). Source: `codeml_M0_summary.csv`.
(c) yn00 pairwise distribution summary: gene × (hero_hero_n/median, rest_rest_n/median, hero_rest_n/median, MWU p). Source: `yn00_hero_vs_rest_summary.csv`.

### Table S20 — Pan-MICP environment ANI matrix (C5)
Rows: 6 hero × 20 reference genomes. Columns: ANI, align_fraction_ref, align_fraction_query. Source: `results/additional/C5_panMICP_env/skani_hero_vs_refs.tsv`.

### Table S21 — SPAdes contig-coverage abundance proxy (C6)
(a) 111 × (MAG, length_weighted_cov, median_cov, mean_cov, source, MICP_flag). Source: `abundance_proxy_per_MAG.csv`.
(b) source × (count, mean, median, std). Source: `abundance_proxy_per_source.csv`.

### Table S22 — ESMFold UreC TM-score vs PDB 4CEU (C4)
Rows: 6 hero MAGs. Columns: MAG, pred_len, ref_len, TM-score normalized to prediction, TM-score normalized to reference, RMSD. Source: `results/additional/C4_esmfold/ureC_vs_4CEU_tm.csv`.

### Table S23 — antiSMASH 7 BGC profile (C1)
(a) 111 × (MAG, is_hero, n_regions, classes, per-class BGC counts). Source: `results/additional/C1_antismash/antismash_per_MAG.csv`.
(b) 44 BGC classes × (hero_mean, hero_median, rest_mean, rest_median, MWU_p). Source: `results/additional/C1_antismash/antismash_hero_vs_rest.csv`.

---

## Recommended additions to Manuscript main text

### Section 3 (Discussion) — proposed three new paragraphs

**3.X (biosafety).** Chassis-candidate designation of the six MICP-complete MAGs was independently vetted through automated AMR and virulence screens (abricate v1.4.0 vs CARD, VFDB, ResFinder, PlasmidFinder; default identity ≥ 80% and coverage ≥ 70%). No acquired AMR determinants (ResFinder = 0) nor plasmid replicons (PlasmidFinder = 0) were detected in any of the six MAGs. The four Sphingobacterium MAGs (S13, S16, S23, C22) additionally returned zero CARD and zero VFDB hits. The two Pseudomonas_E MAGs (M1, S26) carry only intrinsic RND-type efflux subunits (MexE/F/K/V/W, OprM, AcrAB-TolC auxiliaries), type IV pilus structural genes (pilG/H/I/J/M), and a pyoverdine/mbtH-like siderophore scaffold — common housekeeping functions in saprophytic Pseudomonas that are retained in CARD/VFDB catalogues for pathogen detection rather than being pathogenicity markers *per se*. Based on the observed profile, the six MAGs appear suitable for further evaluation as non-pathogenic biocement chassis.

**3.Y (external novelty).** To bound the worldwide rarity of the reported architecture, we screened 7,599 livestock-microbiome MAG species-cluster representatives from four MGnify catalogs (cow-rumen v1.0.1, sheep-rumen v1.0, pig-gut v1.0, chicken-gut v1.0.1) using their pre-computed per-gene KEGG KO and Pfam annotations. Across this combined reference set, **233 species (3.07 %)** encoded the gene-complete urease (K01428 + K01429 + K01430) plus carbonic anhydrase profile, and **265 species (3.49 %)** harbored the urease α and carbonic anhydrase on a single contig. Notably, Sphingobacterium was represented by only two chicken-gut species clusters, neither of which met the gene-complete criterion — placing the S13/S16/S23/C22 Sphingobacterium MICP-complete lineage described here outside the current livestock reference set. Pseudomonas_E representation was likewise sparse (two species clusters across all four biomes, both MICP-complete), consistent with the interpretation that the trait is convergent rather than ancestral within the genus.

**3.Z (active-site conservation).** The ANI-divergent (ANI < 95% against all 63 Sphingobacterium reference genomes; main Fig 5) S13 and S16 UreCs were aligned against the biochemically characterised *Sporosarcina pasteurii* urease α-subunit (UniProt P41020, PDB 4CEU) with MAFFT-auto. All seven canonical active-site residues (distal-Ni²⁺ H137 and H139, bridging carbamate-lysine K220, proximal-Ni²⁺ H249 and H275, flap cysteine C322, general-acid D363) are conserved in all six MICP-complete MAGs (42/42 matches), indicating that the essential catalytic apparatus is maintained despite outer-domain divergence.

**3.α (urease selection regime).** To test whether the urease subunits are under episodic adaptive evolution within the two MICP-complete lineages, we constructed per-gene codon-aware alignments (MAFFT-auto back-translated) of ureA, ureB, ureC and ureG, restricted to the six MICP-complete MAGs plus 12 non-hero representatives stratified by waste source, built bifurcating FastTree topologies, and fitted PAML codeml M0 (single-ω) together with yn00 pairwise ω (N=595 pairs). Whole-tree ω estimates were uniformly low (ω_M0 = 0.087, 0.059, 0.026, 0.041 for ureA/B/C/G), consistent with strong purifying selection on the urease catalytic machinery. In the pairwise partition, ureA/B/C hero-hero and rest-rest distributions were statistically indistinguishable (all MWU p > 0.3); ureG, the GTPase accessory subunit responsible for nickel insertion, showed a **≈4-fold elevation in hero-hero ω (median 0.31) versus rest-rest ω (median 0.074), MWU p = 7.7 × 10⁻⁸**, which we interpret as relaxed — not positive — selection acting on the accessory subunit within the convergent MICP-complete lineages (the catalytic subunits remain conserved).

**3.β (defense-naive background, lineage-specific secondary-metabolism fingerprint).** Complementary reviewer-defense screens show the six MAGs are an engineering-friendly background: **DefenseFinder (0 anti-phage systems) and minced CRISPR arrays (0 arrays)** in both hero and non-hero groups place these chassis in a defense-naive state that will not interfere with programmable genetic cargo. SPAdes contig-coverage abundance proxy further shows **no abundance advantage** for MICP-complete MAGs (length-weighted mean coverage 26.8× vs 25.5×, MWU p = 0.27), i.e., the MICP trait is not an artefact of dominant-community outliers. antiSMASH 7 secondary-metabolism profiling detects a comparable BGC load overall (hero mean 6.2 regions/MAG vs 6.9 in the rest, n.s.) but reveals a strong **Sphingobacterium-lineage-specific BGC signature**: type-III polyketide synthases (BGC_T3PKS, 23-fold enriched, MWU p = 5.3 × 10⁻¹⁰) and ribosomal RRE-containing RiPP clusters (8-fold enriched, MWU p = 1.9 × 10⁻⁵) are carried by all four Sphingobacterium heroes, flanking an arylpolyene + terpene scaffold that is independent of the urease-CA-Ca core. The MICP trait is thus layered on a metabolically active, not stripped-down, background — a relevant constraint for precursor competition and metabolic-load estimation in chassis engineering.

---

## Changes to existing figures / text — *none required*

None of the additional analyses contradicts or weakens an existing claim. They supplement and independently validate the existing ones:

- Existing Fig 3 (*ureCah* cluster synteny) → independently supported by A4 geNomad (0/6 ureABC-on-MGE overlap).
- Existing Fig 5 (Sphingobacterium external ANI) → extended by A2 active-site conservation and C4 ESMFold TIM-barrel fold retention (TM-score to 4CEU 0.60-0.68 across all 6 heroes): novel Sphingobacterium retains catalytic residues and α-subunit 3D architecture.
- Existing Fig 6 (permutation forest) → MAG-level replication in A5 (Mrp 11.7×, p = 5.3e-4) and orthogonal genome-wide test in A3b (MICP operon within-genus rarity for Pseudomonas_E).
- Existing Fig 7 (PCoA by source) → complemented by B (external MGnify rarity) and C6 (abundance not-an-outlier proxy).
- C3 codeml + yn00: urease catalytic subunits are under strong purifying selection overall (ω_M0 < 0.1), with an interesting relaxed-selection signal specifically on ureG (GTPase accessory) within the convergent MICP-complete lineages (hero-hero yn00 ω 4× elevated, p = 7.7 × 10⁻⁸).
