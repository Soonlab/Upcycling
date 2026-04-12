# Figure, Table and Supplementary Legends (revised)

All figures are provided as 300-dpi PNG and vector PDF files at `/data/data/Upcycling/research/` (main figures) and `/data/data/Upcycling/research/extra/` and `/data/data/Upcycling/research/revision/` (supplementary figures). Sample codes: **C** cattle, **M** swine, **S** sheep, **V** poultry. The six "hero-lineage" MAGs (C22, M1, S13, S16, S23, S26) are consistently highlighted in red/bold across all figures.

---

## Main figures

### Figure 1 | Phylogenomic distribution of MICP-related genes across 111 livestock-waste MAGs.
*(`Fig1_Phylogeny_MICP_heatmap.png`, `.pdf`)*

Maximum-likelihood phylogenomic tree of 111 MAGs inferred from the GTDB-Tk bac120 concatenated alignment (IQ-TREE v2.3, ModelFinder, 1,000 ultrafast bootstrap replicates; midpoint-rooted for display). The left colour strip encodes GTDB-Tk genus assignment (legend: top ten genera with sample counts). MAG identifiers and, where available, GTDB species annotations follow the colour strip; hero-lineage MAGs are highlighted in red bold. The right-hand heat map shows presence (dark green) or absence (white) of the eight MICP-associated genes *ureA, ureB, ureC, ureD, ureE, ureF, ureG* and *cah*. The complete 8/8 module is restricted to four *Sphingobacterium* MAGs (C22, S13, S16, S23) and two *Pseudomonas*\_E MAGs (M1, S26); these taxa are distributed across two, non-monophyletic, functionally convergent lineages. Scale bar, substitutions per site.

### Figure 2 | MICP module completeness across genera and gene-level prevalence in hero-lineage MAGs.
*(`Fig2_MICP_completeness_by_genus.png`, `.pdf`)*

**(a)** MICP module score (sum of *ureA–G* + *cah*; maximum = 8) by GTDB-Tk genus, sorted by mean. Boxes show median and inter-quartile range; individual MAGs are overlaid as grey points; hero-lineage MAGs are superimposed as enlarged red circles. **(b)** Per-gene prevalence of the eight MICP genes in the hero-lineage set (red, n = 6) versus the remaining 105 MAGs (grey). Hero-lineage MAGs maintain 100 % prevalence of every *ure* gene and *cah*, whereas the non-hero background retains the same genes at 55–85 %.

### Figure 3 | Gene-order conservation of the *ureABCDEFG* operon across hero-lineage MAGs.
*(`Fig3_ureCah_cluster_synteny.png`, `.pdf`)*

Synteny diagrams of the densest *ure*-containing contig in each hero-lineage Bakta annotation. Each horizontal track represents the ≤ 40 kb window with the greatest number of *ure* genes for one MAG (C22, M1, S13, S16, S23, S26). Arrows are drawn to scale; arrow direction encodes strand; colours encode gene identity (legend: *ureA* blue, *ureB* green, *ureC* red, *ureD* purple, *ureE* brown, *ureF* pink, *ureG* cyan, *cah* orange, other grey). All 7 *ure* genes are recovered on a single contig within a 5.9–28.6 kb window in the three best-assembled hero MAGs (M1, S13, S16); in C22, S23 and S26 a subset of *ure* genes is retained on the main contig (4–5 / 7) with the remainder on secondary contigs, consistent with assembly fragmentation (Table 1). No transposase, integrase, prophage or relaxase gene is detected in ±15 kb flanking any hero cluster (Supplementary Table S3).

### Figure 4 | DRAM metabolic-module heat map and hero-vs-rest comparison of MICP-critical modules.
*(`extra/Fig4_DRAM_metabolism_heatmap.png`, `.pdf`; `extra/Fig4b_DRAM_HeroVsRest.png`, `.pdf`)*

**(a)** Module-completeness heat map (DRAM v1.5 `distill/product.tsv`) for 111 MAGs × 35 curated KEGG/custom modules spanning central carbon metabolism, nitrogen metabolism, stress response, vitamin (cobalamin / B12) biosynthesis, and CAZyme categories. Row labels give MAG identifier and GTDB genus; hero-lineage rows are highlighted in red bold. Colour intensity encodes fractional module completeness (0 to 1). Housekeeping modules (glycolysis, TCA, fatty-acid biosynthesis) are comparably complete in hero and non-hero MAGs. **(b)** Mean module completeness of MICP-critical modules (urease, carbonic anhydrase, nitrogen metabolism, Na⁺/H⁺ antiport, cobalamin biosynthesis, CAZymes/carbohydrate metabolism) in hero MAGs (red, n = 6) versus the remaining MAGs (grey, n = 105).

### Figure 5 | Novel-species delineation of S13 and S16 within *Sphingobacterium*.
*(`extra/Fig_T2d_Novelty_overview.png`, `.pdf`; `extra/novel_species/Fig_NovelSp_AAI.png`, `.pdf`)*

**(a)** Novelty screen. Each MAG is one point; the y-axis is the GTDB closest-reference ANI, and the dashed red line marks the 95 % species cutoff. Hero-lineage MAGs are enlarged and red-labelled. Twenty-one of 111 MAGs fall below 95 %, including S13 and S16 (no species-level ANI available). **(b)** Amino-acid identity of S13 (left) and S16 (right) versus the other five *Sphingobacterium* MAGs, computed from reciprocal-best mmseqs2 easy-search hits. Dashed vertical lines mark the 95 % (species) and 70 % (genus) AAI thresholds. Maximum AAI = 93.15 % (S13 vs V3) and 93.49 % (S16 vs S23); both well below the species threshold.

### Figure 6 | Permutation-tested hero-vs-rest enrichment of trait modules.
*(`revision/Fig_Permutation_forest.png`, `.pdf`)*

Forest-style plot of fold-change enrichment (log-scale x-axis) for trait subcategories with fold change > 1, ranked by magnitude. Points are fold change; error bars are 2,000-iteration bootstrap 95 % CI. Point colour encodes significance level after Benjamini–Hochberg FDR correction on 10,000-iteration one-sided permutation p-values: **red** q < 0.05, **orange** q < 0.10, **grey** n.s. Modules significantly enriched in the hero lineage include Mrp complex (FC 10.85), CBM (9.78), oxidative-stress defence (4.76), glycoside hydrolase (4.66), quorum sensing (2.13) and Na⁺/H⁺ antiporter (2.30). Complete statistics in Supplementary Table S2 (`Hero_vs_Rest_permutation_stats.csv`).

---

## Supplementary figures

### Figure S1 | Biofilm / EPS gene modules across the 111-MAG panel.
*(`extra/Fig_T1a_Biofilm_EPS_heatmap.png`, `.pdf`)*
Keyword-based EPS/biofilm heat map (per 10³ CDS) covering pel, psl, pga, cellulose (bcs), alginate (alg), curli (csg), colanic (wca), capsule (wz*/cps), adhesin/autotransporter, and quorum-sensing subcategories. Hero-lineage rows red bold.

### Figure S2 | Ammonia-handling and nitrogen-assimilation modules.
*(`extra/Fig_T1b_Ammonia_detox_heatmap.png`, `.pdf`)*
Heat map of urea transport (urt, urea transporter), GS-GOGAT (glnA, gltBD), GDH, AmtB, PII regulators (glnB/K) and nitrate/nitrite reductase pathways.

### Figure S3 | Alkaline and osmotic-stress tolerance modules.
*(`extra/Fig_T2a_Alkaline_stress_heatmap.png`, `.pdf`)*
Per-MAG counts of *nhaA–C*, Mrp complex, Kdp / Trk K⁺ uptake, compatible-solute biosynthesis (*opu*, *pro*, *betA/B*, *ect*), chaperones (*groEL/ES*, *dnaK/J*, *clpB*) and oxidative-stress defences (*katA/B/G*, *sodA/B/C*, *ahpCF*).

### Figure S4 | CAZyme profile: keyword proxy and DRAM/dbCAN validation.
*(`extra/Fig_T2b_CAZyme_heatmap.png`, `.pdf`; `revision/Fig_dbCAN_classes.png`, `revision/Fig_dbCAN_topFamilies.png`)*
Keyword-based CAZy proxy heat map and dbCAN-derived class-level and family-level profiles. Both approaches identify hero-lineage enrichment of glycoside hydrolases and carbohydrate-binding modules.

### Figure S5 | Heavy-metal and antibiotic-resistance gene modules.
*(`extra/Fig_T2c_MetalAMR_heatmap.png`, `.pdf`)*
Heavy-metal efflux (*czc, cop, cus, ars, mer, znt*), metallothioneins, β-lactamases, MDR efflux (*acr, mex,* MATE), aminoglycoside/tetracycline/macrolide determinants.

### Figure S6 | Pairwise ANI within hero lineages.
*(`extra/Fig_T2d_HeroANI.png`, `.pdf`)*
skani whole-genome ANI heat map among the six hero-lineage MAGs.

### Figure S7 | UreC gene tree topology.
*(`revision/ureC_tree/ureC.treefile`; pruned versions for comparison in `revision/ureC_tree/`)*
ML UreC gene tree (n = 46 ureC-encoding MAGs) with ultrafast bootstrap support. Hero-lineage tips marked. The gene tree is topologically incongruent with the GTDB-Tk bac120 species tree (normalised Robinson–Foulds = 0.58; SH test rejecting species tree against the *ureC* alignment, p < 0.001), indicating lineage-specific evolutionary dynamics on urease. Within each hero lineage the topology is consistent with vertical inheritance.

---

## Supplementary tables

### Table S1 | MAG statistics and taxonomy.
(`MAGs_FASTA_files/ace_samples_list.csv`, `research/extra/novel_species/MIGS_lite_Sphingobacterium.csv`, `research/revision/hero_cluster_audit.csv`.) Per-MAG waste source, genome size, N50, scaffold/contig counts, GC, CDS/tRNA/rRNA counts, GTDB classification, closest-reference ANI; expanded MIGS-lite view for the six *Sphingobacterium* MAGs; and per-hero cluster-contiguity audit.

### Table S2 | Trait-module keyword dictionaries, raw counts, and permutation statistics.
(`research/extra/gene_category_counts.csv`, `research/extra/gene_category_per1k_CDS.csv`, `research/revision/Hero_vs_Rest_permutation_stats.csv`.) Keyword definitions for the six trait modules; raw and per-10³-CDS normalised hit counts for every MAG × subcategory; observed fold change, bootstrap 95 % CI, one-sided permutation p-value (10,000 iterations) and Mann–Whitney U p-value with Benjamini–Hochberg q-values at the module level.

### Table S3 | Mobile-element and GC analysis of the hero *ure* cluster.
(`research/extra/HGT_ureCah_cluster.csv`, `research/revision/hero_cluster_audit.csv`, `research/revision/hero_gene_contig_distribution.csv`.) Hero-MAG coordinates of the main *ure* contig, detected mobile-element counts in ±15 kb, regional GC, and gene-to-contig distribution for every *ure* and *cah* copy in every hero MAG.

### Table S4 | Whole-genome and amino-acid identity matrices.
(`research/extra/skani_full_matrix.csv`, `research/extra/skani_triangle.tsv.af`, `research/extra/novel_species/AAI_S13_S16_vs_Sphingobacterium.csv`.) Pairwise skani ANI (and alignment fraction) matrix for all 111 MAGs; reciprocal-best mmseqs2 AAI of S13/S16 against the remaining *Sphingobacterium* MAGs.

### Table S5 | DRAM distillate outputs.
(`/data/pangenome_work/dram_output/distillate/product.tsv`, `metabolism_summary.xlsx`, `genome_stats.tsv`.) DRAM v1.5 `distill` output covering 98 modules across 111 MAGs, including Excel workbook and interactive HTML product report.

### Table S6 | dbCAN CAZyme profile.
(`research/revision/dbCAN_class_counts.csv`, `dbCAN_class_per1k_CDS.csv`, `dbCAN_family_counts.csv`, `dbCAN_hero_vs_rest_class.csv`, `dbCAN_hero_vs_rest_family.csv`.) CAZy class- and family-level counts and hero-vs-rest comparison derived from the DRAM dbCAN annotation column.

### Table S7 | *ureC* gene tree versus species tree congruence.
(`research/revision/ureC_tree/RF_result.txt`, `SHtest.iqtree`, `ureC.treefile`, `species_pruned.tre`, `ureC_pruned.tre`, `ureC_samples.csv`.) Robinson–Foulds distance, Shimodaira–Hasegawa and approximately-unbiased test p-values, and the ML gene and species trees used for comparison.

### Table S8 | Novelty screen and MIGS-lite for candidate novel *Sphingobacterium* species.
(`research/extra/novelty_ANI_screen.csv`, `research/extra/novel_species/MIGS_lite_Sphingobacterium.csv`.) Per-MAG ANI-based novelty classification for the full panel and MIGS-lite parameters for S13 and S16.
