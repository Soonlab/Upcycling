# Figure, Table and Supplementary Legends (revised)

All figures are provided as 300-dpi PNG and vector PDF files at `/data/data/Upcycling/research/` (main figures) and `/data/data/Upcycling/research/extra/` and `/data/data/Upcycling/research/revision/` (supplementary figures). Sample codes: **C** cattle, **M** swine, **S** sheep, **V** poultry. The six "MICP-complete" MAGs (C22, M1, S13, S16, S23, S26) are consistently highlighted in red/bold across all figures.

---

## Main figures

### Figure 1 | Phylogenomic distribution of MICP-related genes across 111 livestock-waste MAGs.
*(`Fig1_Phylogeny_MICP_heatmap.png`, `.pdf`)*

Maximum-likelihood phylogenomic tree of 111 MAGs inferred from the GTDB-Tk bac120 concatenated alignment (IQ-TREE v2.3, ModelFinder, 1,000 ultrafast bootstrap replicates; midpoint-rooted for display). The left colour strip encodes GTDB-Tk genus assignment (legend: top ten genera with sample counts). MAG identifiers and, where available, GTDB species annotations follow the colour strip; MICP-complete MAGs are highlighted in red bold. The right-hand heat map shows presence (dark green) or absence (white) of the eight MICP-associated genes *ureA, ureB, ureC, ureD, ureE, ureF, ureG* and *cah*. The complete 8/8 module is restricted to four *Sphingobacterium* MAGs (C22, S13, S16, S23) and two *Pseudomonas*\_E MAGs (M1, S26); these taxa are distributed across two, non-monophyletic, functionally convergent lineages. Scale bar, substitutions per site.

### Figure 2 | MICP module completeness across genera and gene-level prevalence in MICP-complete MAGs.
*(`Fig2_MICP_completeness_by_genus.png`, `.pdf`)*

**(a)** MICP module score (sum of *ureA–G* + *cah*; maximum = 8) by GTDB-Tk genus, sorted by mean. Boxes show median and inter-quartile range; individual MAGs are overlaid as grey points; MICP-complete MAGs are superimposed as enlarged red circles. **(b)** Per-gene prevalence of the eight MICP genes in the MICP-complete set (red, n = 6) versus the remaining 105 MAGs (grey). MICP-complete MAGs maintain 100 % prevalence of every *ure* gene and *cah*, whereas the non-MICP-complete background retains the same genes at 55–85 %.

### Figure 3 | Gene-order conservation of the *ureABCDEFG* operon across MICP-complete MAGs.
*(`Fig3_ureCah_cluster_synteny.png`, `.pdf`)*

Synteny diagrams of the densest *ure*-containing contig in each MICP-complete Bakta annotation. Each horizontal track represents the ≤ 40 kb window with the greatest number of *ure* genes for one MAG (C22, M1, S13, S16, S23, S26). Arrows are drawn to scale; arrow direction encodes strand; colours encode gene identity (legend: *ureA* blue, *ureB* green, *ureC* red, *ureD* purple, *ureE* brown, *ureF* pink, *ureG* cyan, *cah* orange, other grey). All 7 *ure* genes are recovered on a single contig within a 5.9–28.6 kb window in the three best-assembled MICP-complete MAGs (M1, S13, S16); in C22, S23 and S26 a subset of *ure* genes is retained on the main contig (4–5 / 7) with the remainder on secondary contigs, consistent with assembly fragmentation (Table 1). No transposase, integrase, prophage or relaxase gene is detected in ±15 kb flanking any MICP-complete cluster (Supplementary Table S3).

### Figure 4 | DRAM metabolic-module heat map and MICP-complete-vs-rest comparison of MICP-critical modules.
*(`extra/Fig4_DRAM_metabolism_heatmap.png`, `.pdf`; `extra/Fig4b_DRAM_HeroVsRest.png`, `.pdf`)*

**(a)** Module-completeness heat map (DRAM v1.5 `distill/product.tsv`) for 111 MAGs × 35 curated KEGG/custom modules spanning central carbon metabolism, nitrogen metabolism, stress response, vitamin (cobalamin / B12) biosynthesis, and CAZyme categories. Row labels give MAG identifier and GTDB genus; MICP-complete rows are highlighted in red bold. Colour intensity encodes fractional module completeness (0 to 1). Housekeeping modules (glycolysis, TCA, fatty-acid biosynthesis) are comparably complete in MICP-complete and non-MICP-complete MAGs. **(b)** Mean module completeness of MICP-critical modules (urease, carbonic anhydrase, nitrogen metabolism, Na⁺/H⁺ antiport, cobalamin biosynthesis, CAZymes/carbohydrate metabolism) in MICP-complete MAGs (red, n = 6) versus the remaining MAGs (grey, n = 105).

### Figure 5 | Novel-species delineation of S13 and S16 within *Sphingobacterium*.
*(`extra/Fig_T2d_Novelty_overview.png`, `.pdf`; `extra/novel_species/Fig_NovelSp_AAI.png`, `.pdf`; `revision/Fig_ext_Sphingo_ANI.png`, `.pdf`)*

**(a)** Novelty screen within the 111-MAG panel. Each MAG is one point; the y-axis is the GTDB closest-reference ANI, and the dashed red line marks the 95 % species cutoff. MICP-complete MAGs are enlarged and red-labelled. Twenty-one of 111 MAGs fall below 95 %, including S13 and S16 (no species-level ANI available). **(b)** Amino-acid identity of S13 (left) and S16 (right) versus the other five *Sphingobacterium* MAGs, computed from reciprocal-best mmseqs2 easy-search hits. Dashed vertical lines mark the 95 % (species) and 70 % (genus) AAI thresholds. Maximum AAI = 93.15 % (S13 vs V3) and 93.49 % (S16 vs S23); both well below the species threshold. **(c)** External validation — pairwise skani ANI of the six study *Sphingobacterium* MAGs against the 63 RefSeq reference *Sphingobacterium* genomes (NCBI Datasets v16). Each panel ranks RefSeq genomes by descending ANI; the dashed red line marks the 95 % species cutoff, and bar colour encodes classification (red = species match, blue = sub-species). S13 (max ANI 94.57 % to *S. detergens*) and S16 (max ANI 93.85 % to *S. multivorum*) remain below the cutoff against the full public genus catalogue, while S23 (98.96 %) and C22 (99.16 %) validate the pipeline by reproducing their GTDB-Tk species assignments.

### Figure 6 | Permutation-tested MICP-complete-vs-rest enrichment of trait modules.
*(`revision/Fig_Permutation_forest.png`, `.pdf`)*

Forest-style plot of fold-change enrichment (log-scale x-axis) for trait subcategories with fold change > 1, ranked by magnitude. Points are fold change; error bars are 2,000-iteration bootstrap 95 % CI. Point colour encodes significance level after Benjamini–Hochberg FDR correction on 10,000-iteration one-sided permutation p-values: **red** q < 0.05, **orange** q < 0.10, **grey** n.s. Modules significantly enriched in the MICP-complete lineage include Mrp complex (FC 10.85), CBM (9.78), oxidative-stress defence (4.76), glycoside hydrolase (4.66), quorum sensing (2.13) and Na⁺/H⁺ antiporter (2.30). Complete statistics in Supplementary Table S2 (`Hero_vs_Rest_permutation_stats.csv`).

### Figure 7 | Pan-genome and trait-module PCoA by waste source.
*(`revision/Fig_PCoA_source_genus.png`, `.pdf`; `revision/Fig_PCoA_trait_source.png`, `.pdf`)*

**(a)** Jaccard-based Principal Coordinates Analysis of Panaroo gene presence/absence (filtered to 9,668 genes at 5–95 % prevalence). Points are coloured by waste source (red cattle, blue swine, green sheep, orange poultry); MICP-complete MAGs are highlighted with an open red ring and label. PERMANOVA identifies GTDB genus as the dominant structuring variable (pseudo-F = 8.21, p = 0.001) but only weak differentiation by waste source (pseudo-F = 1.25, p = 0.11). **(b)** Z-score–standardised Euclidean PCoA of per-10³-CDS-normalised counts across 38 trait-module subcategories (z-scoring prevents a single high-abundance module from dominating the first principal coordinate). Trait distribution is significantly structured by waste source (PERMANOVA pseudo-F = 2.71, p = 0.001). MICP-complete MAGs consistently project onto the "functionally extreme" tail of the ordination regardless of origin, supporting their interpretation as environmentally specialised lineages rather than source-specific outliers. Pairwise PERMANOVA and post-hoc q-values are listed in Supplementary Table S9 (`PCoA_pairwise_PERMANOVA.csv`).

---

## Supplementary figures

### Figure S1 | Biofilm / EPS gene modules (genus-aggregated).
*(`Figure_S1.png`, `.pdf`)*
Genus-aggregated heat map (mean gene hits per 10³ CDS, log₁₀-transformed for readability) of ten EPS/biofilm subcategories (pel, psl, pga, cellulose/bcs, alginate/alg, curli/csg, colanic/wca, capsule/wz\*, adhesin/autotransporter, quorum sensing) across the nine GTDB-Tk genera with ≥ 2 MAGs. Genera containing MICP-complete lineage members (*Sphingobacterium*, *Pseudomonas*\_E) are placed at the top and rendered in red bold.

### Figure S2 | Ammonia-handling and nitrogen-assimilation modules (genus-aggregated).
*(`Figure_S2.png`, `.pdf`)*
Genus-aggregated heat map (log₁₀-scaled) covering urea transport (urtABC, urea transporter), GS–GOGAT (glnA, gltBD), GDH (gdhA), AmtB, PII regulators (glnB/K) and nitrate/nitrite reductase pathways.

### Figure S3 | Alkaline and osmotic-stress tolerance modules (genus-aggregated).
*(`Figure_S3.png`, `.pdf`)*
Genus-aggregated heat map (log₁₀-scaled) of Na⁺/H⁺ antiporter (nhaA–C), Mrp multi-subunit antiporter complex, Kdp / Trk K⁺ uptake, compatible-solute biosynthesis (opu, pro, betA/B, ect), chaperones (groEL/ES, dnaK/J, clpB) and oxidative-stress defences (katA/B/G, sodA/B/C, ahpCF).

### Figure S4 | CAZyme profile: keyword proxy and DRAM/dbCAN validation.
*(`Figure_S4a_keyword.png`, `Figure_S4b_dbCAN_classes.png`, `Figure_S4c_dbCAN_families.png`)*
Keyword-based CAZy proxy heat map and independent DRAM/dbCAN-derived class-level and family-level profiles. Both approaches identify MICP-complete-lineage enrichment of glycoside hydrolases and carbohydrate-binding modules.

### Figure S5 | Heavy-metal and antibiotic-resistance gene modules (genus-aggregated).
*(`Figure_S5.png`, `.pdf`)*
Genus-aggregated heat map (log₁₀-scaled) of heavy-metal efflux (czc, cop, cus, ars, mer, znt), metallothioneins, β-lactamases, broad MDR efflux (acr, mex, MATE) and aminoglycoside/tetracycline/macrolide determinants.

### Figure S6 | Pairwise ANI within MICP-complete lineages.
*(`extra/Fig_T2d_MICP-completeANI.png`, `.pdf`)*
skani whole-genome ANI heat map among the six MICP-complete MAGs.

### Figure S7 | UreC gene tree topology.
*(`revision/ureC_tree/ureC.treefile`; pruned versions for comparison in `revision/ureC_tree/`)*
ML UreC gene tree (n = 46 ureC-encoding MAGs) with ultrafast bootstrap support. MICP-complete tips marked. The gene tree is topologically incongruent with the GTDB-Tk bac120 species tree (normalised Robinson–Foulds = 0.58; SH test rejecting species tree against the *ureC* alignment, p < 0.001), indicating lineage-specific evolutionary dynamics on urease. Within each MICP-complete lineage the topology is consistent with vertical inheritance.

---

## Supplementary tables

### Table S1 | MAG statistics and taxonomy.
(`MAGs_FASTA_files/ace_samples_list.csv`, `research/extra/novel_species/MIGS_lite_Sphingobacterium.csv`, `research/revision/MICP-complete_cluster_audit.csv`.) Per-MAG waste source, genome size, N50, scaffold/contig counts, GC, CDS/tRNA/rRNA counts, GTDB classification, closest-reference ANI; expanded MIGS-lite view for the six *Sphingobacterium* MAGs; and per-MICP-complete cluster-contiguity audit.

### Table S2 | Trait-module keyword dictionaries, raw counts, and permutation statistics.
(`research/extra/gene_category_counts.csv`, `research/extra/gene_category_per1k_CDS.csv`, `research/revision/Hero_vs_Rest_permutation_stats.csv`.) Keyword definitions for the six trait modules; raw and per-10³-CDS normalised hit counts for every MAG × subcategory; observed fold change, bootstrap 95 % CI, one-sided permutation p-value (10,000 iterations) and Mann–Whitney U p-value with Benjamini–Hochberg q-values at the module level.

### Table S3 | Mobile-element and GC analysis of the MICP-complete *ure* cluster.
(`research/extra/HGT_ureCah_cluster.csv`, `research/revision/MICP-complete_cluster_audit.csv`, `research/revision/MICP-complete_gene_contig_distribution.csv`.) MICP-complete-MAG coordinates of the main *ure* contig, detected mobile-element counts in ±15 kb, regional GC, and gene-to-contig distribution for every *ure* and *cah* copy in every MICP-complete MAG.

### Table S4 | Whole-genome and amino-acid identity matrices.
(`research/extra/skani_full_matrix.csv`, `research/extra/skani_triangle.tsv.af`, `research/extra/novel_species/AAI_S13_S16_vs_Sphingobacterium.csv`.) Pairwise skani ANI (and alignment fraction) matrix for all 111 MAGs; reciprocal-best mmseqs2 AAI of S13/S16 against the remaining *Sphingobacterium* MAGs.

### Table S5 | DRAM distillate outputs.
(`/data/pangenome_work/dram_output/distillate/product.tsv`, `metabolism_summary.xlsx`, `genome_stats.tsv`.) DRAM v1.5 `distill` output covering 98 modules across 111 MAGs, including Excel workbook and interactive HTML product report.

### Table S6 | dbCAN CAZyme profile.
(`research/revision/dbCAN_class_counts.csv`, `dbCAN_class_per1k_CDS.csv`, `dbCAN_family_counts.csv`, `dbCAN_MICP-complete_vs_rest_class.csv`, `dbCAN_MICP-complete_vs_rest_family.csv`.) CAZy class- and family-level counts and MICP-complete-vs-rest comparison derived from the DRAM dbCAN annotation column.

### Table S7 | *ureC* gene tree versus species tree congruence.
(`research/revision/ureC_tree/RF_result.txt`, `SHtest.iqtree`, `ureC.treefile`, `species_pruned.tre`, `ureC_pruned.tre`, `ureC_samples.csv`.) Robinson–Foulds distance, Shimodaira–Hasegawa and approximately-unbiased test p-values, and the ML gene and species trees used for comparison.

### Table S8 | Novelty screen and MIGS-lite for candidate novel *Sphingobacterium* species.
(`research/extra/novelty_ANI_screen.csv`, `research/extra/novel_species/MIGS_lite_Sphingobacterium.csv`.) Per-MAG ANI-based novelty classification for the full panel and MIGS-lite parameters for S13 and S16.

### Table S9 | Pan-genome and trait-module PCoA with PERMANOVA.
(`research/revision/PCoA_panaroo_coords.csv`, `PCoA_PERMANOVA.csv`, `PCoA_pairwise_PERMANOVA.csv`.) PC1–PC3 coordinates, waste-source and genus metadata per MAG; global and pairwise post-hoc PERMANOVA pseudo-F, p-values and BH-FDR q-values for the pan-genome and trait-module ordinations shown in Fig. 7.

### Table S10 | External *Sphingobacterium* ANI comparison.
(`research/revision/ANI_ext_sphingo_matrix.csv`, `ANI_ext_sphingo_novelty.csv`.) Full pairwise skani ANI matrix for the six study *Sphingobacterium* MAGs against 63 RefSeq reference *Sphingobacterium* genomes, and the per-MAG nearest-reference summary used to confirm novelty of S13 and S16.
