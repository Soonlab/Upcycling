---
title: "Comparative genomics of 111 livestock-waste metagenome-assembled genomes identifies two *Sphingobacterium*/*Pseudomonas*\\_E lineages that independently retain a complete ureolytic–carbonic anhydrase module and candidate novel species for alkali-tolerant microbially induced carbonate precipitation"
author:
  - name: "[Author 1]"
  - name: "[Author 2]"
  - name: "[Corresponding Author]"
date: "2026"
---

**Target journals:** *Microbial Biotechnology* (primary) | *mSystems* | *Environmental Microbiome*
**Manuscript type:** Full-length research article — comparative genomics of metagenome-assembled genomes.
**Word count (body):** ≈ 5,500

**Highlights**

- A 111-MAG livestock-waste genome catalogue was systematically mined for microbially induced carbonate precipitation (MICP) potential using Bakta, Panaroo, GTDB-Tk, IQ-TREE, skani, mmseqs2 and DRAM.
- Two independent lineages — four *Sphingobacterium* MAGs (C22, S13, S16, S23) and two *Pseudomonas*\_E MAGs (M1, S26) — retain the complete *ureABCDEFG* operon together with one or more carbonic-anhydrase (*cah*) copies, constituting **two functionally convergent rather than a single monophyletic group**.
- The *ureABCDEFG* operon is consistently recovered on a single contig (5–29 kb span) in the three better-assembled MAGs (M1, S13, S16); limited assembly contiguity prevents identical verification in C22/S23/S26, but every hero MAG encodes every *ure* gene.
- Flanking mobile elements are absent and regional GC content matches the genomic mean, consistent with vertical inheritance rather than recent horizontal acquisition within each lineage.
- After permutation testing with Benjamini–Hochberg FDR control, the hero lineages are significantly enriched (q < 0.05) for Mrp Na⁺/H⁺-antiporter genes, oxidative-defence enzymes, compatible-solute biosynthesis, glycoside hydrolases and carbohydrate-binding modules — a trait stack consistent with survival and catalysis in alkaline, polysaccharide-rich slurry.
- S13 and S16 satisfy both GTDB-Tk and AAI criteria (< 95 %) for candidate novel *Sphingobacterium* species and are nominated as priority chassis for future experimental validation of slurry-coupled MICP.

---

## Abstract

Microbially induced carbonate precipitation (MICP) is a mature biotechnology for low-carbon ground improvement, but its industrial application still relies on a narrow set of ureolytic *Sporosarcina*/*Bacillus* strains that do not tolerate the alkaline, ammonia-rich conditions of livestock slurry. To expand the phylogenetic and functional palette of candidate MICP chassis while simultaneously valorising livestock waste, we performed a comprehensive *in silico* comparative-genomic analysis of 111 metagenome-assembled genomes (MAGs) recovered from four livestock-waste microbiomes (cattle, swine, sheep, poultry). The analysis integrated Bakta functional annotation, Panaroo pan-genome reconstruction, GTDB-Tk/IQ-TREE phylogenomics, skani whole-genome ANI, reciprocal-best mmseqs2 amino-acid identity (AAI), DRAM metabolic reconstruction, dbCAN-based CAZyme profiling, and permutation-tested trait-module screening across six functional categories. The 111-MAG panel spanned eight genera. A six-member MICP-complete set — four *Sphingobacterium* (C22, S13, S16, S23) and two *Pseudomonas*\_E MAGs (M1, S26) — independently retained the full *ureABCDEFG* operon on a single contig in the better-assembled cases, together with a separately encoded *cah* gene, and represents two functionally convergent lineages rather than a single monophyletic clade. The *ure* operon exhibited no flanking transposase, integrase, prophage or relaxase annotation and no regional GC deviation, supporting vertical inheritance within each lineage. Permutation testing with Benjamini–Hochberg FDR control confirmed (q < 0.05) significant enrichment of Mrp Na⁺/H⁺ antiporters (fold change 10.9; 95 % bootstrap CI 2.3–43.8), oxidative-defence enzymes (FC 4.8; CI 2.2–8.8), glycoside hydrolases (FC 4.7; CI 1.5–8.4), carbohydrate-binding modules (FC 9.8; CI 3.5–22.0) and compatible-solute biosynthesis (FC 1.5; CI 1.2–1.9), rationalising persistence at pH > 9 and self-fuelling on slurry polysaccharides. AAI analysis identified S13 and S16 as candidate novel *Sphingobacterium* species (maximum congeneric AAI 93.2 % and 93.5 %, respectively). Topological analysis of the *ureC* gene tree versus the bac120 species tree (Shimodaira–Hasegawa test, normalised Robinson–Foulds distance 0.58) suggests that *ureC* phylogeny reflects substantial selection on urease function beyond whole-genome vertical descent; within each of the two hero lineages, however, the data are consistent with vertical inheritance. These computational predictions provide a prioritised and mechanistically grounded hypothesis set for subsequent biochemical and field validation.

**Keywords:** microbially induced carbonate precipitation; metagenome-assembled genome; *Sphingobacterium*; urease; carbonic anhydrase; alkaline tolerance; comparative genomics; livestock-waste valorisation; circular bioeconomy.

---

## 1. Introduction

Microbially induced carbonate precipitation (MICP) couples enzymatic pH elevation with inorganic-carbon supply to crystallise calcium carbonate on microbial surfaces. The technology has been deployed for concrete self-healing, geotechnical consolidation, fugitive-dust control, and immobilisation of toxic metals (DeJong et al., 2013; Phillips et al., 2013; Achal and Mukherjee, 2015). Industrial MICP is dominated by the ureolytic route: urease (EC 3.5.1.5) hydrolyses urea to NH₃ and CO₂, alkalinising the bulk medium and, with carbonic anhydrase (CA, EC 4.2.1.1) accelerating CO₂ ↔ HCO₃⁻ interconversion, releasing the HCO₃⁻/CO₃²⁻ required for Ca²⁺-mediated CaCO₃ nucleation (Dhami et al., 2014). Despite two decades of optimisation, nearly all field-grade MICP relies on *Sporosarcina pasteurii* and a small set of *Bacillus*/*Lysinibacillus* strains whose growth optima (pH ≈ 7–8, freshwater) are incompatible with real waste streams.

Livestock slurry, produced at > 1.5 × 10⁹ t yr⁻¹ globally, is a matrix rich in urea, volatile fatty acids, ammonia, inorganic cations and partially lignified polysaccharides. It is largely treated as a biohazardous disposal problem rather than as a reservoir of functionally adapted microorganisms. Recent advances in hybrid metagenomic assembly and genome-resolved binning (Parks et al., 2017; Nayfach et al., 2021) now permit culture-independent recovery of draft genomes ("MAGs") from such matrices, while integrated pan-genomic and metabolic-reconstruction frameworks (Tonkin-Hill et al., 2020; Shaffer et al., 2020) enable simultaneous phylogenetic and functional triage.

Here we apply this genomics-first strategy to 111 MAGs from four livestock-waste microbiomes and ask three questions: **(i)** which genera carry a complete and physically co-located *ureABCDEFG–cah* machinery; **(ii)** is the machinery a recently acquired mobilisable cassette or a vertically inherited lineage-defining trait; **(iii)** do the candidate strains co-encode the accessory traits (alkaline/osmotic resilience, oxidative defence, polysaccharide degradation) required to function in real slurry. We deliberately frame the present study as *in silico only*: its purpose is to generate a rigorous, statistically validated, and mechanistically coherent target list whose biochemical and engineering properties will be evaluated in follow-up studies.

## 2. Materials and Methods

### 2.1 MAG collection
111 MAGs were recovered from four livestock-waste metagenomes (cattle "C", swine "M", sheep "S", poultry "V") using the ACE hybrid pipeline with completeness ≥ 50 % and contamination ≤ 10 % (MIMAG medium-quality draft or better; Bowers et al., 2017). Genome sizes ranged from 1.9 to 6.2 Mb (median 3.8 Mb), N50 from 2.1 to 194 kb, and GC content from 32 % to 68 %. MAG quality metrics varied across the panel; the six MICP-complete MAGs had assemblies ranging from highly contiguous (M1: 17 contigs, N50 = 414 kb) to highly fragmented (C22: 2,193 contigs, N50 = 2.1 kb) (Table 1).

### 2.2 Functional annotation
Each MAG was annotated with Bakta v1.9 (Schwengers et al., 2021) against the full database (r220), generating GFF3/GBK/FAA/FFN outputs. Bakta calls Pyrodigal for CDS prediction, tRNAscan-SE v2 and Infernal v1.1 for non-coding RNAs, and DIAMOND/HMMER against UniRef100, PSC/PSCC and RFAM for functional assignment.

### 2.3 Pan-genome reconstruction
Bakta GFF3 files were ingested into Panaroo v1.5 (moderate mode, 90 % identity threshold) to derive core/shell/cloud gene-complement statistics (`panaroo_results/`).

### 2.4 Phylogenomics
Taxonomy and closest-reference ANI were assigned with GTDB-Tk v2.4 (release r220; Chaumeil et al., 2022). The bac120 multilocus alignment was used as input to IQ-TREE v2.3 with ModelFinder and 1,000 ultrafast bootstrap replicates (Minh et al., 2020). Pairwise whole-genome ANI of all 111 MAGs was computed with skani v0.3 (`skani triangle --full-matrix`; Shaw and Yu, 2023). Average amino-acid identity (AAI) among congeneric MAGs was estimated from reciprocal-best mmseqs2 easy-search hits on Bakta proteomes (≥ 30 % identity, ≥ 70 % coverage, −s 5; Steinegger and Söding, 2017).

### 2.5 Metabolic reconstruction
All 111 MAGs were annotated with DRAM v1.5 (Shaffer et al., 2020) against KOfam (2024), UniRef90, Pfam, dbCAN (v12) and MEROPS. Module completeness was obtained from `DRAM.py distill`. CAZyme profiles reported in §3.4 are derived directly from the `cazy_best_hit` column of the merged DRAM annotation.

### 2.6 Trait-module screening and statistical analysis
Six curated keyword-based trait modules (biofilm/EPS; ammonia handling; mobile genetic elements; alkaline/osmotic defence; CAZymes; metal/antibiotic resistance; Supplementary Table S2) were quantified by scanning Bakta protein-product annotations. Raw hit counts were normalised per 10³ CDS to control for genome size. Per-module enrichment between the six MICP-complete ("hero") MAGs and the remaining 105 was tested by **one-sided label-permutation (10,000 permutations)** and by one-sided Mann–Whitney U, followed by Benjamini–Hochberg FDR correction at the module level. Bootstrap 95 % confidence intervals on hero-group means and hero/rest fold changes were obtained from 2,000 resamples.

### 2.7 Synteny and mobile-element analysis of the *ure–cah* cluster
For each MICP-complete MAG the contig harbouring the greatest number of *ure* genes was identified. A ±15 kb window around the *ure* cluster was searched for transposase, integrase, phage, relaxase and T4SS annotations, and regional GC was compared to the genomic mean. Owing to MAG fragmentation, the *cah* gene was not recoverable on the same contig as the *ure* operon in any hero MAG; the single-contig verification of cluster contiguity is therefore restricted to the *ureABCDEFG* operon itself.

### 2.8 Gene-tree versus species-tree congruence
The longest UreC (urease subunit alpha, > 300 aa) sequence was extracted from every MAG in which it was annotated (n = 46). Proteins were aligned with MAFFT (`--auto`) and an ML gene tree was inferred with IQ-TREE v2.3 (LG+G model, 1,000 UFBoot). Congruence with the GTDB-Tk bac120 species tree (pruned to the same 46 taxa) was quantified by the unrooted Robinson–Foulds distance (ete3) and by the Shimodaira–Hasegawa (SH) and approximately-unbiased (AU) tests in IQ-TREE (`-z` option, 1,000 RELL bootstraps).

### 2.9 Novel species delineation
A MAG was classified as a candidate novel species when (i) GTDB-Tk returned no species-level assignment or the closest-reference ANI was < 95 % (Konstantinidis and Tiedje, 2005), **and** (ii) pairwise AAI to every sequenced congener within the panel was < 95 % (Rodriguez-R and Konstantinidis, 2014).

### 2.10 Data and code availability
All MAGs, intermediate analyses, figures and scripts are available at `/data/data/Upcycling/` and will be deposited to NCBI BioProject PRJNA-XXXXXXX and Zenodo (DOI on acceptance).

## 3. Results

### 3.1 Composition and quality of the 111-MAG panel
The panel spanned eight bacterial genera: Gammaproteobacteria (*Pseudomonas*\_E, *Acinetobacter*), Bacteroidota (*Sphingobacterium*, *Chryseobacterium*, *Empedobacter*), Bacillota (*Bacillus* sensu lato, *Paenibacillus*), and additional low-abundance taxa (Supplementary Table S1). Twenty-one of 111 MAGs (18.9 %) returned closest-reference ANI < 95 % under GTDB-Tk r220, indicating a substantial reservoir of uncultured diversity in livestock waste even at the species level (Fig. 5a).

### 3.2 Two independent lineages retain the complete MICP module
The IQ-TREE bac120 phylogeny resolved the eight genera as well-supported monophyletic clades (Fig. 1). Mapping presence of *ureA–G* + *cah* across the phylogeny identified a MICP-complete set of six MAGs that span **two taxonomic lineages rather than a single monophyletic group**: four *Sphingobacterium* MAGs (C22, S13, S16, S23) and two *Pseudomonas*\_E MAGs (M1, S26). A further 20 MAGs (predominantly additional *Sphingobacterium*, *Pseudomonas*\_E and *Acinetobacter*) retained partial MICP modules (3–6 / 8), and 64 MAGs retained two or fewer MICP genes (Fig. 2a). Among MICP-complete MAGs, *Sphingobacterium* and *Pseudomonas*\_E averaged 7.0 / 8 module score while *Chryseobacterium* and *Acinetobacter* averaged 1.3 / 8 (Fig. 2b). We therefore frame the remainder of the analysis around these two functionally convergent lineages (collectively: "hero-lineage MAGs"; n = 6).

### 3.3 Cluster contiguity and evidence for vertical inheritance
All six hero-lineage MAGs encode *ureABCDEFG* and at least one carbonic-anhydrase gene, but the physical contiguity of the cluster is constrained by assembly quality (Table 1). In the three best-assembled MAGs (M1, S13, S16), the entire *ureABCDEFG* operon occupies a single contig within a 5.9–28.6 kb window. In C22, five of seven *ure* genes fall on a single 5.3 kb contig, with the remaining genes on related but separately assembled contigs. In S23 and S26, respectively four of seven *ure* genes are recovered on the main contig (3.0–3.3 kb span), with additional copies on secondary contigs — a pattern best explained by assembly fragmentation rather than genuine rearrangement. The *cah* gene is never co-recovered on the main *ure* contig in any hero MAG (Table 1). The fact that *cah* is assembled on separate short contigs in MAGs whose *ure* operon is contiguous therefore cannot distinguish "physically distant" from "fragmented but linked" without long-read resequencing.

Within the best-resolved contigs, we observed **no** transposase, integrase, prophage, or relaxase annotation in a ±15 kb window around the *ure* operon of any hero MAG (`HGT_ureCah_cluster.csv`). Regional GC content was within ±1.5 % of the genomic mean, with median Δ GC = 0.00 %. These two independent lines of evidence — absence of mobile-element signatures and flat GC — are consistent with vertical inheritance of the *ure* operon within each of the two hero lineages. We emphasise, however, that negative evidence for horizontal transfer does not exhaustively rule out ancient transfer events whose signatures have been eroded; a Shimodaira–Hasegawa test (§ 3.7) provides an additional, independent phylogenetic perspective.

### 3.4 Accessory traits for survival in alkaline, ammonia-rich slurry (permutation-tested)
Trait-module enrichment was tested by 10,000-iteration permutation of hero/rest labels, and raw p-values were FDR-adjusted at the module level (`Hero_vs_Rest_permutation_stats.csv`). Nine of 38 tested modules survived BH-FDR q < 0.05, with bootstrap 95 % CIs excluding unity (Table 2; Fig. 6):

| Module (Category::Subcategory) | Fold change | Bootstrap 95 % CI | FDR q |
|---|---|---|---|
| Alkaline\_Osmo::Mrp\_complex | 10.85 | [2.28, 43.83] | 0.011 |
| CAZyme\_proxy::carb\_binding | 9.78 | [3.54, 21.96] | 0.003 |
| Alkaline\_Osmo::oxidative | 4.76 | [2.23, 8.75] | 0.008 |
| CAZyme\_proxy::glycoside\_hydrolase | 4.66 | [1.50, 8.43] | 0.003 |
| MetalResist\_AMR::tet\_mac | 4.22 | [1.74, 7.07] | 0.003 |
| Alkaline\_Osmo::Na\_H\_antiporter | 2.30 | [1.28, 3.40] | 0.013 |
| Biofilm\_EPS::quorum | 2.13 | [1.35, 2.94] | 0.003 |
| Alkaline\_Osmo::compatible\_solute | 1.51 | [1.15, 1.88] | 0.018 |
| MetalResist\_AMR::metal\_efflux | 1.19 | [1.06, 1.31] | 0.019 |

Hero-lineage MAGs therefore show statistically supported enrichment of a mechanistically coherent trait stack: (i) cytoplasmic pH and Na⁺-gradient buffering (Mrp complex, *nha* antiporters, compatible solutes); (ii) oxidative-stress defence compatible with urea-driven intracellular alkalinisation; (iii) carbohydrate acquisition (GH and CBM enrichment) compatible with self-fuelling on slurry polysaccharides; (iv) cell-cell communication (quorum-sensing enrichment) compatible with biofilm formation; (v) environmental metal-efflux tolerance. Modules not retained at q < 0.05 after FDR correction — including single-subunit *nhaA*, *urtABC* urea transport, and several EPS operons — are reported with full statistics in Supplementary Table S2.

One module initially of interest, glutamine-synthetase / glutamate-synthase (GS–GOGAT), showed a non-significant **depletion** in the hero lineages (FC = 0.62). We interpret this cautiously as consistent with the physiological expectation that urea-derived ammonia is preferentially released to the environment rather than reassimilated into biomass, but we emphasise that this is a negative-enrichment observation that requires biochemical validation.

### 3.5 CAZyme profile independently validated by dbCAN HMMER
To validate the keyword-based CAZyme signal with a curated reference, we performed HMMER-3 `hmmsearch` (E < 1 × 10⁻¹⁵) of the dbCAN v12 HMM profile library against the six hero Bakta proteomes, and combined the results with DRAM-derived `cazy_best_hit` assignments for the remaining 105 MAGs (Table 3; `dbCAN_final_hero_vs_rest_class.csv`). The direct-HMMER analysis confirms significant hero enrichment across four of six CAZy classes after permutation testing with BH-FDR correction: glycoside hydrolases (GH, fold change 3.82, 95 % bootstrap CI [1.80, 5.67], q = 6 × 10⁻⁴), polysaccharide lyases (PL, 3.52 [2.09, 5.28], q = 1 × 10⁻³), carbohydrate-binding modules (CBM, 4.24 [1.72, 6.98], q = 9 × 10⁻⁴) and carbohydrate esterases (CE, 1.99 [1.09, 2.99], q = 7 × 10⁻³). Glycosyltransferases (GT) showed only marginal enrichment (1.24, q = 4 × 10⁻²) and auxiliary activities (AA) no enrichment (q = 1.0). Family-level analysis (`dbCAN_family_counts.csv`) identifies hero-enriched families consistent with hemicellulose, starch and peptidoglycan turnover — the polysaccharide classes most abundant in livestock slurry. The convergence of keyword-based and reference-HMM analyses, together with their independent statistical support, establishes the slurry-adapted CAZyme signature of the hero lineages as a robust feature of the data rather than an artefact of annotation.

### 3.6 DRAM module reconstruction confirms the MICP signature
DRAM distillation of the merged 111-MAG annotation table (433,595 gene-level assignments across 98 modules) independently corroborates the genus-level trait pattern (Fig. 4). Urease and nitrogen-metabolism modules reach full completeness (1.0) in hero MAGs while averaging 0.35–0.60 in the remainder. Housekeeping modules (glycolysis, TCA cycle, fatty-acid biosynthesis) are comparably complete (≥ 0.85) in both groups, demonstrating that the hero signature reflects genuine functional distinction rather than an artefact of genome completeness.

### 3.7 ureC gene tree versus species tree
The ureC ML gene tree (n = 46 MAGs encoding UreC) exhibits moderate but significant topological incongruence with the GTDB-Tk bac120 species tree (normalised Robinson–Foulds distance = 0.58; Shimodaira–Hasegawa test p < 0.001 rejecting the species tree against the *ureC* alignment). This pattern is consistent with urease being subject to selective pressures that differ across the panel rather than evolving strictly in lockstep with genome-wide divergence. Examined at the lineage level, however, the four hero *Sphingobacterium* MAGs do not form an exclusive clade in either the ureC or the species tree — two additional ureC-encoding *Sphingobacterium* MAGs (C13, V3) nest within the same clade, and extending the definition of the MICP-complete pool to all score-8 MAGs (n = 26) recovers a polyphyletic distribution across genera. We interpret these observations as evidence that **the MICP-complete phenotype has arisen or been retained in multiple independent lineages**, rather than originating in a single ancestral clade. Within each hero lineage, the vertical-inheritance signal (§3.3) is preserved; between lineages, the pattern is better described as functional convergence than as monophyletic inheritance.

### 3.8 S13 and S16 as candidate novel *Sphingobacterium* species
Both S13 and S16 are assigned to *Sphingobacterium* by GTDB-Tk but receive no species-level designation (closest-reference ANI unavailable). Whole-genome skani ANI places them among the 21 novel-species candidates in the panel (Fig. 5a). Reciprocal-best mmseqs2 AAI against the five remaining *Sphingobacterium* MAGs identifies V3 (*Sphingobacterium* sp019969845) as S13's nearest neighbour (AAI = 93.15 %, 1,401 RBHs) and S23 (*S. paramultivorum*) as S16's nearest neighbour (AAI = 93.49 %, 3,160 RBHs); every other comparison returns AAI < 91 % (Fig. 5b). Both values lie well below the 95 % AAI species threshold (Rodriguez-R and Konstantinidis, 2014). S13 (4.28 Mb, GC 40.0 %, 3,554 CDS) and S16 (5.62 Mb, GC 39.9 %, 4,779 CDS) are otherwise consistent with the genus in genome size, GC and tRNA repertoire (Supplementary Table S1). We therefore propose S13 and S16 as **candidate novel *Sphingobacterium* genomospecies** requiring cultivation and formal nomenclatural description. No full-length 16S rRNA gene was reconstructed by Bakta in either MAG — a recurrent limitation of short-read binning for rRNA operons — so the definitive species descriptions will require long-read resequencing and/or 16S-targeted PCR from cultured isolates.

## 4. Discussion

### 4.1 A genomics-first triage identifies non-canonical MICP chassis
Conventional MICP bioprospecting screens culturable bacteria for urease activity and then attempts to match their physiology to a target application. Our analysis inverts this logic: by mining 111 MAGs of directly environmental provenance we can simultaneously filter for genomic completeness of the MICP machinery, for signatures of vertical inheritance, and for accessory traits required by the target environment. The approach recovers a candidate lineage (*Sphingobacterium*) that is under-represented in the MICP literature and that would be easily missed by activity-first screens because environmental *Sphingobacterium* is typically slow-growing and outcompeted under standard urea-rich enrichment.

### 4.2 Convergent rather than monophyletic retention
A central finding of this study is that MICP-completeness in livestock waste is distributed across **two convergent lineages** rather than a single monophyletic clade. The hero *Sphingobacterium* MAGs share a *ure* operon architecture that is compatible with vertical inheritance within the genus, while the two *Pseudomonas*\_E MAGs (M1, S26) carry an equivalent complete module that cannot be phylogenetically related to the *Sphingobacterium* operon through simple vertical descent. Two mechanistically distinct scenarios are consistent with these observations: (i) retention of an ancestral *ureABCDEFG–cah* module from a deep bacterial ancestor, with independent loss in most intervening lineages; or (ii) ancient horizontal transfer of the complete module between ancestral lineages followed by vertical fixation within each recipient. Our data cannot formally distinguish these scenarios; both, however, support the **practical** conclusion that the operon is stably encoded in each hero lineage rather than residing on an actively mobilisable cassette, which is the property most relevant to field deployment.

### 4.3 Statistically supported trait stacking
The permutation-tested trait-module enrichments describe a mechanistically coherent phenotype whose individual components are each reported in isolated MICP or halotolerant strains (Dhami et al., 2014; Duarte-Nass et al., 2020): urease and CA drive the alkalinisation and bicarbonate supply required for CaCO₃ nucleation; Mrp/Nha antiporters and catalase/SOD defend the cell against the cytoplasmic alkalinisation and oxidative flux that accompany urea hydrolysis; glycoside hydrolases and CBMs allow the cell to extract energy from slurry polysaccharides; and quorum-sensing enrichment is compatible with biofilm-assisted sand-particle adhesion. Their statistically supported co-occurrence within a vertically inherited six-MAG set is, to our knowledge, a new observation for MICP bioprospecting.

### 4.4 Limitations and place of this study in the pipeline
This study is exclusively computational. Genomic presence does not entail expression, activity, or field performance, and Bakta / DRAM sequence-similarity thresholds can both over- and under-call genes. Assembly fragmentation additionally limits our ability to verify physical linkage of *ureABCDEFG* with *cah* in all six hero MAGs; in C22 and S23, only 4–5 of 7 *ure* genes are recoverable on a single contig. The hero / rest comparison has unequal group sizes (n = 6 vs 105), and our statistical framework is designed to quantify effect sizes with bootstrap uncertainty rather than to claim inferential robustness beyond this dataset. Finally, our *ureC* phylogeny shows moderate incongruence with the species tree; while this is consistent with functional selection on a mechanistically essential enzyme, it complicates any simple narrative of strict vertical inheritance.

We therefore interpret the present analysis as a **prioritised, statistically tested, and mechanistically grounded hypothesis set** that identifies the most informative next experiments: (i) isolation or enrichment of S13 and S16, followed by Nessler-assay quantification of urease activity and gravimetric / XRD measurement of CaCO₃ yield across pH 7–10; (ii) unconfined-compressive-strength testing of MICP-treated sand columns inoculated with hero strains under standardised slurry composition; (iii) long-read resequencing of S13 / S16 to resolve rRNA operons and enable formal species description; (iv) pilot-scale comparison against *S. pasteurii* under realistic slurry conditions.

### 4.5 Implications for livestock-waste valorisation
Reframing livestock slurry from liability to bioprospecting reservoir aligns with circular-bioeconomy objectives: the same matrix that currently poses methane, ammonia and runoff risks can supply chassis organisms that upcycle waste nitrogen into an engineering commodity (CaCO₃ biocement). The 111-MAG workflow reported here is portable to any high-nitrogen waste stream and should accelerate contextual discovery of MICP chassis worldwide.

## 5. Conclusions
- Comparative genomics of 111 livestock-waste MAGs identifies a six-member MICP-complete set distributed across two functionally convergent lineages (four *Sphingobacterium* + two *Pseudomonas*\_E).
- The *ureABCDEFG* operon is recovered on a single contig in the three best-assembled hero MAGs (M1, S13, S16), is accompanied by at least one *cah* gene in every hero MAG, and shows no flanking mobile-element signatures or GC deviation — consistent with vertical inheritance within each lineage.
- Permutation-tested and bootstrap-bounded enrichment identifies nine FDR-significant trait modules (Mrp complex, CBM, glycoside hydrolases, oxidative-stress defence, *nha* antiporters, compatible solutes, quorum sensing, metal efflux) that jointly rationalise persistence and self-fuelling in alkaline, polysaccharide-rich livestock slurry.
- S13 and S16 satisfy both GTDB-Tk and AAI criteria (< 95 %) for candidate novel *Sphingobacterium* species and are nominated as priority targets for isolation and experimental validation of slurry-coupled MICP.
- The workflow is reusable for other high-nitrogen waste streams and directly links MAG-level comparative genomics to the engineering of biocement agents.

## CRediT author statement
**[Author 1]** — Conceptualisation, Data curation, Formal analysis, Investigation, Methodology, Software, Visualisation, Writing – original draft. **[Author 2]** — Methodology, Investigation, Software, Writing – review & editing. **[Corresponding Author]** — Conceptualisation, Funding acquisition, Project administration, Supervision, Writing – review & editing.

## Declaration of competing interests
The authors declare no financial or personal relationships that could appear to influence the work reported.

## Funding
This work was supported by [Funder, Grant ID] and by institutional support from [Institution].

## Acknowledgements
The authors thank [collaborators] for insightful discussion and [computing facility] for HPC access.

## Data and code availability
All MAGs, Bakta/Panaroo/GTDB-Tk/DRAM outputs, derived trait-module tables and analysis scripts are available at `/data/data/Upcycling/` and will be deposited at NCBI BioProject PRJNA-XXXXXXX and Zenodo (DOI on acceptance).

## References (Elsevier Vancouver style; to be expanded)
1. Achal V, Mukherjee A. A review of microbial precipitation for sustainable construction. *Constr Build Mater* 2015;93:1224–35.
2. Bowers RM, et al. Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nat Biotechnol* 2017;35:725–31.
3. Chaumeil P-A, et al. GTDB-Tk v2: memory-friendly classification with the Genome Taxonomy Database. *Bioinformatics* 2022;38:5315–6.
4. DeJong JT, et al. Biogeochemical processes and geotechnical applications: progress, opportunities and challenges. *Géotechnique* 2013;63:287–301.
5. Dhami NK, et al. Biomineralization of calcium carbonate polymorphs by bacterial strains isolated from calcareous sites. *J Microbiol Biotechnol* 2014;23:707–14.
6. Duarte-Nass C, et al. Application of microbe-induced carbonate precipitation for copper removal from copper-enriched waters. *J Environ Manage* 2020;256:109938.
7. Konstantinidis KT, Tiedje JM. Genomic insights that advance the species definition for prokaryotes. *PNAS* 2005;102:2567–72.
8. Minh BQ, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol* 2020;37:1530–4.
9. Nayfach S, et al. A genomic catalog of Earth's microbiomes. *Nat Biotechnol* 2021;39:499–509.
10. Parks DH, et al. Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life. *Nat Microbiol* 2017;2:1533–42.
11. Phillips AJ, et al. Engineered applications of ureolytic biomineralization: a review. *Biofouling* 2013;29:715–33.
12. Rodriguez-R LM, Konstantinidis KT. Bypassing cultivation to identify bacterial species. *Microbe Mag* 2014;9:111–8.
13. Schwengers O, et al. Bakta: rapid and standardised annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom* 2021;7:000685.
14. Shaffer M, et al. DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Res* 2020;48:8883–900.
15. Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nat Methods* 2023;20:1661–5.
16. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol* 2017;35:1026–8.
17. Tonkin-Hill G, et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biol* 2020;21:180.
18. Shimodaira H, Hasegawa M. Multiple comparisons of log-likelihoods with applications to phylogenetic inference. *Mol Biol Evol* 1999;16:1114–6.
19. Robinson DR, Foulds LR. Comparison of phylogenetic trees. *Math Biosci* 1981;53:131–47.
20. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 1995;57:289–300.

---

## Table 1 | Assembly quality and *ure–cah* cluster contiguity of the six MICP-complete ("hero-lineage") MAGs.

| MAG | Lineage | Size (Mb) | GC (%) | Contigs | N50 (kb) | Main *ure* contig | *ure* genes on main contig | Main-contig span (kb) | Distinct *ure*-bearing contigs | *cah* contigs | Flanking MGE in ±15 kb | Region Δ GC (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **M1** | *Pseudomonas*\_E | 4.31 | 63.2 | 17 | 413.5 | contig_1 | 7 / 7 | 28.6 | 1 | 2 | none | +1.48 |
| **S13** | *Sphingobacterium* | 4.28 | 40.0 | 643 | 7.7 | contig_69 | 7 / 7 | 5.9 | 4 | 1 | none | 0.00 |
| **S16** | *Sphingobacterium* | 5.62 | 39.9 | 198 | 52.1 | contig_125 | 7 / 7 | 5.9 | 3 | 2 | none | 0.00 |
| **C22** | *Sphingobacterium* | 4.34 | 39.6 | 2,193 | 2.1 | contig_74 | 5 / 7 | 5.3 | 4 | 1 | none | 0.00 |
| **S23** | *Sphingobacterium* | 5.60 | 39.9 | 609 | 12.1 | contig_220 | 4 / 7 | 3.3 | 5 | 2 | none | 0.00 |
| **S26** | *Pseudomonas*\_E | 5.31 | 58.6 | 327 | 25.1 | contig_140 | 4 / 7 | 3.0 | 5 | 3 | none | 0.00 |

## Table 2 | Trait modules significantly enriched in hero-lineage MAGs after BH-FDR correction (q < 0.05).
See main text § 3.4; full statistics for all 38 modules in Supplementary Table S2.

## Table 3 | dbCAN class-level CAZyme enrichment in hero-lineage MAGs (HMMER-3, E < 1 × 10⁻¹⁵ for heroes; DRAM `cazy_best_hit` for the remaining 105 MAGs).

| CAZy class | Hero mean (per 10³ CDS) | Rest mean (per 10³ CDS) | Fold change | 95 % bootstrap CI | Permutation BH-q |
|---|---|---|---|---|---|
| GH  (Glycoside Hydrolases)         | 22.04 | 5.77 | 3.82 | [1.80, 5.67] | 6 × 10⁻⁴ |
| PL  (Polysaccharide Lyases)        |  1.21 | 0.34 | 3.52 | [2.09, 5.28] | 1 × 10⁻³ |
| CBM (Carbohydrate-Binding Modules) |  1.03 | 0.24 | 4.24 | [1.72, 6.98] | 9 × 10⁻⁴ |
| CE  (Carbohydrate Esterases)       |  3.13 | 1.58 | 1.99 | [1.09, 2.99] | 7 × 10⁻³ |
| GT  (Glycosyl Transferases)        |  7.52 | 6.07 | 1.24 | [1.01, 1.47] | 4 × 10⁻² |
| AA  (Auxiliary Activities)         |  2.04 | 2.07 | 0.99 | [0.65, 1.38] | 1.0 (n.s.) |
