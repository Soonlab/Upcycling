# Upcycling livestock waste into MICP chassis — comparative genomics of 111 MAGs

Comparative-genomic analysis of 111 metagenome-assembled genomes (MAGs) recovered from four livestock-waste microbiomes (cattle, swine, sheep, poultry). The study identifies two functionally convergent lineages (four *Sphingobacterium* + two *Pseudomonas*\_E MAGs) that retain a complete *ureABCDEFG*–carbonic-anhydrase (*cah*) module, and nominates two candidate novel *Sphingobacterium* species (S13, S16) as priority chassis for alkali-tolerant microbially induced carbonate precipitation (MICP).

> **Status:** In-silico analysis complete. Wet-lab validation is planned but not part of the current release.

---

## Repository layout

```
.
├── figures/
│   ├── main/                          Fig. 1–6 (PNG 300 dpi + PDF vector)
│   └── supplementary/                 Fig. S1–S7
├── results/
│   ├── main/                          Pangenome + GTDB-Tk summaries
│   ├── extra/                         Trait-module scans, ANI, AAI, novelty screen
│   └── revision/                      Cluster audit, permutation stats, dbCAN HMMER, ureC tree
└── scripts/                           All Python / shell code used to produce the above
    ├── 01_main_figures.py             Fig. 1–3 (tree + MICP heatmap, synteny)
    ├── 02_trait_module_scans.py       Bakta keyword scan for six trait categories
    ├── 03_ANI_novelty_skani.py        Whole-genome ANI + novelty screen
    ├── 04_novel_species_AAI.py        S13/S16 AAI vs other Sphingobacterium
    ├── 05_dram_figure.py              DRAM distillate heatmap (Fig. 4)
    ├── 06_dram_distill.sh             DRAM distill pipeline wrapper
    ├── 07_dram_hero_annotate.sh       DRAM annotate for hero MAGs
    └── revision/                      Scripts #1–5 from the methodological revision
        ├── 01_cluster_and_quality.py
        ├── 02_ureC_gene_tree.py
        ├── 02b_hero_topology_check.py
        ├── 03_permutation_stats.py
        ├── 04_dbcan_reanalysis.py
        └── 04b_dbcan_final.py
```


> **Note on manuscript files**
> The manuscript draft, cover letters and figure-legend document are kept privately
> and will be released only upon publication (as a tagged GitHub release + Zenodo DOI).
> Figures and all analysis code / tables remain public here.

## Analysis pipeline

1. **Annotation** — Bakta v1.9 (`bakta_results/`) over the 111 MAG FNA files.
2. **Pan-genome** — Panaroo v1.5 (moderate mode, 90 % identity), producing `gene_presence_absence.csv` and a core-gene alignment.
3. **Phylogenomics** — GTDB-Tk v2.4 `classify_wf` (release r220) + IQ-TREE v2.3 (ModelFinder, 1 000 UFBoot) on the bac120 alignment.
4. **Metabolic reconstruction** — DRAM v1.5 annotate + distill against KOfam, UniRef90, Pfam, dbCAN v12, MEROPS.
5. **Genome-wide ANI** — skani v0.3 `triangle --full-matrix`.
6. **Amino-acid identity** — mmseqs2 easy-search reciprocal-best hits (≥ 30 % id, ≥ 70 % cov) on Bakta proteomes.
7. **Trait-module screening** — six keyword-based modules (biofilm/EPS, ammonia handling, mobile-genetic elements, alkaline/osmotic defence, CAZymes, metal/antibiotic resistance); per-1 k-CDS normalisation.
8. **Statistical testing** — 10 000-iteration one-sided label permutation + 2 000-iteration bootstrap 95 % CI; Mann–Whitney U; Benjamini–Hochberg FDR.
9. **CAZyme validation** — direct HMMER-3 `hmmsearch` (E < 1 × 10⁻¹⁵) against dbCAN v12 on the six hero MAGs; combined with DRAM `cazy_best_hit` for the remaining 105.
10. **Gene-tree vs species-tree** — MAFFT alignment of UreC proteins, IQ-TREE ML gene tree, Robinson–Foulds distance (ete3) and Shimodaira–Hasegawa / AU tests against the bac120 species tree.

## Reproducing the analysis

```bash
# conda env used throughout
conda create -n dram_env python=3.10 -y
conda activate dram_env
pip install pandas numpy matplotlib seaborn biopython scipy ete3 six
conda install -c bioconda mafft iqtree skani mmseqs2 hmmer -y

# example (main figures)
python scripts/01_main_figures.py

# revision pipeline
python scripts/revision/01_cluster_and_quality.py
python scripts/revision/02_ureC_gene_tree.py
python scripts/revision/02b_hero_topology_check.py
python scripts/revision/03_permutation_stats.py
python scripts/revision/04_dbcan_reanalysis.py
python scripts/revision/04b_dbcan_final.py
```

Paths inside the scripts are set to `/data/data/Upcycling/` (original analysis host). For reuse, edit the `BASE` / `BAKTA` / `OUT` variables near the top of each script.

## Key results

| Finding | Evidence | File |
|---|---|---|
| Six MICP-complete MAGs across two convergent lineages (4 *Sphingobacterium* + 2 *Pseudomonas*\_E) | Bakta gene presence + IQ-TREE phylogeny | `results/main/MICP_Pangenome_Final_Summary.csv`, `figures/main/Fig1_*` |
| *ureABCDEFG* operon on a single contig in M1 / S13 / S16 (5.9–28.6 kb span) | GFF3 contig audit | `results/revision/hero_cluster_audit.csv` |
| No transposase / integrase / prophage / relaxase in ± 15 kb of any hero *ure* cluster; Δ GC ≤ 1.5 % | Bakta window scan | `results/extra/HGT_ureCah_cluster.csv` |
| Nine trait modules enriched in hero lineages at BH-FDR q < 0.05 (Mrp 10.9×, CBM 9.8×, GH 4.7×, …) | 10 000-perm test + bootstrap CI | `results/revision/Hero_vs_Rest_permutation_stats.csv`, `figures/main/Fig_Permutation_forest.*` |
| Rigorous dbCAN HMMER confirms GH / PL / CBM / CE hero enrichment (q < 0.01) | hmmsearch E < 1e-15 + DRAM | `results/revision/dbCAN_final_hero_vs_rest_class.csv` |
| S13 and S16 candidate novel *Sphingobacterium* species (max congeneric AAI 93.2 % / 93.5 %) | skani ANI + mmseqs2 AAI | `results/extra/novel_species/AAI_S13_S16_vs_Sphingobacterium.csv` |
| UreC gene tree incongruent with species tree (normRF = 0.58; SH / AU p < 10⁻³⁷) | IQ-TREE gene tree + `-z` test | `results/revision/ureC_tree/` |

## Citation

If you use the code or data in this repository, please cite:

> *\[Author list\]*. **Comparative genomics of 111 livestock-waste metagenome-assembled genomes nominates novel *Sphingobacterium* species with a vertically retained ure–cah module for alkali-tolerant microbially induced carbonate precipitation.** *\[Journal\]*, 2026 (in preparation). Repository: https://github.com/Soonlab/Upcycling

## License

Code: MIT. Figures, tables and manuscript drafts: CC-BY 4.0. Raw MAG sequences are not included in this repository; they will be released at NCBI BioProject PRJNA-XXXXXXX upon publication.

## Contact

Maintainer: [@Soonlab](https://github.com/Soonlab). For questions on reproducing the analysis, open a GitHub issue.
