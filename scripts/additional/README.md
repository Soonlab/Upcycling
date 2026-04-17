# Additional in-silico analyses for Upcycling MICP project

Work performed 2026-04-18 after manuscript was at submission-ready state. Goal: strengthen reviewer defense without wet-lab.

## Layout

```
additional/
├── A1_biosafety/      abricate CARD/VFDB/ResFinder/PlasmidFinder scans
├── A2_structure/      UreC active-site residue conservation (MSA vs S. pasteurii)
├── A3_pseudomonas_ani/
│    ├── skani ANI of M1/S26 vs 146 Pseudomonas_E refs
│    └── A3b:  prodigal + hmmsearch screen of MICP operon architecture
├── A4_genomad/        geNomad plasmid/virus detection
├── A5_alkaliphile/    proteome pI + Mrp/Nha antiporter scan
├── A6_metabolic/      MICP pathway stoichiometry (FBA solver unavailable)
├── A7_grodon/         gRodon2 codon-bias doubling-time prediction
├── B_rarity_screen/   MGnify cow/sheep/pig/chicken catalog rarity screen
├── hmms_shared/       Pfam HMMs (urease + CA classes)
├── figures/           all figures (png + pdf)
├── make_figures.py    figure-generation script
└── REPORT_consolidated.md     full consolidated report
```

## Key results (one-liners)

| Analysis | Result |
|---|---|
| A1 Biosafety | Sphingobacterium heroes 0-0-0-0 across CARD/VFDB/ResFinder/PlasmidFinder. Pseudomonas_E heroes have only intrinsic housekeeping pumps/T4P, no acquired AMR or plasmid replicons. |
| A2 UreC active site | 100% (42/42) residue conservation vs S. pasteurii P41020. |
| A3 Pseudomonas_E ANI | M1 = 97.98% vs GCF_025837155.1; S26 = 97.54% vs P_E helleri — assigned species. |
| A3b Pseudomonas_E MICP | Gene-level MICP near-universal (93.8% UreC); single-contig ureC+CA rarer (36.3%). |
| A5 Alkaliphile | Mrp antiporter 11.7× enrichment in MICP-complete (p = 5.3e-4 MWU). |
| A6 Stoichiometry | 6/6 heroes have urease+CA+Ca-handling; Ca pathway 83.3% hero vs 3.8% rest. |
| A7 gRodon | Median doubling 1.06 h hero vs 1.10 h rest — no growth penalty. |
| B MGnify rarity | 3.07% of 7,599 livestock MAGs are MICP gene-complete. Our 6/6 = 100%, ~30× enrichment. Sphingobacterium in livestock n=2, both MICP-negative — our novel Sphingobacterium lineage is absent from global livestock MAG sets. |
| A4 geNomad | In progress. |

## Conda envs

- `biosafety`: abricate 1.4.0 + AMRFinderPlus 4.2.7
- `mge_tools`: geNomad 1.12
- `grodon`: R + gRodon2 + coRdon + Biostrings
- `carveme`: CarveMe + cobra (FBA skipped — solver licensing)
- `dram_env` (pre-existing): skani, mmseqs, hmmer, mafft, prodigal, biopython, pandas

## Reproduce

```bash
# from repo root
for sh in A1_biosafety/run_abricate_all.sh \
          A3_pseudomonas_ani/A3b_resume_hmm.sh \
          A4_genomad/run_genomad.sh \
          B_rarity_screen/01_download_mgnify_catalogs.sh; do
    bash additional/$sh
done
python additional/A5_alkaliphile/run_alkaliphile_signature.py
python additional/A2_structure/ureC_active_site_conservation.py
Rscript additional/A7_grodon/run_grodon.R
python additional/A6_metabolic/run_stoich_summary.py
python additional/B_rarity_screen/03_screen_mgnify_functional.py
python additional/make_figures.py
```

## Manuscript integration plan

Suggested Revised-Supplementary additions:
- **New Supplementary Figure S11** — A1 biosafety panel boxplots.
- **New Supplementary Figure S12** — A2 UreC active-site heatmap (highlight novel S13/S16 catalytic conservation).
- **New Supplementary Figure S13** — B MGnify livestock MAG rarity bar plot.
- **New Supplementary Figure S14** — A7 growth rate + A5 Mrp antiporter box plots.
- **New Supplementary Tables**
    - S11: per-MAG biosafety hit matrix (CARD/VFDB/ResFinder/PlasmidFinder).
    - S12: UreC per-residue alignment vs *S. pasteurii*.
    - S13: per-species MICP profile across 4 MGnify livestock biomes (7,599 rows).
    - S14: gRodon doubling times + per-MAG pI + Mrp counts.

Main text one-paragraph additions:
- **Discussion 3.4 (new)**: "To bound the external novelty of the MICP-complete lineage, we screened 7,599 livestock-associated MAGs from four MGnify biome catalogs (cow-rumen v1.0.1, sheep-rumen v1.0, pig-gut v1.0, chicken-gut v1.0.1). Only 233/7,599 (3.07%) encoded the gene-complete urease + carbonic anhydrase profile, and only 265/7,599 (3.49%) displayed a single-contig ureC + CA architecture. Critically, Sphingobacterium was represented by only two species clusters in the combined catalog (both chicken-gut), and neither met our gene-complete criterion — placing our S13/S16/S23/C22 Sphingobacterium clade outside the existing livestock MAG reference set."
- **Discussion 3.5 (new)**: "To confirm that the novel Sphingobacterium UreCs retain catalytic potential despite amino-acid divergence (ANI <95% in S13/S16), we aligned the MICP-complete UreCs against the biochemically characterised *Sporosarcina pasteurii* UreC (P41020, PDB 4CEU) and examined seven canonical active-site residues. All 42 residue assignments (6 MAGs × 7 sites: H137, H139, K220, H249, H275, C322, D363) were conserved, indicating that both Ni²⁺ coordination and the catalytic dyad are preserved."
- **Biosafety sidebar (Methods/Discussion)**: "Automated scans of CARD, VFDB, ResFinder, and PlasmidFinder returned zero acquired AMR hits across all six MICP-complete MAGs, zero virulence factors and zero plasmid replicons in the Sphingobacterium MAGs (S13/S16/S23/C22), and only intrinsic efflux/type-IV-pilus genes in M1/S26 — consistent with non-pathogenic soil/waste-associated ecotypes suitable as candidate biocement chassis."
