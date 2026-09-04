# Consolidation design — 2026-09-04

User requirement: main figures <= 5, supplementary figures <= 5, tables <= 5 in total,
manuscript body <= 7,000 words. Output must be internally consistent across
figures, tables and manuscript.

Interpretation recorded: "7000자" is read as **7,000 words** of body text
(Introduction -> Conclusions). 7,000 characters would be ~1,100 words, which is not a
viable research article, and every prior word budget in this project has been in words.

## Main figures: 8 -> 5

| new | working title | panels (source) |
|---|---|---|
| Fig 1 | Phylogenomic distribution and completeness of the MICP module | A bac120 tree + genus strip (old Fig1A) · B ureA-G + cah presence matrix (old Fig1B) · C module score by genus (old Fig2A) · D per-gene prevalence, MICP-complete vs rest (old Fig2B) |
| Fig 2 | Architecture and dosage of the ure-cah cluster | A synteny tracks of the six (old Fig3) · B MICP pathway gene dosage heat map (old S14) · C per-MAG urease/CA vs geNomad MGE cross-check (old S15B) |
| Fig 3 | Catalytic and structural conservation of urease | A active-site 6x7 match matrix (old Fig8A) · B ESMFold TM-score + RMSD (old Fig8B) · C codeml M0 omega per urease gene (old S17C) · D yn00 pairwise omega by pair class (old S17D) |
| Fig 4 | Trait-module enrichment and the alkali-tolerance signature | A permutation forest, all rows with FC>1 (old Fig6) · B MICP-critical trait modules, hero vs rest (old Fig4C) · C alkaliphile signature: Mrp, Nha, proteome pI, acidic fraction (old S11) |
| Fig 5 | Novelty, global rarity and engineering background | A novelty ANI screen across the panel (old Fig5A) · B MGnify livestock rarity (old Fig8C) · C antiSMASH BGC class enrichment (old Fig8D) · D pan-genome PCoA by waste source (old Fig7A) |

## Supplementary figures: 21 -> 5

| new | working title | panels (source) |
|---|---|---|
| Fig S1 | Trait-module landscape across genera | A biofilm/EPS (old S1) · B ammonia/nitrogen (old S2) · C alkaline/osmotic stress (old S3) · D CAZy keyword proxy (old S4A) · E metal/antibiotic (old S5) |
| Fig S2 | DRAM metabolic context and dbCAN CAZyme validation | A DRAM fractional modules by genus (old Fig4A) · B DRAM presence/absence modules by genus (old Fig4B) · C dbCAN class abundance hero vs rest (old S4B) · D dbCAN family detail (old S4C/D) |
| Fig S3 | Biosafety, mobile elements and defence systems | A abricate 4-database hits per MAG (old S8) · B geNomad plasmid/virus contigs per MAG (old S15A) · C DefenseFinder systems per MAG (old S16A) · D minced CRISPR arrays per MAG (old S16B) |
| Fig S4 | Genome relatedness and gene-tree congruence | A pairwise ANI among the six (old S6) · B ureC gene tree with species-tree congruence (old S7) · C ANI against the organism-verified MICP reference panel (old S18) · D MICP feature prevalence within Pseudomonas_E (old S10) |
| Fig S5 | Growth rate, in-situ abundance and codon usage | A predicted minimum doubling time (old S12) · B coverage proxy by group (old S19A) · C coverage proxy by waste source (old S19B) · D genome-wide GC3 (old S17A) · E effective number of codons (old S17B) |

**Dropped outright** (each was an explicitly labelled extended duplicate of a main panel):
old S9 (= Fig 8A transposed), S13 (= Fig 8C), S20 (= Fig 8B), S21 (= Fig 8D).
Their underlying numbers survive in the supplementary tables.

## Tables: 26 -> 5

| new | content | absorbs |
|---|---|---|
| Table 1 (main) | The six MICP-complete MAGs: assembly quality, taxonomy, ure-cah contiguity | old Table 1 |
| Table 2 (main) | Trait modules and dbCAN CAZyme classes enriched after BH-FDR (q<0.05) | old Table 2 + Table 3 |
| Table S1 (workbook) | Per-MAG measurements for all 111 MAGs, one sheet per assay | old S1a-d, S2a-b, S6a-c, S11, S15a-b, S16, S17a, S19a, S21a, S23a |
| Table S2 (workbook) | Comparative statistics and group tests | old S2c, S6d, S8, S9a-c, S13a-b, S14a, S18a-b, S19b-c, S21b, S23b |
| Table S3 (workbook) | Reference panels, matrices and method support | old S3a-b, S4a-b, S5a-c, S7a-f, S10a-b, S12, S14b, S17b, S20, S22 |

## Consistency rules for every deliverable
1. Panel letters run A, B, C ... in reading order on each page; the legend letters match.
2. Every manuscript figure callout resolves to a page that exists; every table callout
   resolves to Table 1, 2, S1, S2 or S3 and names the sheet.
3. No value is hard-coded in a build script; every number is loaded or recomputed from a
   repository source and asserted where a stored summary exists.
4. `st.audit` 0 findings and `st.prose_scan` empty for every page.
5. Colour convention unchanged: coral = MICP-complete, grey = rest, blue = Sphingobacterium,
   orange = Pseudomonas_E, green = present/match, source palette for waste source.
