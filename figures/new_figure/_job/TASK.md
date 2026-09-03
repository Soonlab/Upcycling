# Upcycling figure rebuild — build contract (2026-09-04)

Rebuild every manuscript figure (main Fig 1–8, Suppl Fig S1–S21 = 29 pages) as an editable
vector page, applying the figure method settled during the Rectal_Organoid `svg_supp` work
(2026-08-20/21) and used for the SNUH_KMJ `new_figure` rebuild (2026-09-03/04).
The final deliverable the user receives is the **PNG**; SVG and PDF are kept alongside.

Output root: `/data/data/Upcycling/new_figure/`
Shared style: `new_figure/_style.py` (import as `import _style as st`).
Python: **`/home/soon/miniconda3/envs/dram_env/bin/python`** (matplotlib 3.10.8, pandas 1.5.3,
scipy 1.15, numpy, Biopython 1.87 for Newick trees, openpyxl). The base miniconda python has
no Biopython; the system python3 has nothing.

Build scripts: `new_figure/build_fig{N}.py` (N = 1..8) and `new_figure/build_supS{N}.py`
(N = 1..21). Output stems: `figures/Fig{N}.{svg,pdf,png}` and `figures/Fig_S{N}.{svg,pdf,png}`.
**Numbering is unchanged from the SUBMISSION package** — no consolidation, no renumbering, so
the manuscript and legends (`SUBMISSION/02_Figure_legends.md`) keep their callouts. Supplementary
figures that were several files (S4 a/b/c, S16 defense/crispr, S17 a–d, S19 + by_source) become
ONE composed page each.

Old figures for reference (look at them, but rebuild from the SOURCE, not from the picture):
`SUBMISSION/Figures_main/Figure_{N}.png`, `SUBMISSION/Figures_supplementary/Figure_S{N}*.png`.
Legends: `SUBMISSION/02_Figure_legends.md` (what each panel must show).
Old builders (data recipes only — do not copy their titles/annotations):
`/home/soon/Upcycling_repo/scripts/01_main_figures.py` (Fig 1–3),
`scripts/additional/make_figure_8_multipanel.py` (Fig 8),
`scripts/pptx_builder/build_editable_figures.py` (every panel incl. genus aggregation helpers),
`/data/data/Upcycling/research/revision/*.py` (Fig 5c/6/7 recipes, ureC tree).

## 1. Output format

- One **composed page per figure**. Page width is **180 mm** (`st.PAGE_W_MM`); height is whatever
  the content needs, re-fitted after the content is final (no empty band at the bottom).
- Emit **SVG + PDF + PNG** through `st.save(fig, outdir, stem)` (svg.fonttype none, pdf fonttype 42,
  PNG 200 dpi, SVG font-family rewritten to Arial-first).
- Lay out in millimetres with `fig, ax_mm, text_mm, letter = st.page(height_mm)`;
  `ax_mm(l, t, w, h)` places an Axes by its top-left corner from the page top-left.

## 2. What may appear on the figure

Allowed: axis labels, tick labels, axis titles; legends; value labels bound to a datapoint and stat
columns bound to rows (e.g. a `P` column beside a lollipop row); short block headers separating data
blocks inside one panel; table column headers; direction cues (e.g. "higher in MICP-complete →");
panel letters; a bare identifier (MAG ID, gene name, database name) only when the panel is otherwise
unidentifiable.

Forbidden: **panel titles**; **explanatory sentences**, method descriptions, conclusion statements,
italic stat footnotes; any free-floating narrative text (e.g. "present study 6/6 = 100%",
"42/42 matches", "T3PKS 23×"). `st.prose_scan(fig)` must return an empty list.

## 3. Layout

- Panel letters **A, B, C …** (uppercase, no parentheses), strictly left→right then top→bottom,
  left-aligned on a small number of shared x positions (one per column).
- Uniform type scale from `st.setup()`: body/ticks 7 pt, axis titles 8 pt, panel letters 10 pt bold,
  stat annotations 6.5 pt. Do not override sizes per panel.
- Spines left+bottom only (`st.style_axis`). Card-style boxes → borderless mini tables or stat columns.
  Prefer horizontal lollipops/bars when value labels would collide.
- Long category lists (111 MAGs, 63 RefSeq genomes) run vertically with 7 pt labels; the page grows
  in height rather than shrinking the text.

## 4. Data provenance — the hard rule

- **No data value may be hard-coded.** Every number is loaded or recomputed from a repository source
  (`SUBMISSION/Supplementary_tables/*`, `pangenome_work/*`, `MAGs_FASTA_files/bakta_results/*/*.gff3`,
  `research/**`). Style constants and category names (gene names, MAG IDs, module names) are exempt.
- Derived quantities (fold changes, percentages, medians, MWU P, PCoA coordinates, pseudo-F) are
  **recomputed** from the stored inputs and asserted against the stored value where one exists.
- A number in the old figure with **no source anywhere in the repository is dropped**, and the drop is
  recorded in the journal. Where the old figure and the source disagree, **follow the source** and
  record the discrepancy for the user to adjudicate.
- Known source-fidelity defects to FIX in this rebuild (do not carry over):
  * Fig 3 synteny: draw every track from the Bakta `.gff3` (recipe in `01_main_figures.py`), not from
    inline constants; all 7 `ure` genes must appear on the main contig for M1/S13/S16 per Table S1c.
  * Fig 7B: trait-module PCoA coordinates are NOT in Table S9a (that file holds the pan-genome
    ordination only). Recompute the z-score Euclidean PCoA from Table S2b with the recipe in
    `research/revision/05_PCoA_by_source.py`, recompute PERMANOVA pseudo-F, and assert against
    Table S9b / the legend value (2.71). Record in the journal.
  * Fig 4A: the legend says 111 MAGs × 35 modules but the shipped figure is genus-aggregated.
    Rebuild from Table S5a (`product.tsv`) and choose the representation that is legible at 180 mm;
    record which one was drawn.

## 5. Colour

- Palette lives in `_style.py`. One colour carries **one meaning per page**. Before finishing a figure,
  list every colour used and the meaning it carries in each panel, and confirm no colour carries two
  meanings on the same page. HERO coral = MICP-complete group; REST grey = the other 105; SPHINGO blue /
  PSEUDO orange = the two lineages (use `st.hero_col(mag)`); `st.SOURCE` = waste source; GREEN =
  present/match/complete; `st.seq_cmap()` for completeness/copy-number heat maps.
- Colourbar direction must match the cell encoding.

## 6. Verification, per figure

1. `st.audit(fig)` — text overlaps, text outside the canvas, overlapping Axes. Must be zero.
2. `st.prose_scan(fig)` — must be empty.
3. Read the rendered PNG (Read tool) and confirm: reading order, no clipping, legible at page size,
   every panel identifiable without a title. Tall pages: inspect in vertical slices (crop with PIL),
   a 1,900-px page shrunk to a thumbnail merges small text.
4. Re-fit the page height so there is no large empty band at the bottom.
5. Cross-check the printed numbers against the source file.
6. `_job/check_no_hardcoded.py` must report nothing to review for your scripts (or every flagged
   literal must be a documented style constant).

## 7. Journal

Each build group writes `new_figure/_job/journal_<group>.md` (one file per agent, to avoid write
conflicts) with, per figure: source files used; derived values recomputed and how they were checked;
values dropped for lack of a source; discrepancies found (old figure vs source, legend vs source);
colour–meaning table; final canvas height; anything the user must adjudicate. The coordinator merges
them into `_job/JOURNAL.md`.

## 8. Figure inventory (source → page)

| Page | Content | Primary sources |
|---|---|---|
| Fig1 | bac120 ML tree of 111 MAGs (midpoint-rooted), genus colour strip, MAG labels (species where available), presence/absence heat map of ureA–G + cah, scale bar | `pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile`, `pangenome_work/MICP_Pangenome_Final_Summary.csv`, `Table_S1d_GTDB_Tk_classification.tsv` |
| Fig2 | A module score (0–8) by genus box+points, heroes marked; B per-gene prevalence hero (n=6) vs rest (n=105) | `MICP_Pangenome_Final_Summary.csv`, `Table_S1d` |
| Fig3 | synteny tracks of the densest ure contig per hero (6 tracks, arrows to scale, strand, gene colours, kb axis) | `MAGs_FASTA_files/bakta_results/{C22,M1,S13,S16,S23,S26}/*.gff3`, `Table_S1c`, `Table_S3a/b` |
| Fig4 | A DRAM module-completeness heat map; B MICP-critical module completeness hero vs rest | `Table_S5a_DRAM_product.tsv`, `Table_S1d` |
| Fig5 | A novelty screen (closest-ref ANI per MAG, 95 % line, heroes marked); B AAI of S13/S16 vs other Sphingobacterium (95/70 % lines); C skani ANI of 6 study Sphingobacterium vs 63 RefSeq (ranked bars, 95 % line) | `Table_S8`, `Table_S4b`, `Table_S10a/b` |
| Fig6 | permutation forest: fold change (log x) + bootstrap CI, colour by BH-FDR tier | `Table_S2c` |
| Fig7 | A Panaroo Jaccard PCoA by source, heroes ringed, PERMANOVA F/p; B trait-module z-score PCoA (recompute) | `Table_S9a/b/c`, `Table_S2b`, `research/revision/05_PCoA_by_source.py` |
| Fig8 | A active-site 6×7 match heat map; B TM-score + RMSD bars; C MGnify % per biome; D antiSMASH class means hero vs rest | `Table_S12`, `Table_S22`, `Table_S14a`, `Table_S23b` |
| Fig_S1/S2/S3/S5 | genus-aggregated log10 heat maps of trait subcategories (Biofilm_EPS / Ammonia_N / Alkaline_stress / Metal_AMR) | `Table_S2b` (+ `_genus_aggregate` helper in pptx builder) |
| Fig_S4 | A keyword CAZy proxy heat map; B dbCAN class per-1k-CDS; C dbCAN family counts | `Table_S2b` (CAZyme_proxy cols), `Table_S6b`, `Table_S6c/d` |
| Fig_S6 | skani ANI heat map among the 6 heroes | `Table_S4a` |
| Fig_S7 | ureC ML gene tree (46 tips), bootstrap, heroes marked; RF / SH-AU values as a stat column | `Table_S7c` newick, `Table_S7a`, `Table_S7f` |
| Fig_S8 | 4 DBs × hits per MAG, hero vs rest box + points | `Table_S11` |
| Fig_S9 | active-site heat map, extended (7 sites × 6 MAGs, residue letters) | `Table_S12` |
| Fig_S10 | Pseudomonas_E 146 refs: % with each MICP feature | `Table_S13a/b` |
| Fig_S11 | Mrp, Nha, pI median, pI acidic fraction: hero vs rest box + points | `Table_S15a` |
| Fig_S12 | gRodon doubling time hero vs rest | `Table_S16` |
| Fig_S13 | A % gene-complete per biome; B % single-contig per biome (n per bar) | `Table_S14a` |
| Fig_S14 | gene-copy heat map 6 heroes × stoichiometry columns | `Table_S15b` |
| Fig_S15 | A plasmid/virus contigs hero vs rest; B per-hero ureA/B/C-on-MGE cross-check table | `Table_S17a/b` |
| Fig_S16 | defense systems and CRISPR arrays hero vs rest (paired) | `Table_S18a/b` |
| Fig_S17 | A GC3; B ENC (hero vs rest); C codeml M0 ω per gene; D yn00 pairwise ω hero-hero / rest-rest / hero-rest per gene | `Table_S19a/b/c` (+ raw pairwise in `research/additional/C3_dnds_codon/` if needed for D) |
| Fig_S18 | ANI heat map 6 heroes × reference genomes | `Table_S20` |
| Fig_S19 | A coverage hero vs rest; B coverage by source; C heroes over source distribution | `Table_S21a/b` |
| Fig_S20 | A TM-score per hero; B RMSD per hero | `Table_S22` |
| Fig_S21 | A total BGC regions hero vs rest; B class-level means with MWU P column | `Table_S23a/b` |
