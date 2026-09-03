# Build journal — group main_1_3 (main Fig 1, Fig 2, Fig 3)

Built 2026-09-04 with `/home/soon/miniconda3/envs/dram_env/bin/python` (matplotlib 3.10.8,
pandas 1.5.3, Biopython 1.87). Contract: `_job/TASK.md`.

| Figure | Canvas | audit | prose_scan | check_no_hardcoded |
|---|---|---|---|---|
| Fig1 | 180 × 288.75 mm | PASS (135 text items, 0 overlaps) | PASS | clean |
| Fig2 | 180 × 78 mm | PASS (51 text items, 0 overlaps) | PASS | clean |
| Fig3 | 180 × 129 mm | PASS (116 text items, 0 overlaps) | PASS | 1 literal, justified below |

---

## 🔴 BLOCKER for the user — the MICP presence matrix behind Fig 1 and Fig 2

The shipped Fig 1 and Fig 2 were drawn from `pangenome_work/MICP_Pangenome_Final_Summary.csv`.
That file cannot be used and cannot be reproduced:

1. **It covers 100 of the 111 MAGs.** `C1` and `C10`–`C19` are absent. The old Fig 1 code
   initialised the heat map to zeros and filled only the rows it found, so those eleven MAGs
   were drawn as carrying *no* MICP gene. The old Fig 2b x-axis admits this in passing
   ("Other MAGs (n = 94 of 100 Panaroo-QC MAGs)"), while the manuscript and legend describe
   6 vs 105.
   Cause: `pangenome_work/summarize_pangenome_micp.py` takes `df.columns[14:]` as the sample
   columns, but this Panaroo `gene_presence_absence.csv` has only three metadata columns, so
   the first eleven sample columns were skipped.
2. **The named recipe does not reproduce the file.** Re-running that script's logic over
   `pangenome_work/results/gene_presence_absence.csv` matches only `ureD/E/F/G`; the Panaroo
   `Gene` column contains no `ureA`, `ureB`, `ureC` or `cah` entry at all, so the recipe yields
   zero for those four columns while the shipped file has ones. The file's true provenance is
   unknown.
3. **It contradicts the paper's own hero set.** In that file `S26` scores **1/8** and `S13`
   has `cah = 0`, although both are MICP-complete MAGs. This is why the shipped Fig 2b shows
   MICP-complete prevalence of **83.3 % (5/6)** for every gene while its legend claims 100 %.
4. **No source in the repository yields "the complete 8/8 module is restricted to six MAGs".**
   Counts of MAGs carrying all eight genes:

   | source | MAGs covered | MAGs with 8/8 |
   |---|---|---|
   | `MICP_Pangenome_Final_Summary.csv` (shipped) | 100 | 26 |
   | `Table_S1a_ace_samples_list.csv` | 45 ure-bearing (of 111) | 28 |
   | Bakta GFF3 keyword scan, all 111 MAGs (this rebuild's recipe) | 111 | 32 |

   The six-MAG designation is therefore **not** "score == 8". It is a curated flag, carried
   identically in seven analysis tables (S11, S15a, S17a, S18a, S19a, S21a, S23a), and it
   evidently also depends on single-contig operon architecture (Table S1c). **The Fig 1 legend
   sentence "The complete 8/8 module is restricted to four *Sphingobacterium* MAGs and two
   *Pseudomonas*_E MAGs" is false against every available source and must be rewritten**, or
   the MICP-complete criterion must be stated explicitly in the Methods and reflected in the
   figure. This is a manuscript decision, not a figure decision.

**What was drawn instead.** `Table_S1a_ace_samples_list.csv` — a shipped supplementary table
of record giving per-gene copy counts for the 45 ure-bearing MAGs — binarised, with MAGs
absent from it scoring 0. It is the only per-gene source in which all six MICP-complete MAGs
carry all eight genes, and it covers the full 111-MAG panel once absence is read as "no
detected MICP gene". Both Fig 1 and Fig 2 use it, so the two pages now agree with each other.
Asserted in both scripts: `score[heroes] == 8` for all six.

Consequence the user should expect to see: panel Fig 2A now shows MAGs at 8/8 in
*Chryseobacterium*, *Acinetobacter* and *Pseudomonas*_E as well, and Fig 1's heat map shows
28 fully green rows. That is what the source says.

---

## Fig 1 — phylogeny + MICP gene presence

**Sources**
- `pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile` — topology and branch
  lengths, midpoint-rooted for display (`tree.root_at_midpoint()`), drawn as a rectangular
  phylogram from `tree.depths()` rather than through `Phylo.draw`.
- `Table_S1d_GTDB_Tk_classification.tsv` — genus (`g__`) and species (`s__`) per MAG.
- `Table_S1a_ace_samples_list.csv` — gene presence (see BLOCKER above).
- `Table_S15a_alkaliphile_signature_per_MAG.csv` — the MICP-complete flag.

**Recomputed and asserted**
- tip count == 111 == rows of the GTDB table, and the set of tip MAG ids equals the set of
  MAG ids in the taxonomy table.
- the nine genera given their own strip colour are exactly the nine most frequent genera
  (`gcount.head(9)`), asserted rather than typed; the remaining 8 MAGs pool into "Other
  genera" with the count computed, not written.
- every MICP-complete MAG carries all eight genes in the drawn matrix.
- scale bar length is derived (`1/2/5 × 10^floor(log10(xmax/4))`, largest that fits in a third
  of the tree width) and printed with its own computed value — 0.2 substitutions/site here.

**Dropped**
- Ultrafast bootstrap support. The values are present in the treefile but were not drawn in
  the shipped figure either; at 111 tips over a 50 mm tree they cannot be rendered legibly.
  The legend's mention of 1,000 UFBoot replicates describes the method, not a figure element,
  so nothing on the page claims support.

**Discrepancy (minor)** — the legend says the colour strip encodes "the top ten genera";
there are nine genera with ≥ 2 MAGs plus a pooled "Other genera" class of 8 single-MAG genera.
Drawn as nine + pooled, which is what the data supports.

**Colour meanings**

| colour | meaning (whole page) |
|---|---|
| blue `#1F4E78` | genus *Sphingobacterium* (a MICP-complete lineage) |
| orange `#D9772B` | genus *Pseudomonas*_E (the other MICP-complete lineage) |
| lilac / tan / mauve / teal / grey / olive / brown | the other seven strip genera |
| light grey `#D9D9D9` | pooled "Other genera" |
| green `#2F7D4F` | gene present (panel B); white = absent |
| coral `#C53E1F` | a MICP-complete MAG — tip label text only, used for nothing else |
| black | tree branches |

No genus was given green or coral. Teal (*Comamonas*) is the nearest hue to the presence
green and is clearly lighter and bluer; flagged here for the record.

**Canvas** 180 × 288.75 mm (111 tips × 2.25 mm row pitch + 13 mm top + 26 mm for the scale bar
and legend). The row pitch is set by the 6.5 pt tip labels; anything tighter collides.

---

## Fig 2 — module completeness by genus, and per-gene prevalence

**Sources** `Table_S1a` (presence), `Table_S1d` (genus), `Table_S15a` (group flag).

**Recomputed and asserted**
- hero set from the group flag == `st.HEROES`; group sizes 6 and 105 asserted against the
  arrays actually used for the two prevalence series.
- prevalences are means of the binarised matrix, computed at draw time; every bar carries its
  own computed value label.
- genus counts in the y tick labels are computed from the plotted grouping.

**Panel A was turned horizontal.** With ten genera the rotated x labels collided
(`Chryseobacterium` × `Sphingobacterium`, and eight more). Genus names now run along the y
axis, which also puts the highest-scoring genus at the top.

**Discrepancies against the shipped figure and legend**

| item | legend / old figure | this rebuild (from Table S1a) |
|---|---|---|
| MICP-complete prevalence | legend "100 %"; old figure 83.3 % | 100 % for all eight genes |
| rest prevalence | legend "55–85 %"; old figure 28–35 % | 28–37 % |
| rest group size | legend n = 105; old figure n = 94 | n = 105 |

The legend's "55–85 %" matches neither the old figure nor any source and should be corrected
to the drawn values.

**Colour meanings** coral = the MICP-complete group (rings in A, bars in B); grey = every other
MAG; light grey box fill + black medians and points = distribution summary, no group meaning.

**Canvas** 180 × 78 mm.

---

## Fig 3 — ure cluster synteny on the main contig of each MICP-complete MAG

**Sources** `MAGs_FASTA_files/bakta_results/<MAG>/<MAG>.gff3` (CDS coordinates, strand,
product), `Table_S1c_hero_cluster_audit.csv` (contig of record, gene count, span),
`Table_S3a_HGT_ureCah_cluster.csv` (mobile-element count, Δ GC).

**Fixed from the previous rebuild.** The `figure_master_v1` version drew S13/S16 from an inline
schematic constant and showed only *ureA* + *ureB*. Every arrow here comes from the GFF3, and
all seven *ure* genes appear on a single contig for M1, S13 and S16, as the manuscript states.

**Recomputed and asserted, per MAG**
- number of *distinct* `ure` genes on the contig named in Table S1c == the table's
  `ure_genes_on_main_contig`;
- cluster span (max end − min start, kb) == the table's `cluster_span_kb_main` to < 0.02 kb.

Both assertions pass for all six MAGs: M1 7/28.56 kb, S13 7/5.86, S16 7/5.92, C22 5/5.32,
S23 4/3.27, S26 4/3.03.

**Discrepancies found**
1. **The audit column counts distinct genes, not CDS copies.** On S23 `contig_220` there are
   five `ure` CDS (two *ureE* copies) but four distinct genes; the table records 4. The stat
   column on the figure reads "4 / 7" and both *ureE* arrows are drawn, so the reader sees the
   duplication. Worth one clause in the legend.
2. **The old densest-cluster search disagrees with the table of record for S23.** The search in
   `01_main_figures.py` ties at four distinct `ure` genes between `contig_151` (2.85 kb) and
   `contig_220` (3.27 kb) and returns `contig_151`; Table S1c names `contig_220`. This rebuild
   draws the contig named in the table. The shipped figure drew `contig_151`, so its S23 track
   does not match Table S1c.
3. **`cah` is on a separate contig in all six MAGs** (`cah_on_main_contig = False` throughout
   Table S1c), so no *cah* arrow exists on any track. The gene key therefore omits *cah* —
   a key entry with no mark on the page would be misleading. The legend should say explicitly
   that *cah* lies outside the drawn window.

**Per-track scaling, not a shared axis.** M1's window is 32 kb while the other five are 7–10 kb.
On one shared scale the five short clusters occupy the leftmost quarter of the page and their
gene labels are unreadable. Each track now carries its own kb axis with its own ticks, so
absolute scale stays readable per track; arrows remain to scale within a track.

**Label placement** gene labels are staggered over up to three rows by a greedy placer with a
leader line to the arrow, because M1's seven genes sit in two tight sub-clusters inside a 32 kb
window.

**Flagged literal, justified.** `check_no_hardcoded.py` flags `0.353` in `place_levels`
(`build_fig3.py` L140). It is the typographic constant 1 pt = 0.3528 mm, used to convert a
label's point size into millimetres for the collision estimate. It is not a measurement.

**Colour meanings** one colour per MICP gene (`ureA`…`ureG`), light grey for any other CDS.
No other colour is used on the page; the MAG identifiers are plain black bold and the contig
identifiers grey, carrying no group meaning.

**Canvas** 180 × 129 mm.
