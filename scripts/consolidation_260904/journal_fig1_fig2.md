# Build note — consolidated main Fig 1 and Fig 2

Date 2026-09-04 · builders `/data/data/Upcycling/new_figure/build_v2_fig1.py`,
`/data/data/Upcycling/new_figure/build_v2_fig2.py` ·
output `/data/data/Upcycling/new_figure/figures_v2/` (SVG + PDF + PNG@200 dpi) ·
interpreter `/home/soon/miniconda3/envs/dram_env/bin/python`.

No existing file was modified. `_style.py`, `_micp_presence.py` and every `build_fig*.py` /
`build_supS*.py` are untouched; the two new scripts import `_style as st` and
`_micp_presence.presence()` exactly as the sources did.

## Panel mapping

### Fig 1 — phylogenomic distribution and completeness of the MICP module
| new panel | old source | content |
|---|---|---|
| A | `build_fig1.py` panel A | bac120 ML tree (midpoint-rooted), genus colour strip, MAG id + GTDB species |
| B | `build_fig1.py` panel B | *ureA–G* + *cah* presence matrix, same rows, same row pitch |
| C | `build_fig2.py` panel A | MICP module score (0–8) per GTDB genus, box + jittered points, six MICP-complete overplotted as coral rings |
| D | `build_fig2.py` panel B | per-gene prevalence, MICP-complete (n = 6) vs rest (n = 105) |

A and B stay side by side and row-aligned as in old Fig 1; C and D form a second row.
`presence()` is called **once** for the whole page, so A/B and C/D cannot drift apart.

### Fig 2 — architecture and dosage of the ure–cah cluster
| new panel | old source | content |
|---|---|---|
| A | `build_fig3.py` (whole figure) | six synteny tracks of the *ure* cluster, per-row stat columns (ure genes / MGE / ΔGC) |
| B | `build_supS14.py` | MICP pathway gene-dosage heat map, six MAGs × 16 genes, urease-core and Ca-pathway stat columns |
| C | `build_supS15.py` panel **B only** | per-MAG urease/CA contig vs geNomad MGE cross-check table |

Old S15 panel A (geNomad plasmid/virus contigs per MAG, box plots) was deliberately **not**
carried over — it belongs to the consolidated supplementary page built separately.

## Page geometry and acceptance gate (final, after the 235 mm re-lay-out)

| page | width | height | `st.audit` | `st.prose_scan` |
|---|---|---|---|---|
| `figures_v2/Fig1.*` | 180.0 mm | **227.0 mm** | 194 text items · overlaps 0 · outside 0 · axes overlaps 0 → **PASS** | no sentence-like text → **PASS** |
| `figures_v2/Fig2.*` | 180.0 mm | **226.0 mm** | 304 text items · overlaps 0 · outside 0 · axes overlaps 0 → **PASS** | no sentence-like text → **PASS** |

Both scripts now carry `assert H <= 235.0`, so the ceiling cannot be lost in a later edit.

### Which lever was used

**Fig 1: 362.8 → 227.0 mm — lever (3), the two-column split, plus a smaller tip-label font.**
Levers (1) and (2) were evaluated first and rejected on arithmetic:
* (1) pitch alone. Non-tree overhead (letters, key, C/D row, C's axis label and ring key) is
  ~108 mm even after compression, leaving ≤ 127 mm for 111 tips = 1.15 mm per row, i.e. a
  ~3.2 pt tip label. Not legible; rejected.
* (2) label only the six MICP-complete MAGs. This *would* fit, but it removes 105 MAG
  identifiers from the page. Not used — **no MAG identifier was moved to Table S1**, so the
  legend does **not** need a "remaining identifiers in Table S1 sheet S1.1" sentence.
* (3) used. All 111 tips are kept, with their full `MAG  GTDB species` labels. The tip
  order is carried in two adjacent tree columns at an identical row pitch (2.05 mm) and an
  identical branch-length scale; one scale bar under the left tree applies to both.
  The break is **not arbitrary**: the script searches tip boundaries 50–62 for the one
  crossed by the fewest clades. It picks **tip 58** (only **4** clades crossed, columns of
  58 and 53 tips), which falls at a near-natural break in the topology, so only four
  backbone edges run off a column edge. Tip labels dropped from 6.5 pt to **5.5 pt**
  (`FS_TIP`); this is the *only* type-scale change, and it touches tip labels alone — axis
  titles (8 pt), panel letters (10 pt), body/tick (7 pt) and stat annotations (6.5 pt) are
  unchanged, as is the palette.
  Panel letters **A** and **B** sit over the left column; the right column is a
  continuation of the same two element types and carries no letters, so the letters still
  read A, B, C, D in order. **The legend must say that panel A/B continue in the
  right-hand column of the page.**
  One extra guard was added because `st.audit` cannot see this failure mode: a tip label
  that overran its slot would land on the presence matrix, where there is no other text to
  collide with. The script therefore measures the widest rendered tip label and asserts it
  fits the label column — currently **37.0 mm of the 40.0 mm** slot.

**Fig 2: 255.0 → 226.0 mm — row-height compression of A and B only.**
Synteny track height 13.0 → **11.5 mm** and inter-track gap 6.0 → **4.5 mm** (A);
dosage heat-map height 28 → **22 mm** (B); cross-check table height 38 → **34 mm** (C);
inter-panel gaps trimmed. Panel C was **not** moved beside B: after the colour bar and the
two stat columns, panel B's block already reaches ~152 mm across, leaving under 30 mm —
not enough for a six-column table whose narrowest header (`plasmid / contigs`) needs ~8 mm
per column. Nothing was dropped, no font was changed on this page, and all six MAGs and
all 16 dosage columns remain.

Every number on both pages is loaded or recomputed from the repository sources; all asserts
of the five source scripts were carried over verbatim (tree tip count and identity, group
flag == `st.HEROES`, 5-of-6 full module + S26 `ureB == 0`, group sizes 6/105, genus-strip
membership, per-track `ure_genes_on_main_contig` and `cluster_span_kb_main`, recomputed
`Ca_pathway` == stored flag, S17a/S17b plasmid-count agreement), plus the new
split-completeness assert (58 + 53 = 111), the label-width assert and the page-height assert.

## Discrepancies noticed between a source table and what the old scripts drew

1. **Table S1a vs the Bakta annotations (carried over, now resolved on disk).** The shipped
   45-MAG `Table_S1a_ace_samples_list.csv` counted the Bakta `5_ureB_sRNA` non-coding
   feature as a *ureB* copy and so scored **S26 8/8**; the CDS-only scan gives 7/8. The copy
   of S1a now on disk is the 2026-09-04 re-export (111 MAGs, CDS-only) and the build log
   confirms it **agrees with the scan**. Two MAGs change when non-coding features are
   excluded: **S11 7 → 6** and **S26 8 → 7**. 27 MAGs carry all eight genes CDS-only
   (28 before the exclusion).
2. **"MICP-complete" vs the urease core for S26 — an unresolved cross-table inconsistency.**
   `Table_S15a` labels S26 `MICP_complete`, but `Table_S1c` records
   `ureABC_all_present = False` and `ureABCDEFG_all_present = False` for S26, and
   `Table_S15b` records `urease_core_complete = 0` **and** `Ca_pathway = 0` (the only MAG
   with either). The consolidated pages state the fact in three places — Fig 1B white *ureB*
   cell, Fig 1D *ureB* prevalence 83 % rather than 100 %, Fig 2B `ureB = 0` with
   `urease core = 0` — but the group *name* still says complete. This is a manuscript
   wording issue, not a drawing error.
3. **S23 main contig (carried over).** The densest-cluster search in
   `scripts/01_main_figures.py` picks `contig_151` for S23, while `Table_S1c` names
   `contig_220` as `main_contig`. The rebuild draws the table of record (contig_220).
   On that track five *ure* arrows are drawn but the stat column reads **4 / 7**, because
   S1c counts *distinct* ure genes and *ureE* occurs twice on the contig.
4. **cah is absent from the Fig 2A gene key by design.** `Table_S1c` records
   `cah_on_main_contig = False` for all six MAGs, so no *cah* arrow exists on any track; a
   key entry would have no mark. *cah* is nevertheless one of the eight genes coloured in
   Fig 1B/1D, so the two pages use the gene name with different scopes (genome-wide
   presence vs on-contig synteny).
5. **Empty flag columns in `Table_S17b`.** `urease_on_plasmid`, `urease_on_virus` and
   `CA_on_virus` are empty for every MAG and carry no information, so they are not drawn.
   `CA_on_plasmid` is non-empty for S26 alone and is represented by that MAG's
   `CA_MGE_contamination = 1` (the single coral cell in Fig 2C).
6. **New check added this session.** The old `build_supS15.py` asserted only that S17a and
   S17b agree on `n_plasmid_contigs`. `n_virus_contigs` was verified to agree as well for
   all six MAGs — no discrepancy, but the check is now explicit in the data block.
7. **`pangenome_work/MICP_Pangenome_Final_Summary.csv` remains unusable** (100 of 111 MAGs;
   C1 and C10–C19 missing, previously drawn as all-absent). Not read by either new script.
