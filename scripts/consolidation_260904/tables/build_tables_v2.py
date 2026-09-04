"""Consolidate 23 supplementary tables (52 files) into three workbooks.

Design: /data/data/Upcycling/consolidation_260904/DESIGN.md
  Table S1  one row per MAG, one sheet per assay
  Table S2  group comparisons and statistics
  Table S3  matrices, reference panels, trees and method support

Nothing is edited: every sheet is the shipped file read verbatim. The only additions
are a README sheet per workbook and, in S3, two sheets that carry the content of the
small plain-text tree artefacts plus a manifest of the two artefacts too large or too
unstructured to be a sheet.
"""
import gzip
import shutil
from pathlib import Path

import pandas as pd

SRC = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
OUT = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables_v2")
OUT.mkdir(parents=True, exist_ok=True)

def rd(name, **kw):
    p = SRC / name
    sep = "\t" if p.suffix == ".tsv" else ","
    if p.suffix == ".gz":
        with gzip.open(p, "rt") as fh:
            return pd.read_csv(fh, **kw)
    return pd.read_csv(p, sep=sep, **kw)

# ---------------------------------------------------------------- Table S1
S1 = [
    ("S1.1_MAG_stats_taxonomy",   "Table_S1d_GTDB_Tk_classification.tsv",
     "GTDB-Tk r220 classification, ANI to closest reference and assembly statistics, all 111 MAGs"),
    ("S1.2_MICP_gene_presence",   "Table_S1a_ace_samples_list.csv",
     "Protein-coding copy number of ureA-G and cah per MAG (CDS only; re-exported 2026-09-04)"),
    ("S1.3_trait_counts",         "Table_S2a_trait_module_counts.csv",
     "Raw keyword hits per trait subcategory and total CDS per MAG"),
    ("S1.4_trait_per1kCDS",       "Table_S2b_trait_module_per1kCDS.csv",
     "The same 38 trait subcategories normalised to hits per 10^3 CDS"),
    ("S1.5_dbCAN_class_counts",   "Table_S6a_dbCAN_class_counts.csv",
     "dbCAN v12 HMMER class-level CAZyme counts per MAG"),
    ("S1.6_dbCAN_class_per1kCDS", "Table_S6b_dbCAN_class_per1kCDS.csv",
     "dbCAN class counts normalised to hits per 10^3 CDS"),
    ("S1.7_dbCAN_family_counts",  "Table_S6c_dbCAN_family_counts.csv",
     "dbCAN family-level counts per MAG (re-parsed 2026-09-04 to include the six MICP-complete MAGs)"),
    ("S1.8_novelty_ANI",          "Table_S8_novelty_ANI_screen.csv",
     "ANI to the closest GTDB reference genome and the species-boundary call per MAG"),
    ("S1.9_biosafety",            "Table_S11_biosafety_counts_per_MAG.csv",
     "abricate hits per MAG against CARD, VFDB, ResFinder and PlasmidFinder"),
    ("S1.10_alkaliphile",         "Table_S15a_alkaliphile_signature_per_MAG.csv",
     "Mrp and Nha antiporter copy number, proteome median pI and acidic-pI fraction per MAG"),
    ("S1.11_stoichiometry",       "Table_S15b_stoichiometry_per_MAG.csv",
     "MICP pathway gene dosage per MAG (ureB re-exported CDS-only 2026-09-04)"),
    ("S1.12_growth_gRodon",       "Table_S16_gRodon_growth_rates_per_MAG.csv",
     "gRodon2 predicted minimum doubling time, 85 MAGs passing the ribosomal-anchor filter"),
    ("S1.13_geNomad_MGE",         "Table_S17a_genomad_summary_per_MAG.csv",
     "geNomad plasmid- and virus-flagged contig counts per MAG"),
    ("S1.14_defense_CRISPR",      "Table_S18a_defense_per_MAG.csv",
     "DefenseFinder anti-phage systems and minced CRISPR arrays per MAG"),
    ("S1.15_codon_usage",         "Table_S19a_codon_usage_per_MAG.csv",
     "Genome-wide GC3 and effective number of codons per MAG"),
    ("S1.16_abundance_proxy",     "Table_S21a_abundance_proxy_per_MAG.csv",
     "Length-weighted mean SPAdes contig coverage per MAG (source labels corrected 2026-09-04)"),
    ("S1.17_antiSMASH_BGC",       "Table_S23a_antismash_per_MAG.csv",
     "antiSMASH 7 biosynthetic gene cluster regions per MAG, total and by class"),
    ("S1.18_DRAM_product",        "Table_S5a_DRAM_product.tsv",
     "DRAM v1.5 module completeness per MAG (re-exported from the intact full-panel distillate)"),
    ("S1.19_DRAM_genome_stats",   "Table_S5c_DRAM_genome_stats.tsv",
     "DRAM per-genome assembly and annotation statistics"),
    ("S1.20_PCoA_coordinates",    "Table_S9a_PCoA_coordinates.csv",
     "Pan-genome PCoA coordinates PC1-PC3 with waste source, genus and MICP-complete flag"),
]

# ---------------------------------------------------------------- Table S2
S2 = [
    ("S2.1_permutation_traits",     "Table_S2c_permutation_statistics.csv",
     "10,000-iteration one-sided permutation test with BH-FDR for all 38 trait subcategories"),
    ("S2.2_dbCAN_enrichment",       "Table_S6d_dbCAN_hero_vs_rest.csv",
     "dbCAN class-level enrichment, MICP-complete versus rest, with bootstrap CI and BH-q"),
    ("S2.3_PERMANOVA_global",       "Table_S9b_PERMANOVA_global.csv",
     "PERMANOVA of pan-genome and trait-module distance by waste source and by genus"),
    ("S2.4_PERMANOVA_pairwise",     "Table_S9c_PERMANOVA_pairwise.csv",
     "Pairwise PERMANOVA between waste sources"),
    ("S2.5_PseudomonasE_rarity",    "Table_S13a_pseudomonas_e_MICP_rarity_screen.csv",
     "Pfam screen of 146 NCBI Pseudomonas_E reference genomes for MICP gene content"),
    ("S2.6_PseudomonasE_contig",    "Table_S13b_pseudomonas_e_single_contig.csv",
     "Single-contig ureC + carbonic anhydrase architecture across the same 146 genomes"),
    ("S2.7_MGnify_catalog_summary", "Table_S14a_mgnify_catalog_summary.csv",
     "MICP prevalence per MGnify livestock catalog and pooled over 7,599 species clusters"),
    ("S2.8_MGnify_cluster_profile", "Table_S14b_mgnify_species_cluster_profiles.csv.gz",
     "Per-species-cluster MICP gene profile underlying the catalog summary"),
    ("S2.9_defense_group_test",     "Table_S18b_defense_hero_vs_rest.csv",
     "Group comparison of anti-phage and CRISPR counts, MICP-complete versus rest"),
    ("S2.10_codeml_M0",             "Table_S19b_codeml_M0_summary.csv",
     "PAML codeml M0 omega per urease subunit"),
    ("S2.11_yn00_pairwise",         "Table_S19c_yn00_hero_vs_rest_summary.csv",
     "yn00 pairwise omega by gene and pair class, within-shortlist versus within-rest"),
    ("S2.12_abundance_by_source",   "Table_S21b_abundance_proxy_per_source.csv",
     "Coverage proxy aggregated by waste source (labels corrected 2026-09-04)"),
    ("S2.13_antiSMASH_group",       "Table_S23b_antismash_hero_vs_rest.csv",
     "antiSMASH class means and Mann-Whitney P, MICP-complete versus rest"),
]

# ---------------------------------------------------------------- Table S3
S3 = [
    ("S3.1_cluster_audit",       "Table_S1c_hero_cluster_audit.csv",
     "Contig of record, ure gene count and cluster span for each MICP-complete MAG"),
    ("S3.2_MIGS_lite",           "Table_S1b_MIGS_lite_Sphingobacterium.csv",
     "MIGS-lite descriptors for the candidate novel Sphingobacterium species"),
    ("S3.3_ure_cah_cluster_MGE", "Table_S3a_HGT_ureCah_cluster.csv",
     "Flanking mobile-element annotation and regional GC deviation of the ure-cah cluster"),
    ("S3.4_gene_contig_distrib", "Table_S3b_gene_contig_distribution.csv",
     "Distribution of each MICP gene across the contigs of each MICP-complete MAG"),
    ("S3.5_ANI_matrix_111",      "Table_S4a_skani_ANI_matrix_111MAGs.csv",
     "skani whole-genome ANI matrix across all 111 MAGs"),
    ("S3.6_AAI_S13_S16",         "Table_S4b_AAI_S13_S16.csv",
     "mmseqs2 reciprocal-best-hit AAI of S13 and S16 against the in-panel Sphingobacterium"),
    ("S3.7_ext_Sphingo_ANI",     "Table_S10a_ext_Sphingobacterium_ANI_matrix.csv",
     "ANI of the six study Sphingobacterium MAGs against 63 RefSeq genomes of the genus"),
    ("S3.8_ext_Sphingo_novelty", "Table_S10b_ext_Sphingobacterium_novelty.csv",
     "Closest RefSeq match and novelty call from the 63-genome screen"),
    ("S3.9_UreC_active_site",    "Table_S12_UreC_active_site_residues.csv",
     "Observed residue at each of the seven canonical urease-alpha catalytic positions"),
    ("S3.10_ureC_vs_4CEU_TM",    "Table_S22_ureC_vs_4CEU_tm.csv",
     "ESMFold TM-score and backbone RMSD of each UreC against PDB 4CEU chain C"),
    ("S3.11_MICP_ref_panel_ANI", "Table_S20_skani_hero_vs_refs.tsv",
     "skani ANI against the organism-verified MICP reference panel (panel rebuilt 2026-09-04)"),
    ("S3.12_ureC_tree_samples",  "Table_S7b_ureC_samples.csv",
     "The 46 MAGs contributing a UreC sequence to the gene tree"),
    ("S3.13_ureC_MGE_overlap",   "Table_S17b_ureCah_vs_MGE_overlap.csv",
     "Per-MAG overlap between urease/CA-bearing contigs and geNomad-flagged contigs"),
]

TEXT_ARTEFACTS = [
    ("Table_S7a_RF_distance.txt",      "Robinson-Foulds distance, ureC gene tree versus species tree"),
    ("Table_S7e_hero_topology.txt",    "Placement of the MICP-complete tips within the ureC gene tree"),
    ("Table_S7c_ureC_gene_tree.newick","Maximum-likelihood ureC gene tree, Newick"),
    ("Table_S7d_species_tree_pruned.newick","bac120 species tree pruned to the 46 UreC-bearing MAGs, Newick"),
]

DEPOSITED = [
    ("Table_S5b_DRAM_metabolism_summary.xlsx", 8,
     "Full DRAM v1.5 metabolism summary, eight tool-generated sheets; deposited with the "
     "repository rather than reproduced here because it is a raw tool export"),
    ("Table_S7f_SH_AU_test.iqtree", 1,
     "IQ-TREE report of the SH and AU topology tests; a program report, not a table"),
]


def readme(rows, title):
    return pd.DataFrame({"Sheet": [r[0] for r in rows],
                         "Source file (as shipped 2026-09-04)": [r[1] for r in rows],
                         "Contents": [r[2] for r in rows]})


def write(path, rows, extra=None):
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        readme(rows, path.stem).to_excel(xl, sheet_name="README", index=False)
        for sheet, fname, _ in rows:
            df = rd(fname)
            assert len(sheet) <= 31, sheet
            df.to_excel(xl, sheet_name=sheet, index=False)
            print(f"  {sheet:32s} {df.shape[0]:6d} x {df.shape[1]:3d}  <- {fname}")
        if extra:
            for sheet, df in extra:
                df.to_excel(xl, sheet_name=sheet, index=False)
                print(f"  {sheet:32s} {df.shape[0]:6d} x {df.shape[1]:3d}  <- assembled")
    print(f"wrote {path}")


print("Table S1 - per-MAG measurements")
write(OUT / "Table_S1_per_MAG_measurements.xlsx", S1)

print("\nTable S2 - comparative statistics")
write(OUT / "Table_S2_comparative_statistics.xlsx", S2)

print("\nTable S3 - reference panels and method support")
txt = pd.DataFrame(
    [{"Artefact": n, "Contents": d, "Text": (SRC / n).read_text().strip()}
     for n, d in TEXT_ARTEFACTS])
dep = pd.DataFrame(
    [{"Deposited file": n, "Sheets / parts": k, "Why it is not a sheet here": d}
     for n, k, d in DEPOSITED])
write(OUT / "Table_S3_reference_panels_and_methods.xlsx", S3,
      extra=[("S3.14_tree_artefacts", txt), ("S3.15_deposited_files", dep)])

# also ship the two deposited artefacts alongside the workbooks
for n, _, _ in DEPOSITED:
    shutil.copy2(SRC / n, OUT / n)
    print(f"copied {n}")
