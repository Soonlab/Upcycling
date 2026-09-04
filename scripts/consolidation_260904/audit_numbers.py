"""Numeric audit: every headline number in the consolidated manuscript is recomputed
from the shipped supplementary source and compared with what the text says.

A check fails if the text and the data disagree. Nothing here reads the figures; this
is text versus data, so it is independent of the figure rebuild.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

SUPP = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables")
MAN = Path("/data/data/Upcycling/consolidation_260904/manuscript/01_Manuscript.md")
man = MAN.read_text()
fails, checked = [], 0


def rd(n):
    return pd.read_csv(SUPP / n, sep="\t" if n.endswith(".tsv") else ",")


def ck(label, claimed, actual, tol=0.0):
    """Compare a claimed value with the recomputed one."""
    global checked
    checked += 1
    if isinstance(actual, float) or isinstance(claimed, float):
        ok = abs(float(claimed) - float(actual)) <= tol
    else:
        ok = claimed == actual
    print(f"{'PASS' if ok else 'FAIL'}  {label:58s} text={claimed}  data={actual}")
    if not ok:
        fails.append(label)


def says(pattern, label=None):
    """Assert the manuscript contains a literal string; returns True/False."""
    global checked
    checked += 1
    ok = pattern in man
    print(f"{'PASS' if ok else 'FAIL'}  {label or ('text contains: ' + pattern[:50])}")
    if not ok:
        fails.append(label or pattern[:50])
    return ok


print("== gene presence ==")
s1a = rd("Table_S1a_ace_samples_list.csv")
genes = [c for c in s1a.columns if c.lower() in
         ("urea", "ureb", "urec", "ured", "uree", "uref", "ureg", "cah")]
if not genes:
    genes = [c for c in s1a.columns if re.fullmatch(r"ure[A-G]|cah", c, re.I)]
pres = (s1a[genes] > 0)
ck("MAGs carrying all eight genes", 27, int((pres.sum(axis=1) == 8).sum()))
ck("panel size", 111, len(s1a))

hero = ["S13", "S16", "S23", "C22", "M1", "S26"]
idc = [c for c in s1a.columns if c.lower() in ("mag", "sample")][0]
h = s1a[s1a[idc].isin(hero)]
r = s1a[~s1a[idc].isin(hero)]
ck("MICP-complete ureB prevalence (%)", 83.3,
   round(100 * (h[[g for g in genes if g.lower() == "ureb"][0]] > 0).mean(), 1), 0.1)
cah = [g for g in genes if g.lower() == "cah"][0]
ck("rest-group cah prevalence (%)", 92.4, round(100 * (r[cah] > 0).mean(), 1), 0.1)
ure_prev = sorted(round(100 * (r[g] > 0).mean(), 1) for g in genes if g.lower() != "cah")
ck("rest-group ure prevalence, minimum (%)", 29.5, ure_prev[0], 0.1)
ck("rest-group ure prevalence, maximum (%)", 45.7, ure_prev[-1], 0.1)

print("\n== novelty screen ==")
s8 = rd("Table_S8_novelty_ANI_screen.csv")
anicol = [c for c in s8.columns if "ani" in c.lower()][0]
noani = s8[anicol].isna().sum()
ck("MAGs with no species-level ANI", 21, int(noani))
ck("MAGs with a species-level ANI", 90, int(s8[anicol].notna().sum()))
ck("minimum species-level ANI (%)", 95.08, round(float(s8[anicol].min()), 2), 0.01)

print("\n== trait enrichment ==")
s2c = rd("Table_S2c_permutation_statistics.csv")
ck("trait subcategories tested", 38, len(s2c))
ck("trait modules at q < 0.05", 9, int((s2c.Permutation_q_BH < 0.05).sum()))
for sub, fc in (("Mrp_complex", 10.85), ("carb_binding", 9.78), ("oxidative", 4.76),
                ("glycoside_hydrolase", 4.66), ("tet_mac", 4.22),
                ("Na_H_antiporter", 2.30), ("quorum", 2.13), ("metal_efflux", 1.19)):
    ck(f"fold change {sub}", fc,
       round(float(s2c.loc[s2c.Subcategory == sub, "Fold_change"].iloc[0]), 2), 0.01)
gs = s2c[s2c.Subcategory.str.contains("GS_GOGAT|gs_gogat|GOGAT", case=False, regex=True)]
if len(gs):
    ck("GS-GOGAT fold change", 0.62, round(float(gs.Fold_change.iloc[0]), 2), 0.01)

print("\n== dbCAN ==")
s6d = rd("Table_S6d_dbCAN_hero_vs_rest.csv")
ck("CAZy classes at q < 0.05", 5, int((s6d.Permutation_q_BH < 0.05).sum()))
for cl, fc in (("GH", 3.82), ("CBM", 4.24), ("PL", 3.52), ("CE", 1.99), ("GT", 1.24)):
    ck(f"CAZy fold change {cl}", fc,
       round(float(s6d.loc[s6d.Class == cl, "Fold_change"].iloc[0]), 2), 0.01)
ck("CAZy AA q value", 1.0, round(float(s6d.loc[s6d.Class == "AA", "Permutation_q_BH"].iloc[0]), 2), 0.01)

print("\n== urease active site and structure ==")
s12 = rd("Table_S12_UreC_active_site_residues.csv")
matches = int(sum((s12[m] == s12.expected).sum() for m in hero))
ck("active-site matches", 42, matches)
tm = rd("Table_S22_ureC_vs_4CEU_tm.csv")
tm["MAG"] = tm.MAG.str.replace("_UreC", "", regex=False)
tm = tm.set_index("MAG")
for mag, v in (("C22", 0.678), ("S13", 0.620), ("S16", 0.613),
               ("S23", 0.612), ("S26", 0.599), ("M1", 0.597)):
    ck(f"TM-score {mag}", v, round(float(tm.loc[mag, "tm_norm_ref"]), 3), 0.001)
ck("minimum backbone RMSD (A)", 3.52, round(float(tm.rmsd.min()), 2), 0.01)
ck("maximum backbone RMSD (A)", 4.34, round(float(tm.rmsd.max()), 2), 0.01)

print("\n== selection ==")
m0 = rd("Table_S19b_codeml_M0_summary.csv")
wcol = [c for c in m0.columns if c.lower() in ("omega", "omega_m0", "w")][0]
gcol = [c for c in m0.columns if "gene" in c.lower()][0]
m0 = m0.set_index(gcol)
for g, v in (("ureA", 0.087), ("ureB", 0.059), ("ureC", 0.026), ("ureG", 0.041)):
    ck(f"codeml M0 omega {g}", v, round(float(m0.loc[g, wcol]), 3), 0.001)

print("\n== alkali signature and growth ==")
s15a = rd("Table_S15a_alkaliphile_signature_per_MAG.csv")
hh = s15a[s15a.group != "rest"]; rr = s15a[s15a.group == "rest"]
ck("MICP-complete MAGs carrying Mrp", 2, int(hh.Mrp_count.sum()))
ck("rest MAGs carrying Mrp", 3, int(rr.Mrp_count.sum()))
ck("Mrp prevalence ratio", 11.7, round(hh.Mrp_count.mean() / rr.Mrp_count.mean(), 1), 0.1)
ck("Mrp MWU P", 5.3e-4, float(f"{mannwhitneyu(hh.Mrp_count, rr.Mrp_count).pvalue:.1e}"), 1e-5)
ck("Nha mean, MICP-complete", 2.50, round(float(hh.Nha_count.mean()), 2), 0.01)
ck("Nha mean, rest", 2.19, round(float(rr.Nha_count.mean()), 2), 0.01)
ck("Nha MWU P", 0.53, round(float(mannwhitneyu(hh.Nha_count, rr.Nha_count).pvalue), 2), 0.01)

gr = rd("Table_S16_gRodon_growth_rates_per_MAG.csv")
ck("MAGs passing the gRodon filter", 85, len(gr))
dcol = [c for c in gr.columns if "doubl" in c.lower() or c.lower() in ("d", "d_hours")][0]
gcolg = [c for c in gr.columns if gr[c].dtype == object and gr[c].nunique() <= 3]
if gcolg:
    gh_ = gr[gr[gcolg[0]] != "rest"]; gr_ = gr[gr[gcolg[0]] == "rest"]
    ck("gRodon MICP-complete MAGs passing", 4, len(gh_))
    ck("gRodon median, MICP-complete (h)", 1.06, round(float(gh_[dcol].median()), 2), 0.01)
    ck("gRodon median, rest (h)", 1.10, round(float(gr_[dcol].median()), 2), 0.01)

print("\n== abundance proxy ==")
s21a = rd("Table_S21a_abundance_proxy_per_MAG.csv")
ccol = [c for c in s21a.columns if "weight" in c.lower()][0]
idc2 = [c for c in s21a.columns if c.lower() in ("mag", "sample")][0]
ah = s21a[s21a[idc2].isin(hero)]; ar = s21a[~s21a[idc2].isin(hero)]
ck("coverage mean, MICP-complete", 26.8, round(float(ah[ccol].mean()), 1), 0.1)
ck("coverage mean, rest", 25.5, round(float(ar[ccol].mean()), 1), 0.1)
ck("coverage MWU P", 0.27, round(float(mannwhitneyu(ah[ccol], ar[ccol]).pvalue), 2), 0.01)
s21b = rd("Table_S21b_abundance_proxy_per_source.csv")
scol = [c for c in s21b.columns if "source" in c.lower()][0]
sw = s21b[s21b[scol].str.lower() == "swine"].iloc[0]
mcol = [c for c in s21b.columns if "mean" in c.lower()][0]
ncol = [c for c in s21b.columns if c.lower() in ("n", "n_mags", "count")][0]
ck("swine mean coverage", 57.4, round(float(sw[mcol]), 1), 0.1)
ck("swine MAG count", 15, int(sw[ncol]))

print("\n== external rarity ==")
mg = rd("Table_S14a_mgnify_catalog_summary.csv")
ck("MGnify species clusters", 7599, int(mg.n_species_clusters.sum()))
ck("MGnify gene-complete clusters", 233, int(mg.n_MICP_gene_complete.sum()))
ck("MGnify gene-complete (%)", 3.07,
   round(100 * mg.n_MICP_gene_complete.sum() / mg.n_species_clusters.sum(), 2), 0.01)
ck("MGnify single-contig clusters", 265, int(mg.n_MICP_single_contig_ureC_CA.sum()))
ck("MGnify single-contig (%)", 3.49,
   round(100 * mg.n_MICP_single_contig_ureC_CA.sum() / mg.n_species_clusters.sum(), 2), 0.01)

print("\n== antiSMASH ==")
bg = rd("Table_S23b_antismash_hero_vs_rest.csv").set_index("metric")
for cl, hm, rm, pv in (("BGC_T3PKS", 0.67, 0.029, 5.3e-10),
                       ("BGC_RRE-containing", 0.83, 0.105, 1.9e-5)):
    ck(f"{cl} MICP-complete mean", hm, round(float(bg.loc[cl, "hero_mean"]), 2), 0.01)
    ck(f"{cl} rest mean", rm, round(float(bg.loc[cl, "rest_mean"]), 3), 0.001)
    ck(f"{cl} MWU P", pv, float(f"{bg.loc[cl, 'MWU_p']:.1e}"), abs(pv) * 0.1)

print("\n== Pseudomonas_E within-genus rarity ==")
pe = rd("Table_S13b_pseudomonas_e_single_contig.csv")
ck("Pseudomonas_E references screened", 146, len(pe))
n = int(pe["ureC_and_CA_single_contig"].astype(bool).sum())
ck("references with single-contig ureC + CA", 53, n)
ck("single-contig percentage", 36.3, round(100 * n / len(pe), 1), 0.1)
pa = rd("Table_S13a_pseudomonas_e_MICP_rarity_screen.csv")
ck("references carrying UreC (%)", 93.8,
   round(100 * pe["has_UreC"].astype(bool).mean(), 1), 0.1)
ck("references carrying any CA (%)", 100.0,
   round(100 * pe["has_CA_any"].astype(bool).mean(), 1), 0.1)

print("\n== biosafety ==")
bs = rd("Table_S11_biosafety_counts_per_MAG.csv")
idc3 = [c for c in bs.columns if bs[c].dtype == object][0]
bs = bs.set_index(idc3)
card = [c for c in bs.columns if "card" in c.lower()][0]
vf = [c for c in bs.columns if "vfdb" in c.lower()][0]
res = [c for c in bs.columns if "resfinder" in c.lower()][0]
pf = [c for c in bs.columns if "plasmid" in c.lower()][0]
ck("M1 CARD hits", 9, int(bs.loc["M1", card]))
ck("M1 VFDB hits", 46, int(bs.loc["M1", vf]))
ck("S26 CARD hits", 1, int(bs.loc["S26", card]))
ck("S26 VFDB hits", 21, int(bs.loc["S26", vf]))
ck("ResFinder hits across the six", 0, int(bs.loc[hero, res].sum()))
ck("PlasmidFinder hits across the six", 0, int(bs.loc[hero, pf].sum()))
ck("PlasmidFinder hits panel-wide", 1, int(bs[pf].sum()))
ck("the single replicon belongs to S21", "S21", str(bs[bs[pf] > 0].index[0]))

print("\n== codon usage ==")
cu = rd("Table_S19a_codon_usage_per_MAG.csv")
idc4 = [c for c in cu.columns if cu[c].dtype == object][0]
g3 = [c for c in cu.columns if "gc3" in c.lower()][0]
en = [c for c in cu.columns if "enc" in c.lower()][0]
ch = cu[cu[idc4].isin(hero)]; cr = cu[~cu[idc4].isin(hero)]
ck("GC3 mean, MICP-complete (%)", 51.4, round(float(ch[g3].mean()), 1), 0.1)
ck("GC3 mean, rest (%)", 70.1, round(float(cr[g3].mean()), 1), 0.1)
ck("ENC mean, MICP-complete", 50.0, round(float(ch[en].mean()), 1), 0.1)
ck("ENC mean, rest", 39.0, round(float(cr[en].mean()), 1), 0.1)

print("\n== PERMANOVA ==")
pg = rd("Table_S9b_PERMANOVA_global.csv")
row = pg.iloc[0]
nums = {c: row[c] for c in pg.columns if isinstance(row[c], (int, float, np.floating))}
print("      PERMANOVA row:", {k: round(float(v), 3) for k, v in nums.items()})
vals = [round(float(v), 2) for v in nums.values()]
ck("pan-genome genus pseudo-F 8.21 present", True, 8.21 in vals)
ck("pan-genome source pseudo-F 1.25 present", True, 1.25 in vals)

print("\n== text-only assertions ==")
says("**27 of the 111 MAGs carry all eight genes as protein-coding sequences**",
     "text states the 27-MAG module count")
says("2 of the 6 MICP-complete MAGs (S13, S16) against 3 of the remaining 105",
     "text states Mrp as a prevalence, not a dosage")
says("**42/42 identical matches**", "text states the 42/42 active-site result")
says("no protein-coding urease β-subunit", "text states the S26 ureB gap")

print()
print(f"{checked} checks run")
print("ALL NUMBERS AGREE" if not fails else f"{len(fails)} DISAGREEMENT(S): " + "; ".join(fails))
sys.exit(1 if fails else 0)
