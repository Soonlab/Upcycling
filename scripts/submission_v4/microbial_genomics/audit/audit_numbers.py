"""Numeric audit - SUBMISSION_v4 copy (2026-09-23).  Every headline number in the manuscript is recomputed
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
import os
MAN_DIR = Path(os.environ.get("UPCYCLING_MAN_DIR", Path(__file__).resolve().parent.parent))
MAN = MAN_DIR / "01_Manuscript.md"
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


print("\n== data provenance (SUBMISSION_v4, 2026-09-23) ==")
prov = pd.read_csv(MAN_DIR / "data_provenance_260923/MAG_provenance.tsv", sep="\t")
ck("MAGs in the deposited metadata", 111, len(prov))
oc = prov.origin.value_counts()
ck("MAGs from manure", 44, int(oc["Manure"]))
ck("MAGs from manured soil", 38, int(oc["Manured soil"]))
ck("MAGs from soil", 29, int(oc["Soil"]))
pm = prov.groupby(prov.MAG.str[0]).binning_method.agg(lambda s: "/".join(sorted(set(s))))
ck("prefix C = COMEbin", "COMEbin", pm["C"])
ck("prefix M = MaxBin/MetaBAT", "MaxBin/MetaBAT", pm["M"])
ck("prefix S = SemiBin2", "SemiBin2", pm["S"])
ck("prefix V = AVAMB", "AVAMB", pm["V"])
po = prov.set_index("MAG").origin
ck("M1 and S26 from manure", True, bool((po[["M1", "S26"]] == "Manure").all()))
ck("S13, S16, S23, C22 from manured soil", True, bool((po[["S13", "S16", "S23", "C22"]] == "Manured soil").all()))
pe_ = prov.set_index("MAG").experiment
ck("the four Sphingobacterium candidates share one soil series", 1, int(pe_[["S13", "S16", "S23", "C22"]].nunique()))
ps = prov.set_index("MAG").original_sample
ck("S13 and S16 from the same sample", True, bool(ps["S13"] == ps["S16"]))
ck("that sample is a day-28 sample", 28, int(prov.set_index("MAG").sampling_day["S13"]))
cov_h = prov[prov.MAG.isin(hero)].genome_coverage_x; cov_r = prov[~prov.MAG.isin(hero)].genome_coverage_x
ck("read-mapping coverage median, MICP-complete (x)", 9.5, round(float(cov_h.median()), 1), 0.05)
ck("read-mapping coverage median, rest (x)", 17.6, round(float(cov_r.median()), 1), 0.05)
ck("read-mapping coverage MWU P", 0.12, round(float(mannwhitneyu(cov_h, cov_r, alternative="two-sided").pvalue), 2), 0.005)
s1d = rd("Table_S1d_GTDB_Tk_classification.tsv")
phy = s1d.classification.str.extract(r"p__([^;]*)")[0].value_counts()
ck("Pseudomonadota MAGs", 99, int(phy["Pseudomonadota"]))
ck("Bacteroidota MAGs", 12, int(phy["Bacteroidota"]))
gen = s1d.classification.str.extract(r"g__([^;]*)")[0].replace("", np.nan).dropna()
ck("genera in the panel", 17, int(gen.nunique()))
aai = rd("Table_S4b_AAI_S13_S16.csv")
ck("S13-S16 AAI (%)", 82.5, round(float(aai[(aai.Query == "S13") & (aai.Target == "S16")].AAI.iloc[0]), 1), 0.05)
says("PRJNA1231077", "text cites the source BioProject")
says("10.5281/zenodo.15309541", "text cites the source Zenodo record")
# MGen copy: numbered references -> find each paper's number in the list and require that number in a citation bracket
import re as _re
_refs = man.split("## References")[1].split("## Table 1")[0]
_pre = man.split("## References")[0]
def _cited(doi):
    m = _re.search(r"^(\d+)\. .*" + _re.escape(doi), _refs, _re.M)
    if not m:
        return False
    n = int(m.group(1))
    for g in _re.findall(r"\[(\d+(?:\s*[,–]\s*\d+)*)\]", _pre):
        for part in g.split(","):
            a, _, b = part.strip().partition("–")
            if int(a) <= n <= int(b or a):
                return True
    return False
checked += 2
ok_dp, ok_cs = _cited("10.1016/j.dib.2025.111748"), _cited("10.1093/femsec/fiad148")
print(("PASS" if ok_dp else "FAIL") + "  text cites the data paper (Pérez-Valera and Elhottová 2025, numbered)")
print(("PASS" if ok_cs else "FAIL") + "  text cites the companion study (Sardar et al. 2023, numbered)")
if not (ok_dp and ok_cs):
    fails.append("data-paper citations")
for bad in ("swine", "poultry", "MGnify", "ACE hybrid", "PRJNA-XXXXXXX", "livestock source", "waste source"):
    checked += 1
    ok = bad not in man.split("## References")[0]
    print(f"{'PASS' if ok else 'FAIL'}  old framing absent: {bad!r}")
    if not ok:
        fails.append("old framing: " + bad)

print("\n== single-contig ure operons ==")
# the shipped S3b covers the six candidates only; the panel-wide count comes from the definition audit
# (rule C: Bakta GFF3, ureA-G on one shared contig, evaluated over all 111 MAGs)
da = open("/data/data/Upcycling/SUBMISSION/_revision_260904/DEFINITION_AUDIT.md").read()
mC = re.search(r"^\| C \| Bakta GFF3: \*ureA–G\* on one shared contig \| (\d+) \|", da, re.M)
ck("MAGs with all seven ure genes on one contig (definition audit, rule C)", 26, int(mC.group(1)) if mC else -1)

print("\n== novelty of S13 / S16 (revised text) ==")
aai = rd("Table_S4b_AAI_S13_S16.csv")
ck("S13 closest congeneric AAI (%)", 93.15, round(float(aai[aai.Query == "S13"].AAI.max()), 2), 0.01)
ck("S16 closest congeneric AAI (%)", 93.49, round(float(aai[aai.Query == "S16"].AAI.max()), 2), 0.01)
ext = rd("Table_S10b_ext_Sphingobacterium_novelty.csv").set_index("MAG")
for mag, v in (("S13", 94.57), ("S16", 93.85), ("S23", 98.96), ("C22", 99.16)):
    ck(f"{mag} max ANI to RefSeq Sphingobacterium (%)", v, round(float(ext.loc[mag, "Nearest_ANI"]), 2), 0.01)
ck("S13 nearest RefSeq species", "detergens", ext.loc["S13", "Nearest_organism"].split()[-1])
ck("S16 nearest RefSeq species", "multivorum", ext.loc["S16", "Nearest_organism"].split()[-1])
ck("unresolved MAGs (%)", 18.9, round(100 * 21 / 111, 1), 0.05)

print("\n== gene-tree congruence ==")
rf = open(SUPP / "Table_S7a_RF_distance.txt").read()
ck("normalised RF distance", 0.58, round(float(re.search(r"normalized_RF = ([\d.]+)", rf).group(1)), 2), 0.005)
iq = open(SUPP / "Table_S7f_SH_AU_test.iqtree").read()
rows = [[t for t in l.split() if t not in "+-"] for l in iq.splitlines() if re.match(r"\s*2\s+-\d", l)]
# columns after dropping the +/- markers: tree, logL, deltaL, bp-RELL, p-KH, p-SH, p-WKH, p-WSH, c-ELW, p-AU
ck("SH test rejects the species-tree topology (p-SH < 0.001)", True, bool(rows) and float(rows[0][5]) < 0.001)

print("\n== yn00 ureG partition (revised text) ==")
yn = rd("Table_S19c_yn00_hero_vs_rest_summary.csv").set_index("gene")
ck("ureG within-candidate median omega", 0.31, round(float(yn.loc["ureG", "hero_hero_median"]), 2), 0.005)
ck("ureG within-rest median omega", 0.074, round(float(yn.loc["ureG", "rest_rest_median"]), 3), 0.0005)
ck("ureG MWU P", 7.7e-8, float(f"{yn.loc['ureG', 'MWU_hh_vs_rr_p']:.1e}"), 1e-9)
for g in ("ureA", "ureB", "ureC"):
    ck(f"{g} MWU P not significant", True, bool(yn.loc[g, "MWU_hh_vs_rr_p"] > 0.05))

print("\n== gRodon P (revised text) ==")
gr2 = rd("Table_S16_gRodon_growth_rates_per_MAG.csv")
gh2 = gr2[gr2.group != "rest"]; gr_2 = gr2[gr2.group == "rest"]
ck("gRodon MWU P", 0.58, round(float(mannwhitneyu(gh2.d_hours, gr_2.d_hours).pvalue), 2), 0.01)

print("\n== text-only assertions (revised wording) ==")
says("27 MAGs encoded all seven *ure* genes plus at least one protein-coding *cah*", "text states the 27-MAG module count")
says("2 of 6 vs. 3 of 105; 33.3% vs. 2.9%", "text states Mrp as a prevalence, not a dosage")
says("with 42/42 matches", "text states the 42/42 active-site result")
says("S26 lacked protein-coding *ureB*", "text states the S26 ureB gap")

print()
print(f"{checked} checks run")
print("ALL NUMBERS AGREE" if not fails else f"{len(fails)} DISAGREEMENT(S): " + "; ".join(fails))
sys.exit(1 if fails else 0)
