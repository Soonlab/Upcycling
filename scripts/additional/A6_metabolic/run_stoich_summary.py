#!/usr/bin/env python3
"""A6 (revised): Stoichiometric summary of MICP pathway from DRAM + Bakta annotations.

Full FBA via CarveMe requires proprietary solver (CPLEX/Gurobi). Instead we:
  - Enumerate MICP reaction stoichiometry (urease, CA, Ca transporters)
  - Report gene-level presence and dosage for each step in 6 MICP-complete MAGs
  - Contrast with 105 other MAGs (group means)

Pathway (textbook, Mobley 1995 / De Muynck 2010):
  1. urea   + H2O  → NH2COOH + NH3   (urease, EC 3.5.1.5)
  2. NH2COOH + H2O → NH3 + H2CO3      (spontaneous / carbamate hydrolysis)
  3. H2CO3  ⇌ H+ + HCO3−  ⇌ 2H+ + CO32-  (carbonic anhydrase, EC 4.2.1.1)
  4. 2NH3 + 2H2O → 2NH4+ + 2OH−      (pH increase)
  5. Ca2+ + CO32- → CaCO3↓           (cation binding + precipitation)
Net: CO(NH2)2 + 2H2O + Ca2+ → 2NH4+ + CaCO3↓ (+ 1OH−)

Substrate cost per 1 mol CaCO3: 1 urea + 2 H2O + 1 Ca2+
Product: 2 NH4+ + 1 CaCO3 + 1 H+ released to MHCO3-/H2CO3 pool (alkalinizing)
"""
import os, re, glob
import pandas as pd

BAKTA = "/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
OUT = "/data/data/Upcycling/research/additional/A6_metabolic"
os.makedirs(OUT, exist_ok=True)
HERO = {"S13", "S16", "S23", "C22", "M1", "S26"}

# gene patterns (gene name OR product text)
PATTERNS = {
    "ureA":   re.compile(r"\bureA\b|urease.*gamma|gamma.*urease", re.I),
    "ureB":   re.compile(r"\bureB\b|urease.*beta|beta.*urease(?!.*gamma)", re.I),
    "ureC":   re.compile(r"\bureC\b|urease.*alpha|alpha.*urease", re.I),
    "ureD_H": re.compile(r"\bure[DH]\b|urease accessory.*D|urease accessory.*H", re.I),
    "ureE":   re.compile(r"\bureE\b", re.I),
    "ureF":   re.compile(r"\bureF\b", re.I),
    "ureG":   re.compile(r"\bureG\b", re.I),
    "cah_alphaCA": re.compile(r"\bcah\b|alpha.carbonic anhydrase|alpha-class carbonic anhydrase", re.I),
    "canA_gammaCA": re.compile(r"\bcan[AB]\b|gamma.carbonic anhydrase|gamma-class carbonic anhydrase", re.I),
    "cynT_betaCA": re.compile(r"\bcynT\b|beta.carbonic anhydrase|beta-class carbonic anhydrase", re.I),
    "CA_generic": re.compile(r"carbonic anhydrase", re.I),
    "Ca_transporter": re.compile(r"Ca2\+ uniport|calcium.ATPase|Ca2\+.transport|cax\b|chaA", re.I),
    "Ca_ATPase": re.compile(r"\bP-type.*calcium|PMC1|PMR1|calcium.transporting ATPase", re.I),
    "CO3_transporter": re.compile(r"carbonate.transport|HCO3|bicarbonate.transport|bic[AB]", re.I),
    "Na_H_antiporter_Mrp": re.compile(r"mrp[A-G]|mnh[A-G]|multi(ple)? resistance and pH", re.I),
    "K_transport": re.compile(r"kdp[ABCDEF]|TrkA|KtrA|KtrB", re.I),
}

rows = []
for mdir in sorted(glob.glob(f"{BAKTA}/*")):
    if not os.path.isdir(mdir): continue
    mag = os.path.basename(mdir)
    tsv = f"{mdir}/{mag}.tsv"
    if not os.path.exists(tsv): continue
    try:
        df = pd.read_csv(tsv, sep="\t", comment="#", header=None,
                         names=["contig","type","start","end","strand","locus","gene","product","dbxrefs"],
                         quoting=3)
    except Exception:
        continue
    text = (df.gene.fillna("") + " || " + df["product"].fillna("")).tolist()
    row = {"MAG": mag, "group": "MICP_complete" if mag in HERO else "rest"}
    for k, pat in PATTERNS.items():
        row[k] = sum(1 for t in text if pat.search(t))
    # Completeness flags
    row["urease_core_complete"] = int(row["ureA"]>=1 and row["ureB"]>=1 and row["ureC"]>=1)
    row["urease_acc_complete"] = int(row["ureD_H"]>=1 and row["ureE"]>=1 and row["ureF"]>=1 and row["ureG"]>=1)
    row["any_CA"] = int(row["CA_generic"] >= 1)
    row["MICP_stoich_complete"] = int(row["urease_core_complete"] and row["any_CA"])
    row["Ca_pathway"] = int(row["Ca_transporter"]+row["Ca_ATPase"] >= 1)
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/stoichiometry_per_MAG.csv", index=False)

print(f"[A6] MICP stoichiometry scan across {len(df)} MAGs\n")

# group comparison
summary = df.groupby("group").agg(
    n_MAGs=("MAG","count"),
    urease_complete_pct=("urease_core_complete", lambda x: round(100*x.mean(),1)),
    CA_pct=("any_CA", lambda x: round(100*x.mean(),1)),
    MICP_stoich_complete_pct=("MICP_stoich_complete", lambda x: round(100*x.mean(),1)),
    Ca_pathway_pct=("Ca_pathway", lambda x: round(100*x.mean(),1)),
    Mrp_mean_copy=("Na_H_antiporter_Mrp","mean"),
    cah_alphaCA_mean=("cah_alphaCA","mean"),
    ureC_mean=("ureC","mean"),
)
print(summary.to_string())
summary.to_csv(f"{OUT}/stoichiometry_group_summary.csv")

# per-hero detail
print("\n[A6] MICP-complete MAGs: per-gene dosage:")
cols = ["MAG","ureA","ureB","ureC","ureD_H","ureE","ureF","ureG",
        "cah_alphaCA","canA_gammaCA","cynT_betaCA",
        "Ca_transporter","Ca_ATPase","Na_H_antiporter_Mrp","MICP_stoich_complete"]
print(df[df.group=="MICP_complete"][cols].to_string(index=False))

# Net stoichiometry block
print("""
[A6] MICP pathway stoichiometry (net):
     CO(NH2)2 + 2 H2O + Ca2+  ->  2 NH4+ + CaCO3 (precipitate) + OH-
     Substrate: 1 urea + 2 water + 1 Ca2+
     Products : 2 NH4+ + 1 CaCO3 + 1 net OH- (alkalinizing)
     Per-MAG complete pathway:""")
for _, r in df[df.group=="MICP_complete"].iterrows():
    ok = "YES" if r.MICP_stoich_complete else "partial"
    print(f"   {r.MAG}: urease={r.urease_core_complete} CA={r.any_CA} Ca={r.Ca_pathway} Mrp={r.Na_H_antiporter_Mrp}  => {ok}")

print("\n[A6] DONE")
