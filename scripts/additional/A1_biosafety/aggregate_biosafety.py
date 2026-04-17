#!/usr/bin/env python3
"""A1 biosafety aggregation: summarize abricate hits across databases per MAG."""
import os, glob
import pandas as pd

WD = "/data/data/Upcycling/research/additional/A1_biosafety"
HERO = {"S13","S16","S23","C22","M1","S26"}
DBS = ["card","vfdb","resfinder","plasmidfinder"]

# read combined files, if plasmidfinder missing that's ok
per_mag = {}
for db in DBS:
    comb = f"{WD}/combined/{db}_all.tsv"
    if not os.path.exists(comb):
        continue
    try:
        df = pd.read_csv(comb, sep="\t", low_memory=False)
    except Exception:
        continue
    if "#FILE" not in df.columns:
        continue
    df["MAG"] = df["#FILE"].str.extract(r"bakta_results/([^/]+)/")
    counts = df.groupby("MAG").size()
    per_mag[db] = counts

agg = pd.DataFrame(per_mag).fillna(0).astype(int)
# ensure all MAGs present
bakta_mags = sorted([os.path.basename(d) for d in glob.glob("/data/data/Upcycling/MAGs_FASTA_files/bakta_results/*") if os.path.isdir(d)])
agg = agg.reindex(bakta_mags).fillna(0).astype(int)
agg["group"] = ["MICP_complete" if m in HERO else "rest" for m in agg.index]
agg.to_csv(f"{WD}/biosafety_counts_per_MAG.csv")

print(f"\n[A1] biosafety aggregate across {len(agg)} MAGs:")
print(agg.groupby("group").agg(
    card_mean=("card","mean") if "card" in agg.columns else ("group","count"),
    card_max=("card","max") if "card" in agg.columns else ("group","count"),
    vfdb_mean=("vfdb","mean") if "vfdb" in agg.columns else ("group","count"),
    vfdb_max=("vfdb","max") if "vfdb" in agg.columns else ("group","count"),
    resfinder_mean=("resfinder","mean") if "resfinder" in agg.columns else ("group","count"),
    resfinder_max=("resfinder","max") if "resfinder" in agg.columns else ("group","count"),
    plas_mean=("plasmidfinder","mean") if "plasmidfinder" in agg.columns else ("group","count"),
).to_string())

print("\n[A1] MICP-complete MAGs detail:")
print(agg.loc[[m for m in HERO if m in agg.index]].to_string())

# Per-hero gene-level detail for CARD (top AMR classes)
print("\n[A1] MICP-complete — CARD (AMR) hits per MAG:")
card_df = pd.read_csv(f"{WD}/combined/card_all.tsv", sep="\t", low_memory=False)
card_df["MAG"] = card_df["#FILE"].str.extract(r"bakta_results/([^/]+)/")
for mag in HERO:
    sub = card_df[card_df.MAG == mag]
    if len(sub) == 0:
        print(f"  {mag}: 0 CARD hits")
        continue
    gene_list = sub["GENE"].value_counts().head(10)
    print(f"  {mag}: {len(sub)} CARD hits — top genes:")
    for g, n in gene_list.items():
        prod = sub[sub.GENE==g]["PRODUCT"].iloc[0] if "PRODUCT" in sub.columns else ""
        print(f"     {g} (n={n}): {prod[:80]}")

# VFDB by hero
print("\n[A1] MICP-complete — VFDB (virulence) hits per MAG:")
vf_df = pd.read_csv(f"{WD}/combined/vfdb_all.tsv", sep="\t", low_memory=False)
vf_df["MAG"] = vf_df["#FILE"].str.extract(r"bakta_results/([^/]+)/")
for mag in HERO:
    sub = vf_df[vf_df.MAG == mag]
    cats = sub["GENE"].value_counts().head(5) if len(sub)>0 else pd.Series([], dtype=int)
    print(f"  {mag}: {len(sub)} VFDB hits — top: {', '.join(cats.index[:5])}")

# ResFinder
print("\n[A1] MICP-complete — ResFinder (acquired AMR) per MAG:")
rf_df = pd.read_csv(f"{WD}/combined/resfinder_all.tsv", sep="\t", low_memory=False)
rf_df["MAG"] = rf_df["#FILE"].str.extract(r"bakta_results/([^/]+)/")
for mag in HERO:
    sub = rf_df[rf_df.MAG == mag]
    genes = sub["GENE"].tolist()[:5] if len(sub)>0 else []
    print(f"  {mag}: {len(sub)} ResFinder hits: {', '.join(genes)}")

print("\n[A1] summary file: biosafety_counts_per_MAG.csv")
