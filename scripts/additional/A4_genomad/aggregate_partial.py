#!/usr/bin/env python3
"""Aggregate A4 geNomad results across completed MAGs (handles partial run)."""
import os, glob
import pandas as pd

WD = "/data/data/Upcycling/research/additional/A4_genomad/results"
OUT = "/data/data/Upcycling/research/additional/A4_genomad"
HERO = {"S13","S16","S23","C22","M1","S26"}

rows = []
for d in sorted(os.listdir(WD)):
    sdir = os.path.join(WD, d, f"{d}_summary")
    if not os.path.isdir(sdir): continue
    pfn = f"{sdir}/{d}_plasmid_summary.tsv"
    vfn = f"{sdir}/{d}_virus_summary.tsv"
    if not (os.path.exists(pfn) and os.path.exists(vfn)): continue
    try:
        np_ = pd.read_csv(pfn, sep="\t").shape[0]
        nv_ = pd.read_csv(vfn, sep="\t").shape[0]
    except Exception:
        continue
    rows.append({
        "MAG": d,
        "group": "MICP_complete" if d in HERO else "rest",
        "n_plasmid_contigs": np_,
        "n_virus_contigs": nv_,
    })

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/genomad_summary.csv", index=False)
print(f"[A4 agg] {len(df)} MAGs aggregated")
print(df.groupby("group").agg(n=("MAG","count"),
                              mean_plas=("n_plasmid_contigs","mean"),
                              med_plas=("n_plasmid_contigs","median"),
                              max_plas=("n_plasmid_contigs","max"),
                              mean_vir=("n_virus_contigs","mean")))
print("\n[A4 agg] MICP-complete MAGs:")
print(df[df.group=="MICP_complete"][["MAG","n_plasmid_contigs","n_virus_contigs"]].to_string(index=False))
