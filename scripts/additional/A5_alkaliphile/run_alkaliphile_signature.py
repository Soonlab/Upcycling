#!/usr/bin/env python3
"""A5: Alkaliphile signature quantification.

For each of 111 MAGs:
  - Proteome isoelectric point (pI) distribution  -> mean/median, acidic fraction (<5)
  - Mrp / nha* / antiporter gene counts from Bakta annotation
  - GC content, total CDS
Compare MICP-complete (6) vs rest (105) with Mann-Whitney U.
"""
from __future__ import annotations
import os, re, glob, json
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from scipy.stats import mannwhitneyu

BAKTA = "/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
OUTDIR = "/data/data/Upcycling/research/additional/A5_alkaliphile"
os.makedirs(OUTDIR, exist_ok=True)

HERO = {"S13", "S16", "S23", "C22", "M1", "S26"}

# Alkaliphile / pH-adaptation gene name patterns (Bakta 'Name' or 'Product' column)
ANTI_PATTERNS = {
    "Mrp": re.compile(r"\b(mrpA|mrpB|mrpC|mrpD|mrpE|mrpF|mrpG|mnh[A-G]|multi(ple)? resistance and pH adaptation)\b", re.I),
    "Nha":  re.compile(r"\b(nhaA|nhaB|nhaC|nhaD|na\+/h\+ antiporter)\b", re.I),
    "ChaA": re.compile(r"\bchaA\b", re.I),
    "NapA": re.compile(r"\bnapA\b.*antiport", re.I),
    "NhaP": re.compile(r"\bnhaP\b", re.I),
    "Kef":  re.compile(r"\bkef[ABCG]\b", re.I),
    "Mdh":  re.compile(r"\bmdh\b.*antiport", re.I),
}

def pi_stats(faa_path: str) -> dict:
    pIs, lens, acidic_counts = [], [], 0
    for rec in SeqIO.parse(faa_path, "fasta"):
        seq = str(rec.seq).replace("*", "").replace("X", "")
        if len(seq) < 20:
            continue
        try:
            pi = IsoelectricPoint(seq).pi()
        except Exception:
            continue
        pIs.append(pi)
        lens.append(len(seq))
        if pi < 5.0:
            acidic_counts += 1
    if not pIs:
        return {}
    arr = np.asarray(pIs)
    return {
        "n_cds": len(arr),
        "pI_mean": float(arr.mean()),
        "pI_median": float(np.median(arr)),
        "pI_std": float(arr.std()),
        "pI_acidic_frac": acidic_counts / len(arr),     # fraction with pI < 5
        "pI_basic_frac":  float((arr > 9).sum() / len(arr)),
        "pI_skew_acidic_over_basic": float((arr < 5).sum() / max((arr > 9).sum(), 1)),
        "cds_mean_len": float(np.mean(lens)),
    }

def antiporter_counts(tsv_path: str) -> dict:
    counts = {k: 0 for k in ANTI_PATTERNS}
    try:
        df = pd.read_csv(tsv_path, sep="\t", comment="#", header=None,
                         names=["contig","type","start","end","strand","locus","gene","product","dbxrefs"])
    except Exception:
        return counts
    text = (df["gene"].fillna("") + " || " + df["product"].fillna("")).tolist()
    for k, pat in ANTI_PATTERNS.items():
        counts[k] = sum(1 for t in text if pat.search(t))
    return counts

def main():
    rows = []
    faas = sorted(glob.glob(f"{BAKTA}/*/*.faa"))
    # exclude hypotheticals.faa
    faas = [f for f in faas if "hypotheticals" not in os.path.basename(f)]
    print(f"[A5] processing {len(faas)} MAG proteomes")
    for faa in faas:
        mag = os.path.basename(faa).replace(".faa", "")
        tsv = os.path.join(os.path.dirname(faa), f"{mag}.tsv")
        row = {"MAG": mag, "group": "MICP_complete" if mag in HERO else "rest"}
        row.update(pi_stats(faa))
        row.update({f"{k}_count": v for k, v in antiporter_counts(tsv).items()})
        row["antiporter_total"] = sum(antiporter_counts(tsv).values())
        rows.append(row)
        if len(rows) % 20 == 0:
            print(f"  ... {len(rows)}/{len(faas)}")
    df = pd.DataFrame(rows).sort_values("MAG")
    df.to_csv(f"{OUTDIR}/alkaliphile_signature_per_MAG.csv", index=False)

    # group-wise stats
    hero = df[df.group == "MICP_complete"]
    rest = df[df.group == "rest"]
    metrics = ["pI_mean","pI_median","pI_acidic_frac","pI_basic_frac",
               "pI_skew_acidic_over_basic","Mrp_count","Nha_count","antiporter_total"]
    stat_rows = []
    for m in metrics:
        h, r = hero[m].dropna().values, rest[m].dropna().values
        if len(h) == 0 or len(r) == 0:
            continue
        try:
            u, p = mannwhitneyu(h, r, alternative="two-sided")
        except Exception:
            u, p = np.nan, np.nan
        stat_rows.append({
            "metric": m,
            "hero_n": len(h), "rest_n": len(r),
            "hero_mean": float(np.mean(h)), "rest_mean": float(np.mean(r)),
            "hero_median": float(np.median(h)), "rest_median": float(np.median(r)),
            "MWU_U": float(u), "MWU_p": float(p),
            "fold_change_mean": float(np.mean(h)/max(np.mean(r),1e-9)),
        })
    pd.DataFrame(stat_rows).to_csv(f"{OUTDIR}/alkaliphile_MICP_vs_rest_stats.csv", index=False)
    print("\n[A5] group comparison:")
    print(pd.DataFrame(stat_rows).to_string(index=False))
    print(f"\n[A5] outputs: {OUTDIR}/alkaliphile_signature_per_MAG.csv")

if __name__ == "__main__":
    main()
