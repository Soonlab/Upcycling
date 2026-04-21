#!/usr/bin/env python3
"""
C3 v3 - restrict codeml M0 + yn00 to tractable subset.

Strategy:
  - Per gene, take all 6 heroes + up to 12 non-hero representatives (stratified per source).
  - Fit codeml M0 on that subset (16-20 seqs -> tractable).
  - Run yn00 pairwise on same subset; split pairs into hero-hero, hero-rest, rest-rest.
  - Report median omega per bucket + MWU(hero-hero vs rest-rest).
Reads full codon MSAs produced by run_dnds_v2.py; re-runs codeml and adds yn00.
"""
import os, re, subprocess, random
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy.stats import mannwhitneyu

random.seed(0)
BAKTA = Path("/data/data/Upcycling/MAGs_FASTA_files/bakta_results")
OUT = Path("/data/data/Upcycling/research/additional/C3_dnds_codon")
HEROES = ["S13","S16","S23","C22","M1","S26"]
PREFIX_SOURCE = {"C":"cattle","S":"swine","M":"sheep","V":"poultry"}

TARGETS = ["ureA","ureB","ureC","ureG"]

def subset_ids(all_ids, max_nonhero=12):
    heroes = [i for i in all_ids if i in HEROES]
    nonheroes = [i for i in all_ids if i not in HEROES]
    by_src = defaultdict(list)
    for i in nonheroes:
        by_src[PREFIX_SOURCE.get(i[0], "other")].append(i)
    chosen = []
    k_per = max(1, max_nonhero // max(1, len(by_src)))
    for src, ids in by_src.items():
        random.shuffle(ids)
        chosen.extend(ids[:k_per])
    chosen = chosen[:max_nonhero]
    return heroes + chosen

def fasta_subset(in_path, out_path, keep_ids):
    keep = set(keep_ids)
    recs = [r for r in SeqIO.parse(str(in_path), "fasta") if r.id in keep]
    with open(out_path, "w") as fh:
        for r in recs:
            fh.write(f">{r.id}\n{str(r.seq)}\n")
    return [r.id for r in recs]

def write_phy(fasta, phy):
    recs = list(SeqIO.parse(str(fasta), "fasta"))
    L = len(recs[0].seq)
    with open(phy, "w") as fh:
        fh.write(f" {len(recs)} {L}\n")
        for r in recs:
            fh.write(f"{r.id[:30].ljust(31)}{str(r.seq)}\n")
    return len(recs), L

def codeml_M0(label, sub_fasta, work):
    work.mkdir(parents=True, exist_ok=True)
    phy = work/f"{label}.phy"
    n, L = write_phy(sub_fasta, phy)
    if n < 4:
        return None
    tree = work/f"{label}.tre"
    ids = [r.id for r in SeqIO.parse(str(sub_fasta),"fasta")]
    tree.write_text("(" + ",".join(ids) + ");")
    ctl = work/"codeml.ctl"
    ctl.write_text(
        f"seqfile = {phy.name}\n"
        f"treefile = {tree.name}\n"
        f"outfile = {label}.M0.out\n"
        "runmode = 0\nseqtype = 1\nCodonFreq = 2\nmodel = 0\nNSsites = 0\nicode = 0\n"
        "fix_kappa = 0\nkappa = 2\nfix_omega = 0\nomega = 0.4\nclock = 0\n"
    )
    r = subprocess.run("codeml codeml.ctl", shell=True, cwd=str(work),
                       capture_output=True, text=True, timeout=1800)
    out_f = work/f"{label}.M0.out"
    if not out_f.exists():
        return None
    t = out_f.read_text()
    g_omega = re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", t)
    g_dN = re.search(r"tree length for dN:\s+([\d.]+)", t)
    g_dS = re.search(r"tree length for dS:\s+([\d.]+)", t)
    g_ln = re.search(r"lnL\(.*?\):\s*([-\d.]+)", t)
    g_k = re.search(r"kappa \(ts/tv\)\s*=\s*([\d.]+)", t)
    return {
        "gene": label, "n_seqs": n, "aln_len": L,
        "omega_M0": float(g_omega.group(1)) if g_omega else None,
        "tree_dN": float(g_dN.group(1)) if g_dN else None,
        "tree_dS": float(g_dS.group(1)) if g_dS else None,
        "lnL": float(g_ln.group(1)) if g_ln else None,
        "kappa": float(g_k.group(1)) if g_k else None,
    }

def yn00(label, sub_fasta, work):
    work.mkdir(parents=True, exist_ok=True)
    phy = work/f"{label}.yn00.phy"
    n, L = write_phy(sub_fasta, phy)
    if n < 2:
        return []
    ctl = work/"yn00.ctl"
    ctl.write_text(
        f"seqfile = {phy.name}\n"
        f"outfile = {label}.yn00.out\n"
        "verbose = 1\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n"
    )
    r = subprocess.run("yn00 yn00.ctl", shell=True, cwd=str(work),
                       capture_output=True, text=True, timeout=1200)
    out_f = work/f"{label}.yn00.out"
    if not out_f.exists():
        return []
    t = out_f.read_text()
    rows = []
    m_block = re.search(r"seq\. seq\..*?\n(.*?)\n\s*\n", t, flags=re.S)
    if not m_block:
        return rows
    block = m_block.group(1)
    for line in block.splitlines():
        parts = line.split()
        if len(parts) < 10:
            continue
        try:
            a, b = parts[0], parts[1]
            idx_a = int(a); idx_b = int(b)
            S = float(parts[2]); N = float(parts[3])
            t_ = float(parts[4]); kappa = float(parts[5]); omega = float(parts[6])
            dN = float(parts[7].rstrip("+-")); dN_SE = parts[8] if len(parts) > 8 else ""
            dS_field_idx = 9
            if "+-" in parts[8]:
                dS_field_idx = 9
            try:
                dS = float(parts[dS_field_idx].rstrip("+-"))
            except Exception:
                dS = None
            rows.append({"i": idx_a, "j": idx_b, "omega": omega, "dN": dN, "dS": dS})
        except Exception:
            continue
    return rows

rows_M0 = []
rows_yn = []

for label in TARGETS:
    codon_aln = OUT/"per_gene"/f"{label}.codon.aln.fasta"
    if not codon_aln.exists() or codon_aln.stat().st_size == 0:
        continue
    all_ids = [r.id for r in SeqIO.parse(str(codon_aln),"fasta")]
    keep = subset_ids(all_ids, max_nonhero=12)
    sub_fa = OUT/"per_gene"/f"{label}.codon.sub.fasta"
    present = fasta_subset(codon_aln, sub_fa, keep)
    print(f"[{label}] subset n={len(present)} (heroes present: {[x for x in present if x in HEROES]})", flush=True)
    w = OUT/"codeml_v3"/label
    try:
        res = codeml_M0(label, sub_fa, w)
        if res:
            rows_M0.append(res)
            print(f"[codeml v3] {label}: omega={res['omega_M0']} dN={res['tree_dN']} dS={res['tree_dS']}", flush=True)
    except subprocess.TimeoutExpired:
        print(f"[codeml v3] {label} TIMEOUT", flush=True)
    w2 = OUT/"yn00"/label
    try:
        pairs = yn00(label, sub_fa, w2)
        for p in pairs:
            i, j = p["i"] - 1, p["j"] - 1
            a = present[i] if i < len(present) else None
            b = present[j] if j < len(present) else None
            if a and b:
                rows_yn.append({"gene": label, "a": a, "b": b, **p})
        print(f"[yn00] {label}: {len(pairs)} pairs", flush=True)
    except subprocess.TimeoutExpired:
        print(f"[yn00] {label} TIMEOUT", flush=True)

pd.DataFrame(rows_M0).to_csv(OUT/"codeml_M0_summary.csv", index=False)
print(f"wrote codeml_M0_summary.csv: {len(rows_M0)} genes")

yn_df = pd.DataFrame(rows_yn)
yn_df.to_csv(OUT/"yn00_pairwise.csv", index=False)
print(f"wrote yn00_pairwise.csv: {yn_df.shape}")

if not yn_df.empty:
    yn_df = yn_df[(yn_df["omega"] < 99) & (yn_df["omega"] > 0) & (yn_df["dS"] > 0.01)]
    yn_df["bucket"] = yn_df.apply(
        lambda r: "hero_hero" if r.a in HEROES and r.b in HEROES
        else ("hero_rest" if r.a in HEROES or r.b in HEROES else "rest_rest"),
        axis=1,
    )
    summary = []
    for gene in TARGETS:
        gdf = yn_df[yn_df.gene == gene]
        if gdf.empty:
            continue
        hh = gdf[gdf.bucket == "hero_hero"]["omega"]
        rr = gdf[gdf.bucket == "rest_rest"]["omega"]
        hr = gdf[gdf.bucket == "hero_rest"]["omega"]
        row = {
            "gene": gene,
            "hero_hero_n": len(hh), "hero_hero_median": float(hh.median()) if len(hh) else None,
            "rest_rest_n": len(rr), "rest_rest_median": float(rr.median()) if len(rr) else None,
            "hero_rest_n": len(hr), "hero_rest_median": float(hr.median()) if len(hr) else None,
        }
        if len(hh) >= 2 and len(rr) >= 2:
            u, p = mannwhitneyu(hh, rr, alternative="two-sided")
            row["MWU_hh_vs_rr_p"] = float(p)
        summary.append(row)
    pd.DataFrame(summary).to_csv(OUT/"yn00_hero_vs_rest_summary.csv", index=False)
    print("wrote yn00_hero_vs_rest_summary.csv")
    print(pd.DataFrame(summary).to_string(index=False))
