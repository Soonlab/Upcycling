#!/usr/bin/env python3
"""A4 cross-check: is the ureCah cluster on any contig flagged by geNomad as plasmid/virus?

For each hero MAG:
  1. Find all contigs carrying ure[ABCDEFGH] or cah genes (from Bakta TSV)
  2. Compare with geNomad plasmid_summary and virus_summary contig lists
  3. Report any overlap (any overlap = cluster on an MGE-predicted contig = problem for the vertical-inheritance claim)
"""
import os, re
import pandas as pd

BAKTA = "/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
A4 = "/data/data/Upcycling/research/additional/A4_genomad/results"
OUT = "/data/data/Upcycling/research/additional/A4_genomad"
HERO = ["S13","S16","S23","C22","M1","S26"]

def urease_contigs(mag):
    """Return contigs carrying urease *core* subunits (ureA/B/C) only — strict."""
    tsv = f"{BAKTA}/{mag}/{mag}.tsv"
    if not os.path.exists(tsv): return set()
    df = pd.read_csv(tsv, sep="\t", comment="#", header=None,
                     names=["contig","type","start","end","strand","locus","gene","product","dbxrefs"])
    text = df["gene"].fillna("") + " || " + df["product"].fillna("")
    # Strict: ureA/B/C subunit only (drop urease accessories, urea transporters,
    # urea amidolyase, 5'-ureB-sRNA, CA elsewhere on genome)
    mask = text.str.contains(
        r"\b(ureA|ureB|ureC)\b|urease subunit alpha|urease subunit beta|urease subunit gamma|urease.*subunit",
        case=False, na=False, regex=True)
    # exclude sRNA rows (type == ncRNA)
    if "type" in df.columns:
        mask = mask & (df["type"].astype(str).str.lower() == "cds")
    return set(df[mask]["contig"].unique())

def ca_contigs(mag):
    """Return contigs carrying any carbonic anhydrase (all classes)."""
    tsv = f"{BAKTA}/{mag}/{mag}.tsv"
    if not os.path.exists(tsv): return set()
    df = pd.read_csv(tsv, sep="\t", comment="#", header=None,
                     names=["contig","type","start","end","strand","locus","gene","product","dbxrefs"])
    text = df["gene"].fillna("") + " || " + df["product"].fillna("")
    mask = text.str.contains(r"carbonic anhydrase|\bcah\b|\bcan[AB]\b|\bcynT\b",
                             case=False, na=False, regex=True)
    return set(df[mask]["contig"].unique())

def mge_contigs(mag):
    sdir = f"{A4}/{mag}/{mag}_summary"
    plas = set()
    vir  = set()
    if os.path.exists(f"{sdir}/{mag}_plasmid_summary.tsv"):
        d = pd.read_csv(f"{sdir}/{mag}_plasmid_summary.tsv", sep="\t")
        if "seq_name" in d.columns and len(d):
            plas = set(d["seq_name"])
    if os.path.exists(f"{sdir}/{mag}_virus_summary.tsv"):
        d = pd.read_csv(f"{sdir}/{mag}_virus_summary.tsv", sep="\t")
        if "seq_name" in d.columns and len(d):
            # geNomad virus seq_name may include coordinates like contig|NNN-MMM
            vir = set(s.split("|")[0] for s in d["seq_name"])
    return plas, vir

print(f"{'MAG':<5} {'ureABC_contigs':<35} {'CA_contigs':<30} {'n_pl':<5} {'n_vir':<5} {'ure∩pl':<12} {'ure∩vir':<9} {'CA∩pl':<12}")
print("-"*120)
rows = []
for mag in HERO:
    uc = urease_contigs(mag)     # ureA/B/C core only
    cc = ca_contigs(mag)          # CA genes (on any contig)
    if not os.path.exists(f"{A4}/{mag}/{mag}_summary"):
        print(f"{mag:<5} {', '.join(sorted(uc))[:33]:<35} {', '.join(sorted(cc))[:28]:<30} geNomad pending")
        continue
    pl, vr = mge_contigs(mag)
    ov_upl = uc & pl
    ov_uvr = uc & vr
    ov_cpl = cc & pl
    ov_cvr = cc & vr
    print(f"{mag:<5} {', '.join(sorted(uc))[:33]:<35} {', '.join(sorted(cc))[:28]:<30} "
          f"{len(pl):<5} {len(vr):<5} {','.join(sorted(ov_upl)) or '—':<12} "
          f"{','.join(sorted(ov_uvr)) or '—':<9} {','.join(sorted(ov_cpl)) or '—':<12}")
    rows.append({
        "MAG": mag,
        "urease_core_contigs": ",".join(sorted(uc)),
        "CA_contigs": ",".join(sorted(cc)),
        "n_plasmid_contigs": len(pl),
        "n_virus_contigs": len(vr),
        "urease_on_plasmid": ",".join(sorted(ov_upl)),
        "urease_on_virus": ",".join(sorted(ov_uvr)),
        "CA_on_plasmid": ",".join(sorted(ov_cpl)),
        "CA_on_virus": ",".join(sorted(ov_cvr)),
        "urease_core_MGE_contamination": int(len(ov_upl)+len(ov_uvr) > 0),
        "CA_MGE_contamination": int(len(ov_cpl)+len(ov_cvr) > 0),
    })

if rows:
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/ureCah_vs_MGE_overlap.csv", index=False)
    print(f"\n→ {OUT}/ureCah_vs_MGE_overlap.csv")
    n_ure = df.urease_core_MGE_contamination.sum()
    n_ca = df.CA_MGE_contamination.sum()
    print(f"\n[A4] Strict urease core (ureA/B/C) MGE overlap: {n_ure}/{len(df)}")
    print(f"[A4] CA (any class) MGE overlap:                {n_ca}/{len(df)}")
    if n_ure == 0:
        print("[A4]  ⇒ Urease core verticals / no MGE contamination → manuscript claim holds.")
