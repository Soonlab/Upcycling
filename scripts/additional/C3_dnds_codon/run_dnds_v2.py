#!/usr/bin/env python3
"""
C3 v2 - use GFF3 (consistent across all MAGs) instead of TSV (4 were clobbered by CarveMe).
Same output files (overwrites v1 outputs).
"""
import os, re, subprocess, json
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

BAKTA_DIR = Path("/data/data/Upcycling/MAGs_FASTA_files/bakta_results")
OUT = Path("/data/data/Upcycling/research/additional/C3_dnds_codon")
OUT.mkdir(parents=True, exist_ok=True)
(OUT/"per_gene").mkdir(exist_ok=True)
(OUT/"codeml").mkdir(exist_ok=True)
HEROES = {"S13","S16","S23","C22","M1","S26"}

TARGETS = {
    "ureA": [r"urease subunit gamma", r"\bureA\b"],
    "ureB": [r"urease subunit beta", r"\bureB\b"],
    "ureC": [r"urease subunit alpha", r"\bureC\b"],
    "ureG": [r"urease accessory protein UreG", r"\bureG\b"],
}

def parse_gff3_cds(mag):
    gff = BAKTA_DIR/mag/f"{mag}.gff3"
    ffn = BAKTA_DIR/mag/f"{mag}.ffn"
    if not (gff.exists() and ffn.exists()):
        return
    seqs = {r.id.split()[0]: str(r.seq).upper() for r in SeqIO.parse(str(ffn), "fasta")}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = {}
            for kv in parts[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = v.strip()
            lt = attrs.get("locus_tag") or attrs.get("ID", "")
            product = attrs.get("product", "")
            if not product:
                continue
            for label, patterns in TARGETS.items():
                if any(re.search(p, product, re.I) for p in patterns):
                    seq = seqs.get(lt)
                    if seq and len(seq) % 3 == 0 and len(seq) >= 200:
                        yield label, lt, seq
                    break

gene_seqs = defaultdict(dict)
mags = sorted([d.name for d in BAKTA_DIR.iterdir() if d.is_dir()])
for mag in mags:
    for label, lt, seq in parse_gff3_cds(mag):
        existing = gene_seqs[label].get(mag)
        if existing is None or len(seq) > len(existing[1]):
            gene_seqs[label][mag] = (lt, seq)

for label in TARGETS:
    d = gene_seqs.get(label, {})
    nt_fa = OUT/"per_gene"/f"{label}.nt.fasta"
    aa_fa = OUT/"per_gene"/f"{label}.aa.fasta"
    with open(nt_fa, "w") as fnt, open(aa_fa, "w") as faa:
        for mag, (lt, nt) in sorted(d.items()):
            try:
                aa = str(Seq(nt).translate(table=11, to_stop=True))
            except Exception:
                continue
            if "*" in aa or len(aa) < 60:
                continue
            fnt.write(f">{mag}\n{nt}\n")
            faa.write(f">{mag}\n{aa}\n")
    n = sum(1 for _ in open(nt_fa)) // 2
    print(f"{label}: n={n} MAGs have sequence", flush=True)

def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw)

def msa_codon(label):
    aa_fa = OUT/"per_gene"/f"{label}.aa.fasta"
    nt_fa = OUT/"per_gene"/f"{label}.nt.fasta"
    aa_aln = OUT/"per_gene"/f"{label}.aa.aln.fasta"
    codon_aln = OUT/"per_gene"/f"{label}.codon.aln.fasta"
    if aa_fa.stat().st_size == 0:
        return False
    sh(f"mafft --auto --thread 4 {aa_fa} > {aa_aln} 2>/dev/null")
    aa_seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(aa_aln), "fasta")}
    nt_seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(nt_fa), "fasta")}
    with open(codon_aln, "w") as fh:
        for sid, aa in aa_seqs.items():
            nt = nt_seqs.get(sid)
            if nt is None:
                continue
            out = []
            i = 0
            for ch in aa:
                if ch == "-":
                    out.append("---")
                else:
                    out.append(nt[i:i+3])
                    i += 3
            fh.write(f">{sid}\n{''.join(out)}\n")
    return True

for label in TARGETS:
    try:
        ok = msa_codon(label)
        print(f"[msa] {label} {'OK' if ok else 'EMPTY'}", flush=True)
    except Exception as e:
        print(f"[msa] {label} FAIL: {e}", flush=True)

codeml_summary = []
for label in TARGETS:
    codon_aln = OUT/"per_gene"/f"{label}.codon.aln.fasta"
    if not codon_aln.exists() or codon_aln.stat().st_size == 0:
        continue
    work = OUT/"codeml"/label
    work.mkdir(exist_ok=True)
    seqs = list(SeqIO.parse(str(codon_aln), "fasta"))
    if len(seqs) < 4:
        continue
    aln_len = len(seqs[0].seq)
    phy = work/f"{label}.phy"
    with open(phy, "w") as fh:
        fh.write(f" {len(seqs)} {aln_len}\n")
        for r in seqs:
            name = r.id[:30].ljust(31)
            fh.write(f"{name}{str(r.seq)}\n")
    tree = work/f"{label}.tre"
    star = "(" + ",".join(r.id for r in seqs) + ");"
    tree.write_text(star)
    ctl = work/"codeml.ctl"
    ctl.write_text(
        f"seqfile = {phy.name}\n"
        f"treefile = {tree.name}\n"
        f"outfile = {label}.M0.out\n"
        "runmode = 0\nseqtype = 1\nCodonFreq = 2\nmodel = 0\nNSsites = 0\nicode = 0\n"
        "fix_kappa = 0\nkappa = 2\nfix_omega = 0\nomega = 0.4\nclock = 0\nverbose = 1\n"
    )
    try:
        proc = subprocess.run("codeml codeml.ctl",
                              shell=True, cwd=str(work),
                              capture_output=True, text=True, timeout=1800)
        out_file = work/f"{label}.M0.out"
        out_text = out_file.read_text() if out_file.exists() else ""
        m_omega = re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", out_text)
        m_dN = re.search(r"tree length for dN:\s+([\d.]+)", out_text)
        m_dS = re.search(r"tree length for dS:\s+([\d.]+)", out_text)
        m_lnL = re.search(r"lnL\(.*?\):\s*([-\d.]+)", out_text)
        m_kappa = re.search(r"kappa \(ts/tv\)\s*=\s*([\d.]+)", out_text)
        row = {
            "gene": label,
            "n_seqs": len(seqs),
            "aln_len": aln_len,
            "omega_M0": float(m_omega.group(1)) if m_omega else None,
            "tree_dN": float(m_dN.group(1)) if m_dN else None,
            "tree_dS": float(m_dS.group(1)) if m_dS else None,
            "lnL": float(m_lnL.group(1)) if m_lnL else None,
            "kappa": float(m_kappa.group(1)) if m_kappa else None,
        }
        codeml_summary.append(row)
        print(f"[codeml] {label}: omega={row['omega_M0']} dN={row['tree_dN']} dS={row['tree_dS']}", flush=True)
    except subprocess.TimeoutExpired:
        print(f"[codeml] {label} TIMEOUT", flush=True)
    except Exception as e:
        print(f"[codeml] {label} FAIL: {e}", flush=True)

pd.DataFrame(codeml_summary).to_csv(OUT/"codeml_M0_summary.csv", index=False)
print(f"wrote {OUT/'codeml_M0_summary.csv'}")

def codon_usage_for_mag(mag):
    ffn = BAKTA_DIR/mag/f"{mag}.ffn"
    if not ffn.exists():
        return None
    counts = Counter()
    total = 0
    gc3 = 0
    gc3_total = 0
    for r in SeqIO.parse(str(ffn), "fasta"):
        s = str(r.seq).upper()
        if len(s) % 3 != 0:
            continue
        for i in range(0, len(s) - 2, 3):
            codon = s[i:i+3]
            if "N" in codon:
                continue
            counts[codon] += 1
            total += 1
            if codon[2] in "GC":
                gc3 += 1
            gc3_total += 1
    if total == 0:
        return None
    aa_codons = defaultdict(list)
    for codon, count in counts.items():
        try:
            aa = str(Seq(codon).translate(table=11))
        except Exception:
            continue
        if aa == "*":
            continue
        aa_codons[aa].append((codon, count))
    Fk = defaultdict(list)
    for aa, cs in aa_codons.items():
        n = sum(c for _, c in cs)
        if n < 2:
            continue
        k = len(cs)
        if k < 2:
            continue
        s = sum((c/n)**2 for _, c in cs)
        F = (n*s - 1)/(n - 1)
        Fk[k].append(F)
    Fbar = {k: (np.mean(v) if v else 0) for k, v in Fk.items()}
    enc = 2
    for k in [2, 3, 4, 6]:
        f = Fbar.get(k, 0)
        if f > 0:
            n_aa = {2: 9, 3: 1, 4: 5, 6: 3}.get(k, 0)
            enc += n_aa / f
    return {
        "MAG": mag,
        "n_codons": total,
        "GC3_pct": 100*gc3/gc3_total if gc3_total else None,
        "ENC": enc,
        "is_hero": mag in HEROES,
    }

cu_rows = [r for mag in mags if (r := codon_usage_for_mag(mag))]
cu_df = pd.DataFrame(cu_rows)
cu_df.to_csv(OUT/"codon_usage_per_MAG.csv", index=False)
print(f"wrote {OUT/'codon_usage_per_MAG.csv'}: {cu_df.shape}")

from scipy.stats import mannwhitneyu
hero_df = cu_df[cu_df.is_hero]
rest_df = cu_df[~cu_df.is_hero]
summary = []
for metric in ["GC3_pct", "ENC"]:
    h = hero_df[metric].dropna()
    r = rest_df[metric].dropna()
    u, p = mannwhitneyu(h, r, alternative="two-sided")
    summary.append({"metric": metric, "hero_mean": float(h.mean()),
                    "rest_mean": float(r.mean()), "MWU_p": float(p)})
pd.DataFrame(summary).to_csv(OUT/"codon_usage_hero_vs_rest.csv", index=False)
print("done C3 v2")
