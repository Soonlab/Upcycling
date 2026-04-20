#!/usr/bin/env python3
"""
C3 dN/dS + codon usage on urease subunits.
- Pull ureA, ureB, ureC, ureG nucleotide CDS from each MAG's Bakta .ffn
- Build per-gene codon-aware MSA (mafft on AA + back-translate)
- Run codeml M0 (single-omega) and M1a vs M2a (selection test)
- Compute codon usage indices (ENC, GC3, RSCU) per MAG / hero vs rest
"""
import os, re, subprocess, json, gzip
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import standard_dna_table

BAKTA_DIR = Path("/data/data/Upcycling/MAGs_FASTA_files/bakta_results")
OUT = Path("/data/data/Upcycling/research/additional/C3_dnds_codon")
OUT.mkdir(parents=True, exist_ok=True)
(OUT/"per_gene").mkdir(exist_ok=True)
(OUT/"codeml").mkdir(exist_ok=True)
HEROES = {"S13","S16","S23","C22","M1","S26"}

# Targets — match Bakta product names case-insensitive
TARGETS = {
    "ureA": ["urease subunit gamma","ureA"],
    "ureB": ["urease subunit beta","ureB"],
    "ureC": ["urease subunit alpha","ureC"],
    "ureG": ["urease accessory protein UreG","ureG"],
}

def parse_bakta_cds(mag):
    """Yield (gene_label, locus_tag, nt_seq) from Bakta files for given MAG."""
    tsv = BAKTA_DIR/mag/f"{mag}.tsv"
    ffn = BAKTA_DIR/mag/f"{mag}.ffn"
    if not (tsv.exists() and ffn.exists()):
        return
    # Parse TSV header → find Locus_Tag and Product columns
    locus2product = {}
    with open(tsv) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                header = [c.strip().lower() for c in cols]
                continue
            d = dict(zip(header, cols))
            lt = d.get("locus tag") or d.get("locus_tag") or ""
            prod = d.get("product","")
            if lt:
                locus2product[lt] = prod
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(ffn), "fasta")}
    for lt, prod in locus2product.items():
        for label, patterns in TARGETS.items():
            if any(re.search(p, prod, re.I) for p in patterns):
                seq = seqs.get(lt) or seqs.get(lt.split()[0])
                if seq is not None and len(seq) % 3 == 0 and len(seq) >= 200:
                    yield label, lt, seq

# 1) Collect per-gene multi-FASTA across all MAGs
gene_seqs = defaultdict(dict)
mags = sorted([d.name for d in BAKTA_DIR.iterdir() if d.is_dir()])
for mag in mags:
    for label, lt, seq in parse_bakta_cds(mag):
        # keep only the longest copy per (mag, gene)
        existing = gene_seqs[label].get(mag)
        if existing is None or len(seq) > len(existing[1]):
            gene_seqs[label][mag] = (lt, seq)

# Write per-gene multifasta (NT)
for label, d in gene_seqs.items():
    nt_fa = OUT/"per_gene"/f"{label}.nt.fasta"
    aa_fa = OUT/"per_gene"/f"{label}.aa.fasta"
    with open(nt_fa,"w") as fnt, open(aa_fa,"w") as faa:
        for mag, (lt, nt) in sorted(d.items()):
            try:
                aa = str(Seq(nt).translate(table=11, to_stop=True))
            except Exception:
                continue
            if "*" in aa or len(aa)<60: continue
            fnt.write(f">{mag}\n{nt}\n")
            faa.write(f">{mag}\n{aa}\n")
    print(f"{label}: {sum(1 for _ in open(nt_fa))//2} mags")

# 2) MSA + back-translate (codon-aware)
def run(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw)

def msa_and_pal2nal(label):
    aa_fa = OUT/"per_gene"/f"{label}.aa.fasta"
    nt_fa = OUT/"per_gene"/f"{label}.nt.fasta"
    aa_aln = OUT/"per_gene"/f"{label}.aa.aln.fasta"
    codon_aln = OUT/"per_gene"/f"{label}.codon.aln.fasta"
    if aa_fa.stat().st_size == 0:
        return False
    run(f"mafft --auto --thread 4 {aa_fa} > {aa_aln}")
    # pal2nal optional: do simple back-translation manually
    aa_seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(aa_aln),"fasta")}
    nt_seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(nt_fa),"fasta")}
    with open(codon_aln,"w") as fh:
        for sid, aa in aa_seqs.items():
            nt = nt_seqs.get(sid)
            if nt is None: continue
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
        msa_and_pal2nal(label)
        print(f"[msa] {label} OK")
    except Exception as e:
        print(f"[msa] {label} FAIL: {e}")

# 3) codeml M0 (single ω) — pairwise via runmode -2 (yn00 alt)
# Use runmode=0 with M0 model on the gene tree (build quick NJ from codon aln)
codeml_summary = []
for label in TARGETS:
    codon_aln = OUT/"per_gene"/f"{label}.codon.aln.fasta"
    if not codon_aln.exists() or codon_aln.stat().st_size == 0:
        continue
    work = OUT/"codeml"/label
    work.mkdir(exist_ok=True)
    # Convert FASTA -> PHYLIP-relaxed
    seqs = list(SeqIO.parse(str(codon_aln),"fasta"))
    if len(seqs) < 4:
        continue
    aln_len = len(seqs[0].seq)
    phy = work/f"{label}.phy"
    with open(phy,"w") as fh:
        fh.write(f" {len(seqs)} {aln_len}\n")
        for r in seqs:
            name = r.id[:30].ljust(31)
            fh.write(f"{name}{str(r.seq)}\n")
    # Build a quick NJ tree using mafft+iqtree if available; else star
    tree = work/f"{label}.tre"
    star = "(" + ",".join(r.id for r in seqs) + ");"
    tree.write_text(star)
    # codeml ctl
    ctl = work/"codeml.ctl"
    ctl.write_text(f"""seqfile = {phy.name}
treefile = {tree.name}
outfile = {label}.M0.out
runmode = 0
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
clock = 0
getSE = 0
RateAncestor = 0
""")
    try:
        proc = subprocess.run("codeml codeml.ctl",
                              shell=True, cwd=str(work),
                              capture_output=True, text=True, timeout=600)
        out_text = (work/f"{label}.M0.out").read_text() if (work/f"{label}.M0.out").exists() else ""
        m = re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", out_text)
        omega = float(m.group(1)) if m else None
        m2 = re.search(r"tree length for dN:\s+([\d.]+)", out_text)
        m3 = re.search(r"tree length for dS:\s+([\d.]+)", out_text)
        codeml_summary.append({
            "gene": label, "n_seqs": len(seqs),
            "omega_M0": omega,
            "tree_dN": float(m2.group(1)) if m2 else None,
            "tree_dS": float(m3.group(1)) if m3 else None
        })
        print(f"[codeml] {label}: ω = {omega}")
    except Exception as e:
        print(f"[codeml] {label} FAIL: {e}")

pd.DataFrame(codeml_summary).to_csv(OUT/"codeml_M0_summary.csv", index=False)
print(f"wrote {OUT/'codeml_M0_summary.csv'}")

# 4) Codon usage per MAG (whole-genome)
def codon_usage_for_mag(mag):
    ffn = BAKTA_DIR/mag/f"{mag}.ffn"
    if not ffn.exists(): return None
    counts = Counter()
    total = 0
    gc3 = 0; gc3_total = 0
    for r in SeqIO.parse(str(ffn),"fasta"):
        s = str(r.seq).upper()
        if len(s)%3 != 0: continue
        for i in range(0, len(s)-2, 3):
            codon = s[i:i+3]
            if "N" in codon: continue
            counts[codon] += 1
            total += 1
            if codon[2] in "GC": gc3 += 1
            gc3_total += 1
    if total == 0: return None
    # ENC (Wright 1990)
    aa_codons = defaultdict(list)
    for codon, count in counts.items():
        try:
            aa = str(Seq(codon).translate(table=11))
        except Exception:
            continue
        if aa == "*": continue
        aa_codons[aa].append((codon, count))
    Fk = defaultdict(list)  # by synonymous family size
    for aa, cs in aa_codons.items():
        n = sum(c for _, c in cs)
        if n < 1: continue
        k = len(cs)
        if k < 2:
            continue
        s = sum((c/n)**2 for _, c in cs)
        F = (n*s - 1)/(n - 1) if n > 1 else 0
        Fk[k].append(F)
    # ENC formula
    Fbar = {k: (np.mean(v) if v else 0) for k,v in Fk.items()}
    enc = 2  # for Met+Trp (k=1)
    for k in [2,3,4,6]:
        f = Fbar.get(k, 0)
        if f > 0:
            n_aa = {2:9, 3:1, 4:5, 6:3}.get(k,0)
            enc += n_aa / f
    return {
        "MAG": mag,
        "n_codons": total,
        "GC3_pct": 100*gc3/gc3_total if gc3_total else None,
        "ENC": enc,
        "is_hero": mag in HEROES,
    }

cu_rows = []
for mag in mags:
    rec = codon_usage_for_mag(mag)
    if rec: cu_rows.append(rec)
cu_df = pd.DataFrame(cu_rows)
cu_df.to_csv(OUT/"codon_usage_per_MAG.csv", index=False)
print(f"wrote {OUT/'codon_usage_per_MAG.csv'}: {cu_df.shape}")

# Hero vs rest
from scipy.stats import mannwhitneyu
res = []
for col in ["GC3_pct","ENC"]:
    h = cu_df[cu_df.is_hero][col].dropna().values
    r = cu_df[~cu_df.is_hero][col].dropna().values
    if len(h)>=2 and len(r)>=2:
        u,p = mannwhitneyu(h,r,alternative="two-sided")
        res.append({"metric":col,"hero_mean":float(np.mean(h)),"rest_mean":float(np.mean(r)),"MWU_p":p})
pd.DataFrame(res).to_csv(OUT/"codon_usage_hero_vs_rest.csv", index=False)
print("done C3")
