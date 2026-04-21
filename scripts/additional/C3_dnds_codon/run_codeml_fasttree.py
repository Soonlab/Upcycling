#!/usr/bin/env python3
"""Build proper bifurcating tree with FastTree (avoids MAXNSONS polytomy error) then run codeml M0."""
import os, re, subprocess
from pathlib import Path
from Bio import SeqIO
import pandas as pd

OUT = Path("/data/data/Upcycling/research/additional/C3_dnds_codon")
TARGETS = ["ureA","ureB","ureC","ureG"]

def sh(cmd, cwd=None, timeout=1800):
    return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, timeout=timeout)

def write_phy(fasta, phy):
    recs = list(SeqIO.parse(str(fasta), "fasta"))
    L = len(recs[0].seq)
    with open(phy, "w") as fh:
        fh.write(f" {len(recs)} {L}\n")
        for r in recs:
            fh.write(f"{r.id[:30].ljust(31)}{str(r.seq)}\n")
    return len(recs), L

rows = []
for label in TARGETS:
    sub_fa = OUT/"per_gene"/f"{label}.codon.sub.fasta"
    if not sub_fa.exists():
        continue
    work = OUT/"codeml_v4"/label
    work.mkdir(parents=True, exist_ok=True)
    aa_for_tree = work/f"{label}.aa.fasta"
    # Build AA from codon alignment for FastTree
    recs = list(SeqIO.parse(str(sub_fa), "fasta"))
    with open(aa_for_tree, "w") as fh:
        from Bio.Seq import Seq
        for r in recs:
            s = str(r.seq).replace("-", "N")
            try:
                aa = str(Seq(s).translate(table=11))
            except Exception:
                continue
            fh.write(f">{r.id}\n{aa}\n")
    tree_nwk = work/f"{label}.nwk"
    r_ft = sh(f"FastTree -quiet -nosupport -gamma {aa_for_tree} > {tree_nwk}", cwd=str(work))
    if not tree_nwk.exists() or tree_nwk.stat().st_size < 10:
        print(f"[{label}] FastTree FAIL: {r_ft.stderr[:200]}")
        continue
    # Strip branch-support values, keep topology + lengths
    tr = tree_nwk.read_text().strip()
    # Copy sub alignment + write PHY + write ctl
    phy = work/f"{label}.phy"
    n, L = write_phy(sub_fa, phy)
    tre = work/f"{label}.tre"
    tre.write_text(tr + "\n")
    ctl = work/"codeml.ctl"
    ctl.write_text(
        f"seqfile = {phy.name}\n"
        f"treefile = {tre.name}\n"
        f"outfile = {label}.M0.out\n"
        "runmode = 0\nseqtype = 1\nCodonFreq = 2\nmodel = 0\nNSsites = 0\nicode = 0\n"
        "fix_kappa = 0\nkappa = 2\nfix_omega = 0\nomega = 0.4\nclock = 0\ngetSE = 0\nverbose = 0\n"
    )
    r_cm = sh("codeml codeml.ctl", cwd=str(work), timeout=1800)
    out_f = work/f"{label}.M0.out"
    if not out_f.exists():
        print(f"[{label}] codeml produced no output: {r_cm.stderr[:200]}")
        continue
    t = out_f.read_text()
    g_omega = re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", t)
    g_dN = re.search(r"tree length for dN:\s+([\d.]+)", t)
    g_dS = re.search(r"tree length for dS:\s+([\d.]+)", t)
    g_ln = re.search(r"lnL\(.*?\):\s*([-\d.]+)", t)
    g_k = re.search(r"kappa \(ts/tv\)\s*=\s*([\d.]+)", t)
    rows.append({
        "gene": label, "n_seqs": n, "aln_len": L,
        "omega_M0": float(g_omega.group(1)) if g_omega else None,
        "tree_dN": float(g_dN.group(1)) if g_dN else None,
        "tree_dS": float(g_dS.group(1)) if g_dS else None,
        "lnL": float(g_ln.group(1)) if g_ln else None,
        "kappa": float(g_k.group(1)) if g_k else None,
    })
    print(f"[{label}] omega={rows[-1]['omega_M0']} dN={rows[-1]['tree_dN']} dS={rows[-1]['tree_dS']}", flush=True)

pd.DataFrame(rows).to_csv(OUT/"codeml_M0_summary.csv", index=False)
print(f"wrote codeml_M0_summary.csv: {len(rows)} genes")
print(pd.DataFrame(rows).to_string(index=False))
