#!/usr/bin/env python3
"""
C4 ESMFold structure prediction for hero UreC + comparison to PDB 4CEU
- Pulls UreC AA seq for 6 hero MAGs (already aligned to P41020 in A2)
- Predicts 3D structure with ESMFold (RTX 5090)
- Aligns to PDB 4CEU chain A with TM-align → RMSD/TM-score
"""
import os, sys, re
from pathlib import Path
import urllib.request
import torch

OUT = Path("/data/data/Upcycling/research/additional/C4_esmfold")
PDB_DIR = OUT/"pdb"
PDB_DIR.mkdir(parents=True, exist_ok=True)

A2_DIR = Path("/data/data/Upcycling/research/additional/A2_structure")
HEROES = ["S13","S16","S23","C22","M1","S26"]

# Try to find UreC AA seqs from A2 output
ureC_aa = {}
candidate = A2_DIR/"UreC_aligned.faa"
if candidate.exists():
    from Bio import SeqIO
    for r in SeqIO.parse(str(candidate),"fasta"):
        # strip alignment gaps
        seq = str(r.seq).replace("-","").replace("*","")
        # Match MAG name
        for h in HEROES:
            if h in r.id:
                ureC_aa[h] = seq
                break

# Fallback: pull from Bakta
if len(ureC_aa) < len(HEROES):
    BAKTA = Path("/data/data/Upcycling/MAGs_FASTA_files/bakta_results")
    from Bio import SeqIO
    for h in HEROES:
        if h in ureC_aa: continue
        tsv = BAKTA/h/f"{h}.tsv"
        faa = BAKTA/h/f"{h}.faa"
        if not (tsv.exists() and faa.exists()): continue
        # Find ureC locus tag
        target = None
        with open(tsv) as fh:
            header = None
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                cols = line.rstrip("\n").split("\t")
                if header is None:
                    header = [c.strip().lower() for c in cols]; continue
                d = dict(zip(header,cols))
                if re.search("urease subunit alpha|ureC", d.get("product",""), re.I):
                    target = d.get("locus tag") or d.get("locus_tag")
                    break
        if not target: continue
        for r in SeqIO.parse(str(faa),"fasta"):
            if r.id == target or r.id.startswith(target):
                ureC_aa[h] = str(r.seq).replace("*","")
                break

# Save FASTA
fa = OUT/"UreC_hero.fasta"
with open(fa,"w") as fh:
    for h in HEROES:
        if h in ureC_aa:
            fh.write(f">{h}_UreC\n{ureC_aa[h]}\n")
print(f"writing {len(ureC_aa)} sequences to {fa}")

# ESMFold predict
import esm
print("loading ESMFold ...")
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
model.set_chunk_size(64)

for h, seq in ureC_aa.items():
    pdb_path = PDB_DIR/f"{h}_UreC.pdb"
    if pdb_path.exists():
        print(f"[skip] {h}")
        continue
    print(f"[fold] {h} ({len(seq)} aa)")
    with torch.no_grad():
        out = model.infer_pdb(seq)
    pdb_path.write_text(out)
    # Mean pLDDT
    plddts = []
    for line in out.splitlines():
        if line.startswith("ATOM") and " CA " in line:
            try:
                plddts.append(float(line[60:66]))
            except Exception:
                pass
    print(f"   {h}: mean pLDDT = {sum(plddts)/len(plddts):.2f}")

# Download reference 4CEU
ref_pdb = PDB_DIR/"4CEU.pdb"
if not ref_pdb.exists():
    print("downloading 4CEU.pdb ...")
    urllib.request.urlretrieve("https://files.rcsb.org/download/4CEU.pdb", ref_pdb)

# TMalign vs reference
import subprocess, json
tm_summary = []
for h in HEROES:
    pred = PDB_DIR/f"{h}_UreC.pdb"
    if not pred.exists(): continue
    try:
        proc = subprocess.run(["TMalign", str(pred), str(ref_pdb)],
                              capture_output=True, text=True, timeout=600)
        out = proc.stdout
        rmsd_m = re.search(r"RMSD=\s*([\d.]+)", out)
        tm_m   = re.search(r"TM-score=\s*([\d.]+).*\(if normalized by length of Chain_2", out)
        seqid_m = re.search(r"Seq_ID=n_identical/n_aligned=\s*([\d.]+)", out)
        tm_summary.append({
            "MAG": h,
            "RMSD_to_4CEU": float(rmsd_m.group(1)) if rmsd_m else None,
            "TM_score_norm_4CEU": float(tm_m.group(1)) if tm_m else None,
            "Seq_ID": float(seqid_m.group(1)) if seqid_m else None,
        })
    except Exception as e:
        print(f"[TMalign] {h} FAIL: {e}")

import pandas as pd
pd.DataFrame(tm_summary).to_csv(OUT/"esmfold_vs_4CEU.csv", index=False)
print("done C4")
