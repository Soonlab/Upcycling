#!/usr/bin/env python3
"""C4 v2: fold hero UreC with HuggingFace EsmForProteinFolding, compute TM-score vs PDB 4CEU chain C."""
import os, io
from pathlib import Path
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import is_aa
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
import pandas as pd
import numpy as np

ROOT = Path("/data/data/Upcycling/research/additional/C4_esmfold")
OUT_PDB = ROOT/"pdb"
OUT_PDB.mkdir(parents=True, exist_ok=True)
REF_PDB = ROOT/"4CEU.pdb"
FASTA = ROOT/"UreC_hero.fasta"

print("loading ESMFold (HF, CPU fp32)...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
model.trunk.set_chunk_size(64)
model.eval()

def fold(seq):
    ids = tok([seq], return_tensors="pt", add_special_tokens=False)["input_ids"]
    with torch.no_grad():
        out = model(ids)
    return model.output_to_pdb(out)[0]

# Parse 4CEU alpha chain (chain C is UreC - the "alpha" subunit in S. pasteurii urease)
ref_struct = get_structure(str(REF_PDB))
ref_chain_C = None
for model_ in ref_struct:
    for ch in model_:
        if ch.id == "C":
            ref_chain_C = ch
            break
    if ref_chain_C:
        break
if ref_chain_C is None:
    # fallback: pick the longest chain
    longest = max((ch for model_ in ref_struct for ch in model_),
                  key=lambda c: sum(1 for r in c if is_aa(r, standard=True)))
    ref_chain_C = longest
ref_coords, ref_seq = get_residue_data(ref_chain_C)
print(f"reference chain: len={len(ref_seq)}", flush=True)

rows = []
for rec in SeqIO.parse(str(FASTA), "fasta"):
    sid = rec.id
    seq = str(rec.seq).rstrip("*").replace("-", "")
    pdb_path = OUT_PDB/f"{sid}.pdb"
    if not pdb_path.exists():
        print(f"[fold] {sid} len={len(seq)}", flush=True)
        pdb_str = fold(seq)
        pdb_path.write_text(pdb_str)
    # parse pred
    pred_struct = get_structure(str(pdb_path))
    pred_chain = next(ch for model_ in pred_struct for ch in model_)
    pred_coords, pred_seq = get_residue_data(pred_chain)
    res = tm_align(pred_coords, ref_coords, pred_seq, ref_seq)
    rows.append({
        "MAG": sid,
        "pred_len": len(pred_seq),
        "ref_len": len(ref_seq),
        "tm_norm_pred": float(res.tm_norm_chain1),
        "tm_norm_ref": float(res.tm_norm_chain2),
        "rmsd": float(res.rmsd),
    })
    print(f"[{sid}] TM(pred)={rows[-1]['tm_norm_pred']:.3f}  TM(ref)={rows[-1]['tm_norm_ref']:.3f}  RMSD={rows[-1]['rmsd']:.2f}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(ROOT/"ureC_vs_4CEU_tm.csv", index=False)
print(df.to_string(index=False))
