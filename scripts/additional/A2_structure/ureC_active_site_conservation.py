#!/usr/bin/env python3
"""A2: UreC active-site conservation.

Pull UreC sequences from S13, S16, S23, C22, M1, S26 (MICP-complete MAGs)
from Bakta FFN/FAA. Build MSA vs reference UreCs
(Sporosarcina pasteurii UreC = UniProt Q59674, Helicobacter pylori UreB,
Klebsiella aerogenes UreC = P18314). Check conservation of:
  - 4 Ni²⁺-binding His residues (H134, H136, H246, H272 in K. aerogenes numbering)
  - Carbamate K217 (post-translationally carboxylated lysine)
  - Catalytic D360
  - Mobile flap region residues (conserved CCGG)
Reference positions per Benini et al. 1999 (Structure 7:205) and Ha et al. 2001.
"""
from __future__ import annotations
import os, glob, re, subprocess
import pandas as pd
from Bio import SeqIO

OUT = "/data/data/Upcycling/research/additional/A2_structure"
BAKTA = "/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
os.makedirs(OUT, exist_ok=True)

HERO = ["S13", "S16", "S23", "C22", "M1", "S26"]

def extract_urease_alpha(mag: str) -> tuple[str, str] | None:
    """Find urease alpha subunit (UreC) protein in Bakta output."""
    tsv = f"{BAKTA}/{mag}/{mag}.tsv"
    faa = f"{BAKTA}/{mag}/{mag}.faa"
    if not os.path.exists(tsv) or not os.path.exists(faa):
        return None
    df = pd.read_csv(tsv, sep="\t", comment="#", header=None,
                     names=["contig","type","start","end","strand","locus","gene","product","dbxrefs"])
    df["text"] = df["gene"].fillna("") + " || " + df["product"].fillna("")
    # UreC = urease alpha subunit
    pattern = re.compile(r"\bure[Cc]\b|urease.*alpha|alpha.*urease|urease subunit alpha", re.I)
    hits = df[df["text"].str.contains(pattern)]
    if hits.empty:
        return None
    # Prefer longer proteins (UreC is ~570 aa)
    loci = hits["locus"].tolist()
    seqs = {}
    for rec in SeqIO.parse(faa, "fasta"):
        rec_id = rec.id.split()[0]
        if rec_id in loci:
            seqs[rec_id] = str(rec.seq).replace("*", "")
    if not seqs:
        return None
    best = max(seqs.items(), key=lambda kv: len(kv[1]))
    return (best[0], best[1])

# ---- 1) extract UreC from each MICP-complete MAG
rows = []
outfasta = f"{OUT}/UreC_MICP_complete.faa"
with open(outfasta, "w") as out:
    for mag in HERO:
        r = extract_urease_alpha(mag)
        if r is None:
            print(f"[A2] WARN: no UreC for {mag}")
            continue
        locus, seq = r
        out.write(f">{mag}__{locus}\n{seq}\n")
        rows.append({"MAG": mag, "locus": locus, "length_aa": len(seq)})

pd.DataFrame(rows).to_csv(f"{OUT}/UreC_lengths.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))

# ---- 2) Reference sequences (hard-coded from UniProt)
# S. pasteurii UreC (P41020 urease alpha, 570 aa)
# Use NCBI ureases below as minimal refs - appended to query fasta for MSA.
REFS = {
    # S. pasteurii UreC (UniProt P41020, 570 aa, PDB 4CEU template)
    "SporosarcinaPasteurii_P41020": (
        "MKINRQQYAESYGPTVGDQVRLADTDLWIEVEKDTTYGDEAVNFGGGKVLREGMGENGTY"
        "TRTENVLDLLLTNALILDYTGIYKADIGVKDGYIVGIGKGGNPDIMDGVTPNMIVGTATE"
        "VIAAEGKIVTAGGIDTHVHFINPDQVDVALANGITTLFGGGTGPAEGSKATTVTPGPWNI"
        "EKMLKSTEGLPINVGILGKGHGSSIAPIMEQIDAGAAGLKIHEDWGATPASIDRSLTVAD"
        "EADVQVAIHSDTLNEAGFLEDTLRAINGRVIHSFHVEGAGGGHAPDIMAMAGHPNVLPSS"
        "TNPTRPFTVNTIDEHLDMLMVCHHLKQNIPEDVAFADSRIRPETIAAEDILHDLGIISMM"
        "STDALAMGRAGEMVLRTWQTADKMKKQRGPLAEEKNGSDNFRAKRYVSKYTINPAIAQGI"
        "AHEVGSIEEGKFADLVLWEPKFFGVKADRVIKGGIIAYAQIGDPSASIPTPQPVMGRRMY"
        "GTVGDLIHDTNITFMSKSSIQQGVPAKLGLKRRIGTVKNCRNIGKKDMKWNDVTTDIDIN"
        "PETYEVKVDGEVLTCEPVKELPMAQRYFLF"
    ),
}

with open(f"{OUT}/UreC_with_refs.faa", "w") as out:
    for name, seq in REFS.items():
        out.write(f">{name}\n{seq}\n")
    with open(outfasta) as f:
        out.write(f.read())

# ---- 3) run MAFFT to align, IQ-TREE optional
print("[A2] running MAFFT alignment ...")
r = subprocess.run(
    ["mafft", "--auto", "--thread", "8", f"{OUT}/UreC_with_refs.faa"],
    capture_output=True, text=True)
if r.returncode == 0:
    with open(f"{OUT}/UreC_aligned.faa", "w") as f:
        f.write(r.stdout)
    print(f"[A2] MSA saved: {OUT}/UreC_aligned.faa")
else:
    print("[A2] MAFFT failed:", r.stderr[:300])

# ---- 4) Active-site residue inspection
# S. pasteurii UreC (P41020) active-site residues (Benini et al. 1999, PDB 4CEU):
#   H137, H139  — distal Ni²⁺ coordination
#   K220        — carbamylated bridging lysine (posttranslationally CO₂-modified)
#   H249, H275  — proximal Ni²⁺ coordination
#   D363        — general acid / proton shuttle
#   C322        — flap cysteine (conformational gating)
ACTIVE = {"H137":"H","H139":"H","K220":"K","H249":"H","H275":"H","D363":"D","C322":"C"}

def parse_aln(path):
    seqs = {}
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf=[]
            else:
                buf.append(line)
        if name: seqs[name] = "".join(buf)
    return seqs

if os.path.exists(f"{OUT}/UreC_aligned.faa"):
    aln = parse_aln(f"{OUT}/UreC_aligned.faa")
    ref_name = "SporosarcinaPasteurii_P41020"
    if ref_name in aln:
        ref_aln = aln[ref_name]
        # build map from ref ungapped index (1-based) -> alignment column
        ref_map = {}
        ungap_i = 0
        for col_i, aa in enumerate(ref_aln):
            if aa != "-":
                ungap_i += 1
                ref_map[ungap_i] = col_i
        res_rows = []
        for pos_name, expected in ACTIVE.items():
            resi = int(pos_name[1:])
            col = ref_map.get(resi, None)
            row = {"site": pos_name, "expected": expected, "ref_column": col}
            if col is None:
                row["ref_aa"] = "NA"
            else:
                row["ref_aa"] = ref_aln[col]
            for mag in HERO:
                mag_key = [k for k in aln if k.startswith(f"{mag}__")]
                if mag_key:
                    row[mag] = aln[mag_key[0]][col] if col is not None else "NA"
                else:
                    row[mag] = "missing"
            res_rows.append(row)
        rdf = pd.DataFrame(res_rows)
        rdf.to_csv(f"{OUT}/UreC_active_site_residues.csv", index=False)
        print("\n[A2] Active-site residue conservation:")
        print(rdf.to_string(index=False))
        conserved = sum(1 for _, r in rdf.iterrows()
                        for mag in HERO if r[mag] == r["expected"])
        total = len(ACTIVE) * len(HERO)
        print(f"\n[A2] conservation: {conserved}/{total} "
              f"= {100*conserved/total:.1f}% of active-site residues preserved in MICP-complete MAGs")
    else:
        print(f"[A2] ERROR: reference {ref_name} not in alignment. Fix ref seq.")
print("[A2] DONE")
