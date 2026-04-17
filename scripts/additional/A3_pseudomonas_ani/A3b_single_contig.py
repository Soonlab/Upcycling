#!/usr/bin/env python3
"""A3b addendum: per-contig analysis of Pfam hits to quantify
single-contig ureC+CA architecture in 146 Pseudomonas_E refs.
"""
import os, glob, re
import pandas as pd
from collections import defaultdict

WD = "/data/data/Upcycling/research/additional/A3_pseudomonas_ani"

PF_CLASSES = {
    "PF00449": "UreC", "PF00699": "UreB",
    "PF00484": "CA_beta", "PF00194": "CA_alpha", "PF00988": "CA_gamma",
}

rows = []
for dom_path in sorted(glob.glob(f"{WD}/hmm_out/*.dom")):
    acc = os.path.basename(dom_path).replace(".dom","")
    # contig -> set of Pfam classes
    contig2classes = defaultdict(set)
    try:
        with open(dom_path) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split()
                if len(parts) < 22: continue
                pf_full = parts[1]   # e.g. PF00449.27
                pf = pf_full.split(".")[0]
                if pf not in PF_CLASSES: continue
                orf_name = parts[3]  # e.g. NZ_LR134290.1_788
                # contig name = orf_name with trailing "_<int>" stripped
                m = re.match(r"^(.+)_\d+$", orf_name)
                if not m: continue
                contig = m.group(1)
                contig2classes[contig].add(PF_CLASSES[pf])
    except Exception:
        continue

    # per-reference summary
    all_ureC_contigs = [c for c, classes in contig2classes.items() if "UreC" in classes]
    all_CA_contigs = [c for c, classes in contig2classes.items() if any(x in classes for x in ["CA_alpha","CA_beta","CA_gamma"])]
    # co-localized ureC + CA on same contig:
    ureC_CA_same_contig = [c for c, classes in contig2classes.items()
                           if "UreC" in classes and any(x in classes for x in ["CA_alpha","CA_beta","CA_gamma"])]
    # full urease operon same contig: ureC + ureB
    ureC_B_same_contig = [c for c, classes in contig2classes.items()
                          if "UreC" in classes and "UreB" in classes]
    rows.append({
        "accession": acc,
        "has_UreC": int(len(all_ureC_contigs)>0),
        "has_CA_any": int(len(all_CA_contigs)>0),
        "n_contigs_with_ureC": len(all_ureC_contigs),
        "n_contigs_with_CA": len(all_CA_contigs),
        "ureC_and_ureB_single_contig": int(len(ureC_B_same_contig)>0),
        "ureC_and_CA_single_contig": int(len(ureC_CA_same_contig)>0),
    })

df = pd.DataFrame(rows)
df.to_csv(f"{WD}/pseudomonas_e_single_contig.csv", index=False)

n = len(df)
print(f"\n[A3b-contig] n = {n} Pseudomonas_E reference genomes")
print(f"  has UreC:                     {df.has_UreC.sum()}/{n} ({100*df.has_UreC.mean():.1f}%)")
print(f"  has any CA:                   {df.has_CA_any.sum()}/{n} ({100*df.has_CA_any.mean():.1f}%)")
print(f"  ureC + ureB same contig:      {df.ureC_and_ureB_single_contig.sum()}/{n} ({100*df.ureC_and_ureB_single_contig.mean():.1f}%)")
print(f"  ureC + CA same contig:        {df.ureC_and_CA_single_contig.sum()}/{n} ({100*df.ureC_and_CA_single_contig.mean():.1f}%)")

# Note: most Pseudomonas_E refs are complete chromosomal assemblies,
# so "same contig" ≈ "same chromosome" = near universal if genes exist.
# Real rarity criterion is "tight operon" (within N genes) — not captured here.
