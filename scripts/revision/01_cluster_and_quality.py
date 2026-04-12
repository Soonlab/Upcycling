#!/usr/bin/env python3
"""
Revision step #3 — Rigorous audit of:
  (a) MAG quality for the hero clade (size, N50, #contigs, CDS, rRNA/tRNA)
  (b) Exact contig-level distribution of the ureA/B/C/D/E/F/G and cah genes
      in each hero MAG, to test the claim that the cluster is "single-contig".
"""
import os, re, gzip, glob
from collections import defaultdict
import pandas as pd
import numpy as np

BASE="/data/data/Upcycling"
BAKTA=f"{BASE}/MAGs_FASTA_files/bakta_results"
OUT=f"{BASE}/research/revision"
os.makedirs(OUT, exist_ok=True)

HERO=["C22","M1","S13","S16","S23","S26"]

GENE_KEYS = {
    "ureA": [r"urease subunit gamma"],
    "ureB": [r"urease subunit beta(?! small RNA)"],
    "ureC": [r"urease subunit alpha"],
    "ureD": [r"\bureD\b", r"urease accessory protein UreD"],
    "ureE": [r"\bureE\b", r"urease accessory protein UreE"],
    "ureF": [r"\bureF\b", r"urease accessory protein UreF"],
    "ureG": [r"\bureG\b", r"urease accessory protein UreG"],
    "cah":  [r"carbonic anhydrase"],
}

def parse_gff(path):
    rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS": continue
            attrs = dict(kv.split("=",1) for kv in p[8].split(";") if "=" in kv)
            rows.append({"contig":p[0],"start":int(p[3]),"end":int(p[4]),"strand":p[6],
                         "gene":attrs.get("gene",""),
                         "product":attrs.get("product","")})
    return pd.DataFrame(rows)

def classify_gene(product, gene):
    prod = (product or ""); g = (gene or "")
    for key, patterns in GENE_KEYS.items():
        for pat in patterns:
            if re.search(pat, prod, re.I) or re.search(pat, g, re.I):
                return key
    return None

def genome_stats(samp):
    fna=f"{BAKTA}/{samp}/{samp}.fna"
    lens=[]; total=0; gc=0
    with open(fna) as fh:
        seq=""
        for line in fh:
            if line.startswith(">"):
                if seq:
                    lens.append(len(seq)); total+=len(seq)
                    gc+=seq.upper().count("G")+seq.upper().count("C")
                seq=""
            else: seq+=line.strip()
        if seq:
            lens.append(len(seq)); total+=len(seq)
            gc+=seq.upper().count("G")+seq.upper().count("C")
    lens.sort(reverse=True); half=total/2; cum=0; n50=0
    for L in lens:
        cum+=L
        if cum>=half: n50=L; break
    return total, n50, len(lens), gc/total*100 if total else np.nan

# -----------------------------------------------------------------------------
# Per-hero summary
# -----------------------------------------------------------------------------
report_rows=[]
contig_dist=[]
for samp in HERO:
    gff=f"{BAKTA}/{samp}/{samp}.gff3"
    df = parse_gff(gff)
    df["cls"] = [classify_gene(p,g) for p,g in zip(df["product"], df["gene"])]
    annotated = df[df.cls.notna()].copy()
    # Deduplicate multi-copy annotations: keep the longest representative per contig×gene
    annotated["len"] = annotated["end"]-annotated["start"]
    # Distinct contigs per gene
    per_gene = defaultdict(set)
    for _,r in annotated.iterrows():
        per_gene[r.cls].add(r.contig)
    # genome stats
    size, n50, ncontigs, gc = genome_stats(samp)
    # cluster metrics
    all_ure_contigs = set()
    for g in ["ureA","ureB","ureC","ureD","ureE","ureF","ureG"]:
        all_ure_contigs |= per_gene[g]
    cah_contigs = per_gene["cah"]
    # find contig containing most ure genes
    ure_count_per_contig = defaultdict(int)
    for g in ["ureA","ureB","ureC","ureD","ureE","ureF","ureG"]:
        for c in per_gene[g]: ure_count_per_contig[c]+=1
    if ure_count_per_contig:
        main_contig, max_ure = max(ure_count_per_contig.items(), key=lambda x:x[1])
    else:
        main_contig, max_ure = None, 0
    # cluster span on main contig
    if main_contig:
        main_genes = annotated[(annotated.contig==main_contig) &
                               (annotated.cls.isin(["ureA","ureB","ureC","ureD","ureE","ureF","ureG"]))]
        span_kb = (main_genes["end"].max()-main_genes["start"].min())/1000 if len(main_genes) else np.nan
        # does main contig also carry cah?
        cah_on_main = main_contig in cah_contigs
    else:
        span_kb = np.nan; cah_on_main = False

    report_rows.append({
        "MAG": samp,
        "Genome_size_Mb": round(size/1e6, 3),
        "GC_pct": round(gc, 2),
        "Contigs": ncontigs,
        "N50_kb": round(n50/1000, 1),
        "ureABC_all_present": all(len(per_gene[g])>0 for g in ["ureA","ureB","ureC"]),
        "ureABCDEFG_all_present": all(len(per_gene[g])>0 for g in ["ureA","ureB","ureC","ureD","ureE","ureF","ureG"]),
        "cah_present": len(per_gene["cah"])>0,
        "n_contigs_with_any_ure": len(all_ure_contigs),
        "main_contig": main_contig,
        "ure_genes_on_main_contig": max_ure,
        "cluster_span_kb_main": round(span_kb,2) if not np.isnan(span_kb) else "",
        "cah_on_main_contig": cah_on_main,
        "n_contigs_with_cah": len(cah_contigs),
    })
    # per-gene contig list
    for g in ["ureA","ureB","ureC","ureD","ureE","ureF","ureG","cah"]:
        contig_dist.append({"MAG":samp,"Gene":g,
                            "Contigs":", ".join(sorted(per_gene[g])) or "ABSENT",
                            "n_copies":len(per_gene[g])})

rep = pd.DataFrame(report_rows)
rep.to_csv(f"{OUT}/hero_cluster_audit.csv", index=False)
print("=== Hero cluster contiguity audit ===")
print(rep.to_string(index=False))

pd.DataFrame(contig_dist).to_csv(f"{OUT}/hero_gene_contig_distribution.csv", index=False)
print("\nSaved per-gene contig distribution.")

# -----------------------------------------------------------------------------
# Sentence-ready factual summary for manuscript
# -----------------------------------------------------------------------------
lines=["### Cluster contiguity: factual summary"]
for _,r in rep.iterrows():
    lines.append(
        f"- **{r.MAG}** ({r.Genome_size_Mb} Mb, N50 {r.N50_kb} kb, {r.Contigs} contigs): "
        f"main contig ({r.main_contig}) carries {int(r.ure_genes_on_main_contig)}/7 "
        f"*ure* genes over {r.cluster_span_kb_main} kb; "
        f"*cah* {'co-located' if r.cah_on_main_contig else f'on a separate contig (n={r.n_contigs_with_cah})'}; "
        f"*ure* genes split across {int(r.n_contigs_with_any_ure)} contigs."
    )
with open(f"{OUT}/cluster_factual_summary.md","w") as fh:
    fh.write("\n".join(lines)+"\n")
print("\n".join(lines))
