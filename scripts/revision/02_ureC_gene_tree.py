#!/usr/bin/env python3
"""
Revision step #5 — ureC gene tree vs species tree.

Extract ureC (urease alpha subunit) from every MAG that encodes it, align,
build ML tree, and quantify topological congruence with the bac120 species
tree using Robinson-Foulds distance and the Shimodaira-Hasegawa (SH) test.
"""
import os, re, subprocess, shutil
from collections import defaultdict
import pandas as pd

BASE="/data/data/Upcycling"
BAKTA=f"{BASE}/MAGs_FASTA_files/bakta_results"
OUT=f"{BASE}/research/revision/ureC_tree"
os.makedirs(OUT, exist_ok=True)

ENV="/home/soon/miniconda3/envs/dram_env/bin"
MAFFT=f"{ENV}/mafft" if os.path.exists(f"{ENV}/mafft") else shutil.which("mafft")
IQTREE=shutil.which("iqtree") or shutil.which("iqtree2")

# Install if missing
if not MAFFT or not IQTREE:
    subprocess.run(["/home/soon/miniconda3/bin/conda","install","-n","dram_env",
                    "-c","bioconda","mafft","iqtree","-y"], check=False,
                   capture_output=True)
    MAFFT=f"{ENV}/mafft"
    IQTREE=f"{ENV}/iqtree" if os.path.exists(f"{ENV}/iqtree") else f"{ENV}/iqtree2"

print(f"mafft={MAFFT}, iqtree={IQTREE}")

# -----------------------------------------------------------------------------
# 1. Extract ureC protein from each MAG (pick the longest copy, >=300 aa)
# -----------------------------------------------------------------------------
def read_fasta(path):
    seqs={}
    with open(path) as fh:
        name=None; buf=[]
        for line in fh:
            if line.startswith(">"):
                if name: seqs[name]="".join(buf)
                name=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
        if name: seqs[name]="".join(buf)
    return seqs

def extract_ureC(samp):
    gff=f"{BAKTA}/{samp}/{samp}.gff3"
    faa=f"{BAKTA}/{samp}/{samp}.faa"
    if not (os.path.exists(gff) and os.path.exists(faa)): return None
    ids=[]
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p=line.rstrip().split("\t")
            if len(p)<9 or p[2]!="CDS": continue
            if "Urease subunit alpha" not in p[8]: continue
            m=re.search(r"ID=([^;]+)",p[8])
            if m: ids.append(m.group(1))
    if not ids: return None
    seqs=read_fasta(faa)
    best=None
    for i in ids:
        s=seqs.get(i,"")
        if len(s)>=300 and (best is None or len(s)>len(best[1])):
            best=(i,s)
    return best

samples = sorted(os.listdir(BAKTA))
fa=f"{OUT}/ureC.faa"
found={}
with open(fa,"w") as fh:
    for s in samples:
        res=extract_ureC(s)
        if res:
            fh.write(f">{s}\n{res[1]}\n")
            found[s]=len(res[1])
print(f"[info] ureC extracted from {len(found)} MAGs")
pd.DataFrame([(k,v) for k,v in found.items()], columns=["MAG","ureC_len"]).to_csv(f"{OUT}/ureC_samples.csv", index=False)

# -----------------------------------------------------------------------------
# 2. Align with mafft and build ML tree with IQ-TREE
# -----------------------------------------------------------------------------
aln=f"{OUT}/ureC.aln.faa"
subprocess.run(f"{MAFFT} --auto {fa} > {aln} 2> {OUT}/mafft.log", shell=True, check=True)
# Remove old outputs
for ext in [".treefile",".log",".iqtree",".bionj",".ckp.gz",".mldist",".model.gz",".splits.nex",".contree"]:
    p=f"{aln}{ext}"
    if os.path.exists(p): os.remove(p)
subprocess.run([IQTREE, "-s", aln, "-m", "LG+G", "-B", "1000", "-T", "4",
                "--quiet", "--prefix", f"{OUT}/ureC"], check=True)
gene_tree=f"{OUT}/ureC.treefile"
print(f"[ok] gene tree: {gene_tree}")

# -----------------------------------------------------------------------------
# 3. Prune the species tree to the same taxa and compute RF distance
# -----------------------------------------------------------------------------
species_tree_src=f"{BASE}/pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile"
from Bio import Phylo
sp = Phylo.read(species_tree_src, "newick")
# leaf names in species tree include species suffix; strip
id_from = lambda s: s.split("_s__")[0]
for t in sp.get_terminals():
    t.name = id_from(t.name)
# prune to MAGs present in gene tree
gt = Phylo.read(gene_tree, "newick")
gene_leaves = set(t.name for t in gt.get_terminals())
sp_leaves   = set(t.name for t in sp.get_terminals())
common = gene_leaves & sp_leaves
print(f"[info] common taxa: {len(common)} / gene {len(gene_leaves)} / species {len(sp_leaves)}")

# write pruned species tree newick
def prune_tree(tree, keep):
    tree = tree
    removals = [t for t in tree.get_terminals() if t.name not in keep]
    for t in removals:
        tree.prune(t)
    return tree

sp_pruned = Phylo.read(species_tree_src, "newick")
for t in sp_pruned.get_terminals(): t.name = id_from(t.name)
sp_pruned = prune_tree(sp_pruned, common)
gt_pruned = prune_tree(gt, common)
Phylo.write(sp_pruned, f"{OUT}/species_pruned.tre", "newick")
Phylo.write(gt_pruned, f"{OUT}/ureC_pruned.tre", "newick")

# RF via ete3 if available
try:
    from ete3 import Tree
    st = Tree(f"{OUT}/species_pruned.tre", format=1)
    gt_e = Tree(f"{OUT}/ureC_pruned.tre", format=1)
    res = st.compare(gt_e, unrooted=True)
    rf = res.get("rf", "NA")
    max_rf = res.get("max_rf", "NA")
    norm = res.get("norm_rf", "NA")
    print(f"[RF] ureC vs species: RF={rf}/{max_rf}  normRF={norm}")
    with open(f"{OUT}/RF_result.txt","w") as fh:
        fh.write(f"ureC vs species (n={len(common)} taxa)\n")
        fh.write(f"RF = {rf}\nmax_RF = {max_rf}\nnormalized_RF = {norm}\n")
        fh.write("Interpretation: normalized_RF close to 0 = congruent (vertical inheritance)\n")
        fh.write("                 normalized_RF close to 1 = incongruent (HGT likely)\n")
except ImportError:
    # Fallback: install
    subprocess.run(["/home/soon/miniconda3/envs/dram_env/bin/pip","install","ete3","six","-q"], check=False)
    print("Please re-run for RF.")

# -----------------------------------------------------------------------------
# 4. SH test (IQ-TREE -z)
# -----------------------------------------------------------------------------
# Concatenate gene_tree and species_tree into a single file and evaluate against ureC alignment
trees_file=f"{OUT}/candidate_trees.tre"
with open(trees_file,"w") as out, open(f"{OUT}/ureC_pruned.tre") as g, open(f"{OUT}/species_pruned.tre") as s:
    out.write(g.read().strip()+"\n")
    out.write(s.read().strip()+"\n")

# Also subset the alignment to common taxa
aln_subset=f"{OUT}/ureC_common.aln.faa"
with open(aln_subset,"w") as out:
    keep=common
    name=None; buf=[]
    with open(aln) as fh:
        for line in fh:
            if line.startswith(">"):
                if name and name in keep:
                    out.write(f">{name}\n{''.join(buf)}\n")
                name=line[1:].strip(); buf=[]
            else: buf.append(line.strip())
        if name and name in keep:
            out.write(f">{name}\n{''.join(buf)}\n")

# clean prior
for ext in [".iqtree",".log",".ckp.gz"]:
    p=f"{aln_subset}{ext}"
    if os.path.exists(p): os.remove(p)
try:
    subprocess.run([IQTREE,"-s",aln_subset,"-z",trees_file,"-m","LG+G","-n","0",
                    "-zb","1000","-zw","-au","-T","4","--quiet",
                    "--prefix",f"{OUT}/SHtest"], check=True, timeout=600)
    print(f"[ok] SH test: {OUT}/SHtest.iqtree")
except Exception as e:
    print(f"[warn] SH test failed: {e}")
