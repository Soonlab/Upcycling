#!/usr/bin/env python3
"""
(a) Formal novelty description for S13 and S16.

Steps:
  1. Extract 16S rRNA (ffn) for every MAG -> per-sample 16S FASTA
  2. Compute pairwise AAI between S13/S16 and all other Sphingobacterium MAGs
     using reciprocal best hits on Bakta FAA  (fast: mmseqs2 easy-search)
  3. Summarise genome quality & MIGS-lite table
  4. Phylogenetic placement of S13/S16 within Sphingobacterium MAGs using
     skani ANI already computed.
"""
import os, re, gzip, subprocess, shutil, tempfile
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE="/data/data/Upcycling"
BAKTA=f"{BASE}/MAGs_FASTA_files/bakta_results"
OUT=f"{BASE}/research/extra/novel_species"
os.makedirs(OUT, exist_ok=True)

# taxonomy
g = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
g["Genus"]   = g["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "Unclassified")
g["Species"] = g["classification"].apply(lambda s: (re.search(r"s__([^;]*)",str(s)) or [None,""])[1])
tax = g.set_index("user_genome")
sphing_mags = tax.index[tax["Genus"]=="Sphingobacterium"].tolist()
print(f"[info] {len(sphing_mags)} Sphingobacterium MAGs: {sphing_mags}")

# -----------------------------------------------------------------------------
# 1) 16S extraction from Bakta GFF+FFN
# -----------------------------------------------------------------------------
def extract_16S(samp):
    gff = f"{BAKTA}/{samp}/{samp}.gff3"
    ffn = f"{BAKTA}/{samp}/{samp}.ffn"
    if not (os.path.exists(gff) and os.path.exists(ffn)): return []
    ids=[]
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p=line.rstrip().split("\t")
            if len(p)<9 or p[2]!="rRNA": continue
            if "16S ribosomal RNA" not in p[8] and "16S rRNA" not in p[8]: continue
            m=re.search(r"ID=([^;]+)",p[8]);
            if m: ids.append(m.group(1))
    if not ids: return []
    # extract sequences
    seqs=[]
    with open(ffn) as fh:
        name=None; buf=[]
        for line in fh:
            if line.startswith(">"):
                if name in ids and buf: seqs.append((name,"".join(buf)))
                name=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
        if name in ids and buf: seqs.append((name,"".join(buf)))
    return seqs

rrn_out=f"{OUT}/16S_all_sphingobacterium.fasta"
with open(rrn_out,"w") as fh:
    for s in sphing_mags:
        seqs=extract_16S(s)
        for i,(h,seq) in enumerate(seqs):
            fh.write(f">{s}_16S_{i+1} [{tax.loc[s,'Species']}]\n{seq}\n")
print(f"[ok] {rrn_out}")

# -----------------------------------------------------------------------------
# 2) AAI via reciprocal best-hit with mmseqs2 (if available) else BLAST proxy
#    fallback: use proteome-wide average identity via mmseqs easy-search
# -----------------------------------------------------------------------------
def have(cmd):
    from shutil import which; return which(cmd) is not None
mmseqs = shutil.which("mmseqs") or "/home/soon/miniconda3/envs/dram_env/bin/mmseqs"
has_mmseqs = os.path.exists(mmseqs)
print(f"[info] mmseqs2 available: {has_mmseqs} ({mmseqs})")

# Build one combined faa index for all Sphingobacterium + small functional search
tmp=tempfile.mkdtemp(prefix="aai_")
faa_paths={}
for s in sphing_mags:
    fa=f"{BAKTA}/{s}/{s}.faa"
    if os.path.exists(fa):
        faa_paths[s]=fa

def compute_aai(q, t, tmp):
    """Approx AAI: mmseqs easy-search q vs t, min 30% id 70% cov, mean pident of RBH."""
    if not has_mmseqs: return np.nan, 0
    qf, tf = faa_paths.get(q), faa_paths.get(t)
    if not qf or not tf: return np.nan, 0
    res1 = os.path.join(tmp, f"{q}_vs_{t}.m8")
    res2 = os.path.join(tmp, f"{t}_vs_{q}.m8")
    sub = subprocess.run([mmseqs,"easy-search",qf,tf,res1,os.path.join(tmp,"tmp1"),
                          "--min-seq-id","0.3","-c","0.7","--cov-mode","0",
                          "--format-output","query,target,pident,alnlen,evalue",
                          "-s","5","--threads","4","-v","0"],
                         capture_output=True, text=True)
    subprocess.run([mmseqs,"easy-search",tf,qf,res2,os.path.join(tmp,"tmp2"),
                    "--min-seq-id","0.3","-c","0.7","--cov-mode","0",
                    "--format-output","query,target,pident,alnlen,evalue",
                    "-s","5","--threads","4","-v","0"], capture_output=True, text=True)
    if not (os.path.exists(res1) and os.path.exists(res2)): return np.nan, 0
    # best per query
    def best(path):
        d={}
        with open(path) as fh:
            for l in fh:
                p=l.rstrip().split("\t")
                if len(p)<3: continue
                q_,t_,pid = p[0],p[1],float(p[2])
                if q_ not in d or pid>d[q_][1]:
                    d[q_]=(t_,pid)
        return d
    b1=best(res1); b2=best(res2)
    rbh=[]
    for q_,(t_,pid) in b1.items():
        if t_ in b2 and b2[t_][0]==q_:
            rbh.append(pid)
    if not rbh: return np.nan, 0
    return float(np.mean(rbh))*100 if np.mean(rbh)<1.5 else float(np.mean(rbh)), len(rbh)

targets = sphing_mags
rows=[]
for q in ["S13","S16"]:
    if q not in faa_paths:
        print(f"[warn] {q} has no faa"); continue
    for t in targets:
        if t==q: continue
        aai, n = compute_aai(q,t,tmp)
        rows.append({"Query":q,"Target":t,"Target_species":tax.loc[t,"Species"],
                     "AAI":round(aai,2) if not np.isnan(aai) else np.nan,"RBH":n})
aai_df=pd.DataFrame(rows).sort_values(["Query","AAI"],ascending=[True,False])
aai_df.to_csv(f"{OUT}/AAI_S13_S16_vs_Sphingobacterium.csv",index=False)
print(aai_df.to_string(index=False))

# -----------------------------------------------------------------------------
# 3) MIGS-lite table
# -----------------------------------------------------------------------------
def genome_stats(samp):
    fna=f"{BAKTA}/{samp}/{samp}.fna"
    size=0; gc=0; contigs=0; lens=[]
    with open(fna) as fh:
        seq=""
        for line in fh:
            if line.startswith(">"):
                if seq: lens.append(len(seq)); size+=len(seq); gc+=seq.upper().count("G")+seq.upper().count("C")
                seq=""; contigs+=1
            else: seq+=line.strip()
        if seq: lens.append(len(seq)); size+=len(seq); gc+=seq.upper().count("G")+seq.upper().count("C")
    lens.sort(reverse=True); half=size/2; cum=0; n50=0
    for L in lens:
        cum+=L
        if cum>=half: n50=L; break
    # from bakta gff counts
    gff=f"{BAKTA}/{samp}/{samp}.gff3"
    nCDS=ntRNA=nrRNA=0
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p=line.rstrip().split("\t")
            if len(p)<3: continue
            if p[2]=="CDS": nCDS+=1
            elif p[2]=="tRNA": ntRNA+=1
            elif p[2]=="rRNA": nrRNA+=1
    return {"Sample":samp,"Genome_size_Mb":round(size/1e6,3),
            "GC_pct":round(gc/size*100,2) if size else np.nan,
            "Contigs":contigs,"N50_kb":round(n50/1000,1),
            "CDS":nCDS,"tRNA":ntRNA,"rRNA":nrRNA,
            "Genus":tax.loc[samp,"Genus"] if samp in tax.index else "",
            "GTDB_species":tax.loc[samp,"Species"] if samp in tax.index else "",
            "ANI_to_closest":g.set_index("user_genome").loc[samp,"closest_genome_ani"] if samp in g.user_genome.values else ""}

migs=[]
for s in sphing_mags: migs.append(genome_stats(s))
migs_df=pd.DataFrame(migs).sort_values("Sample")
migs_df.to_csv(f"{OUT}/MIGS_lite_Sphingobacterium.csv",index=False)
print("[ok] MIGS table")
print(migs_df.to_string(index=False))

# -----------------------------------------------------------------------------
# 4) Simple visualisation — AAI ranking for S13 and S16
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1,2, figsize=(12,5), sharex=False)
for ax, q in zip(axes, ["S13","S16"]):
    sub = aai_df[aai_df.Query==q].dropna(subset=["AAI"]).copy()
    sub["label"] = sub.apply(lambda r: f"{r.Target} ({r.Target_species or 'sp.'})", axis=1)
    sub = sub.sort_values("AAI", ascending=True).tail(15)
    colors = ["#c0392b" if s in ["S13","S16"] else "#3498db" for s in sub["Target"]]
    ax.barh(sub["label"], sub["AAI"], color=colors, edgecolor="black", lw=0.4)
    ax.axvline(95, ls="--", color="red", lw=0.8, label="95% species")
    ax.axvline(70, ls=":",  color="gray", lw=0.8, label="70% genus")
    ax.set_xlabel("AAI (%)")
    ax.set_xlim(60, 100)
    ax.set_title(f"{q} — AAI vs other Sphingobacterium MAGs", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.tick_params(axis="y", labelsize=7)
fig.suptitle("Amino-acid identity support for novel Sphingobacterium species", fontsize=11)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_NovelSp_AAI.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_NovelSp_AAI.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_NovelSp_AAI")

shutil.rmtree(tmp, ignore_errors=True)
