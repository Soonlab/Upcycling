#!/usr/bin/env python3
"""
Final dbCAN analysis combining:
  - Hero (n=6): direct HMMER-3 hmmsearch on Bakta FAA against dbCAN.hmm (E<1e-15)
  - Non-hero (n=105): DRAM cazy_best_hit column from existing annotations.tsv

Produces a unified per-MAG CAZy class/family table and rigorous hero-vs-rest
comparison using the same permutation framework as the trait-module analysis.
"""
import os, re, glob
from collections import defaultdict
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import false_discovery_control

BASE="/data/data/Upcycling"
OUT =f"{BASE}/research/revision"
RNG = np.random.default_rng(42)

HERO = ["C22","M1","S13","S16","S23","S26"]
CLS  = ["GH","GT","PL","CE","AA","CBM"]

# --- (1) parse hero hmmsearch tables ---------------------------------------
def parse_tbl(path):
    """Best HMM per protein (lowest E)."""
    best = {}  # protein -> (hmm_name, evalue, score)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p = line.split()
            target, hmm = p[0], p[2]
            evalue = float(p[4]); score = float(p[5])
            if target not in best or evalue < best[target][1]:
                best[target] = (hmm, evalue, score)
    return best

def classify_hmm(name):
    """HMM names like GH13.hmm, GT4.hmm, CBM50.hmm; strip suffix, return (class, family)."""
    n = name.replace(".hmm","")
    m = re.match(r"(GH|GT|PL|CE|AA|CBM)(\d+.*)?", n)
    if not m: return None, None
    cls = m.group(1); fam = n
    return cls, fam

hero_rows = []
hero_fam  = defaultdict(lambda: defaultdict(int))
for h in HERO:
    tbl = f"{OUT}/dbcan_direct/{h}.tbl"
    best = parse_tbl(tbl)
    by_cls = defaultdict(int)
    for prot,(hmm,ev,sc) in best.items():
        cls, fam = classify_hmm(hmm)
        if cls:
            by_cls[cls] += 1
            hero_fam[h][fam] += 1
    row = {"fasta":h, **{c: by_cls.get(c,0) for c in CLS}}
    row["Total"] = sum(by_cls.values())
    hero_rows.append(row)
hero_df = pd.DataFrame(hero_rows).set_index("fasta")
print("=== Hero dbCAN class counts (direct hmmsearch, E<1e-15) ===")
print(hero_df)

# --- (2) parse non-hero from existing DRAM all_annotations.tsv ------------
ann = pd.read_csv("/data/pangenome_work/dram_output/all_annotations.tsv", sep="\t",
                  usecols=["fasta","cazy_best_hit","cazy_ids"])
ann = ann[~ann["fasta"].isin(HERO)]  # only non-hero
def family_best(x):
    if pd.isna(x): return None
    m = re.search(r"\b((?:GH|GT|PL|CE|AA|CBM)\d+)", str(x))
    return m.group(1) if m else None
def class_best(x):
    f = family_best(x)
    if not f: return None
    m = re.match(r"(GH|GT|PL|CE|AA|CBM)", f); return m.group(1) if m else None

ann["cls"] = ann.apply(lambda r: class_best(r["cazy_best_hit"]) or class_best(r["cazy_ids"]), axis=1)
ann["fam"] = ann.apply(lambda r: family_best(r["cazy_best_hit"]) or family_best(r["cazy_ids"]), axis=1)
rest_ann = ann[ann["cls"].notna()]
rest_cls = rest_ann.groupby(["fasta","cls"]).size().unstack(fill_value=0)
for c in CLS:
    if c not in rest_cls.columns: rest_cls[c]=0
rest_cls = rest_cls[CLS]; rest_cls["Total"]=rest_cls.sum(axis=1)
print(f"\nNon-hero CAZy summary (n={len(rest_cls)} MAGs)")
print(rest_cls.describe().round(1))

# Rest family counts
rest_fam = rest_ann.groupby(["fasta","fam"]).size().unstack(fill_value=0)

# --- (3) combine + normalise by CDS counts ---------------------------------
cds = pd.read_csv(f"{BASE}/research/extra/gene_category_counts.csv").set_index("Sample")["CDS_total"]
combined = pd.concat([hero_df, rest_cls])
combined = combined.loc[combined.index.isin(cds.index)]
norm = combined[CLS].div(cds.loc[combined.index], axis=0) * 1000.0
norm.to_csv(f"{OUT}/dbCAN_final_class_per1kCDS.csv")
combined.to_csv(f"{OUT}/dbCAN_final_class_counts.csv")

# --- (4) Hero vs Rest statistics (permutation) -----------------------------
hero_mask = norm.index.isin(HERO)
rows = []
for c in CLS:
    x = norm[c].values
    hm = x[hero_mask].mean(); rm = x[~hero_mask].mean()
    fc = hm/rm if rm>0 else np.nan
    # one-sided permutation
    n_perm = 10000
    obs = hm - rm
    if obs <= 0:
        p_perm = 1.0
    else:
        idx = np.arange(len(x)); k=hero_mask.sum(); hits=0
        for _ in range(n_perm):
            RNG.shuffle(idx)
            m = np.zeros(len(x),dtype=bool); m[idx[:k]]=True
            if (x[m].mean()-x[~m].mean()) >= obs: hits+=1
        p_perm = (hits+1)/(n_perm+1)
    # bootstrap CI
    h_idx=np.where(hero_mask)[0]; r_idx=np.where(~hero_mask)[0]
    fcs=[]
    for _ in range(2000):
        hs=RNG.choice(h_idx,size=len(h_idx),replace=True)
        rs=RNG.choice(r_idx,size=len(r_idx),replace=True)
        h=x[hs].mean(); r=x[rs].mean()
        fcs.append(h/r if r>0 else np.nan)
    fcs=np.array(fcs)
    flo,fhi = np.nanpercentile(fcs,2.5), np.nanpercentile(fcs,97.5)
    try: mw = stats.mannwhitneyu(x[hero_mask], x[~hero_mask], alternative="greater").pvalue
    except: mw=np.nan
    rows.append({"Class":c,"Hero_mean_per1kCDS":round(hm,3),"Rest_mean_per1kCDS":round(rm,3),
                 "Fold_change":round(fc,2) if not np.isnan(fc) else np.nan,
                 "Bootstrap_CI95":f"[{flo:.2f}, {fhi:.2f}]",
                 "Permutation_p":round(p_perm,4),
                 "MannWhitney_greater_p":round(mw,4)})
cls_stats=pd.DataFrame(rows)
cls_stats["Permutation_q_BH"]=false_discovery_control(cls_stats["Permutation_p"].values).round(4)
cls_stats.to_csv(f"{OUT}/dbCAN_final_hero_vs_rest_class.csv", index=False)
print("\n=== dbCAN class-level hero vs rest (DRAM+direct hmmsearch combined) ===")
print(cls_stats.to_string(index=False))

# --- (5) Figure ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8,5))
x=np.arange(len(CLS)); w=0.38
ax.bar(x-w/2, cls_stats["Hero_mean_per1kCDS"], w,
       label=f"Hero (n={hero_mask.sum()})", color="#c0392b", edgecolor="black", lw=0.4)
ax.bar(x+w/2, cls_stats["Rest_mean_per1kCDS"], w,
       label=f"Rest (n={(~hero_mask).sum()})", color="#7f8c8d", edgecolor="black", lw=0.4)
ax.set_xticks(x); ax.set_xticklabels(CLS)
ax.set_ylabel("CAZymes per 1k CDS")
ax.set_title("dbCAN class-level profile (direct hmmsearch hero + DRAM rest)", fontsize=11)
# Annotate q
for i,row in cls_stats.iterrows():
    q=row["Permutation_q_BH"]
    if q<0.05: mk="***"
    elif q<0.10: mk="**"
    elif q<0.20: mk="*"
    else: mk=""
    if mk:
        y_top=max(row["Hero_mean_per1kCDS"],row["Rest_mean_per1kCDS"])*1.08
        ax.text(i, y_top, mk, ha="center", fontsize=11)
ax.legend(frameon=False)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_dbCAN_final_classes.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_dbCAN_final_classes.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_dbCAN_final_classes")
