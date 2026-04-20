#!/usr/bin/env python3
"""
C6 In-situ abundance proxy
- Raw FASTQ unavailable → use SPAdes contig header `cov_X` as relative depth
- Per MAG: length-weighted mean coverage (= within-sample dominance proxy)
- Per source group (cattle/swine/sheep/poultry): distribution
- Hero vs rest stats
NOTE: Coverage values are NOT cross-sample comparable in absolute terms,
      because each MAG comes from its own assembly. Use cautiously.
"""
import gzip, re, os
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

MAG_DIR = Path("/data/data/Upcycling/MAGs_FASTA_files")
OUT = Path("/data/data/Upcycling/research/additional/C6_abundance_proxy")
OUT.mkdir(parents=True, exist_ok=True)

HEROES = {"S13","S16","S23","C22","M1","S26"}

# Source mapping: prefix → source
def source_of(name):
    p = name[0].upper()
    return {"C":"cattle","S":"swine","M":"sheep","V":"poultry"}.get(p, "unknown")

rows = []
for fa_gz in sorted(MAG_DIR.glob("*.fasta.gz")):
    mag = fa_gz.stem.replace(".fasta","")
    contigs = []
    with gzip.open(fa_gz,"rt") as fh:
        cur_len = 0; cur_cov = None
        for line in fh:
            if line.startswith(">"):
                if cur_cov is not None:
                    contigs.append((cur_len, cur_cov))
                m = re.search(r"length_(\d+)_cov_([\d.]+)", line)
                if m:
                    cur_len = int(m.group(1)); cur_cov = float(m.group(2))
                else:
                    cur_len = 0; cur_cov = None
            else:
                if cur_cov is None:
                    cur_len += len(line.strip())
        if cur_cov is not None:
            contigs.append((cur_len, cur_cov))
    if not contigs:
        continue
    lens = np.array([c[0] for c in contigs])
    covs = np.array([c[1] for c in contigs])
    total_len = lens.sum()
    weighted_cov = float((lens*covs).sum() / total_len) if total_len else None
    rows.append({
        "MAG": mag,
        "source": source_of(mag),
        "is_hero": mag in HEROES,
        "n_contigs": len(contigs),
        "total_length_bp": int(total_len),
        "median_cov": float(np.median(covs)),
        "mean_cov": float(covs.mean()),
        "length_weighted_cov": weighted_cov,
        "max_cov": float(covs.max()),
        "min_cov": float(covs.min()),
    })

df = pd.DataFrame(rows)
df.to_csv(OUT/"abundance_proxy_per_MAG.csv", index=False)
print(f"wrote {OUT/'abundance_proxy_per_MAG.csv'}: {df.shape}")

# Hero vs rest
res = []
for col in ["length_weighted_cov","median_cov","mean_cov"]:
    h = df[df.is_hero][col].dropna().values
    r = df[~df.is_hero][col].dropna().values
    if len(h)>=2 and len(r)>=2:
        u,p = mannwhitneyu(h,r,alternative="two-sided")
        res.append({"metric":col,"hero_mean":float(np.mean(h)),"rest_mean":float(np.mean(r)),
                    "hero_median":float(np.median(h)),"rest_median":float(np.median(r)),
                    "MWU_p":p})
pd.DataFrame(res).to_csv(OUT/"abundance_proxy_hero_vs_rest.csv", index=False)
print("done C6 hero vs rest")

# Per source breakdown
src = df.groupby("source")["length_weighted_cov"].agg(["count","mean","median","std"]).reset_index()
src.to_csv(OUT/"abundance_proxy_per_source.csv", index=False)
print(src)
