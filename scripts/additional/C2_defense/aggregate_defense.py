#!/usr/bin/env python3
"""Aggregate DefenseFinder + minced into per-MAG counts and hero-vs-rest."""
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu

OUT = Path("/data/data/Upcycling/research/additional/C2_defense")
DF = OUT/"defensefinder"
CR = OUT/"crispr_minced"
HEROES = {"S13","S16","S23","C22","M1","S26"}

rows = []
for d in sorted(DF.iterdir()):
    if not d.is_dir(): continue
    mag = d.name
    sys_tsv = d/"defense_finder_systems.tsv"
    n_sys = 0
    sys_types = []
    if sys_tsv.exists():
        try:
            sdf = pd.read_csv(sys_tsv, sep="\t")
            n_sys = len(sdf)
            if "type" in sdf.columns:
                sys_types = sdf["type"].astype(str).tolist()
            elif "subtype" in sdf.columns:
                sys_types = sdf["subtype"].astype(str).tolist()
        except Exception:
            pass
    # CRISPR
    cr_gff = CR/mag/f"{mag}.minced.gff"
    n_arrays = 0
    if cr_gff.exists():
        with open(cr_gff) as fh:
            n_arrays = sum(1 for line in fh if "\tCRISPR\t" in line)
    rec = {
        "MAG": mag,
        "is_hero": mag in HEROES,
        "n_defense_systems": n_sys,
        "n_crispr_arrays": n_arrays,
        "defense_types": ";".join(sorted(set(sys_types)))
    }
    for t in sys_types:
        key = f"DF_{t.replace(' ','_')}"
        rec[key] = rec.get(key, 0) + 1
    rows.append(rec)

df = pd.DataFrame(rows).fillna(0)
df.to_csv(OUT/"defense_per_MAG.csv", index=False)
print(f"wrote {OUT/'defense_per_MAG.csv'}: {df.shape}")

hero = df[df.is_hero]; rest = df[~df.is_hero]
stat_rows = []
for col in ["n_defense_systems","n_crispr_arrays"] + [c for c in df.columns if c.startswith("DF_")]:
    h = hero[col].astype(float).values
    r = rest[col].astype(float).values
    if len(h)<2 or len(r)<2: continue
    try:
        u,p = mannwhitneyu(h,r,alternative="two-sided")
    except Exception:
        u,p = float("nan"),float("nan")
    stat_rows.append({"metric":col,"hero_mean":h.mean(),"rest_mean":r.mean(),"MWU_p":p})
pd.DataFrame(stat_rows).sort_values("MWU_p").to_csv(OUT/"defense_hero_vs_rest.csv", index=False)
print(f"wrote {OUT/'defense_hero_vs_rest.csv'}")
