#!/usr/bin/env python3
"""Aggregate antiSMASH per-MAG outputs into one BGC table + hero-vs-rest stats."""
import json, re, os, sys, glob
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu

PERMAG = Path("/data/data/Upcycling/research/additional/C1_antismash/per_mag")
OUT = Path("/data/data/Upcycling/research/additional/C1_antismash")
HEROES = {"S13","S16","S23","C22","M1","S26"}

rows = []
for d in sorted(PERMAG.iterdir()):
    if not d.is_dir():
        continue
    mag = d.name
    rec = {"MAG": mag, "is_hero": mag in HEROES, "n_regions": 0}
    cls_counts = {}
    # antiSMASH 7 writes a JSON per input record
    for jf in d.glob("*.json"):
        try:
            with open(jf) as fh:
                data = json.load(fh)
        except Exception:
            continue
        records = data.get("records", []) if isinstance(data, dict) else []
        for r in records:
            for area in r.get("areas", []):
                rec["n_regions"] += 1
                products = area.get("products", []) or []
                for p in products:
                    cls_counts[p] = cls_counts.get(p, 0) + 1
    rec["classes"] = ";".join(sorted(cls_counts)) if cls_counts else ""
    rec.update({f"BGC_{k}": v for k,v in cls_counts.items()})
    rows.append(rec)

df = pd.DataFrame(rows).fillna(0)
out_csv = OUT / "antismash_per_MAG.csv"
df.to_csv(out_csv, index=False)
print(f"wrote {out_csv}: {df.shape}")

# Hero vs rest stats
hero = df[df.is_hero]
rest = df[~df.is_hero]
stat_rows = []
for col in [c for c in df.columns if c.startswith("BGC_")] + ["n_regions"]:
    h = hero[col].astype(float).values
    r = rest[col].astype(float).values
    if len(h)<2 or len(r)<2:
        continue
    try:
        u,p = mannwhitneyu(h, r, alternative="two-sided")
    except Exception:
        u,p = float("nan"), float("nan")
    stat_rows.append({
        "metric": col,
        "hero_mean": h.mean(), "hero_median": float(pd.Series(h).median()),
        "rest_mean": r.mean(), "rest_median": float(pd.Series(r).median()),
        "MWU_p": p
    })
sdf = pd.DataFrame(stat_rows).sort_values("MWU_p")
sdf.to_csv(OUT/"antismash_hero_vs_rest.csv", index=False)
print(f"wrote {OUT/'antismash_hero_vs_rest.csv'}")
print(sdf.head(15).to_string())
