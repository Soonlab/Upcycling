#!/usr/bin/env python3
"""
Revision step #4 — rigorous CAZyme analysis from DRAM dbCAN hits.

Replaces the keyword-based CAZy proxy with proper CAZy-family counts
(GH, GT, PL, CE, AA, CBM) extracted from DRAM's all_annotations.tsv.
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE="/data/data/Upcycling"
OUT=f"{BASE}/research/revision"
os.makedirs(OUT, exist_ok=True)

ann = pd.read_csv("/data/pangenome_work/dram_output/all_annotations.tsv",
                  sep="\t", usecols=["fasta","cazy_ids","cazy_hits","cazy_best_hit"])
print(f"[info] loaded {len(ann)} annotation rows")
ann = ann[ann["cazy_best_hit"].notna() | ann["cazy_ids"].notna()]
print(f"[info] CAZy-annotated rows: {len(ann)}")

# Canonical CAZy class letters
CLS = ["GH","GT","PL","CE","AA","CBM"]
def family(x):
    if pd.isna(x): return None
    m = re.search(r"\b(GH|GT|PL|CE|AA|CBM)\d+", str(x))
    return m.group(1) if m else None
def family_full(x):
    if pd.isna(x): return None
    m = re.search(r"\b((?:GH|GT|PL|CE|AA|CBM)\d+)", str(x))
    return m.group(1) if m else None

ann["cls"]    = ann.apply(lambda r: family(r["cazy_best_hit"]) or family(r["cazy_ids"]), axis=1)
ann["family"] = ann.apply(lambda r: family_full(r["cazy_best_hit"]) or family_full(r["cazy_ids"]), axis=1)
ann = ann[ann["cls"].notna()]
print(f"[info] rows with parsed CAZy class: {len(ann)}")

# Per-MAG class counts
cls_cnt = ann.groupby(["fasta","cls"]).size().unstack(fill_value=0)
for c in CLS:
    if c not in cls_cnt.columns: cls_cnt[c]=0
cls_cnt = cls_cnt[CLS]
cls_cnt["Total_CAZymes"] = cls_cnt.sum(axis=1)

# Per-MAG family counts
fam_cnt = ann.groupby(["fasta","family"]).size().unstack(fill_value=0)

# Normalise per-MAG by CDS count (from genome_category_counts)
counts = pd.read_csv(f"{BASE}/research/extra/gene_category_counts.csv")
cds = counts.set_index("Sample")["CDS_total"]
cls_norm = cls_cnt.div(cds, axis=0).multiply(1000.0).round(3)
cls_norm["CDS_total"] = cds
cls_cnt.to_csv(f"{OUT}/dbCAN_class_counts.csv")
cls_norm.to_csv(f"{OUT}/dbCAN_class_per1k_CDS.csv")
fam_cnt.to_csv(f"{OUT}/dbCAN_family_counts.csv")
print("[ok] dbCAN class/family tables written")

# Load hero labels
counts["Hero"] = counts["Hero"].astype(bool)
hero_list = counts.loc[counts["Hero"],"Sample"].tolist()
print(f"Hero set: {hero_list}")

# Hero-vs-rest class comparison
hero_cls = cls_norm.loc[cls_norm.index.isin(hero_list), CLS].mean()
rest_cls = cls_norm.loc[~cls_norm.index.isin(hero_list), CLS].mean()
cmp = pd.DataFrame({"Hero_mean_per1kCDS":hero_cls, "Rest_mean_per1kCDS":rest_cls,
                    "Fold_change": hero_cls/rest_cls.replace(0,np.nan)}).round(3)
cmp.to_csv(f"{OUT}/dbCAN_hero_vs_rest_class.csv")
print("\n=== dbCAN class-level Hero vs Rest ===")
print(cmp)

# Per-family Hero vs Rest (top enriched)
fam_norm = fam_cnt.div(cds, axis=0).multiply(1000.0)
hero_fam = fam_norm.loc[fam_norm.index.isin(hero_list)].mean()
rest_fam = fam_norm.loc[~fam_norm.index.isin(hero_list)].mean()
fam_cmp = pd.DataFrame({"Hero_mean":hero_fam, "Rest_mean":rest_fam,
                        "Fold_change": hero_fam / rest_fam.replace(0,np.nan),
                        "Hero_present_in_N": (fam_cnt.loc[fam_cnt.index.isin(hero_list)]>0).sum()})
fam_cmp = fam_cmp.sort_values("Fold_change", ascending=False)
fam_cmp.to_csv(f"{OUT}/dbCAN_hero_vs_rest_family.csv")
print(f"\nTotal CAZy families detected: {len(fam_cmp)}")
print("Top 15 hero-enriched families:")
print(fam_cmp.head(15).round(3))

# ---- Figure: class-level bar chart ---------------------------------------
fig, ax = plt.subplots(figsize=(8,5))
x=np.arange(len(CLS)); w=0.4
ax.bar(x-w/2, cmp["Hero_mean_per1kCDS"], w, label=f"Hero (n={len(hero_list)})",
       color="#c0392b", edgecolor="black", lw=0.4)
ax.bar(x+w/2, cmp["Rest_mean_per1kCDS"], w, label=f"Rest (n={len(cls_norm)-len(hero_list)})",
       color="#7f8c8d", edgecolor="black", lw=0.4)
ax.set_xticks(x); ax.set_xticklabels(CLS)
ax.set_ylabel("CAZymes per 1k CDS")
ax.set_title("Rigorous dbCAN class-level profile (DRAM-based)", fontsize=11)
ax.legend(frameon=False)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_dbCAN_classes.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_dbCAN_classes.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_dbCAN_classes")

# ---- Figure: top enriched families ---------------------------------------
top = fam_cmp.dropna().query("Hero_present_in_N>=3").head(20)
fig, ax = plt.subplots(figsize=(8, max(4, len(top)*0.28)))
ax.barh(top.index[::-1], top["Fold_change"].values[::-1],
        color="#c0392b", edgecolor="black", lw=0.4)
ax.axvline(1.0, ls="--", color="gray", lw=0.8)
ax.set_xlabel("Fold change (hero / rest)")
ax.set_title("Top hero-enriched CAZy families (dbCAN; ≥3 hero MAGs)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_dbCAN_topFamilies.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_dbCAN_topFamilies.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_dbCAN_topFamilies")
