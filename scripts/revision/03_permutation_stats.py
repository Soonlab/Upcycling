#!/usr/bin/env python3
"""
Revision step #2 — Permutation tests and bootstrap 95% CIs for every
trait-module enrichment reported in the manuscript.

For each (category, subcategory):
  observed fold change  = mean(Hero) / mean(Rest) using per-1k-CDS counts
  permutation p-value   = fraction of 10,000 random label shuffles giving
                          a fold change >= observed (one-sided, right tail)
  bootstrap 95% CI for hero mean, rest mean and fold change
  Wilcoxon rank-sum two-sided p-value as a parametric-free sanity check
"""
import os, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE="/data/data/Upcycling"
OUT=f"{BASE}/research/revision"
N_PERM = 10000
N_BOOT = 2000
RNG = np.random.default_rng(42)

counts = pd.read_csv(f"{BASE}/research/extra/gene_category_counts.csv")
cds    = counts["CDS_total"].values
subcols = [c for c in counts.columns if "::" in c]
# per-1k-CDS matrix
norm = counts.copy()
for c in subcols: norm[c] = counts[c] / counts["CDS_total"] * 1000.0
hero_mask = counts["Hero"].astype(bool).values

print(f"[info] N MAGs = {len(counts)} (hero n={hero_mask.sum()}, rest n={(~hero_mask).sum()})")
print(f"[info] Testing {len(subcols)} trait subcategories")

def fold_change(x, mask):
    h = x[mask].mean(); r = x[~mask].mean()
    return h / r if r > 0 else np.nan

def perm_p(x, mask, n):
    """one-sided right-tail permutation p of mean(hero) - mean(rest)."""
    obs = x[mask].mean() - x[~mask].mean()
    if obs <= 0: return 1.0  # right-tail only meaningful if hero > rest
    idx = np.arange(len(x))
    hits = 0
    k = mask.sum()
    for _ in range(n):
        RNG.shuffle(idx)
        new_mask = np.zeros(len(x), dtype=bool); new_mask[idx[:k]] = True
        d = x[new_mask].mean() - x[~new_mask].mean()
        if d >= obs: hits += 1
    return (hits+1)/(n+1)

def boot_ci(x, mask, n):
    h_idx = np.where(mask)[0]; r_idx = np.where(~mask)[0]
    hm, rm, fc = [], [], []
    for _ in range(n):
        hs = RNG.choice(h_idx, size=len(h_idx), replace=True)
        rs = RNG.choice(r_idx, size=len(r_idx), replace=True)
        h = x[hs].mean(); r = x[rs].mean()
        hm.append(h); rm.append(r); fc.append(h/r if r>0 else np.nan)
    hm=np.array(hm); rm=np.array(rm); fc=np.array(fc)
    return (np.percentile(hm,2.5), np.percentile(hm,97.5),
            np.percentile(rm,2.5), np.percentile(rm,97.5),
            np.nanpercentile(fc,2.5), np.nanpercentile(fc,97.5))

rows=[]
for col in subcols:
    x = norm[col].values
    cat, sub = col.split("::")
    h_mean = x[hero_mask].mean(); r_mean = x[~hero_mask].mean()
    fc = fold_change(x, hero_mask)
    # Wilcoxon rank-sum (two-sided)
    try:
        w = stats.mannwhitneyu(x[hero_mask], x[~hero_mask], alternative="greater")
        mwu_p = w.pvalue
    except ValueError:
        mwu_p = np.nan
    perm = perm_p(x, hero_mask, N_PERM) if h_mean>0 else 1.0
    hlo,hhi,rlo,rhi,flo,fhi = boot_ci(x, hero_mask, N_BOOT)
    rows.append({"Category":cat,"Subcategory":sub,
                 "Hero_mean_per1kCDS":round(h_mean,3),
                 "Rest_mean_per1kCDS":round(r_mean,3),
                 "Fold_change":round(fc,3) if not np.isnan(fc) else np.nan,
                 "Hero_CI95":f"[{hlo:.3f}, {hhi:.3f}]",
                 "Rest_CI95":f"[{rlo:.3f}, {rhi:.3f}]",
                 "Fold_change_CI95":f"[{flo:.2f}, {fhi:.2f}]" if not np.isnan(flo) else "NA",
                 "Permutation_p": round(perm,4),
                 "MannWhitney_greater_p": round(mwu_p,4) if not np.isnan(mwu_p) else np.nan})

res = pd.DataFrame(rows)
# Benjamini-Hochberg FDR on permutation p
from scipy.stats import false_discovery_control
res["Permutation_q_BH"] = false_discovery_control(res["Permutation_p"].values).round(4)
res["MannWhitney_q_BH"] = false_discovery_control(res["MannWhitney_greater_p"].fillna(1).values).round(4)
res = res.sort_values(["Permutation_q_BH","Fold_change"], ascending=[True,False])
res.to_csv(f"{OUT}/Hero_vs_Rest_permutation_stats.csv", index=False)
print("\n=== Top 12 significantly enriched modules (BH-FDR q<0.1) ===")
sig = res[res["Permutation_q_BH"]<0.1].head(12)
print(sig.to_string(index=False))
print(f"\n[info] modules with q<0.05: {(res['Permutation_q_BH']<0.05).sum()}")
print(f"[info] modules with q<0.10: {(res['Permutation_q_BH']<0.10).sum()}")

# ---- forest-style plot -----------------------------------------------------
plot_df = res[res["Fold_change"]>1].dropna(subset=["Fold_change"]).copy()
plot_df["label"] = plot_df["Category"]+": "+plot_df["Subcategory"]
plot_df = plot_df.sort_values("Fold_change", ascending=True).tail(25)
# extract CI numerics
ci = plot_df["Fold_change_CI95"].str.strip("[]").str.split(", ", expand=True).astype(float)
plot_df["lo"] = ci[0]; plot_df["hi"] = ci[1]

fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df)*0.32)))
y = np.arange(len(plot_df))
colors = ["#c0392b" if q<0.05 else "#e67e22" if q<0.10 else "#95a5a6"
          for q in plot_df["Permutation_q_BH"]]
ax.errorbar(plot_df["Fold_change"], y,
            xerr=[plot_df["Fold_change"]-plot_df["lo"], plot_df["hi"]-plot_df["Fold_change"]],
            fmt="o", color="black", ecolor="gray", elinewidth=0.6, markersize=0, capsize=2)
ax.scatter(plot_df["Fold_change"], y, s=60, c=colors, edgecolor="black", lw=0.4, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(plot_df["label"], fontsize=7)
ax.axvline(1.0, ls="--", color="black", lw=0.6)
ax.set_xscale("log")
ax.set_xlabel("Fold change (hero / rest) — log scale, 95% bootstrap CI")
ax.set_title("Hero-vs-Rest enrichment with permutation test (BH-FDR)",
             fontsize=10)
# legend
from matplotlib.patches import Patch
handles=[Patch(color="#c0392b",label="q < 0.05"),
         Patch(color="#e67e22",label="q < 0.10"),
         Patch(color="#95a5a6",label="n.s.")]
ax.legend(handles=handles, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_Permutation_forest.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_Permutation_forest.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_Permutation_forest")
