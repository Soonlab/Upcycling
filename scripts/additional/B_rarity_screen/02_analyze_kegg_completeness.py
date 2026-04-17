#!/usr/bin/env python3
"""B step 2: screen MGnify livestock MAG catalogs for MICP-complete profile.

Input: kegg_completeness.tar.gz per catalog (already downloaded).
These archives contain KEGG module completeness per pangenome / species cluster.

Look for genomes/species clusters that have:
  - M00029 (urease, ureABCDEFG accessory) — or we proxy via K01428 (ureC)
  - Carbonic anhydrase K01672/K01673/K18246

And compute:
  - How many species clusters have complete urease module + CA
  - Rarity within each livestock biome
"""
from __future__ import annotations
import os, glob, tarfile, re
import pandas as pd

ROOT = "/data/data/Upcycling/research/additional/B_rarity_screen/mgnify"
OUT = "/data/data/Upcycling/research/additional/B_rarity_screen"

# KEGG orthologs of interest
KO_INTEREST = {
    "K01428": "ureC_urease_alpha",
    "K01429": "ureB_urease_beta",
    "K01430": "ureA_urease_gamma",
    "K03187": "ureE",
    "K03188": "ureF",
    "K03189": "ureG",
    "K03190": "ureD",
    "K01672": "cah_alpha_CA",   # alpha carbonic anhydrase
    "K01673": "canA_gamma_CA",  # gamma
    "K18246": "cynT_beta_CA",   # beta
}
MODULE_INTEREST = {
    "M00029": "urease_module",
}

def parse_kegg_archive(tgz_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (ko_per_genome, module_per_genome). File format varies; peek first."""
    genomes = {}
    ko_rows = []
    mod_rows = []
    print(f"[B] opening {tgz_path}")
    with tarfile.open(tgz_path) as tar:
        names = tar.getnames()
        print(f"  {len(names)} files in archive; sample:")
        for n in names[:5]:
            print(f"    {n}")
        # typical format: one TSV per genome with KO/module completeness
        for n in names:
            if not (n.endswith(".tsv") or n.endswith(".txt")): continue
            base = os.path.basename(n).replace(".tsv","").replace(".txt","")
            fh = tar.extractfile(n)
            if fh is None: continue
            try:
                txt = fh.read().decode()
            except Exception:
                continue
            # parse - expect 2-col: KO/module \t completeness
            for line in txt.splitlines():
                parts = line.strip().split("\t")
                if len(parts) < 2: continue
                key = parts[0]
                try:
                    val = float(parts[1])
                except ValueError:
                    continue
                if re.match(r"^K\d{5}$", key):
                    if key in KO_INTEREST:
                        ko_rows.append({"genome": base, "KO": key, "presence": val})
                elif re.match(r"^M\d{5}$", key):
                    if key in MODULE_INTEREST:
                        mod_rows.append({"genome": base, "module": key, "completeness": val})
    return pd.DataFrame(ko_rows), pd.DataFrame(mod_rows)

summaries = []
for cat_dir in sorted(glob.glob(f"{ROOT}/*")):
    cat = os.path.basename(cat_dir)
    tgz = f"{cat_dir}/kegg_completeness.tar.gz"
    if not os.path.exists(tgz):
        print(f"[B] missing {tgz}, skip")
        continue
    ko_df, mod_df = parse_kegg_archive(tgz)
    if ko_df.empty and mod_df.empty:
        print(f"[B] {cat}: no parsed rows — inspect archive format manually")
        continue
    # genome-level summary
    genomes = set(ko_df["genome"]).union(set(mod_df["genome"]))
    micp_complete = []
    for g in genomes:
        kos = set(ko_df[ko_df.genome==g]["KO"])
        has_ureC = "K01428" in kos
        has_ureB = "K01429" in kos
        has_ureA = "K01430" in kos
        has_CA = any(k in kos for k in ["K01672","K01673","K18246"])
        mod_val = mod_df[(mod_df.genome==g) & (mod_df.module=="M00029")]["completeness"]
        urease_mod = float(mod_val.iloc[0]) if len(mod_val) else None
        complete = (has_ureC and has_ureB and has_ureA and has_CA) or \
                   (urease_mod is not None and urease_mod >= 0.75 and has_CA)
        micp_complete.append({
            "genome": g, "catalog": cat,
            "has_ureC": has_ureC, "has_ureB": has_ureB, "has_ureA": has_ureA,
            "has_CA": has_CA, "urease_module_completeness": urease_mod,
            "MICP_complete": complete,
        })
    df = pd.DataFrame(micp_complete)
    df.to_csv(f"{OUT}/{cat}_MICP_screen.csv", index=False)
    summaries.append({
        "catalog": cat,
        "n_genomes": len(df),
        "n_with_ureC": df.has_ureC.sum(),
        "n_with_urease_core": ((df.has_ureC & df.has_ureB & df.has_ureA)).sum(),
        "n_with_CA": df.has_CA.sum(),
        "n_MICP_complete": df.MICP_complete.sum(),
        "pct_MICP_complete": round(100*df.MICP_complete.mean(), 2),
    })
    print(f"[B] {cat}: n={len(df)}  MICP-complete={df.MICP_complete.sum()}  ({100*df.MICP_complete.mean():.2f}%)")

summary_df = pd.DataFrame(summaries)
summary_df.to_csv(f"{OUT}/MICP_rarity_summary_mgnify.csv", index=False)
print("\n[B] Summary across MGnify livestock catalogs:")
print(summary_df.to_string(index=False))
