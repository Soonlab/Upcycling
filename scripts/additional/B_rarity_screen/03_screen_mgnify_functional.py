#!/usr/bin/env python3
"""B step 3: Screen MGnify livestock MAG catalogs for MICP-complete profile
using functional_profiles.tar.gz (per-gene KEGG/Pfam annotations).

Output: per-species MICP completeness flags + rarity summary.
"""
from __future__ import annotations
import os, tarfile, re
import pandas as pd
from collections import defaultdict

ROOT = "/data/data/Upcycling/research/additional/B_rarity_screen/mgnify"
OUT = "/data/data/Upcycling/research/additional/B_rarity_screen"

# KEGG orthologs
UREASE_KOS = {
    "K01428": "ureC", "K01429": "ureB", "K01430": "ureA",
    "K14048": "ureAB_fusion",
    "K03187": "ureE", "K03188": "ureF", "K03189": "ureG",
    "K03190": "ureD",
}
CA_KOS = {
    "K01672": "cah_alpha_CA",
    "K01673": "canA_gamma_CA",
    "K18246": "cynT_beta_CA",
    "K01726": "canA_homolog",  # some variants
}
# Pfam codes (as backup)
UREASE_PFAMS = {"PF00449": "UreC_alpha", "PF00699": "UreB_beta_gamma"}
CA_PFAMS = {"PF00484": "CA_beta", "PF00194": "CA_alpha", "PF00988": "CA_gamma"}

def screen_cluster(lines):
    """Parse a functional_profiles tsv (list of text lines or iterator) for one cluster."""
    kos = set()
    pfams = set()
    contigs = defaultdict(set)  # contig -> set of (KO/Pfam for operon-adjacency test)
    for ln in lines:
        if ln.startswith("#"): continue
        parts = ln.rstrip("\n").split("\t")
        if len(parts) < 8: continue
        contig, gid, start, end, strand, ko, caz, pfam = parts[:8]
        if ko and ko != "-":
            for k in ko.split(","):
                k = k.strip()
                if k: kos.add(k)
                if k in UREASE_KOS or k in CA_KOS:
                    contigs[contig].add(k)
        if pfam and pfam != "-":
            for p in pfam.split(","):
                p = p.strip().split(".")[0]
                if p: pfams.add(p)
                if p in UREASE_PFAMS or p in CA_PFAMS:
                    contigs[contig].add(p)
    # per cluster summary
    has_ureC = "K01428" in kos or "PF00449" in pfams
    has_ureB = "K01429" in kos or "PF00699" in pfams  # beta/gamma
    has_ureA = "K01430" in kos or "K14048" in kos
    has_ureE = "K03187" in kos
    has_ureF = "K03188" in kos
    has_ureG = "K03189" in kos
    has_ureD = "K03190" in kos
    ure_acc = sum([has_ureE, has_ureF, has_ureG, has_ureD])
    has_CA_alpha = "K01672" in kos or "PF00194" in pfams
    has_CA_beta  = "K18246" in kos or "PF00484" in pfams
    has_CA_gamma = "K01673" in kos or "PF00988" in pfams
    has_CA = has_CA_alpha or has_CA_beta or has_CA_gamma
    # single-contig coreness: are ureC + any ureB/A on same contig?
    ureC_operon_contig = 0
    ureC_plus_CA_contig = 0
    for cg, keys in contigs.items():
        has_c = any(k in keys for k in ["K01428","PF00449"])
        has_b = any(k in keys for k in ["K01429","PF00699"])
        has_a = "K01430" in keys or "K14048" in keys
        has_any_ca = any(k in keys for k in ["K01672","K01673","K18246","PF00194","PF00484","PF00988"])
        if has_c and has_b and has_a:
            ureC_operon_contig = 1
        if has_c and has_any_ca:
            ureC_plus_CA_contig = 1
    return {
        "ureC": int(has_ureC), "ureB": int(has_ureB), "ureA": int(has_ureA),
        "urease_accessory_count": ure_acc,
        "CA_any": int(has_CA),
        "CA_alpha": int(has_CA_alpha), "CA_beta": int(has_CA_beta), "CA_gamma": int(has_CA_gamma),
        "urease_core_complete": int(has_ureC and has_ureB and has_ureA),
        "urease_operon_single_contig": ureC_operon_contig,
        "ureC_plus_CA_single_contig": ureC_plus_CA_contig,
        "MICP_gene_complete": int(has_ureC and has_ureB and has_ureA and has_CA),
    }

# load metadata for each catalog
summaries = []
all_genomes = []

for catalog in ["cow-rumen", "sheep-rumen", "pig-gut", "chicken-gut"]:
    cd = f"{ROOT}/{catalog}"
    tgz = f"{cd}/functional_profiles.tar.gz"
    meta_path = f"{cd}/genomes-all_metadata.tsv"
    print(f"[B] === {catalog} ===")

    # metadata for taxonomy
    meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
    if "Genome" in meta.columns:
        # map species-cluster rep to metadata
        meta = meta.set_index("Genome")
    else:
        meta = pd.DataFrame()

    results = []
    with tarfile.open(tgz) as tar:
        members = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(".tsv")]
        for i, m in enumerate(members):
            cluster_rep = os.path.basename(m.name).replace("_clstr.tsv","")
            fh = tar.extractfile(m)
            if fh is None: continue
            try:
                text = fh.read().decode("utf-8", errors="ignore")
            except Exception:
                continue
            screen = screen_cluster(text.splitlines())
            screen["cluster_rep"] = cluster_rep
            screen["catalog"] = catalog
            # taxonomy if available
            if cluster_rep in meta.index:
                r = meta.loc[cluster_rep]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                screen["Lineage"] = r.get("Lineage", None)
                screen["Genome_type"] = r.get("Genome_type", None)
                screen["Completeness"] = r.get("Completeness", None)
                screen["Contamination"] = r.get("Contamination", None)
            results.append(screen)
            if (i+1) % 500 == 0:
                print(f"  ... processed {i+1}/{len(members)} clusters")

    df = pd.DataFrame(results)
    df.to_csv(f"{OUT}/{catalog}_MICP_profile.csv", index=False)

    n = len(df)
    if n == 0:
        continue
    summaries.append({
        "catalog": catalog,
        "n_species_clusters": n,
        "n_with_ureC": int(df.ureC.sum()),
        "n_urease_core": int(df.urease_core_complete.sum()),
        "n_urease_operon_single_contig": int(df.urease_operon_single_contig.sum()),
        "n_with_CA": int(df.CA_any.sum()),
        "n_MICP_gene_complete": int(df.MICP_gene_complete.sum()),
        "n_MICP_single_contig_ureC_CA": int(df.ureC_plus_CA_single_contig.sum()),
        "pct_MICP_gene_complete": round(100*df.MICP_gene_complete.mean(), 3),
        "pct_MICP_single_contig": round(100*df.ureC_plus_CA_single_contig.mean(), 3),
    })
    all_genomes.append(df)

    # print Sphingobacterium + Pseudomonas breakdown
    if "Lineage" in df.columns and df.Lineage.notna().any():
        sph = df[df.Lineage.str.contains("Sphingobacterium", na=False)]
        pse = df[df.Lineage.str.contains("Pseudomonas", na=False)]
        print(f"  total species clusters: {n}")
        print(f"    with ureC:              {df.ureC.sum()} ({100*df.ureC.mean():.2f}%)")
        print(f"    urease core complete:   {df.urease_core_complete.sum()} ({100*df.urease_core_complete.mean():.2f}%)")
        print(f"    MICP gene-complete:     {df.MICP_gene_complete.sum()} ({100*df.MICP_gene_complete.mean():.2f}%)")
        print(f"    single-contig urease+CA:{df.ureC_plus_CA_single_contig.sum()} ({100*df.ureC_plus_CA_single_contig.mean():.3f}%)")
        if len(sph):
            print(f"    Sphingobacterium: n={len(sph)}, MICP complete={sph.MICP_gene_complete.sum()}")
        if len(pse):
            print(f"    Pseudomonas:      n={len(pse)}, MICP complete={pse.MICP_gene_complete.sum()}")

# overall
pd.concat(all_genomes, ignore_index=True).to_csv(f"{OUT}/all_mgnify_MICP_profiles.csv", index=False)
summary_df = pd.DataFrame(summaries)
summary_df.to_csv(f"{OUT}/MICP_rarity_summary_mgnify.csv", index=False)
print("\n[B] Summary across MGnify livestock catalogs:")
print(summary_df.to_string(index=False))

# pooled totals
total_species = summary_df.n_species_clusters.sum()
total_micp = summary_df.n_MICP_gene_complete.sum()
total_single = summary_df.n_MICP_single_contig_ureC_CA.sum()
print(f"\n[B] POOLED (4 livestock biomes): {total_species} species clusters")
print(f"    MICP gene-complete:        {total_micp}/{total_species} = {100*total_micp/total_species:.3f}%")
print(f"    single-contig ureC+CA:     {total_single}/{total_species} = {100*total_single/total_species:.3f}%")
print(f"    i.e. ~1 in {total_species//max(total_single,1)} livestock MAGs has the convergent architecture")
