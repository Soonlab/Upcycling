#!/usr/bin/env python3
"""
Tier 1 + Tier 2 additional analyses for the 111-MAG MICP study.

Outputs a tidy folder of CSVs + figures under /data/data/Upcycling/research/extra/.

Modules:
  T1a  Biofilm / EPS gene profile
  T1b  Ammonia detox & N-assimilation (glnA, gltBD, gdhA, amtB...)
  T1c  HGT signals around the ure-cah cluster (transposase/integrase/IS + GC deviation)
  T2a  Alkaline / osmotic stress tolerance
  T2b  CAZyme proxy (Bakta product keyword scan for GH/GT/PL/CE/CBM-like hits)
  T2c  Heavy-metal & antibiotic resistance keyword scan
"""
import os, re, gzip, glob
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = "/data/data/Upcycling"
BAKTA = f"{BASE}/MAGs_FASTA_files/bakta_results"
OUT  = f"{BASE}/research/extra"
os.makedirs(OUT, exist_ok=True)

samples = sorted(os.listdir(BAKTA))
print(f"[info] {len(samples)} MAGs")

# Load taxonomy
gtdb = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
def g_(s): m=re.search(r"g__([^;]*)",str(s)); return m.group(1) if m and m.group(1) else "Unclassified"
gtdb["Genus"] = gtdb["classification"].apply(g_)
tax = gtdb.set_index("user_genome")["Genus"]

HERO = {"S13","S16","S23","S26","C22","M1"}

# -----------------------------------------------------------------------------
# Gene keyword dictionaries
# -----------------------------------------------------------------------------
CATEGORIES = {
    "Biofilm_EPS": {
        "polysaccharide_intercellular": ["pgaA","pgaB","pgaC","pgaD","ica"],
        "pel":   ["pelA","pelB","pelC","pelD","pelE","pelF","pelG"],
        "psl":   ["pslA","pslB","pslC","pslD","pslE","pslF","pslG","pslH","pslI"],
        "cellulose": ["bcsA","bcsB","bcsC","bcsZ","cellulose synthase"],
        "alginate": ["algA","algD","algE","algG","algK","algL","alginate"],
        "curli":    ["csgA","csgB","csgC","csgD","csgE","csgF","csgG","curli"],
        "colanic":  ["wcaA","wcaB","wcaC","wcaD","wcaE","wcaF"],
        "capsule":  ["wza","wzb","wzc","capB","capC","capD","cpsB","cpsC"],
        "adhesin":  ["fimbri","adhesin","autotransporter adhesin","hemolysin-type calcium"],
        "quorum":   ["luxR","luxS","luxI","agrA","agrB","agrC","agrD","rhlR","lasR","N-acyl-L-homoserine"],
    },
    "Ammonia_N": {
        "urea_transport": ["urtA","urtB","urtC","urtD","urtE","urea transporter","Urea ABC"],
        "GS_GOGAT": ["glnA","gltB","gltD","glutamine synthetase","glutamate synthase"],
        "GDH": ["gdhA","glutamate dehydrogenase"],
        "ammonium_transport": ["amtB","ammonium transport","nrgA"],
        "PII_regulator": ["glnB","glnK","nitrogen regulatory protein P-II"],
        "nitrate_nitrite":["narG","narH","nasA","nasB","nirB","nirD","nitrate reductase","nitrite reductase"],
    },
    "HGT_markers": {
        "transposase": ["transposase","IS element","IS1","IS3","IS110","IS200","IS481","IS630","IS5"],
        "integrase":   ["integrase","recombinase XerC","XerD","tyrosine recombinase"],
        "phage":       ["phage","prophage","tail fiber","portal protein","terminase"],
        "plasmid":     ["plasmid replication","relaxase","mobilization protein","TraD","TraI"],
        "conjugation": ["T4SS","type IV secretion","VirB","VirD"],
    },
    "Alkaline_Osmo": {
        "Na_H_antiporter": ["nhaA","nhaB","nhaC","nhaD","Na(+)/H(+) antiporter","sodium/proton antiporter"],
        "Mrp_complex":     ["mrpA","mrpB","mrpC","mrpD","mrpE","mrpF","mrpG","multicomponent K+:H+"],
        "K_uptake":        ["kdpA","kdpB","kdpC","kdpD","kdpE","trkA","trkH"],
        "compatible_solute": ["opuA","opuB","opuC","opuD","proV","proW","proX","betA","betB",
                              "ectA","ectB","ectC","ectoine","glycine betaine"],
        "HSP_chaperone":   ["groEL","groES","dnaK","dnaJ","grpE","clpB","hsp70"],
        "oxidative":       ["katA","katB","katG","sodA","sodB","sodC","ahpC","ahpF"],
    },
    "CAZyme_proxy": {
        "glycoside_hydrolase": ["glycoside hydrolase","glycosyl hydrolase","GH family","cellulase",
                                "xylanase","chitinase","beta-glucosidase","beta-galactosidase",
                                "alpha-amylase","pullulanase","mannanase","pectinase"],
        "glycosyltransferase": ["glycosyltransferase","glycosyl transferase"],
        "polysaccharide_lyase":["pectate lyase","alginate lyase","polysaccharide lyase"],
        "carb_esterase":       ["carbohydrate esterase","acetyl xylan esterase","feruloyl esterase"],
        "carb_binding":        ["carbohydrate-binding","CBM"],
    },
    "MetalResist_AMR": {
        "metal_efflux": ["czcA","czcB","czcC","czcD","copA","copB","copC","copD","cusA","cusB","cusC",
                         "zntA","zntB","cadA","arsA","arsB","arsC","merA","merT","merP",
                         "heavy metal efflux","heavy metal translocating"],
        "metallothionein":["metallothionein","smtA","bmtA"],
        "beta_lactamase": ["beta-lactamase","metallo-beta-lactamase","penicillin-binding"],
        "efflux_MDR":     ["acrA","acrB","tolC","mexA","mexB","multidrug efflux","MATE efflux"],
        "aminoglycoside": ["aminoglycoside","aph(3')","aac(6')","ant(2'')","streptomycin"],
        "tet_mac":        ["tetracycline resistance","tet(","macrolide","erm(","mph("],
    },
}
ALL_KEYS = [(cat, sub, kw) for cat,d in CATEGORIES.items() for sub,kws in d.items() for kw in kws]

def parse_gff(path):
    rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS": continue
            attrs = dict(kv.split("=",1) for kv in p[8].split(";") if "=" in kv)
            rows.append({"contig":p[0],"start":int(p[3]),"end":int(p[4]),"strand":p[6],
                         "gene":attrs.get("gene",""),"product":attrs.get("product","")})
    return pd.DataFrame(rows)

# -----------------------------------------------------------------------------
# Scan all 111 MAGs — per-category hit counts
# -----------------------------------------------------------------------------
def classify_hits(df):
    """Return dict[(cat,sub)] = count, + list of (cat,sub,idx) for HGT distance calc."""
    prod = df["product"].str.lower().fillna("")
    gene = df["gene"].str.lower().fillna("")
    counts = defaultdict(int)
    tag    = [None]*len(df)
    for cat, subs in CATEGORIES.items():
        for sub, kws in subs.items():
            mask = np.zeros(len(df), dtype=bool)
            for kw in kws:
                kl = kw.lower()
                # word/segment match
                mask |= prod.str.contains(re.escape(kl), regex=True) | (gene == kl)
            counts[(cat,sub)] = int(mask.sum())
            for i,m in enumerate(mask):
                if m and tag[i] is None: tag[i] = (cat,sub)
    return counts, tag

per_sample_rows = []
all_gffs = {}
for samp in samples:
    gff = f"{BAKTA}/{samp}/{samp}.gff3"
    if not os.path.exists(gff): continue
    df = parse_gff(gff); all_gffs[samp] = df
    counts, _ = classify_hits(df)
    row = {"Sample": samp, "Genus": tax.get(samp,"Unclassified"),
           "Hero": samp in HERO, "CDS_total": len(df)}
    for (cat,sub), c in counts.items(): row[f"{cat}::{sub}"] = c
    per_sample_rows.append(row)

mat = pd.DataFrame(per_sample_rows).fillna(0)
mat.to_csv(f"{OUT}/gene_category_counts.csv", index=False)
print(f"[ok] gene_category_counts.csv ({len(mat)} rows, {mat.shape[1]} cols)")

# Per-category totals (presence normalised: gene count / CDS total * 1000)
subcols = [c for c in mat.columns if "::" in c]
norm = mat[subcols].div(mat["CDS_total"], axis=0) * 1000.0
norm["Sample"] = mat["Sample"]; norm["Genus"] = mat["Genus"]; norm["Hero"] = mat["Hero"]
norm.to_csv(f"{OUT}/gene_category_per1k_CDS.csv", index=False)

# -----------------------------------------------------------------------------
# T1c  HGT signals around the ure-cah cluster  (hero samples)
# -----------------------------------------------------------------------------
import math
def ure_cluster_window(df, pad=15000):
    prod = df["product"].str.lower().fillna("")
    gene = df["gene"].str.lower().fillna("")
    ure  = df[prod.str.contains("urease") | gene.str.startswith("ure") |
              prod.str.contains("urease subunit")]
    if ure.empty: return None
    # densest contig
    best=None
    for c,sub in ure.groupby("contig"):
        sub=sub.sort_values("start")
        span=sub["end"].max()-sub["start"].min()
        if len(sub)>=3 and span<40000:
            if best is None or len(sub)>best[3]:
                best=(c, max(0,sub["start"].min()-pad), sub["end"].max()+pad, len(sub))
    return best

def gc_of_region(fna_path, contig, s, e):
    # parse gz fasta
    opener = gzip.open if fna_path.endswith(".gz") else open
    seq=""; keep=False
    with opener(fna_path,"rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep: break
                keep = (line[1:].split()[0]==contig)
            elif keep: seq += line.strip()
    if not seq: return np.nan, np.nan
    region = seq[max(0,s-1):e]
    def gc(x):
        x=x.upper(); n=len(x);
        return (x.count("G")+x.count("C"))/n*100 if n else np.nan
    return gc(seq), gc(region)

hgt_rows=[]
for samp in sorted(HERO):
    df = all_gffs.get(samp)
    if df is None: continue
    win = ure_cluster_window(df)
    if not win: continue
    contig, ws, we, _ = win
    # mobile elements within window
    region = df[(df.contig==contig)&(df.start>=ws)&(df.end<=we)]
    prod_r = region["product"].str.lower().fillna("")
    mob_kw = ["transposase","integrase","recombinase xerc","recombinase xerd",
              "IS element","insertion sequence","phage","prophage","relaxase","mobilization"]
    mob_mask = np.zeros(len(region),dtype=bool)
    for kw in mob_kw: mob_mask |= prod_r.str.contains(re.escape(kw.lower()))
    n_mob = int(mob_mask.sum())
    # GC content
    fna = f"{BAKTA}/{samp}/{samp}.fna"
    genome_gc, region_gc = gc_of_region(fna, contig, ws, we)
    hgt_rows.append({"Sample":samp,"Contig":contig,"Start":ws,"End":we,
                     "Window_kb":(we-ws)/1000,"MobileElements":n_mob,
                     "GenomeGC":round(genome_gc,2),"RegionGC":round(region_gc,2),
                     "DeltaGC":round(region_gc-genome_gc,2)})
hgt_df = pd.DataFrame(hgt_rows)
hgt_df.to_csv(f"{OUT}/HGT_ureCah_cluster.csv", index=False)
print("[ok] HGT_ureCah_cluster.csv")
print(hgt_df.to_string(index=False))

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
def cat_heatmap(cat, fname, title):
    cols = [c for c in mat.columns if c.startswith(f"{cat}::")]
    sub  = mat[["Sample","Genus","Hero"]+cols].copy()
    sub  = sub.sort_values(["Hero","Genus","Sample"], ascending=[False,True,True])
    M    = sub[cols].values
    # normalise to per-1k CDS
    Mn   = M / mat.set_index("Sample").loc[sub["Sample"]]["CDS_total"].values.reshape(-1,1) * 1000.0
    fig_h = max(6, len(sub)*0.14)
    fig, ax = plt.subplots(figsize=(max(6, len(cols)*0.7), fig_h))
    im = ax.imshow(Mn, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.split("::")[1] for c in cols], rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(sub)))
    labels=[]
    for _,r in sub.iterrows():
        lab = f"{r['Sample']} [{r['Genus']}]"
        labels.append(lab)
    ax.set_yticklabels(labels, fontsize=5.5)
    for i,(_,r) in enumerate(sub.iterrows()):
        if r["Hero"]:
            ax.get_yticklabels()[i].set_color("#c0392b")
            ax.get_yticklabels()[i].set_fontweight("bold")
    cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("hits per 1k CDS", fontsize=8)
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{fname}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/{fname}.pdf", bbox_inches="tight")
    plt.close(fig)

cat_heatmap("Biofilm_EPS",    "Fig_T1a_Biofilm_EPS_heatmap",    "T1a  Biofilm / EPS gene modules (per 1k CDS)")
cat_heatmap("Ammonia_N",      "Fig_T1b_Ammonia_detox_heatmap",  "T1b  Ammonia detox & N-assimilation modules")
cat_heatmap("Alkaline_Osmo",  "Fig_T2a_Alkaline_stress_heatmap","T2a  Alkaline / osmotic stress tolerance")
cat_heatmap("CAZyme_proxy",   "Fig_T2b_CAZyme_heatmap",         "T2b  CAZyme keyword profile (proxy)")
cat_heatmap("MetalResist_AMR","Fig_T2c_MetalAMR_heatmap",       "T2c  Heavy-metal & AMR gene profile")

# Hero vs rest summary barplot (mean hits per 1k CDS per subcategory)
def mean_by_group(subcols):
    hero = norm.loc[norm.Hero, subcols].mean()
    rest = norm.loc[~norm.Hero, subcols].mean()
    return hero, rest

summary_rows=[]
for cat, subs in CATEGORIES.items():
    for sub in subs:
        key=f"{cat}::{sub}"
        if key in norm.columns:
            h = norm.loc[norm.Hero,key].mean()
            r = norm.loc[~norm.Hero,key].mean()
            summary_rows.append({"Category":cat,"Subcategory":sub,
                                 "Hero_mean_per1kCDS":round(h,3),
                                 "Rest_mean_per1kCDS":round(r,3),
                                 "Hero_over_Rest": round(h/r,2) if r>0 else np.nan})
summary=pd.DataFrame(summary_rows).sort_values(["Category","Hero_over_Rest"],ascending=[True,False])
summary.to_csv(f"{OUT}/Hero_vs_Rest_summary.csv", index=False)
print("[ok] Hero_vs_Rest_summary.csv")

# Global comparison figure: hero vs rest across categories
fig, ax = plt.subplots(figsize=(12,8))
sub = summary.dropna(subset=["Hero_over_Rest"]).copy()
sub["label"] = sub["Category"]+": "+sub["Subcategory"]
sub = sub.sort_values("Hero_over_Rest", ascending=True).tail(30)
ax.barh(sub["label"], sub["Hero_over_Rest"], color="#c0392b", alpha=0.8, edgecolor="black", lw=0.4)
ax.axvline(1.0, ls="--", color="gray", lw=0.8)
ax.set_xlabel("Hero / Rest  (mean hits per 1k CDS)")
ax.set_title("Top 30 gene modules enriched in Sphingobacterium hero clade", fontsize=10)
ax.tick_params(axis="y", labelsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig_HeroEnrichment_top30.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_HeroEnrichment_top30.pdf", bbox_inches="tight")
plt.close(fig)
print("[ok] Fig_HeroEnrichment_top30")
