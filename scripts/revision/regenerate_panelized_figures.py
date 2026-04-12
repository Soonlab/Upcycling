#!/usr/bin/env python3
"""
Regenerate all main figures as journal-style multi-panel composites with
consistent terminology (MICP-complete lineages; no "hero" wording).

Outputs to /data/data/Upcycling/research/figures_final/
  Figure_1.(png|pdf)   phylogeny + MICP heatmap               (single panel)
  Figure_2.(png|pdf)   MICP completeness by genus + prevalence (panels a, b)
  Figure_3.(png|pdf)   ure-cah synteny across the six MAGs     (single panel)
  Figure_4.(png|pdf)   DRAM metabolism heatmap + MICP bar      (panels a, b)
  Figure_5.(png|pdf)   Novelty + AAI + external ANI            (panels a, b, c)
  Figure_6.(png|pdf)   Permutation forest                      (single panel)
  Figure_7.(png|pdf)   Pan-genome + trait-module PCoA          (panels a, b)
"""
import os, re, gzip
from collections import defaultdict
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, FancyArrow
from matplotlib.gridspec import GridSpec
from Bio import Phylo
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova
from skbio.stats.ordination import pcoa

BASE="/data/data/Upcycling"
OUT =f"{BASE}/research/figures_final"
os.makedirs(OUT, exist_ok=True)

# ---------------- common data ----------------
gtdb = pd.read_csv(f"{BASE}/pangenome_work/gtdbtk_results/gtdbtk.bac120.summary.tsv", sep="\t")
gtdb["Genus"]   = gtdb["classification"].apply(lambda s: (re.search(r"g__([^;]*)",str(s)) or [None,""])[1] or "Unclassified")
gtdb["Species"] = gtdb["classification"].apply(lambda s: (re.search(r"s__([^;]*)",str(s)) or [None,""])[1])
tax = gtdb.set_index("user_genome")[["Genus","Species"]]

LINEAGE = {"C22","M1","S13","S16","S23","S26"}            # MICP-complete set
LABEL_GROUP = "MICP-complete lineages"
LABEL_OTHER = "Other MAGs"
COLOR_HIGHLIGHT = "#c0392b"
COLOR_OTHER     = "#7f8c8d"

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans"
})

# =============================================================================
# Figure 1 — phylogeny + MICP heatmap (rebuild, relabel)
# =============================================================================
def figure1():
    tree_path = f"{BASE}/pangenome_work/gtdbtk_results/align/gtdbtk.bac120.renamed.treefile"
    tree = Phylo.read(tree_path, "newick")
    try: tree.root_at_midpoint()
    except Exception: pass
    leaf_full=[t.name for t in tree.get_terminals()]
    sid=lambda f: f.split("_s__")[0]
    ordered_full=[]
    def collect(cl):
        if cl.is_terminal(): ordered_full.append(cl.name)
        else:
            for c in cl.clades: collect(c)
    collect(tree.root)
    micp = pd.read_csv(f"{BASE}/pangenome_work/MICP_Pangenome_Final_Summary.csv", index_col=0)
    gene_cols=["ureA","ureB","ureC","ureD","ureE","ureF","ureG","cah"]

    genera = sorted(tax["Genus"].unique())
    cmap = plt.get_cmap("tab20", len(genera))
    gcol = {g: cmap(i) for i,g in enumerate(genera)}

    n = len(ordered_full)
    fig = plt.figure(figsize=(12, max(10, n*0.14)))
    gs  = fig.add_gridspec(1, 4, width_ratios=[2.2, 0.1, 0.13, 1.4], wspace=0.02)
    ax_t = fig.add_subplot(gs[0,0]); ax_g = fig.add_subplot(gs[0,1], sharey=ax_t)
    ax_l = fig.add_subplot(gs[0,2], sharey=ax_t); ax_h = fig.add_subplot(gs[0,3], sharey=ax_t)

    y_of = {name:i for i,name in enumerate(ordered_full[::-1])}
    depths = tree.depths(); xmax = max(depths.values()) or 1.0
    def _draw(cl,x0):
        if cl.is_terminal():
            y=y_of[cl.name]; ax_t.plot([x0,depths[cl]],[y,y],color="black",lw=0.5); return y
        ys=[_draw(c,depths[cl]) for c in cl.clades]; ymid=(min(ys)+max(ys))/2
        ax_t.plot([depths[cl],depths[cl]],[min(ys),max(ys)],color="black",lw=0.5)
        if cl is not tree.root: ax_t.plot([x0,depths[cl]],[ymid,ymid],color="black",lw=0.5)
        return ymid
    _draw(tree.root,0)
    ax_t.set_xlim(0,xmax*1.02); ax_t.set_ylim(-0.5,n-0.5); ax_t.invert_yaxis()
    ax_t.set_yticks([]); ax_t.spines["left"].set_visible(False)
    ax_t.set_xlabel("Substitutions / site")

    for full in ordered_full:
        y=y_of[full]; sid_=sid(full)
        g=tax.loc[sid_,"Genus"] if sid_ in tax.index else "Unclassified"
        ax_g.add_patch(Rectangle((0,y-0.5),1,1,color=gcol.get(g,"#cccccc")))
    ax_g.set_xlim(0,1); ax_g.set_xticks([]); ax_g.set_yticks([])
    for s in ax_g.spines.values(): s.set_visible(False)

    ax_l.set_xlim(0,1); ax_l.set_xticks([]); ax_l.set_yticks([])
    for s in ax_l.spines.values(): s.set_visible(False)
    for full in ordered_full:
        sid_=sid(full); y=y_of[full]
        sp=tax.loc[sid_,"Species"] if sid_ in tax.index else ""
        color = COLOR_HIGHLIGHT if sid_ in LINEAGE else "black"
        weight= "bold" if sid_ in LINEAGE else "normal"
        ax_l.text(0,y,f"{sid_}  {sp}" if sp else sid_, va="center", ha="left",
                  fontsize=6.2, color=color, fontweight=weight)

    mat=np.zeros((n,len(gene_cols)))
    for full in ordered_full:
        sid_=sid(full)
        if sid_ in micp.index:
            mat[y_of[full],:] = micp.loc[sid_, gene_cols].values
    ax_h.imshow(mat, aspect="auto", cmap="Greens", vmin=0, vmax=1,
                extent=[-0.5, len(gene_cols)-0.5, n-0.5, -0.5])
    ax_h.set_xticks(range(len(gene_cols)))
    ax_h.set_xticklabels(gene_cols, rotation=45, ha="right", style="italic")
    ax_h.set_yticks([])
    for s in ax_h.spines.values(): s.set_visible(False)

    counts = tax["Genus"].value_counts()
    top = counts.head(10).index.tolist()
    handles=[Patch(facecolor=gcol[g], label=f"{g} (n={counts[g]})") for g in top]
    handles.append(Line2D_patch(COLOR_HIGHLIGHT, "MICP-complete lineage member"))
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5,-0.01))
    fig.suptitle("Figure 1. Phylogenomic distribution of MICP-related genes across 111 livestock-waste MAGs.",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=[0,0.02,1,0.98])
    fig.savefig(f"{OUT}/Figure_1.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_1.pdf", bbox_inches="tight")
    plt.close(fig); print("Figure 1 done")

from matplotlib.patches import Patch
def Line2D_patch(col, lab): return Patch(facecolor="white", edgecolor=col, lw=2, label=lab)

# =============================================================================
# Figure 2 — module score by genus + per-gene prevalence  (a, b)
# =============================================================================
def figure2():
    micp = pd.read_csv(f"{BASE}/pangenome_work/MICP_Pangenome_Final_Summary.csv", index_col=0)
    gene_cols=["ureA","ureB","ureC","ureD","ureE","ureF","ureG","cah"]
    df=micp.copy()
    df["Genus"]=[tax.loc[i,"Genus"] if i in tax.index else "Unclassified" for i in df.index]
    df["Score"]=df[gene_cols].sum(axis=1)
    order=df.groupby("Genus")["Score"].mean().sort_values(ascending=False).index.tolist()

    genera_all = sorted(tax["Genus"].unique())
    cmap = plt.get_cmap("tab20", len(genera_all))
    gcol = {g: cmap(i) for i,g in enumerate(genera_all)}

    fig = plt.figure(figsize=(13, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.28)

    # panel a — boxplot per genus
    ax=fig.add_subplot(gs[0,0])
    data=[df.loc[df.Genus==g,"Score"].values for g in order]
    bp=ax.boxplot(data, positions=range(len(order)), widths=0.6, patch_artist=True,
                  medianprops=dict(color="black"))
    for patch,g in zip(bp["boxes"],order):
        patch.set_facecolor(gcol.get(g,"#bbbbbb")); patch.set_alpha(0.85)
    for i,g in enumerate(order):
        y=df.loc[df.Genus==g,"Score"].values
        x=np.random.normal(i,0.07,size=len(y))
        mask=np.array([idx in LINEAGE for idx in df.loc[df.Genus==g].index])
        ax.scatter(x[~mask], y[~mask], s=14, color="black", alpha=0.5, zorder=3)
        if mask.any():
            ax.scatter(x[mask], y[mask], s=70, color=COLOR_HIGHLIGHT,
                       edgecolor="black", lw=0.6, zorder=4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=40, ha="right", style="italic")
    ax.set_ylabel("MICP module score (ureA–G + cah, max 8)")
    ax.axhline(8, ls="--", lw=0.6, color="gray"); ax.set_ylim(-0.3, 8.8)
    ax.grid(axis="y", ls=":", alpha=0.4)
    ax.text(-0.1, 1.03, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")

    # panel b — prevalence
    ax2=fig.add_subplot(gs[0,1])
    hero_prev = df.loc[df.index.isin(LINEAGE), gene_cols].mean()*100
    rest_prev = df.loc[~df.index.isin(LINEAGE), gene_cols].mean()*100
    xk=np.arange(len(gene_cols)); w=0.38
    ax2.bar(xk-w/2, hero_prev.values, w,
            label=f"{LABEL_GROUP} (n={len(LINEAGE)})",
            color=COLOR_HIGHLIGHT, edgecolor="black", lw=0.5)
    ax2.bar(xk+w/2, rest_prev.values, w,
            label=f"{LABEL_OTHER} (n={len(df)-len(LINEAGE)})",
            color=COLOR_OTHER, edgecolor="black", lw=0.5)
    ax2.set_xticks(xk); ax2.set_xticklabels(gene_cols, style="italic")
    ax2.set_ylabel("Prevalence (% of genomes)")
    ax2.set_ylim(0,105)
    ax2.legend(frameon=False, loc="lower right")
    ax2.grid(axis="y", ls=":", alpha=0.4)
    ax2.text(-0.12, 1.03, "b", transform=ax2.transAxes, fontsize=13, fontweight="bold")

    fig.suptitle("Figure 2. MICP module completeness across genera (a) and per-gene prevalence in MICP-complete lineages vs others (b).",
                 fontsize=10.5, y=1.01)
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_2.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_2.pdf", bbox_inches="tight")
    plt.close(fig); print("Figure 2 done")

# =============================================================================
# Figure 3 — ure-cah synteny across 6 MAGs
# =============================================================================
def figure3():
    BAKTA=f"{BASE}/MAGs_FASTA_files/bakta_results"
    GENE_COLORS={"ureA":"#1f77b4","ureB":"#2ca02c","ureC":"#d62728","ureD":"#9467bd",
                 "ureE":"#8c564b","ureF":"#e377c2","ureG":"#17becf","cah":"#ff7f0e","other":"#bbbbbb"}
    def classify(prod,gene):
        p=(prod or "").lower(); g=(gene or "").lower()
        if "urease subunit alpha" in p: return "ureC"
        if "urease subunit beta"  in p: return "ureB"
        if "urease subunit gamma" in p: return "ureA"
        if g in ["ured","uree","uref","ureg"]: return "ure"+g[3].upper()
        if "ured" in p: return "ureD"
        if "uree" in p: return "ureE"
        if "uref" in p: return "ureF"
        if "ureg" in p: return "ureG"
        if "carbonic anhyd" in p: return "cah"
        return "other"
    def parse_gff(path):
        rows=[]
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                p=line.rstrip().split("\t")
                if len(p)<9 or p[2]!="CDS": continue
                attrs=dict(kv.split("=",1) for kv in p[8].split(";") if "=" in kv)
                rows.append({"contig":p[0],"start":int(p[3]),"end":int(p[4]),"strand":p[6],
                             "gene":attrs.get("gene",""),"product":attrs.get("product","")})
        return pd.DataFrame(rows)
    def find_cluster(df):
        df=df.copy()
        df["cls"]=[classify(p,g) for p,g in zip(df["product"],df["gene"])]
        ure=df[df.cls.str.startswith("ure")]
        best=None
        for c,sub in ure.groupby("contig"):
            sub=sub.sort_values("start")
            if len(sub)<3: continue
            span=sub["end"].max()-sub["start"].min()
            if span>40000: continue
            sc=len(set(sub.cls))
            if best is None or sc>best[0]: best=(sc,c,sub["start"].min(),sub["end"].max())
        if not best: return None
        _,c,s,e=best
        ws=max(0,s-2000); we=e+2000
        region=df[(df.contig==c)&(df.start>=ws)&(df.end<=we)].sort_values("start")
        return c,ws,we,region

    HERO_LIST=["M1","S13","S16","C22","S23","S26"]  # ordered by cluster completeness
    LINEAGE_LBL={"M1":"Pseudomonas_E","S13":"Sphingobacterium","S16":"Sphingobacterium",
                 "C22":"Sphingobacterium","S23":"Sphingobacterium","S26":"Pseudomonas_E"}

    fig, axes = plt.subplots(len(HERO_LIST), 1, figsize=(13, 1.25*len(HERO_LIST)+0.8))
    if len(HERO_LIST)==1: axes=[axes]
    for ax,samp in zip(axes, HERO_LIST):
        gff=f"{BAKTA}/{samp}/{samp}.gff3"
        if not os.path.exists(gff):
            ax.axis("off"); continue
        df=parse_gff(gff); res=find_cluster(df)
        if res is None: ax.axis("off"); continue
        ctg,ws,we,region=res
        off=region["start"].min(); L=region["end"].max()-off
        for _,r in region.iterrows():
            cls=classify(r["product"],r["gene"])
            col=GENE_COLORS.get(cls, GENE_COLORS["other"])
            x0=(r["start"]-off)/1000; x1=(r["end"]-off)/1000
            if r["strand"]=="-": xs,dx = x1, -(x1-x0)
            else: xs,dx = x0, (x1-x0)
            head=min(0.35, abs(dx)*0.4)
            ax.add_patch(FancyArrow(xs, 0, dx, 0, width=0.55, head_width=0.85,
                                    head_length=head, length_includes_head=True,
                                    facecolor=col, edgecolor="black", lw=0.4))
            if cls!="other":
                ax.text((x0+x1)/2, 0.95, cls, ha="center", va="bottom",
                        fontsize=7, style="italic")
        ax.set_xlim(-0.5, L/1000+0.5); ax.set_ylim(-1.2,1.7)
        ax.set_yticks([])
        for s in ["top","right","left"]: ax.spines[s].set_visible(False)
        ax.tick_params(axis="x", labelsize=7)
        ax.set_xlabel("kb" if ax is axes[-1] else "", fontsize=8)
        ax.text(-0.02, 0.5, f"{samp}\n({LINEAGE_LBL[samp]})",
                transform=ax.transAxes, ha="right", va="center",
                fontsize=8.5, fontweight="bold", color=COLOR_HIGHLIGHT)

    handles=[Patch(facecolor=GENE_COLORS[g], edgecolor="black", label=g)
             for g in ["ureA","ureB","ureC","ureD","ureE","ureF","ureG","cah","other"]]
    fig.legend(handles=handles, loc="lower center", ncol=9, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5,-0.02))
    fig.suptitle("Figure 3. Conserved ureABCDEFG operon architecture in MICP-complete MAGs.",
                 fontsize=10.5, y=0.995)
    fig.tight_layout(rect=[0.03,0.03,1,0.97])
    fig.savefig(f"{OUT}/Figure_3.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_3.pdf", bbox_inches="tight")
    plt.close(fig); print("Figure 3 done")

# =============================================================================
# Figure 4 — DRAM metabolism heatmap + MICP-critical bar   (a, b)
# =============================================================================
def figure4():
    prod = pd.read_csv("/data/pangenome_work/dram_output/distillate/product.tsv", sep="\t", index_col=0)
    keep=[c for c in prod.columns if prod[c].std()>0 and (prod[c]>0).sum()>=3]
    curated=[c for c in keep if re.search(
        r"urease|carbonic anhyd|Nitrogen|Na.*H.*antiport|glutamate|glutamine|"
        r"Cobalamin|B12|Urea cycle|Arginine|Carbohydrate|CAZy|biofilm|"
        r"Flagell|Methionine|Serine|Fatty acid|Pyruvate|Glycolysis|TCA", c, re.I)]
    curated=list(dict.fromkeys(curated))[:30]
    sel=prod[curated].copy()
    sel["Genus"]=[tax.loc[i,"Genus"] if i in tax.index else "Unclassified" for i in sel.index]
    sel["Group"]=[i in LINEAGE for i in sel.index]
    # drop MAGs whose DRAM annotation is empty (all-zero rows)
    row_sum = sel[curated].sum(axis=1)
    sel = sel.loc[row_sum > 0]
    sel=sel.sort_values(["Group","Genus"], ascending=[False, True])

    fig=plt.figure(figsize=(17, 10))
    gs=fig.add_gridspec(1, 2, width_ratios=[2.6, 1], wspace=0.25)

    # panel a — heatmap
    ax=fig.add_subplot(gs[0,0])
    mat=sel[curated].values.astype(float)
    im=ax.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks(range(len(curated)))
    ax.set_xticklabels(curated, rotation=60, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(sel)))
    lbls=[f"{i}  [{g}]" for i,g in zip(sel.index, sel["Genus"])]
    ax.set_yticklabels(lbls, fontsize=5.0)
    for i,(_,r) in enumerate(sel.iterrows()):
        if r.Group:
            ax.get_yticklabels()[i].set_color(COLOR_HIGHLIGHT)
            ax.get_yticklabels()[i].set_fontweight("bold")
    cbar=fig.colorbar(im, ax=ax, shrink=0.5, pad=0.01)
    cbar.set_label("Module completeness", fontsize=8)
    ax.text(-0.18, 1.02, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.set_title(f"DRAM metabolic-module completeness ({len(sel)} MAGs × 30 modules)", fontsize=9, loc="left")

    # panel b — Bakta-based trait-module comparison for MICP-critical pathways
    cts=pd.read_csv(f"{BASE}/research/extra/gene_category_counts.csv")
    subcols=[c for c in cts.columns if "::" in c]
    cts_norm=cts.set_index("Sample")[subcols].div(cts.set_index("Sample")["CDS_total"], axis=0)*1000
    cts_norm["Group"]=[i in LINEAGE for i in cts_norm.index]
    # curated list mapping to MICP-critical phenotype
    bakta_keys = [
        ("Alkaline_Osmo::Mrp_complex",   "Mrp Na+/H+ antiporter"),
        ("Alkaline_Osmo::Na_H_antiporter","nhaA–C antiporter"),
        ("Alkaline_Osmo::oxidative",     "Oxidative-stress defence"),
        ("Alkaline_Osmo::compatible_solute","Compatible solute synth."),
        ("CAZyme_proxy::glycoside_hydrolase","Glycoside hydrolases"),
        ("CAZyme_proxy::carb_binding",   "Carbohydrate-binding modules"),
        ("Ammonia_N::urea_transport",    "Urea transport (urtABC)"),
        ("Ammonia_N::GS_GOGAT",          "GS-GOGAT (glnA, gltBD)"),
        ("Biofilm_EPS::quorum",          "Quorum sensing"),
        ("Biofilm_EPS::cellulose",       "Cellulose synthesis (bcs)"),
    ]
    labels=[lbl for k,lbl in bakta_keys if k in cts_norm.columns]
    keys  =[k   for k,lbl in bakta_keys if k in cts_norm.columns]
    hero_mean=cts_norm.loc[cts_norm.Group, keys].mean()
    rest_mean=cts_norm.loc[~cts_norm.Group, keys].mean()
    axb=fig.add_subplot(gs[0,1])
    y=np.arange(len(keys)); w=0.4
    axb.barh(y-w/2, hero_mean.values, w, label=f"{LABEL_GROUP} (n={cts_norm.Group.sum()})",
             color=COLOR_HIGHLIGHT, edgecolor="black", lw=0.4)
    axb.barh(y+w/2, rest_mean.values, w, label=f"{LABEL_OTHER} (n={(~cts_norm.Group).sum()})",
             color=COLOR_OTHER, edgecolor="black", lw=0.4)
    axb.set_yticks(y); axb.set_yticklabels(labels, fontsize=8)
    axb.set_xlabel("Gene hits per 10³ CDS (Bakta keyword scan)")
    axb.set_xscale("log")
    axb.invert_yaxis()
    axb.legend(frameon=False, fontsize=8, loc="lower right")
    axb.grid(axis="x", ls=":", alpha=0.4, which="both")
    axb.text(-0.35, 1.02, "b", transform=axb.transAxes, fontsize=13, fontweight="bold")
    axb.set_title("Trait-module enrichment: MICP-complete vs other MAGs", fontsize=9, loc="left")

    fig.suptitle("Figure 4. DRAM metabolic-module reconstruction across the 111-MAG panel.",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=[0,0,1,0.98])
    fig.savefig(f"{OUT}/Figure_4.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_4.pdf", bbox_inches="tight")
    plt.close(fig); print("Figure 4 done")

# =============================================================================
# Figure 5 — novelty screen + AAI + external ANI  (a, b, c)
# =============================================================================
def figure5():
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[1,1], height_ratios=[1,1.2], hspace=0.35, wspace=0.25)

    # -- panel a: novelty overview --
    nov=pd.read_csv(f"{BASE}/research/extra/novelty_ANI_screen.csv")
    nov_valid=nov.dropna(subset=["ANI"]).sort_values("ANI").reset_index(drop=True)
    nov_valid["MICP_lineage"]=nov_valid["user_genome"].isin(LINEAGE)
    ax=fig.add_subplot(gs[0,0])
    colors=[COLOR_HIGHLIGHT if h else "#7f8c8d" for h in nov_valid.MICP_lineage]
    sizes =[70 if h else 16 for h in nov_valid.MICP_lineage]
    ax.scatter(range(len(nov_valid)), nov_valid["ANI"], c=colors, s=sizes,
               edgecolor="black", lw=0.3, zorder=3)
    ax.axhline(95, ls="--", color="red", lw=0.8, label="95 % species cutoff")
    for i,r in nov_valid.iterrows():
        if r.MICP_lineage:
            ax.annotate(r.user_genome, (i,r.ANI), xytext=(3,4),
                        textcoords="offset points", fontsize=8,
                        color=COLOR_HIGHLIGHT, fontweight="bold")
    ax.set_xticks([]); ax.set_xlabel("MAGs (ranked by ANI)")
    ax.set_ylabel("ANI to closest GTDB reference (%)")
    ax.set_title(f"Within-panel novelty screen — {nov.Novel_sp_candidate.sum()} candidates < 95 % ANI",
                 fontsize=9, loc="left")
    ax.legend(frameon=False, loc="lower right")
    ax.text(-0.1, 1.03, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.grid(axis="y", ls=":", alpha=0.3)

    # -- panel b: AAI S13/S16 vs in-panel Sphingobacterium --
    aai=pd.read_csv(f"{BASE}/research/extra/novel_species/AAI_S13_S16_vs_Sphingobacterium.csv")
    ax=fig.add_subplot(gs[0,1])
    for off,q in enumerate(["S13","S16"]):
        sub=aai[aai.Query==q].dropna(subset=["AAI"]).sort_values("AAI", ascending=True)
        y=np.arange(len(sub)) + off*7
        lbls=[f"{q} → {t} ({s or 'sp.'})" for t,s in zip(sub.Target, sub.Target_species)]
        colors=["#3498db" if qc==q else "#2ecc71" for qc in [q]*len(sub)]
        ax.barh(y, sub.AAI.values, color="#3498db" if q=="S13" else "#16a085",
                edgecolor="black", lw=0.4)
        for yi, lb in zip(y, lbls):
            ax.text(sub.AAI.min()-2, yi, lb, va="center", ha="right", fontsize=7)
    ax.axvline(95, ls="--", color="red", lw=0.8, label="95 % (species)")
    ax.axvline(70, ls=":",  color="gray", lw=0.8, label="70 % (genus)")
    ax.set_xlabel("AAI (%)"); ax.set_xlim(55, 100)
    ax.set_yticks([])
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.set_title("AAI of S13 and S16 vs in-panel Sphingobacterium (RBH, mmseqs2)",
                 fontsize=9, loc="left")
    ax.text(-0.05, 1.03, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")

    # -- panel c: External ANI vs 63 RefSeq --
    ani = pd.read_csv(f"{BASE}/research/revision/ANI_ext_sphingo_matrix.csv", index_col=0)
    ext_nov = pd.read_csv(f"{BASE}/research/revision/ANI_ext_sphingo_novelty.csv")
    our_rows=[n for n in ani.index if n.startswith("OUR_")]
    ref_rows=[n for n in ani.columns if n.startswith("REF_")]
    ax=fig.add_subplot(gs[1,:])
    order_mags=sorted(our_rows, key=lambda x: x.replace("OUR_",""))
    width=0.14
    xs=np.arange(len(our_rows))
    # For each MAG draw max + percentile distribution
    for i,o in enumerate(order_mags):
        mag=o.replace("OUR_","")
        vals=ani.loc[o, ref_rows].dropna().sort_values(ascending=False)
        pos=i
        # box + scatter of all references
        ax.scatter([pos]*len(vals), vals.values, color="#95a5a6", s=6, alpha=0.5, zorder=2)
        # max value
        maxv=vals.max()
        col = COLOR_HIGHLIGHT if maxv<95 else "#3498db"
        ax.scatter(pos, maxv, s=110, color=col, edgecolor="black", lw=0.6, zorder=4)
        nov_row=ext_nov[ext_nov.MAG==mag].iloc[0] if (ext_nov.MAG==mag).any() else None
        if nov_row is not None:
            lab = f"{maxv:.2f}%\n{nov_row['Nearest_organism'].replace('Sphingobacterium ','S. ')}"
            ax.annotate(lab, (pos, maxv), xytext=(8, -8),
                        textcoords="offset points", fontsize=7,
                        color=col, fontweight="bold" if maxv<95 else "normal")
    ax.axhline(95, ls="--", color="red", lw=0.8, label="95 % species cutoff")
    ax.set_xticks(range(len(order_mags)))
    ax.set_xticklabels([m.replace("OUR_","") for m in order_mags], fontweight="bold")
    ax.set_ylabel("ANI vs RefSeq Sphingobacterium (%)")
    ax.set_ylim(72, 101)
    ax.legend(frameon=False, loc="lower left")
    # Custom legend markers
    leg_handles=[
        Patch(facecolor=COLOR_HIGHLIGHT, edgecolor="black", label="< 95 %  novel species candidate"),
        Patch(facecolor="#3498db",      edgecolor="black", label="≥ 95 %  assigned species"),
        Patch(facecolor="#95a5a6",      edgecolor="black", label="per-reference ANI (all 63 RefSeq genomes)"),
    ]
    ax.legend(handles=leg_handles, loc="lower left", frameon=False, fontsize=8)
    ax.set_title("Pairwise skani ANI of 6 study Sphingobacterium MAGs against 63 RefSeq reference genomes",
                 fontsize=9, loc="left")
    ax.text(-0.04, 1.03, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.grid(axis="y", ls=":", alpha=0.3)

    fig.suptitle("Figure 5. Two study MAGs (S13, S16) are candidate novel Sphingobacterium species.",
                 fontsize=11, y=0.995)
    fig.savefig(f"{OUT}/Figure_5.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_5.pdf", bbox_inches="tight")
    plt.close(fig); print("Figure 5 done")

# =============================================================================
# Figure 6 — permutation forest
# =============================================================================
def figure6():
    res=pd.read_csv(f"{BASE}/research/revision/Hero_vs_Rest_permutation_stats.csv")
    plot_df=res[res["Fold_change"]>1].dropna(subset=["Fold_change"]).copy()
    plot_df["label"] = plot_df["Category"]+": "+plot_df["Subcategory"]
    plot_df=plot_df.sort_values("Fold_change", ascending=True).tail(25)
    ci=plot_df["Fold_change_CI95"].str.strip("[]").str.split(", ", expand=True).astype(float)
    plot_df["lo"]=ci[0]; plot_df["hi"]=ci[1]

    fig,ax=plt.subplots(figsize=(11, max(6, len(plot_df)*0.32)))
    y=np.arange(len(plot_df))
    colors=["#c0392b" if q<0.05 else "#e67e22" if q<0.10 else "#95a5a6"
            for q in plot_df["Permutation_q_BH"]]
    ax.errorbar(plot_df["Fold_change"], y,
                xerr=[plot_df["Fold_change"]-plot_df["lo"], plot_df["hi"]-plot_df["Fold_change"]],
                fmt="o", color="black", ecolor="gray", elinewidth=0.6, markersize=0, capsize=2)
    ax.scatter(plot_df["Fold_change"], y, s=70, c=colors, edgecolor="black", lw=0.4, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(plot_df["label"], fontsize=8)
    ax.axvline(1.0, ls="--", color="black", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("Fold change  (MICP-complete lineages / other MAGs;  log scale, 95 % bootstrap CI)")
    handles=[Patch(color="#c0392b", label="q < 0.05"),
             Patch(color="#e67e22", label="q < 0.10"),
             Patch(color="#95a5a6", label="n.s.")]
    ax.legend(handles=handles, frameon=False, loc="lower right", title="BH-FDR")
    fig.suptitle("Figure 6. Permutation-tested enrichment of trait modules in MICP-complete lineages vs others.",
                 fontsize=10.5, y=0.995)
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_6.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_6.pdf", bbox_inches="tight")
    plt.close(fig); print("Figure 6 done")

# =============================================================================
# Figure 7 — PCoA a (pan-genome) + b (trait-module)
# =============================================================================
def figure7():
    # pan-genome
    rtab=pd.read_csv(f"{BASE}/pangenome_work/panaroo_results/gene_presence_absence.Rtab",
                     sep="\t", index_col=0)
    prev=rtab.sum(axis=1)/rtab.shape[1]
    rtab_f=rtab.loc[(prev>=0.05)&(prev<=0.95)]
    mags=rtab_f.columns.tolist()
    dist=pdist(rtab_f.T.values.astype(int), metric="jaccard")
    dm=DistanceMatrix(squareform(dist), ids=mags)
    SRC=lambda m: {"C":"Cattle","M":"Swine","S":"Sheep","V":"Poultry"}.get(m[0],"Other")
    meta=pd.DataFrame({"Source":[SRC(m) for m in mags],
                       "Genus":[tax.loc[m,"Genus"] if m in tax.index else "Unclassified" for m in mags],
                       "MICP":[m in LINEAGE for m in mags]}, index=mags)
    perm=permanova(dm, grouping=meta.loc[mags,"Source"], permutations=999)
    perm_g=permanova(dm, grouping=meta.loc[mags,"Genus"], permutations=999)
    ord_=pcoa(dm, number_of_dimensions=2); var=ord_.proportion_explained*100
    coords=ord_.samples.copy(); coords.index=mags

    # trait
    cts=pd.read_csv(f"{BASE}/research/extra/gene_category_counts.csv").set_index("Sample")
    subcols=[c for c in cts.columns if "::" in c]
    norm=cts[subcols].div(cts["CDS_total"], axis=0)*1000
    mag_t=[m for m in mags if m in norm.index]
    dm_t=DistanceMatrix(squareform(pdist(norm.loc[mag_t].values, metric="euclidean")), ids=mag_t)
    perm_t=permanova(dm_t, grouping=meta.loc[mag_t,"Source"], permutations=999)
    ord_t=pcoa(dm_t, number_of_dimensions=2); var_t=ord_t.proportion_explained*100
    coords_t=ord_t.samples.copy(); coords_t.index=mag_t

    SCOL={"Cattle":"#e74c3c","Swine":"#3498db","Sheep":"#2ecc71","Poultry":"#f39c12"}
    fig,axes=plt.subplots(1,2, figsize=(14,6))

    # panel a
    ax=axes[0]
    for s,col in SCOL.items():
        sub=meta.loc[meta.Source==s]
        ax.scatter(coords.loc[sub.index,"PC1"], coords.loc[sub.index,"PC2"],
                   s=35, color=col, alpha=0.8, edgecolor="black", lw=0.4,
                   label=f"{s} (n={len(sub)})")
    lin_=meta[meta.MICP]
    ax.scatter(coords.loc[lin_.index,"PC1"], coords.loc[lin_.index,"PC2"],
               s=150, facecolor="none", edgecolor="red", lw=1.8, zorder=5)
    for m in lin_.index:
        ax.annotate(m, (coords.loc[m,"PC1"], coords.loc[m,"PC2"]),
                    xytext=(5,5), textcoords="offset points", fontsize=8,
                    color=COLOR_HIGHLIGHT, fontweight="bold")
    ax.set_xlabel(f"PC1 ({var.iloc[0]:.1f} %)"); ax.set_ylabel(f"PC2 ({var.iloc[1]:.1f} %)")
    ax.set_title("a  Pan-genome PCoA (Jaccard on Panaroo presence/absence)\n"
                 f"Source: pseudo-F={perm['test statistic']:.2f}, p={perm['p-value']:.3f}  |  "
                 f"Genus: pseudo-F={perm_g['test statistic']:.2f}, p={perm_g['p-value']:.3f}",
                 fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(ls=":", alpha=0.3)

    # panel b
    ax=axes[1]
    for s,col in SCOL.items():
        sub=meta.loc[mag_t].loc[meta.loc[mag_t,"Source"]==s]
        ax.scatter(coords_t.loc[sub.index,"PC1"], coords_t.loc[sub.index,"PC2"],
                   s=35, color=col, alpha=0.8, edgecolor="black", lw=0.4,
                   label=f"{s} (n={len(sub)})")
    lin_=meta.loc[mag_t][meta.loc[mag_t,"MICP"]]
    ax.scatter(coords_t.loc[lin_.index,"PC1"], coords_t.loc[lin_.index,"PC2"],
               s=150, facecolor="none", edgecolor="red", lw=1.8, zorder=5)
    for m in lin_.index:
        ax.annotate(m, (coords_t.loc[m,"PC1"], coords_t.loc[m,"PC2"]),
                    xytext=(5,5), textcoords="offset points", fontsize=8,
                    color=COLOR_HIGHLIGHT, fontweight="bold")
    ax.set_xlabel(f"PC1 ({var_t.iloc[0]:.1f} %)"); ax.set_ylabel(f"PC2 ({var_t.iloc[1]:.1f} %)")
    ax.set_title("b  Trait-module PCoA (Euclidean on 38 per-1k-CDS module counts)\n"
                 f"Source: pseudo-F={perm_t['test statistic']:.2f}, p={perm_t['p-value']:.3f}",
                 fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(ls=":", alpha=0.3)

    fig.suptitle("Figure 7. Pan-genome is structured by taxonomy; trait-module distribution is structured by waste source.",
                 fontsize=10.5, y=1.01)
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_7.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/Figure_7.pdf", bbox_inches="tight")
    plt.close(fig); print("Figure 7 done")

# =============================================================================
if __name__=="__main__":
    figure1(); figure2(); figure3(); figure4(); figure5(); figure6(); figure7()
    print(f"\nAll figures written to {OUT}")
