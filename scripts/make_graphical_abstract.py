#!/usr/bin/env python3
"""Single-panel graphical abstract summarising the study workflow + findings."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
from matplotlib.lines import Line2D
import matplotlib.image as mpimg
import numpy as np

fig, ax = plt.subplots(figsize=(14, 6.5))
ax.set_xlim(0,14); ax.set_ylim(0,7); ax.axis("off")

# -- Title --
ax.text(7, 6.6, "Comparative genomics of 111 livestock-waste MAGs nominates novel "
                "$\\it{Sphingobacterium}$ species with a vertically inherited\n"
                "$\\it{ure}$–$\\it{cah}$ module as chassis for alkali-tolerant MICP",
        ha="center", va="center", fontsize=11.5, fontweight="bold")

# -- Panel 1: Input (livestock waste) --
box1 = FancyBboxPatch((0.2, 2.2), 2.6, 3.2, boxstyle="round,pad=0.08",
                      facecolor="#fdf2e9", edgecolor="#e67e22", lw=1.5)
ax.add_patch(box1)
ax.text(1.5, 5.0, "Livestock waste\nmetagenomes", ha="center", fontsize=10, fontweight="bold", color="#d35400")
# 4 source icons
sources = [("🐄 Cattle","#e74c3c"), ("🐖 Swine","#3498db"), ("🐑 Sheep","#2ecc71"), ("🐔 Poultry","#f39c12")]
for i,(s,c) in enumerate(sources):
    ax.text(1.5, 4.4 - i*0.45, s, ha="center", fontsize=9, color=c)
ax.text(1.5, 2.5, "111 MAGs", ha="center", fontsize=11, fontweight="bold", color="#2c3e50")

# -- Arrow 1 --
ax.add_patch(FancyArrowPatch((2.9, 3.8), (3.7, 3.8), arrowstyle="->", mutation_scale=20, lw=2, color="#34495e"))

# -- Panel 2: Integrated pipeline --
box2 = FancyBboxPatch((3.8, 2.2), 2.9, 3.2, boxstyle="round,pad=0.08",
                      facecolor="#eaf2f8", edgecolor="#2980b9", lw=1.5)
ax.add_patch(box2)
ax.text(5.25, 5.05, "Integrated in-silico\npipeline", ha="center", fontsize=10, fontweight="bold", color="#1f618d")
tools=["Bakta · Panaroo","GTDB-Tk · IQ-TREE","DRAM · dbCAN-HMMER","skani · mmseqs2 AAI","Permutation + BH-FDR"]
for i,t in enumerate(tools):
    ax.text(5.25, 4.4 - i*0.36, t, ha="center", fontsize=8.5, color="#1a5276")

# -- Arrow 2 --
ax.add_patch(FancyArrowPatch((6.8, 3.8), (7.6, 3.8), arrowstyle="->", mutation_scale=20, lw=2, color="#34495e"))

# -- Panel 3: Key findings --
box3 = FancyBboxPatch((7.7, 2.2), 4.0, 3.2, boxstyle="round,pad=0.08",
                      facecolor="#fdedec", edgecolor="#c0392b", lw=1.5)
ax.add_patch(box3)
ax.text(9.7, 5.05, "Key findings", ha="center", fontsize=10, fontweight="bold", color="#922b21")
findings=[
    "• 2 convergent lineages carry complete",
    "   $\\it{ureABCDEFG}$–$\\it{cah}$ module (n=6)",
    "• Vertical inheritance (no flanking MGE,",
    "   ∆GC ≈ 0 %)",
    "• 9 trait modules enriched (q<0.05)",
    "   Mrp 10.9× · CBM 9.8× · GH 4.7×",
    "• $\\bf{S13, S16: novel Sphingobacterium}$",
    "   $\\bf{species}$ (ANI < 95 % vs public)",
]
for i,t in enumerate(findings):
    ax.text(7.9, 4.55 - i*0.32, t, ha="left", fontsize=8.7,
            color="#922b21" if "novel" in t else "#566573",
            fontweight="bold" if "novel" in t else "normal")

# -- Arrow 3 (down) --
ax.add_patch(FancyArrowPatch((12.0, 3.8), (12.8, 3.8), arrowstyle="->", mutation_scale=20, lw=2, color="#34495e"))

# -- Panel 4: Application --
box4 = FancyBboxPatch((12.9, 2.2), 1.05, 3.2, boxstyle="round,pad=0.08",
                      facecolor="#e8f8f5", edgecolor="#16a085", lw=1.5)
ax.add_patch(box4)
ax.text(13.42, 5.05, "Waste-\ncoupled\nMICP\nbiocement", ha="center", va="center",
        fontsize=9, fontweight="bold", color="#0e6655")
# CaCO3 icon
ax.text(13.42, 3.3, "CaCO$_3$", ha="center", fontsize=10, fontweight="bold", color="#16a085")

# -- Bottom strip: equation --
ax.text(7, 1.5,
        "Urea + H$_2$O  →  2 NH$_3$ + CO$_2$  (urease) · CO$_2$ + H$_2$O ⇌ HCO$_3^-$ + H$^+$  (carbonic anhydrase) · Ca$^{2+}$ + HCO$_3^-$  →  CaCO$_3$ ↓",
        ha="center", fontsize=9, color="#2c3e50", style="italic")

# -- Bottom metadata line --
ax.text(7, 0.8, "Circular bioeconomy · Livestock-waste valorisation · Alkali-tolerant MICP",
        ha="center", fontsize=9, fontweight="bold", color="#7f8c8d")
ax.text(7, 0.35, "Code & data: github.com/Soonlab/Upcycling", ha="center", fontsize=8, color="#95a5a6")

fig.tight_layout()
out="/data/data/Upcycling/research/revision/Graphical_abstract"
fig.savefig(out+".png", dpi=300, bbox_inches="tight")
fig.savefig(out+".pdf", bbox_inches="tight")
plt.close(fig)
print("saved", out)
