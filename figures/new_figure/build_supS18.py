"""Suppl Fig S18 - ANI of the MICP-complete MAGs against the curated MICP reference panel (C5).

Panels
  A  the complete comparison: six MICP-complete MAGs against every verified reference of the
     curated panel.  A cell carries the ANI where skani reported an alignment and a dash where
     it reported none.  The sparsity is the result: the MAGs share no measurable whole-genome
     identity with any canonical ureolytic reference, only with references of their own genus.
  B  the reported pairs, ranked, against the 95 % species boundary.

REFERENCE PANEL - REBUILT 2026-09-04
The first C5 run is superseded.  Its `curated_refs.tsv` accessions did not correspond to the
organisms they were labelled with: of 20 accessions, only GCF_000009045.1 (Bacillus subtilis 168)
and GCF_001043025.1 (Pseudomonas helleri DSM 29165) named the intended taxon.  The rest pointed
at unrelated organisms (the accession labelled Sporosarcina ureae is Escherichia coli; the two
labelled Sphingobacterium are Marinobacter and Enterobacter) or did not exist, Sporosarcina
pasteurii among the latter.  The downloads had been faithful to those accessions, so the fault
was upstream of the download and the whole panel was fictitious.
The panel was rebuilt by resolving each intended TAXON against the NCBI datasets API,
downloading by accession, and verifying the organism NCBI reports for that accession against
the intended taxon before accepting the file
(`research/additional/C5_panMICP_env_v2/`, manifest `reference_manifest.tsv`).
Twenty of twenty references are now organism-verified.  Two taxa carry current names that
differ from the manuscript's: Halomonas pacifica is now Bisbaumannia pacifica and Bacillus
megaterium is now Priestia megaterium, both the same tax_id.  Sphingobacterium sp. 21 has no
genome under that strain designation and was replaced in the panel by the two Sphingobacterium
type strains the manuscript actually compares against.

Sources
  research/additional/C5_panMICP_env_v2/reference_manifest.tsv  accession -> verified organism
  research/additional/C5_panMICP_env_v2/skani_panMICP.matrix    all-vs-all ANI, verified panel
  Table_S20_skani_hero_vs_refs.tsv                              the 6 x 20 pair table
A matrix entry of 0.00 means skani reported no alignment for that pair, not 0 % identity, and
is drawn as absent rather than as a value.

Colour meanings on this page
  blue   = Sphingobacterium lineage MAG (S13, S16, S23, C22), and its cells and points
  orange = Pseudomonas_E lineage MAG (M1, S26), and its cells and points
  grey   = the 95 % species boundary, the no-alignment dashes and the panel rules; none of
           them a lineage
Cell and point colour therefore carries lineage only; ANI magnitude is printed, not encoded,
so that a single hue is never asked to mean two things on one page.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import GREY, TEXT, LIGHT, FS_STAT, FS_BODY

st.setup()
OUT = HERE / "figures"
SUPP = Path(gh.SUPP)
C5 = Path(gh.ADDITIONAL) / "C5_panMICP_env_v2"

SPECIES_ANI = 95.0  # standard prokaryotic species boundary (Konstantinidis & Tiedje 2005)

# ------------------------------------------------------------------ verified reference panel
man = pd.read_csv(C5 / "reference_manifest.tsv", sep="\t")
assert (man.status == "verified").all(), man.loc[man.status != "verified", "accession"].tolist()
n_refs = len(man)


def binomial(organism):
    """Genus + species of the organism NCBI reports, for use as a column label."""
    return " ".join(str(organism).split()[:2])


man["label"] = [binomial(o) for o in man.reported_organism]
# where two accessions share a binomial, the strain designation distinguishes them; it is
# read from the assembly report, not typed in, so the label always names the genome drawn
dupe = man.label.duplicated(keep=False)
man.loc[dupe, "label"] = [f"{lab} {str(strain).replace('type strain: ', '')}"
                          for lab, strain in zip(man.loc[dupe, "label"],
                                                 man.loc[dupe, "strain"])]
assert man.label.is_unique

# ------------------------------------------------------------------ data: skani matrix
lines = (C5 / "skani_panMICP.matrix").read_text().splitlines()
n_declared = int(lines[0])
names, rows = [], []
for ln in lines[1:]:
    parts = ln.split("\t")
    names.append(Path(parts[0]).stem)
    rows.append([float(v) for v in parts[1:]])
M = np.array(rows)
assert M.shape == (n_declared, len(names)) and M.shape[0] == M.shape[1], M.shape
assert np.allclose(M, M.T), "skani matrix is not symmetric"
idx = {nm: i for i, nm in enumerate(names)}
assert sorted(n.replace("HERO_", "") for n in names if n.startswith("HERO_")) == sorted(st.HEROES)

# ------------------------------------------------------------------ data: Table S20, asserted
tab = pd.read_csv(SUPP / "Table_S20_skani_hero_vs_refs.tsv", sep="\t")
# pandas reads the TRUE/FALSE column as bool; normalise so the check below is unambiguous
tab["Alignment_reported"] = tab.Alignment_reported.astype(str).str.upper() == "TRUE"
assert len(tab) == len(st.HEROES) * n_refs, (len(tab), n_refs)
assert set(tab.Ref_accession) == set(man.accession)

ani = pd.DataFrame(np.nan, index=st.HEROES, columns=list(man.label))
acc2label = dict(zip(man.accession, man.label))
for r in tab.itertuples():
    m = M[idx[f"HERO_{r.MAG}"], idx[r.Ref_accession]]
    reported = bool(r.Alignment_reported)
    assert reported == (m > 0), (r.MAG, r.Ref_accession, m, r.Alignment_reported)
    if reported:
        assert abs(m - float(r.ANI)) < 0.01, (r.MAG, r.Ref_accession, m, r.ANI)
        ani.loc[r.MAG, acc2label[r.Ref_accession]] = m

reported = ani.stack().sort_values(ascending=False)
assert len(reported) == int(tab.Alignment_reported.sum())
print(f"  verified references {n_refs} | reported MAG-reference pairs {len(reported)}")

# every reported pair is with a reference of the MAG's own genus - state it as a check,
# derived from the manifest rather than asserted from a list
for (mag, lab), v in reported.items():
    genus_of_ref = lab.split()[0]
    expected = "Sphingobacterium" if st.LINEAGE[mag] == "Sphingobacterium" else "Pseudomonas"
    assert genus_of_ref == expected, (mag, lab)

# ------------------------------------------------------------------ page
H = 132.0
fig, ax_mm, text_mm, letter = st.page(H)

# ---- A: the full 6 x n_refs comparison; column labels sit ABOVE the matrix so that the
#         long organism names do not hang into panel B
axA = ax_mm(26.0, 50.0, 148.0, 26.0)
nr, nc = len(st.HEROES), n_refs
for i, mag in enumerate(st.HEROES):
    for j, lab in enumerate(man.label):
        v = ani.loc[mag, lab]
        if np.isnan(v):
            axA.text(j, i, "\u2013", ha="center", va="center", fontsize=FS_STAT, color=GREY)
        else:
            axA.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                    facecolor=st.hero_col(mag), edgecolor="none"))
            axA.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=FS_STAT,
                     color="white")
axA.set_xlim(-0.5, nc - 0.5)
axA.set_ylim(nr - 0.5, -0.5)
axA.xaxis.set_ticks_position("top")
axA.xaxis.set_label_position("top")
axA.set_xticks(range(nc))
axA.set_xticklabels(man.label, rotation=90, fontsize=FS_STAT, fontstyle="italic",
                    ha="center", va="bottom")
axA.set_yticks(range(nr))
axA.set_yticklabels(st.HEROES)
for tick, mag in zip(axA.get_yticklabels(), st.HEROES):
    tick.set_color(st.hero_col(mag))
axA.tick_params(length=0)
for sp in axA.spines.values():
    sp.set_visible(False)
axA.set_xticks(np.arange(-0.5, nc, 1), minor=True)
axA.set_yticks(np.arange(-0.5, nr, 1), minor=True)
axA.grid(which="minor", color=LIGHT, lw=0.5)

# ---- B: the reported pairs, ranked
axB = ax_mm(58.0, 92.0, 106.0, 30.0)
y = np.arange(len(reported), dtype=float)
XMIN = 88.0
for yy, ((mag, lab), v) in zip(y, reported.items()):
    col = st.hero_col(mag)
    axB.plot([XMIN, v], [yy, yy], lw=1.0, color=col, zorder=2)
    axB.scatter([v], [yy], s=20, color=col, zorder=3)
    axB.text(v + 0.3, yy, f"{v:.2f}", ha="left", va="center", fontsize=FS_STAT, color=TEXT)
axB.axvline(SPECIES_ANI, color=GREY, lw=0.8, ls="--", zorder=1)
axB.text(SPECIES_ANI, 1.04, "95 % species boundary", transform=axB.get_xaxis_transform(),
         ha="center", va="bottom", fontsize=FS_STAT, color=GREY)
axB.set_yticks(y)
axB.set_yticklabels([f"{mag}  vs  {lab}" for (mag, lab), _ in reported.items()],
                    fontsize=FS_STAT)
for tick, ((mag, _lab), _v) in zip(axB.get_yticklabels(), reported.items()):
    tick.set_color(st.hero_col(mag))
axB.set_ylim(len(reported) - 0.5, -0.5)
axB.set_xlim(XMIN, 101.0)
axB.set_xlabel("skani ANI (%)")
axB.grid(axis="x", color=LIGHT, lw=0.5, zorder=0)
axB.set_axisbelow(True)
st.style_axis(axB, left=False)
axB.tick_params(left=False)

letter(13.0, 8.0, "A")
letter(13.0, 84.0, "B")
text_mm(26.0, 8.0, f"MICP-complete MAG  \u00d7  {n_refs} organism-verified reference genomes",
        fontsize=FS_BODY, color=TEXT)
text_mm(26.0, 12.5, "\u2013  no alignment reported by skani", fontsize=FS_STAT, color=GREY)

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S18")
