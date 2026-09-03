"""Suppl Fig S18 - ANI of the MICP-complete MAGs against the curated MICP reference panel (C5).

Panel
  A  nearest curated reference genome for each MICP-complete MAG, with the 95 % species
     threshold.  Only one MAG-reference pair produced an alignment above skani's reporting
     threshold; the other five MAGs are marked as having no alignment reported, which is the
     result rather than a missing series.

REFERENCE-SET INTEGRITY (why this page is sparse)
The build verifies each reference FASTA against the organism its filename claims, by reading
the first sequence header.  Twelve of the fourteen downloaded references contain a DIFFERENT
organism (e.g. the file named Halomonas_elongata holds Bdellovibrio bacteriovorus HD100, and
both files named Sphingobacterium hold Enterobacter and Marinobacter).  Those files are
EXCLUDED here, because drawing them would put false organism names on the figure.  Five
further references in curated_refs.tsv never downloaded at all, Sporosarcina pasteurii - the
canonical MICP organism - among them.  Only Bacillus subtilis 168 and Pseudomonas helleri
survive verification.  The two impossible ANI values in the stored matrix (99.02 % between
"Halomonas elongata" and "Lysinibacillus fusiformis", 96.45 % between "Halomonas pacifica"
and "Sporosarcina ureae") are artefacts of the same mislabelling: each pair is really two
strains of one species.  The consequence for the manuscript is recorded in the journal - the
claim that the Sphingobacterium MAGs stay below 95 % against every reference examined is NOT
supported by C5, because no genuine Sphingobacterium reference was ever in this comparison.
It is supported instead by the independent 63-genome RefSeq screen drawn in Fig. 5c / Table S10,
whose reference set was checked here and is correct.

Sources
  research/additional/C5_panMICP_env/skani_panMICP.matrix  20 x 20 all-vs-all ANI
  research/additional/C5_panMICP_env/curated_refs.tsv      intended reference panel
  research/additional/C5_panMICP_env/refs/*.fna            first header, for identity check
  Table_S20_skani_hero_vs_refs.tsv                         the reported pairs (asserted against)
A matrix entry of 0.00 means skani reported no alignment for that pair, not 0 % identity, so
it is drawn as absent rather than as a value.

Colour meanings on this page
  blue   = Sphingobacterium lineage MAG (S13, S16, S23, C22)
  orange = Pseudomonas_E lineage MAG (M1, S26)
  grey   = the 95 % species threshold and the no-alignment markers, neither of them a lineage
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _style as st
import _grp_supp_hi as gh
from _style import GREY, TEXT, LIGHT, FS_STAT, FS_BODY

st.setup()
OUT = HERE / "figures"
SUPP = Path(gh.SUPP)
C5 = Path(gh.ADDITIONAL) / "C5_panMICP_env"

SPECIES_ANI = 95.0  # standard prokaryotic species boundary, a published constant

# ------------------------------------------------------------------ data: matrix
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

heroes = [n for n in names if n.startswith("HERO_")]
assert sorted(h.replace("HERO_", "") for h in heroes) == sorted(st.HEROES)

# ------------------------------------------------------------------ reference identity check
curated = pd.read_csv(C5 / "curated_refs.tsv", sep="\t", header=None,
                      names=["accession", "organism", "habitat"])


def first_header(path):
    with open(path) as fh:
        return fh.readline().lstrip(">").strip()


def identity_ok(expected_name, header):
    """Does the FASTA header name the organism the filename claims?

    Compares the genus and species tokens of the expected name against the header text,
    ignoring strain designations and the numeric suffixes used in the file names.
    """
    tokens = [t for t in expected_name.split("_") if not t.isdigit()]
    genus = tokens[0].lower()
    species = tokens[1].lower() if len(tokens) > 1 else ""
    h = header.lower()
    if genus not in h:
        return False
    return species in ("", "sp") or species in h


audit = []
for _, row in curated.iterrows():
    fna = C5 / "refs" / f"{row.organism}.fna"
    if not fna.exists():
        audit.append(dict(organism=row.organism, habitat=row.habitat, status="not downloaded",
                          actual=""))
        continue
    header = first_header(fna)
    ok = identity_ok(row.organism, header)
    audit.append(dict(organism=row.organism, habitat=row.habitat,
                      status="verified" if ok else "wrong organism",
                      actual=header))
audit = pd.DataFrame(audit)
verified = set(audit.loc[audit.status == "verified", "organism"])
print(f"  reference panel: {len(audit)} curated | "
      f"{(audit.status == 'verified').sum()} verified | "
      f"{(audit.status == 'wrong organism').sum()} wrong organism | "
      f"{(audit.status == 'not downloaded').sum()} not downloaded")
for _, r in audit[audit.status != "verified"].iterrows():
    print(f"    {r.organism:34s} {r.status:16s} {r.actual[:70]}")
audit.to_csv(HERE / "_job" / "S18_reference_identity_audit.csv", index=False)

# ------------------------------------------------------------------ nearest verified reference
idx = {nm: i for i, nm in enumerate(names)}
records = []
for mag in st.HEROES:
    hi = idx[f"HERO_{mag}"]
    best_ref, best_ani = None, np.nan
    for ref in verified:
        if ref not in idx:
            continue
        ani = M[hi, idx[ref]]
        if ani > 0 and (np.isnan(best_ani) or ani > best_ani):
            best_ref, best_ani = ref, ani
    records.append(dict(MAG=mag, ref=best_ref, ani=best_ani))
res = pd.DataFrame(records)

# assert every reported pair in the shipped table is reproduced by the matrix
tab = pd.read_csv(SUPP / "Table_S20_skani_hero_vs_refs.tsv", sep="\t")
for _, r in tab.iterrows():
    a, b = Path(r.Ref_file).stem, Path(r.Query_file).stem
    assert abs(M[idx[a], idx[b]] - r.ANI) < 0.01, (a, b, M[idx[a], idx[b]], r.ANI)
got = res.dropna(subset=["ani"])
assert len(got) == 1 and got.iloc[0].MAG == "S26", res
assert abs(got.iloc[0].ani - tab.query("Ref_name.str.contains('helleri')").ANI.iloc[0]) < 0.01

# ------------------------------------------------------------------ page
H = 59.0
fig, ax_mm, text_mm, letter = st.page(H)

ax = ax_mm(28.0, 16.0, 118.0, 34.0)
y = np.arange(len(res), dtype=float)
XMIN = 88.0

for yy, r in zip(y, res.itertuples()):
    col = st.hero_col(r.MAG)
    if np.isnan(r.ani):
        ax.text(XMIN + 0.4, yy, "no alignment reported", ha="left", va="center",
                fontsize=FS_STAT, color=GREY)
    else:
        ax.plot([XMIN, r.ani], [yy, yy], lw=1.0, color=col, zorder=2)
        ax.scatter([r.ani], [yy], s=22, color=col, zorder=3)
        ax.text(r.ani + 0.35, yy, f"{r.ani:.2f}  {r.ref.replace('_', ' ')}",
                ha="left", va="center", fontsize=FS_STAT, color=TEXT)

ax.axvline(SPECIES_ANI, color=GREY, lw=0.8, ls="--", zorder=1)
ax.text(SPECIES_ANI, 1.02, "95 % species threshold", transform=ax.get_xaxis_transform(),
        ha="center", va="bottom", fontsize=FS_STAT, color=GREY)
ax.set_yticks(y)
ax.set_yticklabels(res.MAG)
ax.invert_yaxis()
ax.set_xlim(XMIN, 100.5)
ax.set_xlabel("ANI to nearest verified reference genome (%)")
ax.grid(axis="x", color=LIGHT, lw=0.5, zorder=0)
ax.set_axisbelow(True)
st.style_axis(ax, left=False)
ax.tick_params(left=False)
for tick, mag in zip(ax.get_yticklabels(), res.MAG):
    tick.set_color(st.hero_col(mag))

letter(15.0, 8.5, "A")
text_mm(146.0, 9.0, f"verified references: {len(verified)} of {len(audit)} curated",
        fontsize=FS_STAT, color=GREY, ha="right")

st.audit(fig)
st.prose_scan(fig)
st.save(fig, OUT, "Fig_S18")
