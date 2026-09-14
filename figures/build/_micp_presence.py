"""CDS-only presence/absence of the eight MICP genes across the 111 MAGs.

Used by build_fig1.py (panel B) and build_fig2.py (both panels), which must agree.

Why this module exists
----------------------
Both figures previously read `Table_S1a_ace_samples_list.csv`.  That table's per-gene
copy counts come from a keyword scan of the Bakta annotation that matched the gene-name
column over **every feature type**, so it counted the Bakta `5_ureB_sRNA` non-coding RNA
(RFAM RF02514, product "5' ureB small RNA") as a copy of the urease beta subunit.  The
sRNA occurs 58 times across the panel.  For most MAGs it merely inflates the *ureB* copy
count beside a real CDS, but for two MAGs it is the only *ureB* feature there is, so the
table records the gene as present when no protein-coding beta subunit is annotated:

    S26   0 CDS + 2 sRNA  -> Table S1a reports ureB = 2   (score 8/8 -> 7/8)
    S11   0 CDS + 1 sRNA  -> Table S1a reports ureB = 1   (score  7  ->  6)

S26 is one of the six MICP-complete MAGs, so this is the source of the long-standing
disagreement between Table S1a (which scores S26 8/8) and both Table S1c
(`ureABC_all_present = False`) and the shipped pangenome summary (which scores S26 1/8).
The disagreement is an annotation-parsing artefact, not a biological ambiguity: the
beta subunit is absent from the S26 Bakta annotation.

A second, independent defect surfaced while checking the first.  **Table S1a lists only
45 of the 111 MAGs**, but 58 MAGs carry at least one *ure* CDS and 58 carry a carbonic
anhydrase.  Thirteen MAGs with real *ure* genes (C2, C3, C4, C10, C17, C19, M10, M15,
S15, S28, S29, V2, V9 - up to 7 of the 8 genes, in V9) are absent from the table
altogether, as are the many MAGs whose only MICP gene is *cah*.  Reading absence from
Table S1a as "no MICP gene" therefore drew those MAGs as gene-free and understated the
prevalence of every gene in the non-MICP-complete group.

`presence()` therefore rebuilds the matrix from the Bakta annotations of all 111 MAGs,
counting **CDS features only**.  It asserts that the same keyword rule applied to every
feature type reproduces Table S1a's presence/absence for the 45 MAGs that table covers,
so the departure from the shipped figures is attributable to exactly two causes: the
excluded non-coding features, and the 66 MAGs the table never listed.  Only *ureB* is
affected by the CDS restriction; the other seven genes have no non-CDS features anywhere
in the panel.
"""

import glob
import os

import pandas as pd

BAKTA = "/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
SUPP = "/data/data/Upcycling/SUBMISSION/Supplementary_tables"
GENES = ["ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG", "cah"]


def _classify(product, gene):
    """The Table S1a keyword rule: gene name first, then product text.

    Reproduces Table S1a's presence/absence exactly.  Copy counts differ for four MAGs
    (S23, S26, V3, C13) where a feature carries both a gene name and a product naming a
    different subunit, so alpha and gamma swap between the two routes; presence is
    unaffected, and these figures use presence only.
    """
    p = (product or "").lower()
    g = (gene or "").lower()
    for key in GENES[:7]:
        if key.lower() in g:
            return key
    if "urease subunit alpha" in p:
        return "ureC"
    if "urease subunit beta" in p:
        return "ureB"
    if "urease subunit gamma" in p:
        return "ureA"
    if "ured" in p:
        return "ureD"
    if "carbonic anhyd" in p or g == "cah":
        return "cah"
    return None


def _scan():
    """Per-MAG copy counts over all features and over CDS features only."""
    every, cds = {}, {}
    mags = sorted(os.path.basename(d.rstrip("/")) for d in glob.glob(BAKTA + "/*/"))
    for mag in mags:
        a = dict.fromkeys(GENES, 0)
        c = dict.fromkeys(GENES, 0)
        path = os.path.join(BAKTA, mag, mag + ".tsv")
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                key = _classify(f[7], f[6])
                if key is None:
                    continue
                a[key] += 1
                if f[1].lower() == "cds":
                    c[key] += 1
        every[mag] = a
        cds[mag] = c
    return pd.DataFrame(every).T[GENES], pd.DataFrame(cds).T[GENES]


def presence(verbose=True):
    """Return the CDS-only presence matrix (111 x 8, 0/1), asserting its provenance."""
    every, cds = _scan()
    assert len(cds) == 111, len(cds)

    pres_all = (every > 0).astype(int)
    pres_cds = (cds > 0).astype(int)

    # Provenance check against Table S1a.  That table exists in two states: the shipped
    # 45-MAG version whose keyword scan counted every feature type, and the 2026-09-04
    # re-export (CDS-only, all 111 MAGs) written by
    # SUBMISSION/_revision_260904/tables/reexport_tables.py.  Both are pinned here, so this
    # module stays correct whichever version is on disk and cannot silently drift from it.
    s1a = pd.read_csv(os.path.join(SUPP, "Table_S1a_ace_samples_list.csv"), index_col=0)
    covered = [m for m in s1a.index if m in pres_all.index]
    s1a_pres = (s1a.loc[covered, GENES] > 0).astype(int)
    reexported = len(s1a) == 111

    if reexported:
        # the re-export is this scan, so it must agree with the CDS-only matrix exactly
        assert (pres_cds.loc[covered, GENES] == s1a_pres).all().all(), \
            "the CDS-only scan no longer reproduces the re-exported Table S1a"
        omitted = []
        omitted_with_ure = []
    else:
        # the shipped table used the all-feature rule and listed only part of the panel
        assert (pres_all.loc[covered, GENES] == s1a_pres).all().all(), \
            "the keyword rule no longer reproduces the shipped Table S1a"
        omitted = sorted(set(pres_all.index) - set(s1a.index))
        omitted_with_ure = sorted(m for m in omitted
                                  if pres_cds.loc[m, GENES[:7]].sum() > 0)
        assert len(covered) == 45 and len(omitted) == 66, (len(covered), len(omitted))
        assert not any(pres_cds.loc[m].sum() == 8 for m in omitted), omitted_with_ure

    # excluding non-CDS features may only remove ureB, and only where the sRNA was the
    # sole feature; every other gene must be untouched
    changed_genes = [g for g in GENES if (pres_all[g] != pres_cds[g]).any()]
    assert changed_genes == ["ureB"], changed_genes
    assert (pres_cds <= pres_all).all().all()

    if verbose:
        drop = pres_all.index[pres_all.sum(1) != pres_cds.sum(1)].tolist()
        print(f"  MICP presence: CDS-only scan of {len(pres_cds)} Bakta annotations; "
              f"{int((pres_cds.sum(1) == 8).sum())} MAGs carry all eight genes "
              f"({int((pres_all.sum(1) == 8).sum())} before excluding non-coding features)")
        for mag in drop:
            print(f"    {mag}: score {int(pres_all.sum(1)[mag])} -> "
                  f"{int(pres_cds.sum(1)[mag])} (ureB was 5_ureB_sRNA only)")
        if reexported:
            print(f"    Table S1a (re-exported 2026-09-04) covers all {len(covered)} MAGs "
                  f"and agrees with this scan")
        else:
            print(f"    Table S1a covers {len(covered)} MAGs; of the {len(omitted)} it "
                  f"omits, {len(omitted_with_ure)} carry a ure CDS: "
                  f"{', '.join(omitted_with_ure)}")
    return pres_cds
