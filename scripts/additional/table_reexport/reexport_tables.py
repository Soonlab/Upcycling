#!/usr/bin/env python
"""Re-export the defective Upcycling supplementary tables from their true sources.

Run: /home/soon/miniconda3/envs/dram_env/bin/python reexport_tables.py [--check]

Four defects were found during the 2026-09-04 figure rebuild
(`/data/data/Upcycling/new_figure/_job/JOURNAL.md`).  Each is repaired here from the
source of record, with assertions that pin the repair to the shipped values wherever the
shipped values are known to be right:

  S5a/S5b/S5c  the shipped DRAM distillate is the 2026-04-12 hero-only re-run, whose
               `DRAM.py distill` step crashed (FileNotFoundError on the annotation glob)
               after the hero annotate step had already overwritten the output directory.
               The result carries all-zero rows for the six MICP-complete MAGs.  The
               intact 2026-02-09 full-panel distillate survives in a separate directory.

  S11          `aggregate_biosafety.py` writes four database columns, but it was run at
               04:48 while `combined/plasmidfinder_all.tsv` was only written at 04:58, so
               the PlasmidFinder column silently never appeared.  Re-aggregating from the
               raw abricate output restores it.

  S21a/S21b    `run_abundance_proxy.py` maps the MAG prefix with
               {"C":"cattle","S":"swine","M":"sheep","V":"poultry"} - S and M are swapped
               against the sample-code key (C cattle, M swine, S sheep, V poultry) used
               everywhere else, including Table S9a.

  S6c          S6a/S6b/S6d come from the hero-inclusive `04b_dbcan_final.py`, but S6c is
               the older DRAM-only family table covering 105 MAGs.  The six MICP-complete
               MAGs are added from their hmmsearch tables, with family names normalised to
               the bare family so that the two annotation routes agree
               (`GT2_Glycos_transf_2` vs `GT2`).

  S1a          the per-gene copy counts come from a keyword scan that matched the gene-name
               column over **every feature type**, so the Bakta `5_ureB_sRNA` non-coding RNA
               (RFAM RF02514) was counted as a copy of the urease beta subunit.  For S26 and
               S11 the sRNA is the only *ureB* feature there is, so the table records the
               gene as present where no protein-coding beta subunit is annotated.  The table
               also covers only 45 of the 111 MAGs, and 13 of the MAGs it omits carry a real
               *ure* CDS.  Re-exported CDS-only over all 111 MAGs.

  S15b         inherits the same *ureB* defect (it reports two copies for S26).  The *ureB*
               column is re-exported CDS-only, and the two derived flags that are a pure
               function of it (`urease_core_complete`, `MICP_stoich_complete`) are recomputed
               so the row stays internally consistent.

Upstream analysis files carrying the same defect (the A1 aggregate, the C6 outputs) are
refreshed too, so that a future re-export cannot reproduce the bug.  Everything replaced
is copied to `originals/` first.
"""

import argparse
import re
import shutil
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/data/data/Upcycling")
TAB = BASE / "SUBMISSION/Supplementary_tables"
ZIP = BASE / "SUBMISSION/Supplementary_tables.zip"
REV = Path(__file__).resolve().parent
ORIG = REV / "originals"

DRAM_GOOD = Path("/data/pangenome_work/dram_distilled")
A1 = BASE / "research/additional/A1_biosafety"
C6 = BASE / "research/additional/C6_abundance_proxy"
REVDIR = BASE / "research/revision"
DRAM_ANN = Path("/data/pangenome_work/dram_output/all_annotations.tsv")

HEROES = ["S13", "S16", "S23", "C22", "M1", "S26"]
# sample-code key, stated in SUBMISSION/02_Figure_legends.md
SOURCE_BY_PREFIX = {"C": "cattle", "M": "swine", "S": "sheep", "V": "poultry"}
CLS = ["GH", "GT", "PL", "CE", "AA", "CBM"]

CHECK_ONLY = False
LOG = []


def say(msg):
    print(msg)
    LOG.append(msg)


def backup(path):
    """Copy a file to originals/ once, keeping the first version seen."""
    path = Path(path)
    dest = ORIG / path.name
    if not dest.exists():
        ORIG.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    return dest


def write(path, writer):
    """Back up then write, unless running with --check."""
    path = Path(path)
    if CHECK_ONLY:
        say(f"    [check] would rewrite {path}")
        return
    if path.exists():
        backup(path)
    writer(path)
    say(f"    wrote {path}")


# --------------------------------------------------------------- S5a / S5b / S5c
def fix_dram():
    say("\n## Table S5a / S5b / S5c — DRAM distillate")
    ship_a = TAB / "Table_S5a_DRAM_product.tsv"
    good_a = DRAM_GOOD / "product.tsv"
    assert good_a.exists(), f"intact distillate missing: {good_a}"

    ship = pd.read_csv(ship_a, sep="\t", index_col=0)
    good = pd.read_csv(good_a, sep="\t", index_col=0)
    assert list(ship.index) == list(good.index), "row order differs between distillates"
    assert list(ship.columns) == list(good.columns), "columns differ between distillates"
    assert len(good) == 111, f"expected 111 genomes, found {len(good)}"

    num = ship.select_dtypes(include=[np.number]).columns
    zero = ship.loc[HEROES, num].sum(axis=1)
    if (zero > 0).all():
        # re-run against an already-repaired package
        assert ship.equals(good), "S5a is non-zero but does not match the intact distillate"
        say("  already repaired: Table S5a matches the intact distillate; nothing to do")
        return
    assert (zero == 0).all(), f"shipped hero rows are not all-zero: {zero.to_dict()}"
    good_sums = good.loc[HEROES, num].sum(axis=1)
    assert (good_sums > 0).all(), f"intact hero rows are zero too: {good_sums.to_dict()}"

    differing = [m for m in ship.index if not ship.loc[m].equals(good.loc[m])]
    assert sorted(differing) == sorted(HEROES), \
        f"distillates differ outside the six MICP-complete rows: {sorted(differing)}"
    say(f"  asserted: the two distillates differ in exactly the 6 MICP-complete rows, "
        f"which are all-zero in the shipped file")
    say("  intact per-MAG module-completeness sums: "
        + ", ".join(f"{m} {good_sums[m]:.2f}" for m in HEROES))
    say(f"  shape unchanged: {ship.shape} -> {good.shape}")
    write(ship_a, lambda p: shutil.copy2(good_a, p))

    for shipped_name, src_name in (("Table_S5b_DRAM_metabolism_summary.xlsx",
                                    "metabolism_summary.xlsx"),
                                   ("Table_S5c_DRAM_genome_stats.tsv", "genome_stats.tsv")):
        src = DRAM_GOOD / src_name
        if not src.exists():
            say(f"  {shipped_name}: no matching source in {DRAM_GOOD}, left untouched")
            continue
        dst = TAB / shipped_name
        note = ""
        if src_name.endswith(".tsv"):
            old = pd.read_csv(dst, sep="\t", index_col=0)
            new = pd.read_csv(src, sep="\t", index_col=0)
            same_content = old.fillna("NA").astype(str).equals(new.fillna("NA").astype(str))
            note = ("; content was already identical, replaced only so that all three S5 "
                    "files carry one provenance" if same_content else "")
        say(f"  {shipped_name}: replaced from the intact distillate "
            f"({dst.stat().st_size} -> {src.stat().st_size} bytes){note}")
        write(dst, lambda p, s=src: shutil.copy2(s, p))

    # genome_stats is small enough to check content-wise
    gs_new = pd.read_csv(DRAM_GOOD / "genome_stats.tsv", sep="\t")
    say(f"  S5c rows {len(gs_new)}, columns {list(gs_new.columns)}")


# --------------------------------------------------------------------------- S11
def fix_biosafety():
    say("\n## Table S11 — biosafety counts (PlasmidFinder column restored)")
    dbs = ["card", "vfdb", "resfinder", "plasmidfinder"]
    per_db = {}
    for db in dbs:
        comb = A1 / "combined" / f"{db}_all.tsv"
        assert comb.exists(), f"missing abricate output: {comb}"
        df = pd.read_csv(comb, sep="\t", low_memory=False)
        df["MAG"] = df["#FILE"].str.extract(r"bakta_results/([^/]+)/")
        per_db[db] = df.groupby("MAG").size()
        say(f"  {db}: {len(df)} hits over {df['MAG'].nunique()} MAG(s)")

    mags = sorted(d.name for d in (BASE / "MAGs_FASTA_files/bakta_results").iterdir()
                  if d.is_dir())
    assert len(mags) == 111, f"expected 111 Bakta directories, found {len(mags)}"
    agg = pd.DataFrame(per_db).reindex(mags).fillna(0).astype(int)
    agg.index.name = "MAG"
    agg["group"] = ["MICP_complete" if m in HEROES else "rest" for m in agg.index]

    ship = pd.read_csv(TAB / "Table_S11_biosafety_counts_per_MAG.csv", index_col=0)
    assert list(ship.index) == list(agg.index), "MAG order changed"
    for db in ["card", "vfdb", "resfinder"]:
        assert (ship[db].values == agg[db].values).all(), \
            f"{db} does not reproduce the shipped column"
    assert (ship["group"].values == agg["group"].values).all()
    say("  asserted: CARD / VFDB / ResFinder reproduce the shipped values for all 111 MAGs")

    hits = agg[agg.plasmidfinder > 0]
    say(f"  PlasmidFinder: {int(agg.plasmidfinder.sum())} hit(s) in the whole panel"
        + (f" — {', '.join(f'{m} ({n})' for m, n in hits.plasmidfinder.items())}" if len(hits) else ""))
    say(f"  MICP-complete PlasmidFinder hits: {int(agg.loc[HEROES, 'plasmidfinder'].sum())}")
    say(f"  shape {ship.shape} -> {agg.shape} (column added: plasmidfinder)")

    agg = agg[["card", "vfdb", "resfinder", "plasmidfinder", "group"]]
    write(TAB / "Table_S11_biosafety_counts_per_MAG.csv", lambda p: agg.to_csv(p))
    # refresh the upstream analysis file so a re-run cannot reproduce the omission
    write(A1 / "biosafety_counts_per_MAG.csv", lambda p: agg.to_csv(p))


# ------------------------------------------------------------------- S21a / S21b
def fix_abundance():
    say("\n## Table S21a / S21b — waste-source labels (swine and sheep were swapped)")
    a_path = TAB / "Table_S21a_abundance_proxy_per_MAG.csv"
    a = pd.read_csv(a_path)
    fixed = a.copy()
    fixed["source"] = [SOURCE_BY_PREFIX[m[0].upper()] for m in fixed.MAG]

    s9 = pd.read_csv(TAB / "Table_S9a_PCoA_coordinates.csv")[["MAG", "Source"]]
    m = fixed.merge(s9, on="MAG", how="left")
    assert m.Source.notna().all(), "some MAGs are absent from Table S9a"
    assert (m.source.str.lower() == m.Source.str.lower()).all(), \
        "corrected labels still disagree with Table S9a"
    n_changed = int((a.source.values != fixed.source.values).sum())
    say(f"  asserted: corrected labels agree with Table S9a for all {len(m)} MAGs")
    say(f"  rows relabelled: {n_changed}")
    say("  before: " + ", ".join(f"{k} {v}" for k, v in a.source.value_counts().items()))
    say("  after:  " + ", ".join(f"{k} {v}" for k, v in fixed.source.value_counts().items()))

    # the coverage numbers themselves are untouched
    other = [c for c in a.columns if c != "source"]
    assert a[other].equals(fixed[other]), "a column other than `source` changed"

    b_path = TAB / "Table_S21b_abundance_proxy_per_source.csv"
    b_old = pd.read_csv(b_path)
    b_new = (fixed.groupby("source")["length_weighted_cov"]
             .agg(["count", "mean", "median", "std"]).reset_index())
    b_new = b_new.sort_values("source").reset_index(drop=True)
    say("  per-source length-weighted coverage, shipped -> corrected:")
    old_by = b_old.set_index("source")
    new_by = b_new.set_index("source")
    for s in sorted(set(old_by.index) | set(new_by.index)):
        o = old_by.loc[s] if s in old_by.index else None
        n = new_by.loc[s] if s in new_by.index else None
        say(f"    {s:8s} n {int(o['count']) if o is not None else '-':>3} -> "
            f"{int(n['count']) if n is not None else '-':>3}   "
            f"mean {o['mean']:.2f} -> {n['mean']:.2f}   "
            f"median {o['median']:.2f} -> {n['median']:.2f}")
    say("  the shipped 'sheep 57.4x' is the SWINE group (n = 15); corrected sheep mean is "
        f"{new_by.loc['sheep', 'mean']:.2f}x (n = {int(new_by.loc['sheep', 'count'])})")

    hero_src = fixed[fixed.MAG.isin(HEROES)].set_index("MAG").source
    say("  MICP-complete MAG sources after the fix: "
        + ", ".join(f"{m} {hero_src[m]}" for m in HEROES))

    write(a_path, lambda p: fixed.to_csv(p, index=False))
    write(b_path, lambda p: b_new.to_csv(p, index=False))
    # refresh the upstream C6 outputs
    if (C6 / "abundance_proxy_per_MAG.csv").exists():
        write(C6 / "abundance_proxy_per_MAG.csv", lambda p: fixed.to_csv(p, index=False))
    if (C6 / "abundance_proxy_per_source.csv").exists():
        write(C6 / "abundance_proxy_per_source.csv", lambda p: b_new.to_csv(p, index=False))
    say("  NOTE: the bug lives in scripts/additional/C6_abundance_proxy/run_abundance_proxy.py "
        "line 26 — that prefix map must be corrected before the analysis is re-run.")


# --------------------------------------------------------------------------- S6c
def parse_tbl(path):
    """Best HMM per protein (lowest E-value) — the recipe of 04b_dbcan_final.py."""
    best = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            target, hmm, evalue = p[0], p[2], float(p[4])
            if target not in best or evalue < best[target][1]:
                best[target] = (hmm, evalue)
    return best


FAM_RE = re.compile(r"^(GH|GT|PL|CE|AA|CBM)(\d+)")


def bare_family(name):
    """`GT2_Glycos_transf_2.hmm` -> `GT2`; `GH13` -> `GH13`; else None."""
    m = FAM_RE.match(name.replace(".hmm", ""))
    return f"{m.group(1)}{m.group(2)}" if m else None


def fix_dbcan_families():
    say("\n## Table S6c — dbCAN family counts (the six MICP-complete MAGs were missing)")
    ship = pd.read_csv(TAB / "Table_S6c_dbCAN_family_counts.csv", index_col=0)
    say(f"  shipped: {ship.shape[0]} MAGs x {ship.shape[1]} families, "
        f"MICP-complete present: {sum(h in ship.index for h in HEROES)}/6")

    # --- hero side: direct hmmsearch tables
    hero_fam = defaultdict(lambda: defaultdict(int))
    for h in HEROES:
        tbl = REVDIR / "dbcan_direct" / f"{h}.tbl"
        assert tbl.exists(), f"missing hmmsearch table: {tbl}"
        for hmm, _ in parse_tbl(tbl).values():
            fam = bare_family(hmm)
            if fam:
                hero_fam[h][fam] += 1
    hero_df = pd.DataFrame(hero_fam).T.fillna(0).astype(int)

    # --- rest side: DRAM cazy assignments
    assert DRAM_ANN.exists(), f"missing DRAM annotations: {DRAM_ANN}"
    ann = pd.read_csv(DRAM_ANN, sep="\t", usecols=["fasta", "cazy_best_hit", "cazy_ids"],
                      low_memory=False)
    ann = ann[~ann["fasta"].isin(HEROES)]

    def fam_of(row):
        for field in ("cazy_best_hit", "cazy_ids"):
            v = row[field]
            if pd.notna(v):
                m = re.search(r"\b((?:GH|GT|PL|CE|AA|CBM)\d+)", str(v))
                if m:
                    return m.group(1)
        return None

    ann["fam"] = ann.apply(fam_of, axis=1)
    rest = ann[ann.fam.notna()]
    rest_df = rest.groupby(["fasta", "fam"]).size().unstack(fill_value=0)
    say(f"  rest side: {rest_df.shape[0]} MAGs from {DRAM_ANN.name}")

    combined = pd.concat([hero_df, rest_df]).fillna(0).astype(int)
    combined = combined.reindex(sorted(combined.columns), axis=1)
    combined.index.name = "fasta"

    # --- assertion: per-class family sums must reproduce Table S6a for every MAG
    s6a = pd.read_csv(TAB / "Table_S6a_dbCAN_class_counts.csv", index_col=0)
    combined = combined.reindex(s6a.index)
    assert combined.notna().all().all(), "a MAG in Table S6a is missing from the family table"
    combined = combined.astype(int)
    for c in CLS:
        cols = [f for f in combined.columns if FAM_RE.match(f) and FAM_RE.match(f).group(1) == c]
        got = combined[cols].sum(axis=1)
        assert (got.values == s6a[c].values).all(), (
            f"class {c} family sums do not reproduce Table S6a; "
            f"first mismatch: {got[got.values != s6a[c].values].head(3).to_dict()}")
    say("  asserted: per-class family sums reproduce Table S6a for all "
        f"{len(combined)} MAGs and all six CAZy classes")

    # what the bare-family normalisation changed
    hero_gt2 = int(combined.loc[HEROES, "GT2"].mean()) if "GT2" in combined else None
    rest_gt2 = combined.drop(index=HEROES)["GT2"].mean() if "GT2" in combined else None
    if hero_gt2 is not None:
        say(f"  family-name normalisation example — GT2: MICP-complete mean "
            f"{combined.loc[HEROES, 'GT2'].mean():.2f}, rest mean {rest_gt2:.2f} "
            f"(the shipped table had no MICP-complete rows at all)")
    say(f"  shape {ship.shape} -> {combined.shape}")
    write(TAB / "Table_S6c_dbCAN_family_counts.csv", lambda p: combined.to_csv(p))
    # keep a copy of the hero-inclusive family table beside the other `final` artefacts
    if not CHECK_ONLY:
        combined.to_csv(REVDIR / "dbCAN_final_family_counts.csv")
        say(f"    wrote {REVDIR / 'dbCAN_final_family_counts.csv'}")


# --------------------------------------------------------------------------- S1a
def _cds_counts():
    """CDS-only and all-feature MICP gene copy counts, from the shared scan.

    Imported from `new_figure/_micp_presence.py` rather than reimplemented, so the tables
    and the figures cannot drift apart.
    """
    sys.path.insert(0, str(BASE / "new_figure"))
    from _micp_presence import _scan, GENES as SCAN_GENES
    every, cds = _scan()
    return every, cds, list(SCAN_GENES)


def fix_s1a():
    say("\n## Table S1a — per-gene copy counts (CDS-only, all 111 MAGs)")
    ship_p = TAB / "Table_S1a_ace_samples_list.csv"
    # idempotency: once re-exported the live file holds 111 CDS-only rows, so the
    # provenance assertions below must be checked against the preserved original
    orig_p = ORIG / "Table_S1a_ace_samples_list.csv"
    ship = pd.read_csv(orig_p if orig_p.exists() else ship_p, index_col=0)
    every, cds, genes = _cds_counts()
    say(f"    original: {len(ship)} MAGs x {len(ship.columns)} columns"
        f"{' (live file already re-exported)' if orig_p.exists() else ''}")

    covered = [m for m in ship.index if m in every.index]
    assert len(covered) == len(ship), (len(covered), len(ship))

    # the all-feature rule must reproduce the shipped presence/absence for those 45
    ship_pres = (ship.loc[covered, genes] > 0).astype(int)
    all_pres = (every.loc[covered, genes] > 0).astype(int)
    assert (ship_pres == all_pres).all().all(), \
        "the keyword rule no longer reproduces Table S1a presence/absence"

    # only ureB may change when non-coding features are excluded
    changed = [g for g in genes if (every[g] > 0).ne(cds[g] > 0).any()]
    assert changed == ["ureB"], changed
    lost = sorted(every.index[(every[genes] > 0).sum(1) != (cds[genes] > 0).sum(1)])
    assert lost == ["S11", "S26"], lost
    say(f"    non-coding 5_ureB_sRNA counted as ureB in the shipped rule; "
        f"MAGs losing a gene when CDS-only: {', '.join(lost)}")

    out = cds.copy()
    out["Total_Count"] = out[genes].sum(axis=1)
    out.index.name = ship.index.name or "Sample"
    out = out.sort_values("Total_Count", ascending=False)
    n8 = int((out[genes] > 0).sum(1).eq(8).sum())
    n8_before = int((every[genes] > 0).sum(1).eq(8).sum())
    say(f"    re-exported: {len(out)} MAGs x {len(out.columns)} columns; "
        f"{n8} carry all eight genes as CDS ({n8_before} counting non-coding features)")
    for mag in ("S26", "S11"):
        say(f"      {mag}: ureB {int(every.loc[mag, 'ureB'])} -> {int(cds.loc[mag, 'ureB'])}, "
            f"score {int((every.loc[mag, genes] > 0).sum())} -> "
            f"{int((cds.loc[mag, genes] > 0).sum())}")
    say("    note: ureA/ureC copy counts differ from the shipped table for four MAGs "
        "(C13, S23, S26, V3) where a feature carries a gene name and a product naming a "
        "different subunit; presence is unaffected")
    write(ship_p, lambda q: out.to_csv(q))


# -------------------------------------------------------------------------- S15b
def fix_s15b():
    say("\n## Table S15b — urease beta subunit (CDS-only) and its derived flags")
    ship_p = TAB / "Table_S15b_stoichiometry_per_MAG.csv"
    # idempotency: the provenance assertions describe the ORIGINAL shipped table, whose
    # ureB column is the all-feature scan; once re-exported the live file is CDS-only
    orig_p = ORIG / "Table_S15b_stoichiometry_per_MAG.csv"
    df = pd.read_csv(orig_p if orig_p.exists() else ship_p, index_col=0)
    every, cds, _ = _cds_counts()
    say(f"    original: {len(df)} MAGs x {len(df.columns)} columns"
        f"{' (live file already re-exported)' if orig_p.exists() else ''}")

    # the shipped ureB column is exactly the all-feature scan, so the CDS-only value is
    # well defined; ureA and ureG use a broader keyword rule and are left untouched
    assert (df["ureB"] == every.loc[df.index, "ureB"]).all(), \
        "Table S15b ureB no longer matches the all-feature scan"
    for col in ("ureC", "ureE", "ureF"):
        assert (df[col] == every.loc[df.index, col]).all(), col

    # derived flags, verified against the shipped table before being recomputed
    ca = df[["cah_alphaCA", "canA_gammaCA", "cynT_betaCA", "CA_generic"]].sum(1) > 0
    core = (df.ureA > 0) & (df.ureB > 0) & (df.ureC > 0)
    assert (core.astype(int) == df.urease_core_complete).all()
    assert ((core & ca).astype(int) == df.MICP_stoich_complete).all()

    out = df.copy()
    out["ureB"] = cds.loc[df.index, "ureB"].values
    new_core = (out.ureA > 0) & (out.ureB > 0) & (out.ureC > 0)
    out["urease_core_complete"] = new_core.astype(int)
    out["MICP_stoich_complete"] = (new_core & ca).astype(int)

    changed = out.index[(out != df).any(axis=1)].tolist()
    flagged = out.index[(out.urease_core_complete != df.urease_core_complete)
                        | (out.MICP_stoich_complete != df.MICP_stoich_complete)].tolist()
    say(f"    {len(changed)} MAGs have a lower ureB copy count once the sRNA is excluded "
        f"(the feature occurs 58 times panel-wide)")
    say(f"    MAGs whose derived flags change: {', '.join(flagged) if flagged else 'none'}")
    for mag in flagged:
        say(f"      {mag}: ureB {df.loc[mag,'ureB']} -> {out.loc[mag,'ureB']}, "
            f"urease_core_complete {df.loc[mag,'urease_core_complete']} -> "
            f"{out.loc[mag,'urease_core_complete']}, MICP_stoich_complete "
            f"{df.loc[mag,'MICP_stoich_complete']} -> {out.loc[mag,'MICP_stoich_complete']}")
    say("    ureA and ureG use a broader keyword rule than the shared scan "
        "(75 and 17 MAGs differ) and were left as shipped; only the documented "
        "ureB sRNA defect is repaired")
    write(ship_p, lambda q: out.to_csv(q))


# --------------------------------------------------------------------------- zip
def rebuild_zip():
    say("\n## Supplementary_tables.zip")
    files = sorted(p for p in TAB.iterdir() if p.is_file())
    if CHECK_ONLY:
        say(f"    [check] would rebuild with {len(files)} entries")
        return
    if ZIP.exists():
        backup(ZIP)
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for p in files:
            z.write(p, p.name)
    say(f"    wrote {ZIP} with {len(files)} entries "
        f"({ZIP.stat().st_size / 1e6:.2f} MB)")


def main():
    global CHECK_ONLY
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="run every assertion but write nothing")
    CHECK_ONLY = ap.parse_args().check
    say(f"# Table re-export {'(check only)' if CHECK_ONLY else ''}")
    fix_dram()
    fix_biosafety()
    fix_abundance()
    fix_dbcan_families()
    fix_s1a()
    fix_s15b()
    rebuild_zip()
    say("\nAll assertions passed.")
    (REV / "reexport_run.log").write_text("\n".join(LOG) + "\n")


if __name__ == "__main__":
    sys.exit(main())
