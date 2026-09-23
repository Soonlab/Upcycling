#!/usr/bin/env python
"""SUBMISSION_v4 workbook edits (2026-09-23): remove every table that depended on the
false "waste source" variable (the MAG prefixes are binning tools; Perez-Valera and
Elhottova, 2025) and add the deposited provenance / read-mapping coverage.

Table S2  S2P_abundance_proxy  -> S2P_MAG_provenance  (origin, sample, day, assembly type,
                                  binner, CheckM2 quality, read-mapping coverage, SRA /
                                  BioSample accessions, from Zenodo 10.5281/zenodo.15309541)
          S2T_PCoA_coordinates    Source column (binning tool) -> Origin (true sample origin)
Table S3  S3C_PERMANOVA_global    source columns dropped (genus only)
          S3D pairwise-by-source, S3G/S3H MGnify, S3L abundance-by-source  removed
          remaining sheets re-lettered S3A-S3I in the original order; new S3J = read-mapping
          coverage, MICP-complete vs rest (Mann-Whitney U)
Input: the SUBMISSION_v2 workbooks (read-only).  Idempotent: always rebuilds from v2.
"""
import shutil
from pathlib import Path
import numpy as np, pandas as pd, openpyxl
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
V2 = HERE.parent / "SUBMISSION_v2/Supplementary_tables"
OUT = HERE / "Supplementary_tables"
PROV = pd.read_csv(HERE / "data_provenance_260923/MAG_provenance.tsv", sep="\t")
assert len(PROV) == 111
HERO = ["S13", "S16", "S23", "C22", "M1", "S26"]

def copy_sheet(ws_src, wb_dst, title):
    ws = wb_dst.create_sheet(title)
    for row in ws_src.iter_rows(values_only=True):
        ws.append(list(row))
    return ws

def df_to_sheet(wb, title, df):
    ws = wb.create_sheet(title)
    ws.append(list(df.columns))
    for r in df.itertuples(index=False):
        ws.append([None if (isinstance(v, float) and np.isnan(v)) else v for v in r])
    return ws

# ------------------------------------------------------------------ Table S2
src2 = openpyxl.load_workbook(V2 / "Table_S2_per_MAG_measurements.xlsx")
wb2 = openpyxl.Workbook(); wb2.remove(wb2.active)
readme2 = []
for ws in src2.worksheets:
    name = ws.title
    if name == "README":
        rows = list(ws.iter_rows(values_only=True))
        header, body = rows[0], rows[1:]
        continue
    if name == "S2P_abundance_proxy":
        df = PROV.copy()
        df_to_sheet(wb2, "S2P_MAG_provenance", df)
        continue
    if name == "S2T_PCoA_coordinates":
        df = pd.DataFrame(ws.values); df.columns = df.iloc[0]; df = df.iloc[1:].reset_index(drop=True)
        origin = PROV.set_index("MAG").origin
        df["Source"] = df.MAG.map(origin)
        df = df.rename(columns={"Source": "Origin"})
        assert df.Origin.notna().all()
        df_to_sheet(wb2, name, df)
        continue
    copy_sheet(ws, wb2, name)
# README
new_rows = []
for r in body:
    if r[0] == "S2P_abundance_proxy":
        new_rows.append(("S2P_MAG_provenance", "Zenodo 10.5281/zenodo.15309541 MAGs_data.xlsx (Perez-Valera and Elhottova, 2025)",
                         "Sample origin (manure / manured soil / soil), microcosm experiment and day, assembly type, binning tool, CheckM2 quality, read-mapping genome coverage and NCBI SRA / BioSample accessions of every MAG"))
    elif r[0] == "S2T_PCoA_coordinates":
        new_rows.append((r[0], r[1], "Pan-genome PCoA coordinates PC1-PC3 with true sample origin, genus and MICP-complete flag (the former Source column encoded the binning tool and was replaced 2026-09-23)"))
    else:
        new_rows.append(r)
ws = wb2.create_sheet("README", 0)
ws.append(list(header))
for r in new_rows:
    ws.append(list(r))
order = ["README"] + [n for n in [w.title for w in src2.worksheets] if n != "README"]
order = [("S2P_MAG_provenance" if n == "S2P_abundance_proxy" else n) for n in order]
wb2._sheets = [wb2[n] for n in order]
wb2.save(OUT / "Table_S2_per_MAG_measurements.xlsx")

# ------------------------------------------------------------------ Table S3
src3 = openpyxl.load_workbook(V2 / "Table_S3_comparative_statistics.xlsx")
DROP = {"S3D_PERMANOVA_pairwise", "S3G_MGnify_catalog_summary", "S3H_MGnify_cluster_profile", "S3L_abundance_by_source"}
keep = [w.title for w in src3.worksheets if w.title != "README" and w.title not in DROP]
relabel = {old: f"S3{chr(65 + i)}_" + old.split("_", 1)[1] for i, old in enumerate(keep)}
wb3 = openpyxl.Workbook(); wb3.remove(wb3.active)
rows = list(src3["README"].iter_rows(values_only=True)); header3, body3 = rows[0], rows[1:]
ws = wb3.create_sheet("README"); ws.append(list(header3))
for r in body3:
    if r[0] in DROP:
        continue
    new = relabel[r[0]]
    if r[0] == "S3C_PERMANOVA_global":
        ws.append([new, r[1], "PERMANOVA of pan-genome (Jaccard) distance by GTDB-Tk genus, with PC1/PC2 variance (the by-source test was dropped 2026-09-23: the source labels were the binning tools)"])
    else:
        ws.append([new, r[1], r[2]])
ws.append(["S3J_read_coverage_group_test", "Zenodo 10.5281/zenodo.15309541 MAGs_data.xlsx",
           "Read-mapping genome coverage (CoverM, deposited by the data producers), MICP-complete versus rest, two-sided Mann-Whitney U"])
for old in keep:
    if old == "S3C_PERMANOVA_global":
        df = pd.DataFrame(src3[old].values); df.columns = df.iloc[0]; df = df.iloc[1:]
        df = df[["pseudo_F_genus", "p_genus", "PC1_var", "PC2_var"]]
        df_to_sheet(wb3, relabel[old], df)
    else:
        copy_sheet(src3[old], wb3, relabel[old])
h = PROV[PROV.MAG.isin(HERO)].genome_coverage_x; r = PROV[~PROV.MAG.isin(HERO)].genome_coverage_x
cov = pd.DataFrame([
    {"group": "MICP-complete", "n": len(h), "mean_x": round(h.mean(), 2), "median_x": round(h.median(), 2), "min_x": round(h.min(), 2), "max_x": round(h.max(), 2)},
    {"group": "rest", "n": len(r), "mean_x": round(r.mean(), 2), "median_x": round(r.median(), 2), "min_x": round(r.min(), 2), "max_x": round(r.max(), 2)}])
cov["MWU_two_sided_P"] = round(float(mannwhitneyu(h, r, alternative="two-sided").pvalue), 4)
df_to_sheet(wb3, "S3J_read_coverage_group_test", cov)
wb3.save(OUT / "Table_S3_comparative_statistics.xlsx")

# S1 is carried over unchanged
shutil.copy2(V2 / "Table_S1_reference_panels_and_methods.xlsx", OUT / "Table_S1_reference_panels_and_methods.xlsx")
print("relabelled:", relabel)
print(cov.to_string(index=False))
print("S2 sheets:", openpyxl.load_workbook(OUT / "Table_S2_per_MAG_measurements.xlsx", read_only=True).sheetnames)
print("S3 sheets:", openpyxl.load_workbook(OUT / "Table_S3_comparative_statistics.xlsx", read_only=True).sheetnames)
