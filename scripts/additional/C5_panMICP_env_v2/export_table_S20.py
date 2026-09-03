#!/usr/bin/env python
"""Export Table S20 from the accession-verified C5 re-run.

Every reference is named from `reference_manifest.tsv` — the organism NCBI reports for
that accession — never from a filename.  The original table named references by
filename, and 18 of 20 filenames did not describe the genome inside.

skani reports a pair only when the genomes are similar enough to align (roughly
ANI > 80 % with the run's --min-af 30).  A pair with no reported alignment is written
with an empty ANI and alignment_reported = FALSE rather than as 0.00, because 0.00 in
the previous table was read as an ANI value.
"""
import csv
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUT_TSV = Path("/data/data/Upcycling/SUBMISSION/Supplementary_tables/"
               "Table_S20_skani_hero_vs_refs.tsv")

man = pd.read_csv(HERE / "reference_manifest.tsv", sep="\t")
man = man[man.status == "verified"].copy()
acc2org = dict(zip(man.accession, man.reported_organism))
acc2hab = dict(zip(man.accession, man.habitat))
acc2strain = dict(zip(man.accession, man.strain.fillna('')))
acc2intended = dict(zip(man.accession, man.taxon))
acc2old = dict(zip(man.accession, man.old_accession))

d = pd.read_csv(HERE / "skani_hero_vs_refs.tsv", sep="\t")
d["query_mag"] = d.Query_file.str.extract(r"HERO_(\w+)\.fna")[0]
d["ref_stem"] = d.Ref_file.str.replace(r".*/", "", regex=True).str.replace(".fna", "", regex=False)
pairs = {(q, r): row for (q, r), row in
         zip(zip(d.query_mag, d.ref_stem), d.to_dict("records"))}

HEROES = ["S13", "S16", "S23", "C22", "M1", "S26"]
LINEAGE = {"S13": "Sphingobacterium", "S16": "Sphingobacterium",
           "S23": "Sphingobacterium", "C22": "Sphingobacterium",
           "M1": "Pseudomonas_E", "S26": "Pseudomonas_E"}

rows = []
for mag in HEROES:
    for acc in man.accession:
        rec = pairs.get((mag, acc))
        rows.append({
            "MAG": mag,
            "MAG_lineage": LINEAGE[mag],
            "Ref_accession": acc,
            "Ref_organism": acc2org[acc],
            "Ref_strain": acc2strain[acc],
            "Ref_intended_taxon": acc2intended[acc],
            "Ref_habitat": acc2hab[acc],
            "Superseded_accession_in_v1": acc2old[acc],
            "ANI": f"{rec['ANI']:.2f}" if rec else "",
            "Align_fraction_ref": f"{rec['Align_fraction_ref']:.2f}" if rec else "",
            "Align_fraction_query": f"{rec['Align_fraction_query']:.2f}" if rec else "",
            "Alignment_reported": "TRUE" if rec else "FALSE",
        })

with open(OUT_TSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

n_hit = sum(1 for r in rows if r["Alignment_reported"] == "TRUE")
print(f"wrote {OUT_TSV}")
print(f"{len(rows)} rows = {len(HEROES)} MAGs x {len(man)} verified references; "
      f"{n_hit} with a reported alignment")
for r in rows:
    if r["Alignment_reported"] == "TRUE":
        print(f"  {r['MAG']:4s} vs {r['Ref_organism']:45s} ANI {r['ANI']}")
