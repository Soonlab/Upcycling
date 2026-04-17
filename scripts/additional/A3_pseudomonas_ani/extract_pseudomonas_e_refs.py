#!/usr/bin/env python3
"""Extract Pseudomonas_E reference accessions from GTDB-Tk output
for M1 and S26 (both Pseudomonas_E) and produce accession list.
"""
import pandas as pd, re, os

GTDB_TSV = "/data/data/Upcycling/SUBMISSION/Supplementary_tables/Table_S1d_GTDB_Tk_classification.tsv"
OUT = "/data/data/Upcycling/research/additional/A3_pseudomonas_ani"
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(GTDB_TSV, sep="\t")
queries = ["M1", "S26"]
acc_set = set()
per_query = {}

for q in queries:
    row = df[df.user_genome == q].iloc[0]
    accs = []
    # closest genome
    if isinstance(row.closest_genome_reference, str) and row.closest_genome_reference.startswith("GC"):
        accs.append(row.closest_genome_reference)
    if isinstance(row.closest_placement_reference, str) and row.closest_placement_reference.startswith("GC"):
        accs.append(row.closest_placement_reference)
    # other_related_references is a comma-or-semicolon separated list with "acc, species, radius, ANI, AF"
    rel = row["other_related_references(genome_id,species_name,radius,ANI,AF)"]
    if isinstance(rel, str):
        for seg in rel.split(";"):
            m = re.match(r"\s*(GC[AF]_\d+\.\d+)", seg)
            if m:
                accs.append(m.group(1))
    per_query[q] = accs
    acc_set.update(accs)

print(f"[A3] M1 references: {len(per_query['M1'])}")
print(f"[A3] S26 references: {len(per_query['S26'])}")
print(f"[A3] Union (unique): {len(acc_set)}")

with open(f"{OUT}/pseudomonas_e_accessions.txt", "w") as f:
    for a in sorted(acc_set):
        f.write(a + "\n")
print(f"[A3] wrote {OUT}/pseudomonas_e_accessions.txt")

# Also write per-query for reporting
with open(f"{OUT}/pseudomonas_e_refs_per_query.tsv", "w") as f:
    f.write("query\taccession\n")
    for q, accs in per_query.items():
        for a in accs:
            f.write(f"{q}\t{a}\n")
