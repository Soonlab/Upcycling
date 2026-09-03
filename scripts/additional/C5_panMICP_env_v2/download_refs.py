#!/usr/bin/env python
"""Download each resolved C5 reference genome BY ACCESSION and verify its organism.

Files are named by accession (never by an organism label), and the organism recorded
in the assembly data report is compared against the intended taxon before the file is
accepted.  This is the check the original C5 run lacked: there, the accession list
itself pointed at unrelated organisms and the mislabelled downloads were used as-is.
"""
import csv
import json
import re
import subprocess
import sys
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFS = HERE / "refs"
LOGS = HERE / "logs"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

# extra panel member: the strain the manuscript names for the S26 comparison.
EXTRA = [dict(name="Pseudomonas_helleri_DSM29165", habitat="soil_pseudomonas",
              old_accession="GCF_001043025.1", taxon="Pseudomonas helleri",
              accession="GCF_001043025.1", organism="", level="", category="",
              n_candidates="", resolved_by="manuscript_strain")]


def norm(s):
    """Comparable form of an organism name: lowercase alphanumeric tokens."""
    return re.sub(r"[^a-z0-9 ]", " ", s.lower()).split()


# NCBI taxonomy renames: the intended taxon and the reported organism are the same
# tax_id under a new genus name.  Verified against the taxonomy endpoint, not assumed.
SYNONYM = {"halomonas pacifica": "bisbaumannia pacifica",   # tax_id 77098
           "bacillus megaterium": "priestia megaterium"}    # tax_id 1404


def intended_matches(intended_taxon, reported_organism):
    """The reported organism must contain the genus and species of the intended taxon,
    or of its current name where NCBI has renamed the taxon."""
    got = norm(reported_organism)
    names = [intended_taxon]
    syn = SYNONYM.get(" ".join(norm(intended_taxon)[:2]))
    if syn:
        names.append(syn)
    return any(all(w in got for w in norm(n)[:2]) for n in names)


rows = list(csv.DictReader(open(HERE / "resolved_accessions.tsv"), delimiter="\t"))
rows = [r for r in rows if r["accession"]] + EXTRA

manifest = []
for r in rows:
    acc = r["accession"]
    zpath = LOGS / f"{acc}.zip"
    if not zpath.exists():
        cmd = ["curl", "-sS", "-L", "--max-time", "600", "-o", str(zpath),
               f"{API}/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"]
        out = subprocess.run(cmd, capture_output=True, text=True)
        if out.returncode != 0:
            print(f"{acc}\tDOWNLOAD_FAILED\t{out.stderr[:120]}", flush=True)
            manifest.append(dict(r, status="download_failed", reported_organism="",
                                 file="", n_contigs="", total_bp=""))
            continue
        time.sleep(0.3)

    try:
        z = zipfile.ZipFile(zpath)
    except zipfile.BadZipFile:
        head = zpath.read_bytes()[:200].decode("utf-8", "replace")
        print(f"{acc}\tBAD_ZIP\t{head[:120]}", flush=True)
        manifest.append(dict(r, status="bad_zip", reported_organism="", file="",
                             n_contigs="", total_bp=""))
        continue

    # organism as NCBI reports it for this accession
    reported, strain = "", ""
    for line in z.read("ncbi_dataset/data/assembly_data_report.jsonl").decode().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        if d.get("accession") == acc or not reported:
            org = d.get("organism", {})
            reported = org.get("organismName", "")
            strain = (org.get("infraspecificNames", {}) or {}).get("strain", "")
    fna = [n for n in z.namelist() if n.endswith(".fna")]
    if not fna:
        print(f"{acc}\tNO_FASTA", flush=True)
        manifest.append(dict(r, status="no_fasta", reported_organism=reported,
                             strain=strain, file="", n_contigs="", total_bp=""))
        continue

    ok = intended_matches(r["taxon"], reported)
    if not ok:
        print(f"{acc}\tORGANISM_MISMATCH\tintended={r['taxon']!r} reported={reported!r}",
              flush=True)
        manifest.append(dict(r, status="organism_mismatch", reported_organism=reported,
                             strain=strain, file="", n_contigs="", total_bp=""))
        continue

    dest = REFS / f"{acc}.fna"
    dest.write_bytes(z.read(fna[0]))
    n_contigs = total = 0
    with open(dest) as fh:
        for ln in fh:
            if ln.startswith(">"):
                n_contigs += 1
            else:
                total += len(ln.strip())
    print(f"{acc}\tOK\t{reported}\t{n_contigs} contigs\t{total/1e6:.2f} Mb", flush=True)
    manifest.append(dict(r, status="verified", reported_organism=reported,
                         strain=strain, file=dest.name, n_contigs=n_contigs,
                         total_bp=total))

fields = ["name", "habitat", "old_accession", "taxon", "accession", "organism",
          "level", "category", "n_candidates", "resolved_by", "status",
          "reported_organism", "strain", "file", "n_contigs", "total_bp"]
with open(HERE / "reference_manifest.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(manifest)
n_ok = sum(1 for m in manifest if m["status"] == "verified")
print(f"\n{n_ok}/{len(manifest)} references downloaded and organism-verified")
