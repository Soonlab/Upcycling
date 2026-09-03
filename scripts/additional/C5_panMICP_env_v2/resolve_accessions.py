#!/usr/bin/env python
"""Resolve the correct RefSeq accession for each intended C5 reference taxon.

The accessions in the original `curated_refs.tsv` do not correspond to the organisms
they are labelled with (18 of 20 point at unrelated taxa or do not exist), so the
reference panel is rebuilt by resolving each intended TAXON NAME against the NCBI
datasets API.  Preference order:
  1. the RefSeq reference genome for the taxon
  2. the best-assembled current RefSeq genome (Complete Genome > Chromosome > Scaffold > Contig)
The organism returned by NCBI is recorded so that the downloaded file can be verified
against the intended taxon rather than against a filename.
"""
import json
import subprocess
import sys
import time
import urllib.parse

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}

# intended taxon query for each entry of curated_refs.tsv (category names, not data)
TAXON_QUERY = {
    "Sporosarcina_pasteurii": "Sporosarcina pasteurii",
    "Sporosarcina_psychrophila": "Sporosarcina psychrophila",
    "Sporosarcina_ureae": "Sporosarcina ureae",
    "Bacillus_subtilis_168": "Bacillus subtilis subsp. subtilis str. 168",
    "Bacillus_megaterium": "Priestia megaterium",
    "Bacillus_licheniformis": "Bacillus licheniformis",
    "Lysinibacillus_sphaericus": "Lysinibacillus sphaericus",
    "Lysinibacillus_fusiformis": "Lysinibacillus fusiformis",
    "Halomonas_pacifica": "Halomonas pacifica",
    "Halomonas_elongata": "Halomonas elongata",
    "Idiomarina_loihiensis": "Idiomarina loihiensis",
    "Pseudoalteromonas_haloplanktis": "Pseudoalteromonas haloplanktis",
    "Helicobacter_pylori_26695": "Helicobacter pylori 26695",
    "Klebsiella_aerogenes": "Klebsiella aerogenes",
    "Proteus_mirabilis": "Proteus mirabilis",
    "Yersinia_enterocolitica": "Yersinia enterocolitica",
    "Sphingobacterium_sp_21": "Sphingobacterium sp. 21",
    "Sphingobacterium_multivorum": "Sphingobacterium multivorum",
    "Sphingobacterium_paramultivorum": "Sphingobacterium paramultivorum",
    "Pseudomonas_helleri": "Pseudomonas helleri",
}


def fetch(url):
    out = subprocess.run(["curl", "-sS", "--max-time", "90", url],
                         capture_output=True, text=True)
    if out.returncode != 0:
        return None
    try:
        return json.loads(out.stdout)
    except json.JSONDecodeError:
        return None


def query(taxon, reference_only):
    q = urllib.parse.quote(taxon)
    url = (f"{API}/genome/taxon/{q}/dataset_report"
           f"?filters.assembly_source=refseq&filters.assembly_version=current"
           f"&page_size=200")
    if reference_only:
        url += "&filters.reference_only=true"
    return fetch(url)


def pick(taxon):
    for ref_only in (True, False):
        d = query(taxon, ref_only)
        time.sleep(0.4)
        if not d or not d.get("reports"):
            continue
        cands = []
        for r in d["reports"]:
            info = r.get("assembly_info", {})
            if info.get("assembly_status") != "current":
                continue
            cands.append(dict(
                accession=r["accession"],
                organism=r["organism"]["organism_name"],
                level=info.get("assembly_level", "?"),
                category=info.get("refseq_category", "-"),
                rank=(0 if info.get("refseq_category") == "reference genome" else 1,
                      LEVEL_RANK.get(info.get("assembly_level", "?"), 9)),
            ))
        if cands:
            cands.sort(key=lambda c: c["rank"])
            return cands[0], len(cands), ref_only
    return None, 0, None


rows = []
with open("intended_refs.tsv") as fh:
    for line in fh:
        acc_old, name, habitat = line.rstrip("\n").split("\t")
        taxon = TAXON_QUERY[name]
        best, n, ref_only = pick(taxon)
        if best is None:
            rows.append(dict(name=name, habitat=habitat, old_accession=acc_old,
                             taxon=taxon, accession="", organism="", level="",
                             category="", n_candidates=0, resolved_by=""))
            print(f"{name}\tUNRESOLVED", flush=True)
            continue
        rows.append(dict(name=name, habitat=habitat, old_accession=acc_old,
                         taxon=taxon, accession=best["accession"],
                         organism=best["organism"], level=best["level"],
                         category=best["category"], n_candidates=n,
                         resolved_by="reference" if ref_only else "best_current"))
        print(f"{name}\t{best['accession']}\t{best['organism']}\t{best['level']}\t"
              f"{best['category']}\t(of {n})", flush=True)

import csv
with open("resolved_accessions.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"\nwrote resolved_accessions.tsv ({len(rows)} rows)")
