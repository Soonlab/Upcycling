#!/usr/bin/env python3
"""Reconstruct Bakta TSV from .gff3 (for S13, S16 which were clobbered)."""
import re, sys, os

def parse_attrs(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def gff3_to_tsv(gff3_path, tsv_out):
    with open(gff3_path) as fin, open(tsv_out, "w") as fout:
        fout.write("# Annotated with Bakta (reconstructed from .gff3)\n")
        for line in fin:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            contig, source, ftype, start, end, score, strand, phase, attrs_s = parts
            if ftype.lower() == "region": continue
            attrs = parse_attrs(attrs_s)
            locus = attrs.get("locus_tag", attrs.get("ID",""))
            gene = attrs.get("gene","")
            product = attrs.get("product","")
            dbxrefs = attrs.get("Dbxref","")
            # write: contig type start end strand locus gene product dbxrefs
            fout.write("\t".join([
                contig, ftype.lower(), start, end, strand, locus, gene, product, dbxrefs
            ]) + "\n")

for mag in ["S13", "S16"]:
    gff = f"/data/data/Upcycling/MAGs_FASTA_files/bakta_results/{mag}/{mag}.gff3"
    tsv = f"/data/data/Upcycling/MAGs_FASTA_files/bakta_results/{mag}/{mag}.tsv"
    # backup the corrupted file first
    os.rename(tsv, tsv + ".carveme_clobbered.bak")
    gff3_to_tsv(gff, tsv)
    n = sum(1 for _ in open(tsv))
    print(f"{mag}: reconstructed {n} lines into {tsv}")
