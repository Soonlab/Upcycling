#!/usr/bin/env python
"""Apply the reference decisions to the master manuscript ../01_Manuscript.md.

  1. in-text edits listed in TEXT_EDITS (each must match its stated count, or already be applied);
  2. the typed reference list is regenerated from the library metadata (CrossRef by DOI)
     as an alphabetical author-date (Elsevier Harvard) list holding only cited references.

Titles of references that were already correct are taken from the pre-edit manuscript
(sentence case, species names in italics); corrected / new titles are in TITLES.
The pre-edit manuscript is kept once as ../01_Manuscript.pre_refs_260919.md.
Idempotent: re-running on an already updated manuscript changes nothing.
"""
import csv, json, re, shutil, sys, unicodedata
from pathlib import Path

HERE = Path(__file__).resolve().parent
MD = HERE.parent / "01_Manuscript.md"
BACKUP = HERE.parent / "01_Manuscript.pre_refs_260919.md"

# (old, new) - old must occur exactly `count` times in the body
TEXT_EDITS = [
    ("Achal and Mukherjee, 2015; Cheng et al., 2017; Hamdan et al., 2017; Omoregie et al., 2022)",
     "Achal and Mukherjee, 2015; Kumari et al., 2016; Hamdan et al., 2017; Omoregie et al., 2019)", 1),
    ("(Jiménez-Martínez et al., 2022)", "(Hommel et al., 2015)", 1),
    ("whose growth and nutrient optima (pH ≈ 7–8, freshwater) are incompatible with real waste streams (Krawczyk et al., 2021)",
     "whose auxotrophies and reliance on complex laboratory media are poorly matched to real waste streams (Lapierre et al., 2020)", 1),
    ("(Lagesen et al., 2007; Lee et al., 2022)", "(Lagesen et al., 2007; Yuan et al., 2015)", 1),
    ("Parks et al., 2017; Nayfach et al., 2021; Li et al., 2021)", "Parks et al., 2017; Chen et al., 2021; Nayfach et al., 2021)", 1),
    ("(Stegen et al., 2013; Zamanzadeh et al., 2023)", "(Stegen et al., 2013; Gupta et al., 2016)", 1),
    ("by back-translating MAFFT protein MSAs (Suzuki et al., 2022)", "by back-translating MAFFT protein MSAs in Biopython (Cock et al., 2009)", 1),
    ("isoelectric point with Biopython (Cock et al., 2009) and", "isoelectric point with Biopython and", 1),
    ("(Dhami et al., 2014", "(Dhami et al., 2013", 3),
    # 2026-09-19 (2): Stegen 2013 removed (did not support the sentence); tools named in Methods now cited
    ("(Stegen et al., 2013; Gupta et al., 2016)", "(Gupta et al., 2016)", 1),
    ("(r220; Chaumeil et al., 2022)", "(r220; Parks et al., 2020; Chaumeil et al., 2022)", 1),
    ("UniRef90, Pfam, dbCAN v12 and MEROPS", "UniRef90, Pfam, dbCAN v12 (Zheng et al., 2023) and MEROPS", 1),
    ("computed with PAML yn00 over all taxon pairs", "computed with PAML yn00 (Yang and Nielsen, 2000) over all taxon pairs", 1),
]

TITLES = {  # sentence case; only where the pre-edit title was wrong, incomplete or absent
    6: "Microbially-induced carbonate precipitation for immobilization of toxic metals",
    9: "Application of microbe-induced carbonate precipitation for copper removal from copper-enriched waters: challenges to future industrial application",
    13: "A revised model for microbially induced calcite precipitation: improvements and new insights based on recent experiments",
    16: "Revealing nutritional requirements of MICP-relevant *Sporosarcina pasteurii* DSM33 for growth improvement in chemically defined and complex media",
    18: "Reconstructing 16S rRNA genes in metagenomic data",
    19: "Expanded catalog of microbial genes and metagenome-assembled genomes from the pig gut microbiome",
    22: "Biocementation of sand by *Sporosarcina pasteurii* strain and technical-grade cementation reagents through surface percolation treatment method",
    32: "Quantifying community assembly processes and identifying features that impose them",
    37: "Current status of cow dung as a bioresource for sustainable development",
    55: "Evolutionary-scale prediction of atomic-level protein structure with a language model",
}
JOURNAL_FIX = {"Proceedings of the National Academy of Sciences": "Proceedings of the National Academy of Sciences of the United States of America",
               "Journal of the Royal Statistical Society Series B: Statistical Methodology": "Journal of the Royal Statistical Society: Series B (Methodological)"}


def fold(s):
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()


def initials(given):
    out = []
    for part in re.split(r"\s+", given.replace(".", ". ").strip()):
        if not part:
            continue
        out.append("-".join(p[0].upper() + "." for p in part.split("-") if p.strip(".")))
    return "".join(out)


def old_titles(md_text):
    block = md_text.split("\n## References", 1)[1].split("\n---", 1)[0]
    out = {}
    for line in block.splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if not m:
            continue
        t = re.match(r"^.*?\.\s+(.*?)([\.\?])\s+\*[^*]+\*\s+\d{4}", m.group(2))
        if t:
            out[int(m.group(1))] = t.group(1) + ("?" if t.group(2) == "?" else "")
    return out


def main():
    sys.path.insert(0, str(HERE))
    lib = __import__("02_build_library")
    if not BACKUP.exists():
        shutil.copy2(MD, BACKUP)
    text = MD.read_text(encoding="utf-8")
    titles0 = old_titles(BACKUP.read_text(encoding="utf-8"))
    assert len(titles0) >= 55, f"only {len(titles0)} titles parsed from the pre-edit list"

    head, rest = text.split("\n## References", 1)
    _old, tail = rest.split("\n---", 1)

    for old, new, count in TEXT_EDITS:
        if head.count(old) == count:
            head = head.replace(old, new)
        else:
            # already applied (new present), or superseded by a later edit in this list (neither present)
            assert head.count(old) == 0, f"edit does not apply cleanly: {old[:60]!r} (found {head.count(old)}x, expected {count})"

    keys = {int(k): v for k, v in json.loads((HERE / "citekeys.json").read_text()).items()}
    cache = json.loads((HERE / "crossref_by_doi.json").read_text())
    entries, skipped = [], []
    for n, k in keys.items():
        if k["status"] == "pending":
            raise SystemExit(f"ref {n} is still pending - decide it in ref_decisions.tsv first")
        # cited?  "Family ... , Year" anywhere in the body
        if not re.search(re.escape(k["family"]) + r"[^();]{0,40}?,\s*" + k["year"], head):
            skipped.append(n)
            continue
        if k["status"] == "manual":
            d = lib.MANUAL[n]
            fam, giv = d["authors"][0].split(", ")
            entries.append((fold(fam) + " 0", k["year"], f"{fam}, {initials(giv)}, {d['year']}. {d['title']}. {d['publisher']}. {d['url']}"))
            continue
        m = lib.fetch(k["doi"], cache)
        people = []
        for a in lib.authors_of(m):
            fam, _, giv = a.partition(", ")
            people.append(f"{fam}, {initials(giv)}" if giv else fam.rstrip(","))
        if len(people) > 6:
            people = people[:6] + ["et al."]
        title = TITLES.get(n) or titles0.get(n)
        assert title, f"no title for ref {n}"
        journal = lib.clean((m.get("container-title") or [""])[0])
        journal = JOURNAL_FIX.get(journal, journal)
        page = (m.get("page") or m.get("article-number") or "").replace("-", "–")
        vol = m.get("volume", "")
        loc = ", ".join(x for x in (vol, page) if x)
        end = "" if title.endswith("?") else "."
        entries.append((fold(people[0].split(",")[0]) + (" 0" if len(people) == 1 else " 1"), k["year"],
                        f"{', '.join(people)}, {k['year']}. {title}{end} *{journal}* {loc}. https://doi.org/{m['DOI']}"))
    entries.sort(key=lambda e: (e[0], e[1]))

    new_list = "\n## References\n\n" + "\n\n".join(e[2] for e in entries) + "\n\n---"
    out = head + new_list + tail

    # keep the stated body word count honest
    lines = out.split("\n")
    s = next(i for i, l in enumerate(lines) if l.startswith("## 1. Introduction"))
    e = next(i for i, l in enumerate(lines) if l.startswith("## CRediT"))
    w = len(" ".join(lines[s:e]).split())
    out, nsub = re.subn(r"(\*\*Word count \(body Intro→Conclusions\):\*\* )[\d,]+", lambda mm: mm.group(1) + f"{w:,}", out)
    assert nsub == 1, f"word count banner not found exactly once ({nsub})"

    MD.write_text(out, encoding="utf-8")
    print(f"reference list: {len(entries)} entries | not cited, left out of the list: {skipped} | body {w} words")


if __name__ == "__main__":
    main()
