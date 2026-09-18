#!/usr/bin/env python
"""Resolve every entry of the manuscript reference list against CrossRef.

Input : ../01_Manuscript.md  (numbered list under '## References')
Output: crossref_candidates.json  (top-3 CrossRef hits per entry, raw metadata)
        crossref_match_report.tsv (best hit + similarity checks, for human review)

Nothing is accepted automatically here; 02_build_library.py consumes the report
after the DOI column has been reviewed (column `accept`).
"""
import json, re, sys, time, unicodedata, urllib.parse, urllib.request
from difflib import SequenceMatcher
from pathlib import Path

HERE = Path(__file__).resolve().parent
# the numbered list this script parses survives only in the pre-edit copy (the master now holds the regenerated list)
MD = HERE.parent / "01_Manuscript.pre_refs_260919.md"
if not MD.exists():
    MD = HERE.parent / "01_Manuscript.md"
UA = "upcycling-reflib/1.0 (reference verification script)"


def parse_refs(md_text):
    block = md_text.split("\n## References", 1)[1].split("\n---", 1)[0]
    refs = []
    for line in block.splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if not m:
            continue
        n, raw = int(m.group(1)), m.group(2)
        plain = raw.replace("*", "")
        first_author = plain.split(",")[0].split(" et al")[0].split(".")[0]
        surname = first_author.rsplit(" ", 1)[0] if " " in first_author else first_author
        ym = re.search(r"\b(19|20)\d{2}\b(?=;|\.|\s*;)", plain) or re.search(r"\b(19|20)\d{2}\b", plain)
        # title = text between the author block and the journal (first '. ' to the italic journal)
        tm = re.match(r"^(.*?)\.\s+(.*?)[\.\?]\s+\*", raw)
        title = tm.group(2).replace("*", "") if tm else plain
        refs.append(dict(n=n, raw=raw, surname=surname, year=int(ym.group(0)) if ym else None, title=title))
    return refs


def norm(s):
    s = unicodedata.normalize("NFKD", s or "").encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9 ]+", " ", re.sub(r"<[^>]+>", "", s)).split()


def sim(a, b):
    return SequenceMatcher(None, " ".join(norm(a)), " ".join(norm(b))).ratio()


def query(ref):
    q = urllib.parse.quote(ref["raw"].replace("*", ""))
    url = f"https://api.crossref.org/works?query.bibliographic={q}&rows=3"
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.load(r)["message"]["items"]
        except Exception as e:  # noqa
            time.sleep(3 * (attempt + 1))
            err = e
    print(f"  !! ref {ref['n']} lookup failed: {err}", file=sys.stderr)
    return []


def item_year(it):
    for k in ("published-print", "issued", "published-online", "published"):
        dp = it.get(k, {}).get("date-parts", [[None]])
        if dp and dp[0] and dp[0][0]:
            return dp[0][0]
    return None


def main():
    refs = parse_refs(MD.read_text(encoding="utf-8"))
    assert len(refs) > 0, "reference parser returned nothing"
    print(f"parsed {len(refs)} references")
    cand, rows = {}, []
    for ref in refs:
        items = query(ref)
        cand[ref["n"]] = items
        best, best_s = None, -1
        for it in items:
            s = sim(ref["title"], (it.get("title") or [""])[0])
            if s > best_s:
                best, best_s = it, s
        if best is None:
            rows.append([ref["n"], ref["surname"], ref["year"], "", "", "", "", "", "NO_HIT", ref["title"], ""])
            continue
        fam = (best.get("author") or [{}])[0].get("family", "")
        a_ok = " ".join(norm(fam)) == " ".join(norm(ref["surname"]))
        y = item_year(best)
        y_ok = (y == ref["year"])
        flag = "OK" if (best_s >= 0.90 and a_ok and y_ok) else "REVIEW"
        rows.append([ref["n"], ref["surname"], ref["year"], best.get("DOI", ""), f"{best_s:.2f}", fam, y,
                     (best.get("container-title") or [""])[0], flag, ref["title"], (best.get("title") or [""])[0]])
        print(f"  {ref['n']:>2} {flag:6} sim={best_s:.2f} {ref['surname']} {ref['year']} -> {best.get('DOI','')}")
        time.sleep(0.4)
    (HERE / "crossref_candidates.json").write_text(json.dumps(cand, ensure_ascii=False, indent=1), encoding="utf-8")
    hdr = ["n", "ms_first_author", "ms_year", "doi", "title_sim", "cr_first_author", "cr_year", "cr_journal", "flag", "ms_title", "cr_title"]
    with open(HERE / "crossref_match_report.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    n_ok = sum(r[8] == "OK" for r in rows)
    print(f"OK {n_ok} / REVIEW {sum(r[8]=='REVIEW' for r in rows)} / NO_HIT {sum(r[8]=='NO_HIT' for r in rows)}")


if __name__ == "__main__":
    main()
