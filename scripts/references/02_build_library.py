#!/usr/bin/env python
"""Build the EndNote-importable library from ref_decisions.tsv.

Every record is fetched from CrossRef by DOI (full author list, journal, volume,
pages), never typed by hand. Outputs:

  Upcycling_MICP_library.ris      verified references (action keep / manual / replace)
  Upcycling_MICP_candidates.ris   proposed replacements for the unverifiable entries
  citekeys.json                   ref number -> {family, year, status, doi}, used by 03 and 04
"""
import csv, html, json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
CACHE = HERE / "crossref_by_doi.json"
UA = "upcycling-reflib/1.0 (reference verification script)"

MANUAL = {
    42: dict(TY="COMP", authors=["Seemann, Torsten"], year="2020",
             title="ABRicate: mass screening of contigs for antimicrobial resistance and virulence genes",
             publisher="GitHub", url="https://github.com/tseemann/abricate"),
}

# Fields CrossRef does not hold for a DOI; each value checked against PubMed.
OVERRIDES = {
    "10.1093/nar/gkz1035": {"volume": "48", "issue": "D1", "page": "D570-D578", "_year": "2020"},  # PMID 31696235
    "10.1099/mgen.0.000685": {"page": "000685"},  # article number
    # CrossRef holds only the first author for these two; full bylines from PubMed (PMID 23648862, 15849316)
    "10.4014/jmb.1212.11087": {"author": [{"family": "Dhami", "given": "Navdeep Kaur"}, {"family": "Reddy", "given": "M. Sudhakara"},
                                          {"family": "Mukherjee", "given": "Abhijit"}]},
    "10.1093/nar/gki524": {"author": [{"family": "Zhang", "given": "Yang"}, {"family": "Skolnick", "given": "Jeffrey"}]},
    "10.1128/microbe.9.111.1": {"subtitle": []},  # CrossRef "subtitle" is the magazine standfirst, not part of the title
    "10.1029/2011wr011714": {"page": "W07519"},  # WRR article number
    "10.1016/bs.aambs.2015.12.002": {"volume": "94"},  # Adv Appl Microbiol serial volume
}

# CrossRef deposits this record in capitals; particles that str.title() gets wrong.
CAPS_FIX = {"DEJONG": "DeJong", "VAN PAASSEN": "van Paassen", "AL QABANY": "Al Qabany"}


def fix_caps(name):
    if name in CAPS_FIX:
        return CAPS_FIX[name]
    return name.title() if (len(name) > 2 and name.isupper()) else name


def fetch(doi, cache):
    key = doi.lower()
    if key in cache:
        return {**cache[key], **OVERRIDES.get(key, {})}
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/()")
    err = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=40) as r:
                cache[key] = json.load(r)["message"]
                time.sleep(0.3)
                break
        except Exception as e:  # noqa
            err = e
            time.sleep(3 * (attempt + 1))
    if key not in cache:
        raise RuntimeError(f"CrossRef fetch failed for {doi}: {err}")
    return {**cache[key], **OVERRIDES.get(key, {})}


def clean(s):
    s = html.unescape(re.sub(r"<[^>]+>", "", s or ""))
    return re.sub(r"\s+", " ", s).strip()


def year_of(m):
    if m.get("_year"):
        return m["_year"]
    for k in ("published-print", "issued", "published-online"):
        dp = m.get(k, {}).get("date-parts", [[None]])
        if dp and dp[0] and dp[0][0]:
            return str(dp[0][0])
    return ""


def authors_of(m):
    people, corporate = [], []
    for a in m.get("author", []):
        if a.get("family"):
            people.append(f"{fix_caps(a['family'])}, {a.get('given', '')}".rstrip(", "))
        elif a.get("name"):
            corporate.append(a["name"] + ",")  # trailing comma = corporate author in EndNote
    # CrossRef lists a consortium first for nbt.3893; the person authors lead the byline.
    return people + corporate if people else corporate


def ris_record(m, label, note=""):
    L = ["TY  - JOUR"]
    L += [f"AU  - {a}" for a in authors_of(m)]
    L.append(f"PY  - {year_of(m)}")
    title = clean((m.get("title") or [""])[0])
    if m.get("subtitle") and clean(m["subtitle"][0]).lower() not in title.lower():
        title += ": " + clean(m["subtitle"][0])
    L.append(f"TI  - {title}")
    L.append(f"T2  - {clean((m.get('container-title') or [''])[0])}")
    if m.get("short-container-title"):
        L.append(f"J2  - {clean(m['short-container-title'][0])}")
    if m.get("volume"):
        L.append(f"VL  - {m['volume']}")
    if m.get("issue"):
        L.append(f"IS  - {m['issue']}")
    page = m.get("page") or m.get("article-number") or ""
    if page:
        sp, _, ep = page.partition("-")
        L.append(f"SP  - {sp}")
        if ep:
            L.append(f"EP  - {ep}")
    if m.get("ISSN"):
        L.append(f"SN  - {m['ISSN'][0]}")
    L.append(f"DO  - {m['DOI']}")
    L.append(f"UR  - https://doi.org/{m['DOI']}")
    L.append(f"LB  - {label}")
    if note:
        L.append(f"N1  - {note}")
    L.append("ER  - ")
    return "\n".join(L) + "\n"


def ris_manual(d, label):
    L = [f"TY  - {d['TY']}"] + [f"AU  - {a}" for a in d["authors"]]
    L += [f"PY  - {d['year']}", f"TI  - {d['title']}", f"PB  - {d['publisher']}", f"UR  - {d['url']}", f"LB  - {label}", "ER  - "]
    return "\n".join(L) + "\n"


def main():
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    rows = list(csv.DictReader(open(HERE / "ref_decisions.tsv", encoding="utf-8"), delimiter="\t"))
    assert len(rows) > 0, "ref_decisions.tsv is empty"

    lib, cand, keys = [], [], {}
    for r in rows:
        n = int(r["n"])
        label = f"MS-R{n:02d}"
        act = r["action"]
        if act == "manual":
            d = MANUAL[n]
            lib.append(ris_manual(d, label))
            keys[n] = dict(family=d["authors"][0].split(",")[0], year=d["year"], status=act, doi="")
        elif act in ("keep", "replace", "uncited"):
            doi = r["replacement_doi"].split("|")[0] if act == "replace" else r["doi"]
            m = fetch(doi, cache)
            lib.append(ris_record(m, label, r["note"].split(" || was:")[0]))
            keys[n] = dict(family=authors_of(m)[0].split(",")[0], year=year_of(m), status=act, doi=doi)
        elif act == "pending":
            keys[n] = dict(family="", year="", status=act, doi="")
            for i, doi in enumerate(filter(None, r["replacement_doi"].split("|"))):
                cand.append(ris_record(fetch(doi, cache), f"CAND-for-MS-R{n:02d}-{chr(97 + i)}",
                                       f"Proposed replacement for unverifiable manuscript ref {n}"))
        else:
            assert act == "drop", f"unknown action {act!r} for ref {n}"

    CACHE.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
    (HERE / "Upcycling_MICP_library.ris").write_text("\n".join(lib), encoding="utf-8")
    (HERE / "Upcycling_MICP_candidates.ris").write_text("\n".join(cand), encoding="utf-8")
    (HERE / "citekeys.json").write_text(json.dumps(keys, ensure_ascii=False, indent=1), encoding="utf-8")
    from collections import Counter
    print(f"library {len(lib)} records | candidates {len(cand)} | actions {dict(Counter(r['action'] for r in rows))}")
    dup = [k for k, c in Counter((fold_key(v["family"]), v["year"]) for v in keys.values() if v["family"]).items() if c > 1]
    assert not dup, f"two library records share family+year (temporary citations would be ambiguous): {dup}"


def fold_key(s):
    import unicodedata
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()


if __name__ == "__main__":
    main()
