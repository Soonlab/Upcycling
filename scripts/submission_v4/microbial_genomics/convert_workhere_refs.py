"""Convert the references of the authors' working copy (01_Manuscript_MicrobialGenomics_workhere.docx) from Harvard
author-date to Microbial Genomics Vancouver style, editing the .docx in place and leaving every other edit untouched.

Same rules as build_micgen_package.py: square-bracket numbers in order of first citation, ranges only for >= 3
consecutive numbers, "N. Surname AB, ... Title. Journal Year;Vol:Pages. DOI" with full journal names.
Citations before the Introduction (Data Summary) are removed, so numbering starts in the Introduction.
Only the characters of each citation and of each reference entry's author/journal/volume segment are rewritten;
run formatting elsewhere (for example italic species names in titles) is kept.

Usage: python convert_workhere_refs.py IN.docx OUT.docx
"""
import copy
import re
import sys

import docx
from docx.oxml.ns import qn
from docx.text.run import Run

AVR = "AUTHOR VERIFICATION REQUIRED"
NAME = r"[A-Z][A-Za-zÀ-ɏ'’\-]+(?: [A-Z][A-Za-zÀ-ɏ'’\-]+)*"
CITE = re.compile(rf"^({NAME})(?: and ({NAME}))?(?: et al\.)?, ((?:19|20)\d\d)[a-z]?$")
PAREN = re.compile(r"\(([^()]*?, (?:19|20)\d\d[a-z]?(?=[;)])[^()]*)\)")
NARR = re.compile(rf"({NAME}(?: and {NAME}| et al\.)?) \(((?:19|20)\d\d)[a-z]?\)")


def runs_of(p):
    """All text runs of a paragraph in order, including runs inside hyperlinks."""
    return [Run(r, p) for r in p._p.iter(qn("w:r"))]


def ptext(p):
    return "".join(r.text for r in runs_of(p))


def rewrite(p, start, end, new, italic=None):
    """Replace characters [start, end) of the paragraph text with `new`, keeping run boundaries elsewhere."""
    pos = 0
    first = True
    for r in runs_of(p):
        t = r.text
        a, b = pos, pos + len(t)
        pos = b
        if b <= start or a >= end:
            continue
        s, e = max(start, a) - a, min(end, b) - a
        if first:
            r.text = t[:s] + new + t[e:]
            if italic is not None:
                r.italic = italic
            first = False
        else:
            r.text = t[:s] + t[e:]
    assert not first, (start, end, new)


def vanc_authors(auth):
    etal = auth.endswith(", et al.")
    auth = auth[:-len(", et al.")] if etal else auth
    names = re.findall(r"([^,]+), ((?:[A-Z][a-z]?\.-?)+)(?:, |$)", auth)
    assert ", ".join(f"{a}, {b}" for a, b in names) == auth, auth
    return ", ".join(f"{a.strip()} {b.replace('.', '').replace('-', '')}" for a, b in names) + (", et al." if etal else ".")


def fmt(nums):
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else ", ".join(map(str, nums[i:j + 1])))
        i = j + 1
    return "[" + ", ".join(out) + "]"


def main(src, dst):
    d = docx.Document(src)
    P = d.paragraphs
    ri = [i for i, p in enumerate(P) if p.text.strip() == "References"]
    assert len(ri) == 1, ri
    ri = ri[0]
    # reference entries: consecutive author-date paragraphs after the heading
    ents = []
    for p in P[ri + 1:]:
        t = ptext(p).strip()
        if not t:
            continue
        if not re.match(r"^.+?, \d{4}[a-z]?\. ", t) or t.startswith("Table "):
            break
        ents.append(p)
    keyed = {}
    for p in ents:
        t = ptext(p).strip()
        first = t.split(", ")[0]
        year = re.search(r", (\d{4})[a-z]?\. ", t).group(1)
        second = re.match(r"[^,]+, [A-Z.\-]+, ([^,]+), ", t)
        assert (first, year) not in keyed, (first, year)
        keyed[(first, year)] = (p, second.group(1) if second else None)

    order = []

    def num(item):
        m = CITE.match(item.strip())
        if not m:
            return None
        key = (m.group(1), m.group(3))
        assert key in keyed, f"citation without entry: {item}"
        if m.group(2):
            assert keyed[key][1] == m.group(2), f"second author mismatch: {item}"
        if key not in order:
            order.append(key)
        return order.index(key) + 1

    # front matter (title page to Data Summary / Impact Statement): no citations; numbering starts in the Introduction
    ii = [i for i, p in enumerate(P[:ri]) if re.fullmatch(r"(?:1\. )?Introduction", p.text.strip())]
    assert len(ii) == 1, ii
    ii = ii[0]
    n_removed = 0
    for p in P[:ii]:
        t = ptext(p)
        reps = []
        for m in PAREN.finditer(t):
            keep = [x.strip() for x in m.group(1).split(";") if not CITE.match(x.strip())]
            if len(keep) < len(m.group(1).split(";")):
                a, b = m.span()
                if not keep and t[a - 1:a] == " ":
                    a -= 1
                reps.append((a, b, f"({'; '.join(keep)})" if keep else ""))
        for m in NARR.finditer(t):
            if CITE.match(f"{m.group(1)}, {m.group(2)}"):
                reps.append((m.start(), m.end(), m.group(1)))
        for a, b, new in sorted(reps, reverse=True):
            rewrite(p, a, b, new)
        n_removed += len(reps)

    # in-text citations, document order (paragraph by paragraph, left to right)
    n_groups = 0
    for p in P[ii:ri]:
        t = ptext(p)
        found = []
        for m in PAREN.finditer(t):
            found.append((m.start(), m.end(), "paren", m))
        for m in NARR.finditer(t):
            if not any(a <= m.start() < b for a, b, _, _ in found):
                found.append((m.start(), m.end(), "narr", m))
        found.sort()
        reps = []
        for a, b, kind, m in found:
            if kind == "paren":
                items = [x.strip() for x in m.group(1).split(";")]
                nums, keep = [], []
                for it in items:
                    n = num(it)
                    (nums.append(n) if n else keep.append(it))
                if not nums:
                    continue
                reps.append((a, b, (f"({'; '.join(keep)}) " if keep else "") + fmt(nums)))
            else:
                n = num(f"{m.group(1)}, {m.group(2)}")
                if n:
                    reps.append((a, b, f"{m.group(1)} {fmt([n])}"))
        for a, b, new in sorted(reps, reverse=True):
            rewrite(p, a, b, new)
        n_groups += len(reps)
    assert len(order) == len(ents), f"{len(ents) - len(order)} entries never cited"

    # reference entries: rewrite author/year prefix and journal/volume tail, then reorder
    for i, key in enumerate(order, 1):
        p = keyed[key][0]
        t = ptext(p)
        lead = len(t) - len(t.lstrip())
        m = re.match(r"^(?P<auth>.+?), (?P<year>\d{4})[a-z]?\. ", t[lead:])
        auth, year = m.group("auth"), m.group("year")
        tail = re.search(r" (?P<jour>[^.]+?) (?P<vol>\d+)(?P<iss>\(\d+\))?, (?P<pages>[A-Za-z]*\d+(?:–[A-Za-z]*\d+)?)\."
                         r"(?: (?P<doi>https://doi\.org/\S+))?\s*$", t)
        if tail:
            # the journal must be the italic run(s) just before the volume
            new_tail = f" {tail.group('jour')} {year};{tail.group('vol')}{tail.group('iss') or ''}:{tail.group('pages')}."
            if tail.group("doi"):
                new_tail += f" {tail.group('doi')}"
            rewrite(p, tail.start(), len(t), new_tail, italic=False)
        else:
            sm = re.search(r" (?P<url>https?://\S+)\s*$", t)
            assert sm, t
            rewrite(p, sm.start(), len(t), f" {sm.group('url')} ({year}; accessed {AVR})")
        rewrite(p, lead, lead + m.end(), f"{i}. {vanc_authors(auth)} ")
    anchor = P[ri]._p
    for key in reversed(order):
        el = keyed[key][0]._p
        el.getparent().remove(el)
        anchor.addnext(el)
    d.save(dst)
    print(f"front-matter citations removed: {n_removed}; citation groups converted: {n_groups}; references numbered: {len(order)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
