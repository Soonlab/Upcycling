#!/usr/bin/env python
"""Port the authors' revised manuscript (Manuscript_revised.docx, last saved 2026-09-17)
into the markdown master chain of this package.

  Manuscript_revised.docx  --pandoc-->  01_Manuscript.md + 02_Figure_legends.md

What is taken from the revised docx: title, authors, affiliations, highlights, abstract,
keywords, every body section (Introduction ... Conclusions), the figure and table legends,
the ethics / competing-interests / funding / data-availability statements.
What is kept from the previous master: the two markdown tables (Table 1, Table 2 bodies)
and the reference machinery (references/04_apply_to_master.py regenerates the list).
What is added: numbered section headings (the audits and the reference scripts key on
them), a CRediT placeholder, the Elsevier generative-AI declaration, and one parenthetical
that ties the figures' group label ("MICP-complete") to the text's "prioritized candidates".

The pre-port files are kept once as *.pre_revised_260923.*  Idempotent.
"""
import json, re, subprocess, shutil, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "Manuscript_revised.docx"
MD = HERE / "01_Manuscript.md"
LEG = HERE / "02_Figure_legends.md"
BK_MD = HERE / "01_Manuscript.pre_revised_260923.md"
BK_LEG = HERE / "02_Figure_legends.pre_revised_260923.md"
PANDOC = "/home/soon/miniconda3/bin/pandoc"

TOP = {"Authors", "Affiliations", "Highlights", "Abstract", "Introduction", "Materials and Methods",
       "Results", "Discussion", "5. Conclusions", "Figure and Table Legends", "Supplementary figures",
       "Tables", "Ethics and biosafety statement", "Declaration of competing interests", "Funding",
       "Data and code availability", "References"}

# text edits that are not citation changes (those live in references/04_apply_to_master.py)
TEXT_EDITS = [
    # the figures and the main tables label the six-MAG group "MICP-complete"; the revised text calls
    # them "prioritized candidates" and never uses the figure label -> tie the two once, in Methods
    ("were prioritized for detailed comparative analysis based on",
     "were prioritized for detailed comparative analysis (labelled \"MICP-complete\" in the figures and tables) based on", 1),
]

GENAI = ("During the preparation of this work the authors used Claude (Anthropic) to check the internal "
         "consistency of the manuscript, to verify the reference list against CrossRef and PubMed, and to "
         "edit language. After using this tool, the authors reviewed and edited the content as needed and "
         "take full responsibility for the content of the published article. "
         "AUTHOR VERIFICATION REQUIRED: confirm the tool name and the scope of use.")

CREDIT = ("AUTHOR VERIFICATION REQUIRED: assign CRediT roles to the five authors. Role sets carried from the "
          "previous master for reassignment: **[Author]** — Conceptualization, Data curation, Formal analysis, "
          "Investigation, Methodology, Software, Visualization, Writing – original draft. **[Author]** — "
          "Methodology, Investigation, Software, Writing – review & editing. **[Corresponding author]** — "
          "Conceptualization, Funding acquisition, Project administration, Supervision, Writing – review & editing.")


def unescape(s):
    return (s.replace("\\<", "<").replace("\\>", ">").replace("\\~", "~").replace("\\$", "$")
             .replace("\\[", "[").replace("\\]", "]").replace("\\'", "'"))


def paragraphs():
    raw = subprocess.run([PANDOC, str(SRC), "-t", "markdown-smart", "--wrap=none"],
                         capture_output=True, text=True, check=True).stdout
    raw = re.sub(r"^\*\*\\\*\*\s*$", "", raw, flags=re.M)  # page-break runs
    return [unescape(p.strip()) for p in re.split(r"\n\s*\n", raw) if p.strip()]


def sectionise(paras):
    """-> title, ordered list of (heading, [paragraphs])."""
    title, rest = paras[0], paras[1:]
    secs, cur = [], None
    for p in rest:
        m = re.fullmatch(r"\*\*(.+?)\*\*", p)
        h = m.group(1).strip() if m else (re.sub(r"\\", "", p) if re.fullmatch(r"5\\?\. Conclusions", p) else None)
        if h is not None:
            cur = (h, [])
            secs.append(cur)
        else:
            assert cur is not None, f"paragraph before any heading: {p[:60]!r}"
            cur[1].append(p)
    return title, secs


def main():
    if not BK_MD.exists():
        shutil.copy2(MD, BK_MD)
    if not BK_LEG.exists():
        shutil.copy2(LEG, BK_LEG)
    old = BK_MD.read_text(encoding="utf-8")

    title, secs = sectionise(paragraphs())
    heads = [h for h, _ in secs]
    for h in ("Authors", "Affiliations", "Highlights", "Abstract", "Introduction", "Materials and Methods",
              "Results", "Discussion", "5. Conclusions", "Figure and Table Legends", "Ethics and biosafety statement",
              "Declaration of competing interests", "Funding", "Data and code availability", "References"):
        assert h in heads, f"heading not found in the revised docx: {h}"
    get = lambda h: next(ps for hh, ps in secs if hh == h)

    # ---------------------------------------------------------------- front matter
    authors = [a.strip() for a in re.split(r",\s*(?=[A-Z])", get("Authors")[0])]
    affils = [p for p in get("Affiliations") if not p.startswith("^†^")]
    highlights = [p for p in get("Highlights") if p.startswith("- ")]
    abstract = [p for p in get("Abstract") if not p.startswith("**Keywords")]
    keywords = next(p for p in get("Abstract") if p.startswith("**Keywords"))
    assert len(abstract) == 1 and len(highlights) == 5, (len(abstract), len(highlights))
    n_abs = len(abstract[0].split())

    # ---------------------------------------------------------------- body
    def numbered(parent_heading, stop, n):
        i0 = heads.index(parent_heading)
        i1 = heads.index(stop)
        out = []
        k = 0
        for h, ps in secs[i0 + 1:i1]:
            k += 1
            out.append(f"### {n}.{k} {h}\n" + "\n\n".join(ps))
        return out

    body = ["## 1. Introduction\n" + "\n\n".join(get("Introduction")),
            "## 2. Materials and Methods\n" + "\n\n".join(numbered("Materials and Methods", "Results", 2)),
            "## 3. Results\n" + "\n\n".join(numbered("Results", "Discussion", 3)),
            "## 4. Discussion\n" + "\n\n".join(numbered("Discussion", "5. Conclusions", 4)),
            "## 5. Conclusions\n" + "\n\n".join(get("5. Conclusions"))]
    body = "\n\n".join(body)
    for a, b, cnt in TEXT_EDITS:
        if body.count(a) == cnt:
            body = body.replace(a, b)
        else:
            assert body.count(b) >= 1 and body.count(a) == 0, f"edit does not apply: {a[:50]!r}"

    # ---------------------------------------------------------------- back matter
    back = [
        "## CRediT author statement\n" + CREDIT,
        "## Declaration of generative AI and AI-assisted technologies in the writing process\n" + GENAI,
        "## Ethics and biosafety statement\n" + "\n\n".join(get("Ethics and biosafety statement")),
        "## Declaration of competing interests\n" + "\n\n".join(get("Declaration of competing interests")),
        "## Funding\n" + "\n\n".join(get("Funding")),
        "## Data and code availability\n" + "\n\n".join(get("Data and code availability")),
    ]

    # ---------------------------------------------------------------- legends and tables
    items = {}  # "Figure 1" / "Table S1" -> [title, [paras]]
    order = []
    for h in ("Figure and Table Legends", "Supplementary figures", "Tables"):
        cur = None
        for p in get(h):
            m = re.match(r"^\*\*(Figure S?\d+)\.\*\*\s*(.+)$", p) or re.match(r"^(Table S?\d+)\.\s*(.+)$", p)
            if m:
                cur = m.group(1)
                items[cur] = [m.group(2).strip(), []]
                order.append(cur)
            else:
                assert cur, p[:60]
                items[cur][1].append(p)
    expect = [f"Figure {i}" for i in range(1, 6)] + [f"Figure S{i}" for i in range(1, 6)] + \
             ["Table 1", "Table 2", "Table S1", "Table S2", "Table S3"]
    assert order == expect, order

    def leg_block(k):
        return f"### {k} | {items[k][0]}\n" + "\n\n".join(items[k][1])

    legends = ["# Figure and Table Legends (ported from the authors' revised manuscript of 2026-09-17; port 2026-09-23)",
               "Figures are the 2026-09-09 pages in `Figures/`; the supplementary workbooks follow the S1A/S2A/S3A sheet scheme of the revised text.",
               "## Main figures", *[leg_block(f"Figure {i}") for i in range(1, 6)],
               "## Supplementary figures", *[leg_block(f"Figure S{i}") for i in range(1, 6)],
               "## Tables", *[leg_block(k) for k in ("Table 1", "Table 2", "Table S1", "Table S2", "Table S3")]]

    def pipe_table(n):
        m = re.search(rf"^## Table {n} \|.*?\n(.*?)(?=^## Table|\Z)", old, re.S | re.M)
        rows = [l for l in m.group(1).splitlines() if l.startswith("|")]
        assert len(rows) >= 3, f"Table {n} body not found in the previous master"
        return "\n".join(rows)

    tables = [f"## Table {n} | {items[f'Table {n}'][0]}\n" + "\n\n".join(items[f"Table {n}"][1]) + "\n\n" + pipe_table(n)
              for n in (1, 2)]

    # ---------------------------------------------------------------- assemble
    fm = ["---", "title: " + json.dumps(title, ensure_ascii=False), "author:"] + \
         [f"  - {json.dumps(a, ensure_ascii=False)}" for a in authors] + ['date: "2026"', "---"]
    banner = ["<!--",
              "**Target journal:** *Microbiological Research* (Elsevier; Harvard author–date references)",
              "**Manuscript type:** Full-length research article — comparative genomics of metagenome-assembled genomes.",
              "**Display items:** 5 main figures, 2 main tables, 5 supplementary figures, 3 supplementary table workbooks.",
              f"**Word count (body Intro→Conclusions):** 0 | Abstract: {n_abs}",
              "**Source:** Manuscript_revised.docx (authors' revision, saved 2026-09-17), ported by port_revised_260923.py; "
              "the reference list is regenerated by references/04_apply_to_master.py.",
              "-->"]
    front = "\n".join(fm) + "\n\n" + "\n".join(banner) + "\n\n" + "\n\n".join(affils) + \
            "\n\n^†^Corresponding author. E-mail: AUTHOR VERIFICATION REQUIRED\n\n**Highlights**\n\n" + \
            "\n\n".join(highlights) + "\n\n---\n\n## Abstract\n\n" + abstract[0] + "\n\n" + keywords
    out = front + "\n\n" + body + "\n\n" + "\n\n".join(back) + \
          "\n\n## References\n\n(regenerated by references/04_apply_to_master.py)\n\n---\n\n" + "\n\n".join(tables) + "\n"
    MD.write_text(out, encoding="utf-8")
    LEG.write_text("\n\n".join(legends) + "\n", encoding="utf-8")
    print(f"ported: {len(authors)} authors, {len(affils)} affiliations, abstract {n_abs} words, "
          f"{len(items)} legends; tables 1-2 bodies carried from the previous master")

    r = subprocess.run([sys.executable, str(HERE / "references" / "04_apply_to_master.py")], capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode:
        print(r.stderr)
        sys.exit(r.returncode)


if __name__ == "__main__":
    main()
