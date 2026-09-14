# figures/build — builders for the final figure set

The files one level up (`figures/Fig1–5`, `figures/Fig_S1–S5`, `figures/Graphical_abstract`,
PDF vector + PNG) are the current manuscript figures (2026-09-04 consolidation, last revised
2026-09-09). Each page is built by one script here:

| script | output |
|---|---|
| `build_v2_fig1.py` … `build_v2_fig5.py` | `Fig1` … `Fig5` |
| `build_v2_supS1.py` … `build_v2_supS5.py` | `Fig_S1` … `Fig_S5` |
| `build_v2_graphical_abstract.py` | `Graphical_abstract` |

Shared modules: `_style.py` (palette, type scale, mm layout, self-audit), `_micp_presence.py`
(MICP gene presence from the Bakta CDS tables), `_supp_traits.py`, `_grp_supp_hi.py`.

The scripts read the analysis tables from `/data/data/Upcycling/` (original analysis host; see
the `BASE` / `SUPP` variables at the top of each file) and write to `figures_v2/` next to
the scripts. `bash run_all.sh` rebuilds all eleven pages with the `dram_env` interpreter.

Old → new panel numbering for the previous 29-page set is recorded in
`scripts/consolidation_260904/MAPPING.md`; the superseded builders are in the git history
(commit 02de6ab and earlier, `figures/new_figure/`).
