#!/usr/bin/env python3
"""Check hero-clade monophyly and topology in ureC gene tree vs species tree."""
from ete3 import Tree

HERO=["C22","M1","S13","S16","S23","S26"]
OUT="/data/data/Upcycling/research/revision/ureC_tree"

gt = Tree(f"{OUT}/ureC.treefile", format=1)
st = Tree(f"{OUT}/species_pruned.tre", format=1)

def analyze(tree, label):
    leaves = [t.name for t in tree.get_leaves()]
    hero_present = [h for h in HERO if h in leaves]
    print(f"\n== {label} ==")
    print(f"total leaves = {len(leaves)}; hero present = {hero_present}")
    # check hero monophyly
    if len(hero_present)>=2:
        try:
            tree_ = tree.copy()
            tree_.set_outgroup(tree_.get_leaves()[0].name if tree_.get_leaves()[0].name not in hero_present else tree_.get_leaves()[-1].name)
            ca = tree_.get_common_ancestor(hero_present)
            ca_leaves = set(t.name for t in ca.get_leaves())
            is_mono = ca_leaves == set(hero_present)
            print(f"Hero monophyletic: {is_mono}")
            print(f"Common ancestor subtree has {len(ca_leaves)} taxa; {len(ca_leaves - set(hero_present))} non-hero intruders")
            intruders = sorted(ca_leaves - set(hero_present))
            if intruders: print(f"Intruders: {intruders}")
        except Exception as e:
            print(f"error: {e}")
    # check ureC subclade of Sphingobacterium hero
    sph_hero = [h for h in ["C22","S13","S16","S23"] if h in leaves]
    if len(sph_hero)>=2:
        try:
            ca = tree.get_common_ancestor(sph_hero)
            ca_leaves = set(t.name for t in ca.get_leaves())
            mono = ca_leaves == set(sph_hero)
            print(f"Sphingobacterium-hero monophyletic: {mono}")
            if not mono:
                print(f"  intruders: {sorted(ca_leaves - set(sph_hero))}")
        except Exception as e:
            print(f"error: {e}")

analyze(gt, "ureC gene tree")
analyze(st, "species tree (pruned to ureC-bearing MAGs)")

# Write summary
with open(f"{OUT}/hero_topology_summary.txt","w") as fh:
    import io
    from contextlib import redirect_stdout
    buf=io.StringIO()
    with redirect_stdout(buf):
        analyze(gt, "ureC gene tree")
        analyze(st, "species tree (pruned to ureC-bearing MAGs)")
    fh.write(buf.getvalue())
