"""Main Table 2 - trait modules and dbCAN CAZyme classes enriched in the MICP-complete
group after BH-FDR correction. Merges the shipped Table 2 and Table 3 into one.

Every value is read from the supplementary source; nothing is typed by hand.
"""
import pandas as pd

SRC = "/data/data/Upcycling/SUBMISSION/Supplementary_tables"
trait = pd.read_csv(f"{SRC}/Table_S2c_permutation_statistics.csv")
cazy = pd.read_csv(f"{SRC}/Table_S6d_dbCAN_hero_vs_rest.csv")

Q = 0.05
t = trait[trait.Permutation_q_BH < Q].sort_values("Fold_change", ascending=False)
c = cazy[cazy.Permutation_q_BH < Q].sort_values("Fold_change", ascending=False)

NICE = {"carb_binding": "Carbohydrate-binding modules", "tet_mac": "Tetracycline / macrolide determinants",
        "quorum": "Quorum sensing", "glycoside_hydrolase": "Glycoside hydrolases",
        "oxidative": "Oxidative-stress defence", "Mrp_complex": "Mrp Na+/H+ antiporter complex",
        "Na_H_antiporter": "Na+/H+ antiporter (nha)", "compatible_solute": "Compatible-solute biosynthesis",
        "metal_efflux": "Heavy-metal efflux"}
CAZ = {"GH": "Glycoside hydrolases (GH)", "PL": "Polysaccharide lyases (PL)",
       "CBM": "Carbohydrate-binding modules (CBM)", "CE": "Carbohydrate esterases (CE)",
       "GT": "Glycosyl transferases (GT)", "AA": "Auxiliary activities (AA)"}

SUP = str.maketrans("0123456789", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def q(v):
    """Format a q value: plain decimal at or above 0.01, else m x 10^-e with a real
    superscript exponent."""
    if v >= 0.01:
        return f"{v:.3f}"
    m, e = f"{v:.0e}".split("e")
    return f"{m} \u00d7 10\u207b{str(abs(int(e))).translate(SUP)}"

rows = []
for _, r in t.iterrows():
    rows.append((NICE.get(r.Subcategory, r.Subcategory), "keyword trait scan",
                 f"{r.Hero_mean_per1kCDS:.2f}", f"{r.Rest_mean_per1kCDS:.2f}",
                 f"{r.Fold_change:.2f}", r.Fold_change_CI95, q(r.Permutation_q_BH)))
for _, r in c.iterrows():
    rows.append((CAZ.get(r.Class, r.Class), "dbCAN v12 HMMER",
                 f"{r.Hero_mean_per1kCDS:.2f}", f"{r.Rest_mean_per1kCDS:.2f}",
                 f"{r.Fold_change:.2f}", r.Bootstrap_CI95, q(r.Permutation_q_BH)))

hdr = ("| Module or CAZy class | Source | MICP-complete mean (per 10³ CDS) | Rest mean "
       "(per 10³ CDS) | Fold change | 95 % bootstrap CI | Permutation BH-*q* |")
sep = "|---|---|---|---|---|---|---|"
body = "\n".join("| " + " | ".join(str(x) for x in r) + " |" for r in rows)
print(hdr); print(sep); print(body)
print(f"\n<!-- {len(t)} trait subcategories and {len(c)} CAZy classes reach q < {Q} -->")
