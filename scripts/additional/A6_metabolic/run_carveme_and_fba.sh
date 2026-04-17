#!/bin/bash
# A6: CarveMe GEM reconstruction for 6 MICP-complete MAGs + FBA
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate carveme

WD=/data/data/Upcycling/research/additional/A6_metabolic
BAKTA=/data/data/Upcycling/MAGs_FASTA_files/bakta_results
mkdir -p $WD/models
cd $WD

HERO="S13 S16 S23 C22 M1 S26"
echo "[A6] reconstructing GEMs for 6 MICP-complete MAGs"

for m in $HERO; do
    faa=$BAKTA/$m/$m.faa
    sbml=$WD/models/${m}.xml
    if [ -f $sbml ]; then
        echo "[A6] $m model exists, skip"
        continue
    fi
    echo "[A6] --- $m ---"
    carve $faa -o $sbml --solver glpk -v 2>&1 | tail -3 || echo "[A6] $m carve FAILED"
done

ls -la $WD/models/

# FBA on urea hydrolysis + biomass
python - <<'PY'
import os, glob
import pandas as pd
try:
    import cobra
    from cobra.io import read_sbml_model
except ImportError:
    raise SystemExit("cobra not installed")

mdir = "/data/data/Upcycling/research/additional/A6_metabolic/models"
out = "/data/data/Upcycling/research/additional/A6_metabolic"
rows = []
for xml in sorted(glob.glob(f"{mdir}/*.xml")):
    mag = os.path.basename(xml).replace(".xml","")
    try:
        model = read_sbml_model(xml)
    except Exception as e:
        print(f"[A6] {mag} load failed: {e}")
        continue
    # baseline biomass (complete media)
    try:
        model.solver = "glpk"
        # open all exchanges
        for ex in model.exchanges:
            ex.lower_bound = -10
        bio = model.slim_optimize()
    except Exception as e:
        bio = None
    # urea uptake & hydrolysis: try urea exchange
    urea_ids = [ex.id for ex in model.exchanges if "urea" in ex.id.lower() or ex.id.lower()=="ex_urea_e"]
    urea_flux = None
    if urea_ids:
        try:
            with model:
                model.reactions.get_by_id(urea_ids[0]).lower_bound = -5
                urea_flux = model.slim_optimize()
        except Exception:
            pass
    # urease reaction check
    urease_rxn = [r.id for r in model.reactions if "urea" in r.id.lower() or "URE" in r.name]
    rows.append({"MAG": mag,
                 "n_genes": len(model.genes),
                 "n_reactions": len(model.reactions),
                 "n_metabolites": len(model.metabolites),
                 "has_urease": int(len(urease_rxn)>0),
                 "urease_rxn_ids": ";".join(urease_rxn[:5]),
                 "biomass_flux_openmedia": bio,
                 "biomass_with_urea_uptake": urea_flux})
df = pd.DataFrame(rows)
df.to_csv(f"{out}/carveme_FBA_summary.csv", index=False)
print(df.to_string(index=False))
PY

echo "[A6] DONE"
