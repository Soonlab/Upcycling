#!/bin/bash
# Download Pfam HMMs for urease + carbonic anhydrase
set -e
WD=/data/data/Upcycling/research/additional
HMM_DIR=$WD/hmms_shared
mkdir -p $HMM_DIR
cd $HMM_DIR

# Pfam IDs:
# PF00449 Urease_alpha
# PF00699 Urease_beta_gamma (Urease, subunit beta & gamma)
# PF01979 Amidohydro_1 (urease catalytic superfamily)
# PF00484 Pro_CA     (beta-class carbonic anhydrase, prokaryotic)
# PF00194 Carb_anhydrase (alpha-class CA)
# PF01682 Urease_acc (urease accessory)
# PF00988 CA_gamma   (gamma-class CA)

for pf in PF00449 PF00699 PF01979 PF00484 PF00194 PF01682 PF00988; do
    if [ ! -s ${pf}.hmm ]; then
        echo "Downloading $pf..."
        curl -sL -o ${pf}.hmm.gz "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${pf}?annotation=hmm"
        if [ -s ${pf}.hmm.gz ]; then
            gunzip -f ${pf}.hmm.gz
            echo "  $pf size: $(wc -c < ${pf}.hmm)"
        else
            echo "  $pf download failed"
            rm -f ${pf}.hmm.gz
        fi
    fi
done

# Concat into one
cat *.hmm > urease_CA.hmm 2>/dev/null
ls -la
echo "HMM download done"
