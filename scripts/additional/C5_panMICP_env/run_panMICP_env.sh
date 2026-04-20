#!/usr/bin/env bash
# C5 Pan-MICP environment comparison
# Curated MICP-confirmed type strains from soil/marine/karst/freshwater + 6 hero MAGs
# → skani ANI matrix + summary
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate align_env

OUT="/data/data/Upcycling/research/additional/C5_panMICP_env"
REFS="$OUT/refs"
LOG="$OUT/logs"
mkdir -p "$REFS" "$LOG"

# Curated MICP-relevant type strain accessions
# (urease-positive ureolytic strains from diverse environments)
cat > "$OUT/curated_refs.tsv" <<'EOF'
GCF_000189475.1	Sporosarcina_pasteurii	soil_canonical
GCF_000293575.1	Sporosarcina_psychrophila	cold_soil
GCF_000299475.1	Sporosarcina_ureae	soil
GCF_000009045.1	Bacillus_subtilis_168	soil
GCF_002209305.2	Bacillus_megaterium	soil
GCF_001548095.1	Bacillus_licheniformis	soil
GCF_001399775.1	Lysinibacillus_sphaericus	soil
GCF_000691605.1	Lysinibacillus_fusiformis	soil
GCF_900119685.1	Halomonas_pacifica	marine
GCF_000196175.1	Halomonas_elongata	marine_saline
GCF_002288455.1	Idiomarina_loihiensis	marine_deep
GCF_000613485.1	Pseudoalteromonas_haloplanktis	marine
GCF_001077795.1	Helicobacter_pylori_26695	host_acidic
GCF_000196215.1	Klebsiella_aerogenes	clinical_ureolytic
GCF_000648555.1	Proteus_mirabilis	clinical_ureolytic
GCF_000009305.1	Yersinia_enterocolitica	clinical
GCF_000284615.1	Sphingobacterium_sp_21	soil_sphingo
GCF_000620325.1	Sphingobacterium_multivorum	soil_sphingo
GCF_001719105.1	Sphingobacterium_paramultivorum	soil_sphingo
GCF_001043025.1	Pseudomonas_helleri	soil_pseudomonas
EOF

# Download with NCBI datasets CLI
echo "[$(date)] downloading curated reference genomes"
cd "$REFS"
while IFS=$'\t' read -r acc name env; do
    if [[ -f "${name}.fna" ]]; then
        echo "[have] $name"; continue
    fi
    echo "[get ] $acc -> $name"
    datasets download genome accession "$acc" --include genome --filename "${name}.zip" \
        > "$LOG/${name}_download.log" 2>&1 || { echo "[fail] $acc"; continue; }
    unzip -o "${name}.zip" -d "${name}_unzip" > /dev/null 2>&1 || true
    fna=$(find "${name}_unzip" -name "*.fna" | head -1)
    if [[ -n "$fna" ]]; then
        cp "$fna" "${name}.fna"
    fi
    rm -rf "${name}.zip" "${name}_unzip"
done < "$OUT/curated_refs.tsv"

# Stage hero genomes alongside
HEROES=(S13 S16 S23 C22 M1 S26)
for h in "${HEROES[@]}"; do
    if [[ ! -f "$REFS/HERO_${h}.fna" ]]; then
        zcat "/data/data/Upcycling/MAGs_FASTA_files/${h}.fasta.gz" > "$REFS/HERO_${h}.fna"
    fi
done

# Run skani all-vs-all
echo "[$(date)] running skani"
ls "$REFS"/*.fna > "$OUT/genome_list.txt"
# Use dram_env (has skani)
conda deactivate
conda activate dram_env
skani triangle -l "$OUT/genome_list.txt" --full-matrix \
    -o "$OUT/skani_panMICP.matrix" 2> "$LOG/skani.log"
skani dist -q "$REFS"/HERO_*.fna -r "$REFS"/*.fna --min-af 30 \
    -o "$OUT/skani_hero_vs_refs.tsv" 2>> "$LOG/skani.log"

echo "[$(date)] done C5"
