#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/diversity/pi"
mkdir -p "$OUT"

awk -F'\t' 'NR>1 && $2=="Farm" {print $1}' "$SAMPLES_WGS" > "${OUT}/farm.samples"
awk -F'\t' 'NR>1 && $2=="Wild" {print $1}' "$SAMPLES_WGS" > "${OUT}/wild.samples"

# RELAXED_VCF has already been filtered at >=80% call rate across all 20 samples.
# No additional population-specific missingness filter is applied for pi.
for pop in farm wild; do
  keep="${OUT}/${pop}.samples"
  vcftools --gzvcf "$RELAXED_VCF" --keep "$keep" \
    --window-pi 50000 --window-pi-step 10000 --out "${OUT}/${pop}.50kb_10kb"
  vcftools --gzvcf "$RELAXED_VCF" --keep "$keep" \
    --window-pi 5000000 --window-pi-step 5000000 --out "${OUT}/${pop}.5Mb"
done

Rscript "$(dirname "$0")/pi_block_bootstrap.R" \
  "${OUT}/wild.5Mb.windowed.pi" "${OUT}/farm.5Mb.windowed.pi" \
  10000 "${OUT}/pi_block_bootstrap_summary.tsv" "${OUT}/PI_5Mb_blocks.tsv"
