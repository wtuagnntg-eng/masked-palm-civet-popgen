#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/differentiation"
mkdir -p "$OUT"

awk -F'\t' 'NR>1 && $2=="Farm" {print $1}' "$SAMPLES_WGS" > "${OUT}/farm.samples"
awk -F'\t' 'NR>1 && $2=="Wild" {print $1}' "$SAMPLES_WGS" > "${OUT}/wild.samples"

# Use the same less-restrictive autosomal dataset directly. No additional
# population-specific missingness filter is applied for pi or FST.
vcftools --gzvcf "$RELAXED_VCF" --keep "${OUT}/wild.samples" \
  --window-pi 50000 --window-pi-step 10000 --out "${OUT}/wild.50kb_10kb"
vcftools --gzvcf "$RELAXED_VCF" --keep "${OUT}/farm.samples" \
  --window-pi 50000 --window-pi-step 10000 --out "${OUT}/farm.50kb_10kb"
vcftools --gzvcf "$RELAXED_VCF" \
  --weir-fst-pop "${OUT}/farm.samples" --weir-fst-pop "${OUT}/wild.samples" \
  --fst-window-size 50000 --fst-window-step 10000 --out "${OUT}/farm_vs_wild"

python "$(dirname "$0")/identify_differentiated_regions.py" \
  --wild-pi "${OUT}/wild.50kb_10kb.windowed.pi" \
  --farm-pi "${OUT}/farm.50kb_10kb.windowed.pi" \
  --fst "${OUT}/farm_vs_wild.windowed.weir.fst" \
  --min-snps 10 --quantile 0.95 \
  --windows-out "${OUT}/matched_windows.tsv" \
  --candidate-windows-out "${OUT}/highly_differentiated_windows.bed" \
  --regions-out "${OUT}/highly_differentiated_regions.bed" \
  --thresholds-out "${OUT}/thresholds.tsv"
