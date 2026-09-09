#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/roh"
mkdir -p "$OUT"

plink --vcf "$RELAXED_VCF" --double-id --allow-extra-chr --make-bed --out "${OUT}/roh_input"
plink --bfile "${OUT}/roh_input" --allow-extra-chr \
  --homozyg \
  --homozyg-window-snp 50 \
  --homozyg-snp 50 \
  --homozyg-kb 100 \
  --homozyg-density 50 \
  --homozyg-gap 1000 \
  --homozyg-window-het 1 \
  --homozyg-window-missing 5 \
  --homozyg-window-threshold 0.05 \
  --out "${OUT}/civet_roh"

python "$(dirname "$0")/summarize_roh.py" \
  --hom "${OUT}/civet_roh.hom" \
  --metadata "$SAMPLES_WGS" \
  --fai "$REF.fai" \
  --rename "$AUTOSOME_RENAME" \
  --output "${OUT}/roh_summary.tsv"

Rscript "$(dirname "$0")/wilcoxon_individual_metrics.R" \
  "${OUT}/roh_summary.tsv" FROH "${OUT}/FROH_wilcoxon.tsv"
