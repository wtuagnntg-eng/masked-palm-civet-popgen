#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/diversity"
mkdir -p "$OUT"

plink --vcf "$RELAXED_VCF" --double-id --allow-extra-chr --het --out "${OUT}/heterozygosity"
python "$(dirname "$0")/calculate_observed_heterozygosity.py" \
  --het "${OUT}/heterozygosity.het" --metadata "$SAMPLES_WGS" \
  --output "${OUT}/observed_heterozygosity.tsv"
Rscript "$(dirname "$0")/wilcoxon_individual_metrics.R" \
  "${OUT}/observed_heterozygosity.tsv" Ho "${OUT}/Ho_wilcoxon.tsv"
