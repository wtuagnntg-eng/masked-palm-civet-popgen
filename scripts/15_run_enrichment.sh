#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/differentiation"

python "$(dirname "$0")/map_transcripts_to_genes.py" \
  --transcripts "${OUT}/candidate_transcripts.txt" \
  --mapping "$TRANSCRIPT_TO_GENE" \
  --output "${OUT}/candidate_genes.txt"

Rscript "$(dirname "$0")/15_enrichment.R" \
  "${OUT}/candidate_genes.txt" "$GENE2GO" "$GENE2KO" "$KO2PATHWAY" \
  "${OUT}/enrichment" 3
