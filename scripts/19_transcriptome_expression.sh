#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/rnaseq"
mkdir -p "$OUT"

python "$(dirname "$0")/17_build_transcript_lengths.py" --gff "$GFF" \
  --output "${OUT}/transcript_lengths.tsv"
python "$(dirname "$0")/18_calculate_tpm.py" \
  --counts "${OUT}/counts/transcript_counts.txt" \
  --lengths "${OUT}/transcript_lengths.tsv" \
  --transcript-to-gene "$TRANSCRIPT_TO_GENE" \
  --sample-metadata "$SAMPLES_RNA" \
  --output-transcript "${OUT}/transcript_TPM.tsv" \
  --output-gene "${OUT}/gene_TPM.tsv"
Rscript "$(dirname "$0")/20_transcriptome_descriptive_plots.R" \
  "${OUT}/gene_TPM.tsv" "$IMMUNE_PANEL" "$OUT"
