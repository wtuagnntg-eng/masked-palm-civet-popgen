#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../config/config.sh"

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${RESULTS}/rnaseq"

mkdir -p "$OUT"

########################################
# 1. Build gene-level exon annotation
########################################

echo "===== Build gene-level SAF and effective exon lengths ====="

python "${HERE}/17_build_gene_annotation.py" \
  --gff "$GFF" \
  --saf "${OUT}/gene_exons.saf" \
  --lengths "${OUT}/gene_lengths.tsv" \
  --gene-map "${OUT}/gene_map.tsv"

########################################
# 2. RNA-seq mapping and gene-level counting
########################################

echo "===== RNA-seq mapping and featureCounts ====="

bash "${HERE}/16_rnaseq_mapping_counts.sh"

########################################
# 3. Calculate gene-level TPM
########################################

echo "===== Calculate gene-level TPM ====="

python "${HERE}/18_calculate_tpm.py" \
  --counts "${OUT}/counts/gene_counts.txt" \
  --lengths "${OUT}/gene_lengths.tsv" \
  --sample-metadata "$SAMPLES_RNA" \
  --output-all "${OUT}/gene_TPM_all.tsv" \
  --output-filtered "${OUT}/gene_TPM.tsv" \
  --output-counts "${OUT}/gene_counts_filtered.tsv"

########################################
# 4. Transcriptomic descriptive analyses
########################################

echo "===== Transcriptomic descriptive analyses ====="

Rscript "${HERE}/20_transcriptome_descriptive_plots.R" \
  "${OUT}/gene_TPM.tsv" \
  "$IMMUNE_PANEL" \
  "${OUT}/gene_map.tsv" \
  "$OUT"

########################################
# 5. Summary
########################################

echo
echo "Transcriptomic workflow completed."
echo
echo "Main outputs:"
echo "  Gene-level SAF:"
echo "    ${OUT}/gene_exons.saf"
echo
echo "  Effective gene exon lengths:"
echo "    ${OUT}/gene_lengths.tsv"
echo
echo "  Gene ID/name mapping:"
echo "    ${OUT}/gene_map.tsv"
echo
echo "  Gene-level counts:"
echo "    ${OUT}/counts/gene_counts.txt"
echo
echo "  TPM before low-expression filtering:"
echo "    ${OUT}/gene_TPM_all.tsv"
echo
echo "  TPM retained for downstream analyses:"
echo "    ${OUT}/gene_TPM.tsv"
echo
