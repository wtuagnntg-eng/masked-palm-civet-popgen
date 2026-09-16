#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../config/config.sh"

OUT="${RESULTS}/rnaseq"

mkdir -p \
  "$OUT/bam" \
  "$OUT/counts"

########################################
# 1. Build HISAT2 index
########################################

IDX="${PROJECT_DIR}/reference/civet_hisat2"

if [[ ! -f "${IDX}.1.ht2" && ! -f "${IDX}.1.ht2l" ]]; then
    hisat2-build \
      -p "$THREADS" \
      "$REF" \
      "$IDX"
fi

########################################
# 2. RNA-seq mapping
########################################

bams=()

while IFS=$'\t' read -r sample tissue r1 r2; do

    [[ "$sample" == "sample_id" ]] && continue
    [[ -z "$sample" ]] && continue

    echo "===== Mapping ${sample} (${tissue}) ====="

    hisat2 \
      -p "$THREADS" \
      --dta \
      -x "$IDX" \
      -1 "$r1" \
      -2 "$r2" | \
    samtools sort \
      -@ "$THREADS" \
      -o "${OUT}/bam/${sample}.sorted.bam" \
      -

    samtools index \
      -@ "$THREADS" \
      "${OUT}/bam/${sample}.sorted.bam"

    bams+=(
      "${OUT}/bam/${sample}.sorted.bam"
    )

done < "$SAMPLES_RNA"

########################################
# 3. Check SAF annotation
########################################

SAF="${OUT}/gene_exons.saf"

if [[ ! -s "$SAF" ]]; then
    echo "ERROR: gene-level SAF file not found:"
    echo "       $SAF"
    echo
    echo "Run 17_build_gene_annotation.py first."
    exit 1
fi

########################################
# 4. Gene-level read counting
########################################

echo "===== featureCounts gene-level counting ====="

featureCounts \
  -T "$THREADS" \
  -p \
  -F SAF \
  -a "$SAF" \
  -o "${OUT}/counts/gene_counts.txt" \
  "${bams[@]}"

########################################
# 5. Summary
########################################

echo
echo "RNA-seq mapping and gene-level counting completed."
echo
echo "BAM directory:"
echo "  ${OUT}/bam"
echo
echo "Gene-level count matrix:"
echo "  ${OUT}/counts/gene_counts.txt"
