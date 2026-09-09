#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/rnaseq"
mkdir -p "$OUT/bam" "$OUT/counts"

# Build HISAT2 index once if it does not already exist.
IDX="${PROJECT_DIR}/reference/civet_hisat2"
if [[ ! -f "${IDX}.1.ht2" && ! -f "${IDX}.1.ht2l" ]]; then
  hisat2-build -p "$THREADS" "$REF" "$IDX"
fi

bams=()
while IFS=$'\t' read -r sample tissue r1 r2; do
  [[ "$sample" == "sample_id" ]] && continue
  hisat2 -p "$THREADS" --dta -x "$IDX" -1 "$r1" -2 "$r2" | \
    samtools sort -@ "$THREADS" -o "${OUT}/bam/${sample}.sorted.bam" -
  samtools index -@ "$THREADS" "${OUT}/bam/${sample}.sorted.bam"
  bams+=("${OUT}/bam/${sample}.sorted.bam")
done < "$SAMPLES_RNA"

# Paired-end counting of exon features; Parent is used to assign exons to transcripts.
featureCounts -T "$THREADS" -p -t exon -g Parent \
  -a "$GFF" -o "${OUT}/counts/transcript_counts.txt" "${bams[@]}"
