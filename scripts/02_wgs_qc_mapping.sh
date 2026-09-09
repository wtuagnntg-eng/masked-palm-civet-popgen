#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

OUT="${RESULTS}/mapping"
CLEAN="${RESULTS}/clean_reads"
mkdir -p "$OUT" "$CLEAN"

while IFS=$'\t' read -r sample ecotype r1 r2; do
    [[ "$sample" == "sample_id" ]] && continue

    fastp \
      -i "$r1" -I "$r2" \
      -o "${CLEAN}/${sample}.R1.clean.fastq.gz" \
      -O "${CLEAN}/${sample}.R2.clean.fastq.gz" \
      -w "$THREADS" \
      -h "${CLEAN}/${sample}.fastp.html" \
      -j "${CLEAN}/${sample}.fastp.json"

    bwa-mem2 mem -t "$THREADS" "$REF" \
      "${CLEAN}/${sample}.R1.clean.fastq.gz" \
      "${CLEAN}/${sample}.R2.clean.fastq.gz" | \
      samtools sort -@ "$THREADS" -o "${OUT}/${sample}.sorted.bam" -

    picard AddOrReplaceReadGroups \
      I="${OUT}/${sample}.sorted.bam" \
      O="${OUT}/${sample}.rg.bam" \
      RGID="$sample" \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM="$sample"

    picard MarkDuplicates \
      I="${OUT}/${sample}.rg.bam" \
      O="${OUT}/${sample}.dedup.bam" \
      M="${OUT}/${sample}.dedup.metrics.txt" \
      REMOVE_DUPLICATES=true

    samtools index -@ "$THREADS" "${OUT}/${sample}.dedup.bam"
    samtools flagstat -@ "$THREADS" "${OUT}/${sample}.dedup.bam" \
      > "${OUT}/${sample}.flagstat.txt"
done < "$SAMPLES_WGS"
