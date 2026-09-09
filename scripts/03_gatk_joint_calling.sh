#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

BAMDIR="${RESULTS}/mapping"
GVCFDIR="${RESULTS}/gvcf"
VCFDIR="${RESULTS}/variants/raw"
mkdir -p "$GVCFDIR" "$VCFDIR"

# Per-sample GVCF calling.
while IFS=$'\t' read -r sample ecotype r1 r2; do
    [[ "$sample" == "sample_id" ]] && continue
    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
      -R "$REF" \
      -I "${BAMDIR}/${sample}.dedup.bam" \
      -O "${GVCFDIR}/${sample}.g.vcf.gz" \
      --emit-ref-confidence GVCF \
      --min-base-quality-score 30
done < "$SAMPLES_WGS"

# Combine and jointly genotype chromosome by chromosome.
mapfile -t CHROMS < <(cut -f1 "$REF.fai")
for chr in "${CHROMS[@]}"; do
    args=()
    while IFS=$'\t' read -r sample ecotype r1 r2; do
        [[ "$sample" == "sample_id" ]] && continue
        args+=(--variant "${GVCFDIR}/${sample}.g.vcf.gz")
    done < "$SAMPLES_WGS"

    gatk --java-options "-Xmx${JAVA_MEM}" CombineGVCFs \
      -R "$REF" \
      -L "$chr" \
      "${args[@]}" \
      -O "${GVCFDIR}/${chr}.combined.g.vcf.gz"

    gatk --java-options "-Xmx${JAVA_MEM}" GenotypeGVCFs \
      -R "$REF" \
      -V "${GVCFDIR}/${chr}.combined.g.vcf.gz" \
      -L "$chr" \
      -O "${VCFDIR}/${chr}.vcf.gz"

done

# Concatenate chromosome VCFs in the exact reference-FASTA order.
VCF_LIST="${VCFDIR}/chromosome_vcfs.list"
: > "$VCF_LIST"
for chr in "${CHROMS[@]}"; do
    printf '%s\n' "${VCFDIR}/${chr}.vcf.gz" >> "$VCF_LIST"
done

bcftools concat \
  -f "$VCF_LIST" \
  -Oz \
  -o "${VCFDIR}/joint.raw.vcf.gz"

bcftools index -t "${VCFDIR}/joint.raw.vcf.gz"
