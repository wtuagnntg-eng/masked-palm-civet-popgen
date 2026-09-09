#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

RAW="${RESULTS}/variants/raw/joint.raw.vcf.gz"
OUT="${RESULTS}/variants"
mkdir -p "$OUT"

# Build list of retained autosomes from the rename table.
awk 'BEGIN{FS="\t"} $0!~/^#/ {print $1}' "$AUTOSOME_RENAME" > "${OUT}/autosomes.list"

AUTOSOMES=$(paste -sd, "${OUT}/autosomes.list")
bcftools view -r "$AUTOSOMES" "$RAW" -Oz -o "${OUT}/autosomes.raw.vcf.gz"
bcftools index -t "${OUT}/autosomes.raw.vcf.gz"

# Strict dataset: common variants / QC dataset described in Section 2.3.
vcftools --gzvcf "${OUT}/autosomes.raw.vcf.gz" \
  --remove-indels --min-alleles 2 --max-alleles 2 \
  --min-meanDP 5 --max-meanDP 200 --max-missing 0.8 \
  --maf 0.05 --hwe 0.001 \
  --recode --recode-INFO-all --stdout | bgzip -c > "${OUT}/strict.autosomes.vcf.gz"
bcftools index -t "${OUT}/strict.autosomes.vcf.gz"
bcftools annotate --rename-chrs "$AUTOSOME_RENAME" "${OUT}/strict.autosomes.vcf.gz" -Oz -o "$STRICT_VCF"
bcftools index -t "$STRICT_VCF"

# Less restrictive dataset: diversity / ROH / LD / differentiation.
vcftools --gzvcf "${OUT}/autosomes.raw.vcf.gz" \
  --remove-indels --min-alleles 2 --max-alleles 2 \
  --min-meanDP 5 --max-meanDP 200 --max-missing 0.8 \
  --recode --recode-INFO-all --stdout | bgzip -c > "${OUT}/relaxed.autosomes.vcf.gz"
bcftools index -t "${OUT}/relaxed.autosomes.vcf.gz"
bcftools annotate --rename-chrs "$AUTOSOME_RENAME" "${OUT}/relaxed.autosomes.vcf.gz" -Oz -o "$RELAXED_VCF"
bcftools index -t "$RELAXED_VCF"

# QC summaries: compare these values with those reported in the manuscript.
bcftools view -H "$STRICT_VCF" | wc -l > "${OUT}/strict.snp_count.txt"
bcftools stats "$STRICT_VCF" > "${OUT}/strict.bcftools.stats.txt"
bcftools view -H "$RELAXED_VCF" | wc -l > "${OUT}/relaxed.snp_count.txt"
bcftools stats "$RELAXED_VCF" > "${OUT}/relaxed.bcftools.stats.txt"
