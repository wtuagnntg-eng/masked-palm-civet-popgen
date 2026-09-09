#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/ld"
mkdir -p "$OUT"

awk -F'\t' 'NR>1 && $2=="Farm" {print $1}' "$SAMPLES_WGS" > "${OUT}/farm.samples"
awk -F'\t' 'NR>1 && $2=="Wild" {print $1}' "$SAMPLES_WGS" > "${OUT}/wild.samples"

# For LD only, MAF and call rate are recalculated independently in each population.
for pop in farm wild; do
  vcftools --gzvcf "$RELAXED_VCF" --keep "${OUT}/${pop}.samples" \
    --maf 0.05 --max-missing 0.8 \
    --recode --recode-INFO-all --stdout | bgzip -c > "${OUT}/${pop}.ld_input.vcf.gz"
  tabix -f -p vcf "${OUT}/${pop}.ld_input.vcf.gz"
  PopLDdecay -InVCF "${OUT}/${pop}.ld_input.vcf.gz" \
    -OutStat "${OUT}/${pop}.LDdecay.stat.gz" -MaxDist 300
done

Rscript "$(dirname "$0")/summarize_ld.R" \
  "${OUT}/wild.LDdecay.stat.gz" "${OUT}/farm.LDdecay.stat.gz" \
  "${OUT}/LD_0_30kb_summary.tsv" "${OUT}/LD_decay_0_30kb.tsv"
