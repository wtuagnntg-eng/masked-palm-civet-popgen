#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/population_structure"

# Recover the LD-pruned SNP VCF used by PCA and ADMIXTURE.
plink --bfile "${OUT}/structure.pruned" --allow-extra-chr --recode vcf bgz \
  --out "${OUT}/structure.pruned"

python "$(dirname "$0")/06_vcf_to_iupac_fasta.py" \
  -i "${OUT}/structure.pruned.vcf.gz" \
  -o "${OUT}/structure.pruned.iupac.fa"

# SNP-only alignment: ModelFinder + ascertainment-bias correction; 1,000 UFBoot and SH-aLRT.
iqtree3 -s "${OUT}/structure.pruned.iupac.fa" \
  -m MFP+ASC -B 1000 --alrt 1000 -T AUTO \
  --prefix "${OUT}/civet_tree"

# No external outgroup is used. Midpoint rooting, if desired, is for visualization only.
