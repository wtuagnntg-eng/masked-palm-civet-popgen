#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
steps=(
  01_prepare_reference.sh
  02_wgs_qc_mapping.sh
  02b_mapping_qc.sh
  03_gatk_joint_calling.sh
  04_filter_variants.sh
  05_population_structure.sh
  07_phylogeny_iqtree.sh
  08_relatedness_king.sh
  09_diversity_heterozygosity.sh
  10_nucleotide_diversity.sh
  11_ld_decay.sh
  12_roh.sh
  13_genomic_differentiation.sh
  14_overlap_annotation.sh
  15_run_enrichment.sh
  16_rnaseq_mapping_counts.sh
  19_transcriptome_expression.sh
)
for s in "${steps[@]}"; do
  echo "===== ${s} ====="
  bash "${HERE}/scripts/${s}"
done
