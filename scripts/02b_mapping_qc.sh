#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/mapping_qc"
BAMDIR="${RESULTS}/mapping"
mkdir -p "$OUT"
echo -e 'sample_id\tecotype\tmean_depth\tmapping_rate\tproperly_paired_rate\tcov_1x\tcov_5x\tcov_10x' > "${OUT}/mapping_qc.tsv"
while IFS=$'\t' read -r sample ecotype r1 r2; do
  [[ "$sample" == "sample_id" ]] && continue
  bam="${BAMDIR}/${sample}.dedup.bam"
  flag="${OUT}/${sample}.flagstat.txt"
  samtools flagstat -@ "$THREADS" "$bam" > "$flag"
  mapping=$(awk '/ mapped \(/ {gsub(/[()%]/,"",$5); print $5; exit}' "$flag")
  proper=$(awk '/ properly paired \(/ {gsub(/[()%]/,"",$6); print $6; exit}' "$flag")
  read mean c1 c5 c10 < <(samtools depth -aa -@ "$THREADS" "$bam" | awk '{n++;s+=$3;if($3>=1)a++;if($3>=5)b++;if($3>=10)c++} END{if(n) printf "%.8f %.8f %.8f %.8f\n",s/n,a/n,b/n,c/n; else print "NA NA NA NA"}')
  echo -e "${sample}\t${ecotype}\t${mean}\t${mapping}\t${proper}\t${c1}\t${c5}\t${c10}" >> "${OUT}/mapping_qc.tsv"
done < "$SAMPLES_WGS"
