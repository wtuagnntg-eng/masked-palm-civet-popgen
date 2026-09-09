#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/differentiation"

# BED contains transcript ID in column 4. Any physical overlap with a merged
# highly differentiated region is retained; no arbitrary coverage cutoff is used.
awk 'BEGIN{FS=OFS="\t"} $0!~/^#/ && ($3=="mRNA" || $3=="transcript") {
  id=""; n=split($9,a,";");
  for(i=1;i<=n;i++){if(a[i] ~ /^ID=/){sub(/^ID=/,"",a[i]); id=a[i]; break}}
  if(id!="") print $1,$4-1,$5,id,$3,$7
}' "$GFF" > "${OUT}/reference_transcripts.bed"

bedtools intersect -a "${OUT}/reference_transcripts.bed" \
  -b "${OUT}/highly_differentiated_regions.bed" -wa -u \
  > "${OUT}/transcripts_in_differentiated_regions.bed"
cut -f4 "${OUT}/transcripts_in_differentiated_regions.bed" | sort -u \
  > "${OUT}/candidate_transcripts.txt"
