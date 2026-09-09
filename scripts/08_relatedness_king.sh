#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT="${RESULTS}/population_structure"

king -b "${OUT}/structure.pruned.bed" --kinship --prefix "${OUT}/KING"
