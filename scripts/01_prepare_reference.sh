#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

mkdir -p "${PROJECT_DIR}/reference"

bwa-mem2 index "$REF"
samtools faidx "$REF"
gatk CreateSequenceDictionary -R "$REF" -O "${REF%.*}.dict"
