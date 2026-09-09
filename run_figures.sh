#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "${HERE}/config/config.sh"
OUT="${RESULTS}/figures"; mkdir -p "$OUT"
Rscript "${HERE}/scripts/21_plot_figure1.R" "${RESULTS}/mapping_qc/mapping_qc.tsv" "$OUT"
Rscript "${HERE}/scripts/22_plot_figure2.R" "${RESULTS}/population_structure" "$SAMPLES_WGS" "$OUT"
Rscript "${HERE}/scripts/23_plot_figure3.R" "$RESULTS" "$OUT"
Rscript "${HERE}/scripts/24_plot_figure4.R" "${RESULTS}/differentiation" "$OUT"
# Figure 5 and Figure S2 are produced by 20_transcriptome_descriptive_plots.R.
