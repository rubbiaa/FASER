#!/usr/bin/env bash
# Loop over scifi_hits_all_*GeV.root files and run the plotting macro for each.
# Usage:
#   ./run_plots_for_all_energies.sh [input_dir] [pattern]
# Defaults:
#   input_dir = build
#   pattern   = scifi_hits_all_*GeV.root
# Notes:
# - Outputs are saved to build/<TAG>/ by the macro (TAG inferred from filename suffix),
#   unless you modify the macro call to pass an explicit output dir override.
# - Run this from the Display/ directory (where plot_muon_genfit_taubin.C is located).

set -u
INPUT_DIR=${1:-build}
PATTERN=${2:-scifi_hits_all_*GeV.root}

shopt -s nullglob
files=("${INPUT_DIR}"/${PATTERN})
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
  echo "No files matched: ${INPUT_DIR}/${PATTERN}"
  exit 0
fi

echo "Found ${#files[@]} files in ${INPUT_DIR}"
for f in "${files[@]}"; do
  echo "Processing ${f}"
  # Call the ROOT macro; it will infer the output directory (e.g., build/100GeV)
  root -l -b -q "plot_muon_genfit_taubin.C(\"${f}\")" || echo "ROOT macro failed for ${f}"
 done

echo "All done. Check build/<TAG>/ directories for outputs."
