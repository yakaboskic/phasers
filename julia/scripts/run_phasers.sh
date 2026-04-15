#!/bin/bash
# Interactive single-merged-file phasers continuous runner.
# Prompts for the merged-data file and a list of similarity thresholds, then
# runs phasers continuous for each threshold. Outputs land in julia/results/.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
JULIA_DIR="$(dirname "$SCRIPT_DIR")"
RESULTS_DIR="$JULIA_DIR/results"
mkdir -p "$RESULTS_DIR"

read -p "Enter merged data file: " MERGED_FILE
read -p "Enter similarity threshold(s) (e.g., 0.6 0.7 0.8): " -a THRESHOLDS

BASENAME=$(basename "$MERGED_FILE" .tsv.gz)

for THRESHOLD in "${THRESHOLDS[@]}"; do
    echo "Running threshold: $THRESHOLD"
    TAG="${THRESHOLD/./p}"
    OUTPUT_FILE="${RESULTS_DIR}/rs_curve_${BASENAME}_${TAG}.csv"
    LOG_FILE="${RESULTS_DIR}/rs_run_${BASENAME}_${TAG}.log"

    uv run phasers continuous \
        --merged-data-file "$MERGED_FILE" \
        --score-column combined \
        --source-name pigean \
        --similarity-matrix-threshold "$THRESHOLD" \
        --output-file "$OUTPUT_FILE" 2>&1 | tee "$LOG_FILE"
done
