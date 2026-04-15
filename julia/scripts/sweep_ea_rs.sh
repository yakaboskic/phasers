#!/bin/bash
# Sweep similarity thresholds for ontology and semantic merged files,
# capturing EA_RS / RS_full / MLRS metrics for comparison.
#
# Anchor: ontology @ 0.8
# Semantic sweep covers a wide band so we can find the threshold whose
# EA_RS matches the ontology anchor.
#
# Outputs land in julia/results/sweep_logs/

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
JULIA_DIR="$(dirname "$SCRIPT_DIR")"
RESULTS_DIR="$JULIA_DIR/results"
LOG_DIR="$RESULTS_DIR/sweep_logs"
mkdir -p "$LOG_DIR"

ONTOLOGY_FILE="$JULIA_DIR/merged_data.tsv.gz"
SEMANTIC_FILE="$JULIA_DIR/merged_cosine_miniLM.tsv.gz"

ONTOLOGY_THRESHOLD="0.80"
SEMANTIC_THRESHOLDS=(0.30 0.40 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85)

run_phasers() {
    local label=$1
    local file=$2
    local threshold=$3
    local tag="${threshold/./p}"   # 0.80 -> 0p80
    local out="${LOG_DIR}/rs_curve_${label}_${tag}.csv"
    local log="${LOG_DIR}/rs_run_${label}_${tag}.log"

    echo ">>> ${label} @ ${threshold}"
    uv run phasers continuous \
        --merged-data-file "$file" \
        --score-column combined \
        --source-name pigean \
        --similarity-matrix-threshold "$threshold" \
        --output-file "$out" 2>&1 | tee "$log"
}

run_phasers ontology "$ONTOLOGY_FILE" "$ONTOLOGY_THRESHOLD"
for t in "${SEMANTIC_THRESHOLDS[@]}"; do
    run_phasers semantic "$SEMANTIC_FILE" "$t"
done

echo
echo "Sweep complete. Logs and curves are in $LOG_DIR"
