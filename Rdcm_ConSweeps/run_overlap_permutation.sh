#!/usr/bin/env bash
set -euo pipefail

INTRO_DIR="intro_w100s10"
SWEEPS_DIR="shared_quant92"
GENOME_FILE="vitis.chrom.sizes"
PLOT_SCRIPT="plot_permutation_results.R"
PERMUTATIONS=1000

mkdir -p "$OUTPUT_DIR"

if [ ! -f "$PLOT_SCRIPT" ]; then
    echo "Error: R script $PLOT_SCRIPT not found."
    exit 1
fi

for INTRO_FILE in "${INTRO_DIR}"/*_merged_filtered.bed; do
    [ -e "$INTRO_FILE" ] || continue
    
    BASENAME=$(basename "$INTRO_FILE")
    echo "----------------------------------------------------------------"
    echo "Processing: $BASENAME"

    POP_A=$(echo "$BASENAME" | cut -d'_' -f 2)
    POP_B=$(echo "$BASENAME" | cut -d'_' -f 3)

    echo "  -> Identified pair: $POP_A vs $POP_B"

    SWEEP_FILE_1="${SWEEPS_DIR}/shared_sweeps_${POP_A}_vs_${POP_B}.bed"
    SWEEP_FILE_2="${SWEEPS_DIR}/shared_sweeps_${POP_B}_vs_${POP_A}.bed"
    
    TARGET_SWEEP=""

    if [ -f "$SWEEP_FILE_1" ]; then
        TARGET_SWEEP="$SWEEP_FILE_1"
    elif [ -f "$SWEEP_FILE_2" ]; then
        TARGET_SWEEP="$SWEEP_FILE_2"
    else
        echo "  [WARNING] No matching shared sweeps file found for ${POP_A} and ${POP_B}."
        echo "  Checked: $(basename "$SWEEP_FILE_1") and $(basename "$SWEEP_FILE_2")"
        continue
    fi

    echo "  -> Found sweeps file: $(basename "$TARGET_SWEEP")"

    PERM_RESULT_TXT="${OUTPUT_DIR}/perm_dist_${PAIR_NAME}.txt"
    PERM_PLOT_PNG="${OUTPUT_DIR}/perm_plot_${PAIR_NAME}.png"

    echo "  -> Calculating observed overlap..."
    OBSERVED_OVERLAP=$(bedtools intersect -a "$INTRO_FILE" -b "$TARGET_SWEEP" -wo | awk '{s+=$NF} END {print s+0}')
    
    if [ -z "$OBSERVED_OVERLAP" ]; then OBSERVED_OVERLAP=0; fi
    echo "     Observed: $OBSERVED_OVERLAP bp"

    echo "  -> Running $PERMUTATIONS permutations..."
    > "$PERM_RESULT_TXT"
    for i in $(seq 1 "$PERMUTATIONS"); do
            bedtools intersect -a - -b "$TARGET_SWEEP" -wo | \
            awk '{s+=$NF} END {print s+0}' >> "$PERM_RESULT_TXT"
    done

    sed -i 's/^$/0/' "$PERM_RESULT_TXT"

    echo "  -> Generating plot..."
    Rscript "$PLOT_SCRIPT" \
        --results_file "$PERM_RESULT_TXT" \
        --observed_value "$OBSERVED_OVERLAP" \
        --output_plot "$PERM_PLOT_PNG" \
        --title "Intro: ${POP_A}-${POP_B} vs Sweeps"

    echo "  -> Done. Results saved to $OUTPUT_DIR"
done

echo "================================================================"
echo "All pairs processed."
