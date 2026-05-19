#!/bin/bash
set -euo pipefail

TEMPLATE="temp.rdmc.slurm"
SEARCH_DIR="data_convergent_sweeps_quant92"

if [ ! -f "$TEMPLATE" ]; then
    echo "Error: Template file $TEMPLATE not found in current directory."
    exit 1
fi

for FILE in "$SEARCH_DIR"/*_convergent_regions_quant92.tsv; do
    if [ -f "$FILE" ]; then
        BASENAME=$(basename "$FILE")
        
        PREFIX=${BASENAME%_convergent_regions_quant92.tsv}
        
        POP1=$(echo "$PREFIX" | awk -F'_vs_' '{print $1}')
        POP2=$(echo "$PREFIX" | awk -F'_vs_' '{print $2}')
        
        NEW_SLURM="run.slurm.${POP1}_${POP2}.rdmc"
        
        echo "Creating $NEW_SLURM for $POP1 vs $POP2..."
        
        cp "$TEMPLATE" "$NEW_SLURM"
        
        
        sed -i "s/#SBATCH -J rdmcRup/#SBATCH -J rdmc_${POP1}_${POP2}/" "$NEW_SLURM"
        
        sed -i "s/SHARED_SWEEPS_FILE=\"rupestris_vs_riparia_convergent_regions_quant92.tsv\"/SHARED_SWEEPS_FILE=\"${BASENAME}\"/" "$NEW_SLURM"
        
        sed -i "s/POP1=\"rupestris\"/POP1=\"${POP1}\"/" "$NEW_SLURM"
        
        sed -i "s/POP2=\"riparia\"/POP2=\"${POP2}\"/" "$NEW_SLURM"
        
    fi
done

echo "Done generating SLURM scripts."
