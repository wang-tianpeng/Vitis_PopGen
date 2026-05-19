#!/usr/bin/env bash
set -euo pipefail

TARGET_GROUP="rupestris" # The population you are analyzing
TARGET_CHR="VITVarB40-14_v2.0.hap1.chr01" # The chromosome you are simulating
NUM_SIMULATIONS=100 # Number of neutral simulations to run

SAMPLE_SIZE=$(awk -v g="${TARGET_GROUP}" '$2==g' taxa.sample451.group20 | wc -l)

SEQ_LENGTH=$(bcftools view -h doaniana_sample94.mac2.vcf.gz | grep "##contig=<ID=${TARGET_CHR}" | perl -ne '/length=(\d+)/ && print $1')

GRID_POINTS=$(awk -v r="${TARGET_CHR}" '$1==r {print $2}' chrom_grid10k.txt)

NE=20000                  # Effective Population Size (Example value)
MUT_RATE=5.4e-9           # Mutation rate for Vitis (from literature, example)
REC_RATE=6.2e-8           # Recombination rate for Vitis (from literature, example)

SIM_WORKDIR="sweed_simulations_${TARGET_GROUP}/sweed_simulations_${TARGET_GROUP}_${TARGET_CHR}"
MAX_CLR_LIST="${SIM_WORKDIR}/max_clr_values.txt"
mkdir -p "${SIM_WORKDIR}/vcf"
mkdir -p "${SIM_WORKDIR}/sweed_out"
rm -f "${MAX_CLR_LIST}" # Clear previous results
touch "${MAX_CLR_LIST}"

echo "--- Starting Neutral Simulation Workflow ---"
echo "Target Group: ${TARGET_GROUP}"
echo "Target Chromosome: ${TARGET_CHR}"
echo "Sample Size: ${SAMPLE_SIZE} individuals"
echo "Sequence Length: ${SEQ_LENGTH} bp"
echo "Grid Points for SweeD: ${GRID_POINTS}"
echo "-------------------------------------------"


for i in $(seq 1 ${NUM_SIMULATIONS}); do
    echo "*** Running Simulation ${i} of ${NUM_SIMULATIONS} ***"
    SIM_VCF="${SIM_WORKDIR}/vcf/sim_${i}.vcf"
    
    python run_msprime_simulation.py \
        --sample-size "${SAMPLE_SIZE}" \
        --seq-length "${SEQ_LENGTH}" \
        --Ne "${NE}" \
        --mut-rate "${MUT_RATE}" \
        --rec-rate "${REC_RATE}" \
        --output-vcf "${SIM_VCF}" \
        --seed ${i} # Use loop counter as a reproducible seed

    SweeD-P -name "sim_${i}" \
          -input "${SIM_VCF}" \
          -grid "${GRID_POINTS}" \
          -folded
    
    # NR>2 skips the header lines. $2 is the Likelihood column.
    MAX_CLR=$(awk 'NR>2 {print $2}' "SweeD_Report.sim_${i}" | sort -nr | head -n 1)
    echo "${MAX_CLR}" >> "${MAX_CLR_LIST}"
    
    echo "Simulation ${i}: Max CLR = ${MAX_CLR}"

    mv "SweeD_Report.sim_${i}" "${SIM_WORKDIR}/sweed_out/"
    mv "SweeD_Info.sim_${i}" "${SIM_WORKDIR}/sweed_out/"
done

echo "--- Simulation loop finished ---"

P95_LINE=$(echo "0.95 * ${NUM_SIMULATIONS}" | bc | xargs printf "%.0f")

THRESHOLD=$(sort -n "${MAX_CLR_LIST}" | sed -n "${P95_LINE}p")

echo "Step complete."
echo "----------------------------------------------------------------"
echo "The 95% significance threshold for CLR is: ${THRESHOLD}"
echo "----------------------------------------------------------------"
echo "A list of all maximum CLR values is saved in: ${MAX_CLR_LIST}"




# // ...existing code...
echo "A list of all maximum CLR values is saved in: ${MAX_CLR_LIST}"

#### after running this script, collect all the max CLR values from the simulations

echo "--- Consolidating simulation results for all chromosomes ---"

TARGET_GROUP="arizonica" # Change this to the desired group name
# BASE_DIR="sweed_simulations_${TARGET_GROUP}"
BASE_DIR="."
OUTPUT_FILE="${BASE_DIR}/simulations_all_chromosomes_max_${TARGET_GROUP}.txt"

echo -e "MaxCLR\tSpecies\tChromosome" > "${OUTPUT_FILE}"

for chr_dir in $(find "${BASE_DIR}" -mindepth 1 -maxdepth 1 -type d); do
    
    CHR_NAME=$(basename "${chr_dir}" | sed "s/sweed_simulations_${TARGET_GROUP}_//")
    
    MAX_CLR_LIST_FILE="${chr_dir}/max_clr_values.txt"
    
    if [ -f "${MAX_CLR_LIST_FILE}" ]; then
        echo "Processing: ${CHR_NAME}"
        awk -v species="${TARGET_GROUP}" -v chromosome="${CHR_NAME}" '{print $1"\t"species"\t"chromosome}' "${MAX_CLR_LIST_FILE}" >> "${OUTPUT_FILE}"
    else
        echo "Warning: Could not find results file for ${CHR_NAME} at ${MAX_CLR_LIST_FILE}"
    fi
done

echo "Step complete."
echo "All maximum CLR values from all simulations are saved in: ${OUTPUT_FILE}"
