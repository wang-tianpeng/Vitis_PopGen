#!/usr/bin/env bash
set -euo pipefail



# ---# ---# ---# ---# ---# ---# ---# ---# ---
RAW_VCF="vitis_sample451_phased_beagle.group20.taxa.sample451.group8.mac3.vcf.gz"
RAW_VCF="vitis_sample451_phased_beagle.group20.vcf.gz"
# FILTERED_VCF="${RAW_VCF%.vcf.gz}.$(basename "${TAXA_FILE}").mac3.vcf.gz"

TAXA_FILE="taxa.sample451.group8.v2"

cut -f 1 $TAXA_FILE >taxa_keep_samples.fam.v2

#     "${TAXA_FILE}" > taxa_keep_samples.fam.v2

FILTERED_VCF="${RAW_VCF%.vcf.gz}.$(basename "${TAXA_FILE}").mac3.vcf.gz"
if [ ! -s "${FILTERED_VCF}" ]; then
    plink2 \
        --allow-extra-chr \
        --vcf "${RAW_VCF}" \
        --keep taxa_keep_samples.fam.v2 \
        --mac 3 \
        --export vcf bgz \
        --out plink_filtered

    mv plink_filtered.vcf.gz "${FILTERED_VCF}"
fi

# Index
if [ ! -s "${FILTERED_VCF}.tbi" ]; then
    if command -v bcftools >/dev/null 2>&1; then
        bcftools index -t -f "${FILTERED_VCF}"
    else
        tabix -f -p vcf "${FILTERED_VCF}"
    fi
fi
# ---# ---# ---# ---# ---# ---# ---# ---# ---# ---# ---# ---

SWEEPS_DIR="data_convergent_sweeps_quant"
RAW_VCF="vitis_sample451_phased_beagle.group20.taxa.sample451.group8.v2.mac3.vcf.gz"
for SHARED_SWEEPS_FILE_PATH in "${SWEEPS_DIR}"/*_convergent_regions_quant.tsv; do

    if [ ! -f "$SHARED_SWEEPS_FILE_PATH" ]; then
        echo "No files found matching pattern in ${SWEEPS_DIR}"
        continue
    fi

    SHARED_SWEEPS_FILE=$(basename "$SHARED_SWEEPS_FILE_PATH")
    echo "=========================================================="
    echo "Processing: $SHARED_SWEEPS_FILE"
    echo "=========================================================="

    NAME=${SHARED_SWEEPS_FILE%_convergent_regions_quant.tsv}


    OUTPUT_DIR="rdmc_quant/rdmc_analysis_advanced.${NAME}"
    mkdir -p "${OUTPUT_DIR}/sweep_freq_files"
    
    cut -f 3-5 "$SHARED_SWEEPS_FILE_PATH" | tail -n +2 |perl -alne '$F[1]=~s/\.\d+//; $F[2]=~s/\.\d+//; print "$F[0]\t$F[1]\t$F[2]"' >${OUTPUT_DIR}/shared_sweeps_$NAME.bed

    SHARED_SWEEPS_BED="${OUTPUT_DIR}/shared_sweeps_${NAME}.bed"

    # RAW_VCF="${FILTERED_VCF}"

    R_SCRIPT="perform_rdmc_analysis_advanced.R"

    mkdir -p "${OUTPUT_DIR}/sweep_freq_files"
    mkdir -p "${OUTPUT_DIR}/results"

    # POPS=($(cut -f 2 "${TAXA_FILE}" | sort | uniq))
    # done



    # # --- Main Loop ---
    # SWEEP_CHR=VITVarB40-14_v2.0.hap1.chr01
    # SWEEP_START=20325740
    # SWEEP_END=20445571
   POPS=($(cut -f 2 "${TAXA_FILE}" | sort | uniq))
    while read -r SWEEP_CHR SWEEP_START SWEEP_END; do
        [[ "$SWEEP_CHR" =~ ^# ]] || [ -z "$SWEEP_CHR" ] && continue
        
        REGION_ID="${SWEEP_CHR}:${SWEEP_START}-${SWEEP_END}"
        echo "--- Preparing frequency file for region: ${REGION_ID} ---"

        LENGTH=$(awk "BEGIN {print $SWEEP_END - $SWEEP_START}")
        EXTEND_BP=3000
        EXT_START=$(awk -v START="$SWEEP_START" -v EXTBP="$EXTEND_BP" 'BEGIN {print START - EXTBP}')
        # [ "${EXART}" -l] && EXTRT=0
        EXT_END=$(awk -v SWPEND="$SWEEP_END" -v EXTBP="$EXTEND_BP" 'BEGIN {print SWPEND + EXTBP}')
        EXT_REGION="${SWEEP_CHR}:${EXT_START}-${EXT_END}"

        SWEEP_FREQ_FILE="${OUTPUT_DIR}/sweep_freq_files/sweep_${SWEEP_CHR}_start${EXT_START}_end${EXT_END}.tsv"
        
        QUERY_STR="%CHROM\t%POS"
        PASTE_FILES=""
        for POP in "${POPS[@]}"; do
            bcftools view -r "${EXT_REGION}" -S "group_samples/${POP}.samples" "${RAW_VCF}" | \
                bcftools +fill-tags - -- -t AF > "tmp_af_${POP}.vcf"
            QUERY_STR="${QUERY_STR}\t%AF"
            PASTE_FILES="${PASTE_FILES} <(bcftools query -f '%AF\\n' tmp_af_${POP}.vcf)"
        done

        HEADER="chrom\tpos\t$(printf "%s\t" "${POPS[@]}")"
        echo -e "${HEADER%\\t}" > "${SWEEP_FREQ_FILE}"
        eval "paste <(bcftools query -f '%CHROM\\t%POS\\n' tmp_af_${POPS[0]}.vcf) ${PASTE_FILES}" >> "${SWEEP_FREQ_FILE}"
        
        rm tmp_af_*.vcf

    done < "${SHARED_SWEEPS_BED}" 

    ls ${OUTPUT_DIR}/sweep_freq_files/*.tsv |xargs -i perl -i -pe "s/\t$//" {}

done
