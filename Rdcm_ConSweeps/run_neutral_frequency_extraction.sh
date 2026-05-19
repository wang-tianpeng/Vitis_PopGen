#!/usr/bin/env bash
set -euo pipefail



RAW_VCF="vitis_sample451_phased_beagle.group20.vcf.gz"
TAXA_FILE="taxa.sample451.group8"

#     "${TAXA_FILE}" > taxa_keep_samples.fam
cut -f 1 "${TAXA_FILE}" | awk '{print $1}' > taxa_keep_samples.fam

FILTERED_VCF="${RAW_VCF%.vcf.gz}.$(basename "${TAXA_FILE}").mac3.vcf.gz"
if [ ! -s "${FILTERED_VCF}" ]; then
    plink2 \
        --allow-extra-chr \
        --vcf "${RAW_VCF}" \
        --keep "taxa_keep_samples.fam" \
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
RAW_VCF="${FILTERED_VCF}"



RAW_VCF="vitis_sample451_phased_beagle.group20.taxa.sample451.group8.v2.mac3.vcf.gz"

TAXA_FILE="taxa.sample451.group8.v2"
OUTPUT_FILE="neutral_freqs_group8.v2.tsv"
N_SITES=100000
readarray -t POPS < <(cut -f 2 "${TAXA_FILE}" | sort -u)

for POP in "${POPS[@]}"; do
    awk -v g="${POP}" '$2==g {print $1}' "${TAXA_FILE}" > "group_samples/${POP}.samples"
done


echo "Extracting random sites and calculating frequencies..."
bcftools view "${RAW_VCF}" | \
    bcftools query -f '%CHROM\t%POS\n' | \
    shuf -n ${N_SITES} | \
    sort -k1,1V -k2,2n > random_sites.bed

perl -alne '$start=$F[1]-1; print"$F[0]\t$start\t$F[1]"' random_sites.bed > random_sites.bed.2

PASTE_FILES=""

POP="arizonica"

for POP in "${POPS[@]}"; do
    echo "Processing ${POP}..."
    bcftools view -R random_sites.bed.2 -S "group_samples/${POP}.samples" "${RAW_VCF}" | \
        bcftools +fill-tags - -- -t AF > "af_${POP}.vcf"
    
    QUERY_STR="${QUERY_STR}\t%AF"
    PASTE_FILES="${PASTE_FILES} <(bcftools query -f '%AF\\n' af_${POP}.vcf)"
done

echo -e "${HEADER%\\t}" > "${OUTPUT_FILE}"

eval "paste <(bcftools query -f '%CHROM\\t%POS\\n' af_${POPS[0]}.vcf) ${PASTE_FILES}" >> "${OUTPUT_FILE}"

rm *.samples af_*.vcf 
#random_sites.bed
echo "Step complete."
