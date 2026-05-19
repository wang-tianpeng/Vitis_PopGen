#!/usr/bin/env bash
set -euo pipefail

RAW_VCF="vitis_sample451_phased_beagle.group20.vcf.gz"
TAXA_FILE="taxa.sample451.group20"
OUTPUT_FILE="neutral_freqs.tsv"
N_SITES=100000
echo "Step 1: Setting up sample lists..."
mkdir -p "${SAMPLE_DIR}"

for POP in "${POPS[@]}"; do
    awk -v g="${POP}" '$2==g {print $1}' "${TAXA_FILE}" > "${SAMPLE_DIR}/${POP}.samples"
done

echo "Step 2: Extracting random sites..."
    shuf -n ${N_SITES} | \
    sort -k1,1V -k2,2n | perl -alne 'print "$F[0]\t$F[1]\t$F[1]"' > random_sites.bed

echo "Step 3: Calculating allele frequencies for each population..."
#TMP_VCF_DIR=$(mktemp -d -p ".")
for POP in "${POPS[@]}"; do
    echo "  -> Processing ${POP}..."
        bcftools annotate -x INFO/AF - | \
        bcftools +fill-tags - -- -t AF > "af_${POP}.vcf"
    

    PASTE_FILES="${PASTE_FILES} <(bcftools query -f '%AF\\n' af_${POP}.vcf)"
done

echo "Step 4: Assembling the final wide-format frequency file..."
echo -e "${HEADER%\\t}" > "${OUTPUT_FILE}"


echo "Step 5: Cleaning up temporary files..."
rm -r "${TMP_VCF_DIR}" random_sites.bed
echo "Step complete."
echo "Sample lists are stored in '${SAMPLE_DIR}/' for downstream use."


echo "Step 6: Filtering out invariant sites..."
awk '
BEGIN {FS="\t"; OFS="\t"}
NR==1 {print; next}
{
    has_value = 0
    for(i=3; i<=NF; i++) {
        if($i > 0) {
            has_value = 1
            break
        }
    }
    if(has_value == 1) print
}
' "${OUTPUT_FILE}" > "${OUTPUT_FILE%.tsv}.filter.tsv"
