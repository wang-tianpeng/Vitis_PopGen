#!/usr/bin/env bash
set -euo pipefail

RAW_VCF="${RAW_VCF:-vitis_sample451_phased_beagle.group20.vcf.gz}"
TAXA_FILE="${TAXA_FILE:-taxa.sample451.group20}"
SHARED_SWEEPS_BED="${SHARED_SWEEPS_BED:-shared_sweeps.bed}"
GENETIC_MAP="${GENETIC_MAP:-vitis_genetic_map.txt}"
R_SCRIPT="${R_SCRIPT:-perform_rdmc_analysis.R}"
OUTPUT_DIR="${OUTPUT_DIR:-rdmc_analysis}"
POP1="${POP1:-rupestris}"
POP2="${POP2:-riparia}"
MAX_SNPS="${MAX_SNPS:-100000}"
MIN_SNPS="${MIN_SNPS:-1000}"

mkdir -p "${OUTPUT_DIR}"

awk -v g="${POP1}" '$2 == g {print $1}' "${TAXA_FILE}" > "${OUTPUT_DIR}/${POP1}.samples"
awk -v g="${POP2}" '$2 == g {print $1}' "${TAXA_FILE}" > "${OUTPUT_DIR}/${POP2}.samples"
cat "${OUTPUT_DIR}/${POP1}.samples" "${OUTPUT_DIR}/${POP2}.samples" > "${OUTPUT_DIR}/both_pops.samples"

echo "region_id,pop1,pop2,n_snps,status" > "${OUTPUT_DIR}/rdmc_summary_results.csv"

while read -r sweep_chr sweep_start sweep_end _; do
    [[ -z "${sweep_chr}" || "${sweep_chr}" =~ ^# ]] && continue

    region_id="${sweep_chr}:${sweep_start}-${sweep_end}"
    sweep_len=$((sweep_end - sweep_start))
    extend_bp=$((sweep_len / 10))
    ext_start=$((sweep_start - extend_bp))
    ext_end=$((sweep_end + extend_bp))
    (( ext_start < 0 )) && ext_start=0
    ext_region="${sweep_chr}:${ext_start}-${ext_end}"

    temp_dir=$(mktemp -d -p "${OUTPUT_DIR}" "temp_${sweep_chr}_${sweep_start}_XXXX")
    trap 'rm -rf "${temp_dir}"' RETURN

    bcftools view -r "${ext_region}" \
        -S "${OUTPUT_DIR}/both_pops.samples" \
        --min-alleles 2 --max-alleles 2 \
        -Oz -o "${temp_dir}/region.vcf.gz" "${RAW_VCF}"
    bcftools index "${temp_dir}/region.vcf.gz"

    bcftools +fill-tags "${temp_dir}/region.vcf.gz" -Oz -o "${temp_dir}/pop1.af.vcf.gz" -- \
        -S "${OUTPUT_DIR}/${POP1}.samples" -t AF
    bcftools +fill-tags "${temp_dir}/region.vcf.gz" -Oz -o "${temp_dir}/pop2.af.vcf.gz" -- \
        -S "${OUTPUT_DIR}/${POP2}.samples" -t AF

    paste <(bcftools query -f '%CHROM\t%POS\t%AF\n' "${temp_dir}/pop1.af.vcf.gz") \
          <(bcftools query -f '%AF\n' "${temp_dir}/pop2.af.vcf.gz") \
        | awk -v n1="$(wc -l < "${OUTPUT_DIR}/${POP1}.samples")" \
              -v n2="$(wc -l < "${OUTPUT_DIR}/${POP2}.samples")" '
            BEGIN { OFS = "\t" }
            {
                af1 = $3
                af2 = $4
                total_freq = ((af1 * n1 * 2) + (af2 * n2 * 2)) / ((n1 + n2) * 2)
                if (total_freq > 0.05 && total_freq < 0.95) print $1, $2, af1, af2
            }' > "${temp_dir}/filtered_freqs.txt"

    n_snps=$(wc -l < "${temp_dir}/filtered_freqs.txt")
    if (( n_snps > MAX_SNPS )); then
        shuf "${temp_dir}/filtered_freqs.txt" | head -n "${MAX_SNPS}" | sort -k1,1V -k2,2n > "${temp_dir}/final_freqs.txt"
    else
        cp "${temp_dir}/filtered_freqs.txt" "${temp_dir}/final_freqs.txt"
    fi

    awk '
        BEGIN { OFS = "\t" }
        FNR == NR {
            if ($1 ~ /^#/) next
            chr = $1
            pos[chr, count[chr]] = $2
            cm[chr, count[chr]] = $3
            count[chr]++
            next
        }
        {
            chr = $1
            snp = $2
            if (!(chr in count)) { print "NA"; next }
            for (i = 0; i < count[chr] - 1; i++) {
                if (pos[chr, i] <= snp && pos[chr, i + 1] >= snp) {
                    p1 = pos[chr, i]; p2 = pos[chr, i + 1]
                    c1 = cm[chr, i]; c2 = cm[chr, i + 1]
                    print (p1 == p2) ? c1 : c1 + (c2 - c1) * (snp - p1) / (p2 - p1)
                    next
                }
            }
            print (snp < pos[chr, 0]) ? cm[chr, 0] : cm[chr, count[chr] - 1]
        }' "${GENETIC_MAP}" <(cut -f1,2 "${temp_dir}/final_freqs.txt") > "${temp_dir}/genetic_pos.txt"

    cut -f3 "${temp_dir}/final_freqs.txt" > "${temp_dir}/p1.txt"
    cut -f4 "${temp_dir}/final_freqs.txt" > "${temp_dir}/p2.txt"

    if (( n_snps < MIN_SNPS )); then
        echo "${region_id},${POP1},${POP2},${n_snps},too_few_snps" >> "${OUTPUT_DIR}/rdmc_summary_results.csv"
        rm -rf "${temp_dir}"
        continue
    fi

    Rscript "${R_SCRIPT}" \
        --p1_file "${temp_dir}/p1.txt" \
        --p2_file "${temp_dir}/p2.txt" \
        --pos_file "${temp_dir}/genetic_pos.txt" \
        --region_id "${region_id}" \
        --pop1 "${POP1}" \
        --pop2 "${POP2}" \
        --output_file "${OUTPUT_DIR}/rdmc_summary_results.csv"

    rm -rf "${temp_dir}"
done < "${SHARED_SWEEPS_BED}"

rm -f "${OUTPUT_DIR}/${POP1}.samples" "${OUTPUT_DIR}/${POP2}.samples" "${OUTPUT_DIR}/both_pops.samples"
