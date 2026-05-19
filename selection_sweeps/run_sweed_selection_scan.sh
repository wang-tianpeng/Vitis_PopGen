#!/usr/bin/env bash
set -euo pipefail

RAW_VCF="${RAW_VCF:-vitis_sample451_phased_beagle.group20.vcf.gz}"
TAXA_FILE="${TAXA_FILE:-taxa.sample451.group20}"
TARGET_GROUPS="${TARGET_GROUPS:-$(cut -f2 "${TAXA_FILE}" | sort -u)}"
GRID_FILE="${GRID_FILE:-chrom_grid10k.txt}"

if ! command -v SweeD-P >/dev/null 2>&1; then
    echo "SweeD-P was not found in PATH." >&2
    exit 1
fi

if [[ ! -f "${GRID_FILE}" ]]; then
    bcftools view -h "${RAW_VCF}" \
        | grep "contig" \
        | perl -alne 'if (/ID=(.*?),length=(\d+)>/) { $grid = int($2 / 10000); print "$1\t$grid" }' \
        > "${GRID_FILE}"
fi

for target_group in ${TARGET_GROUPS}; do
    workdir="sweed_${target_group}"
    mkdir -p "${workdir}/sample_lists" "${workdir}/vcf_split" "${workdir}/sweed_out"

    awk -v g="${target_group}" '$2 == g {print $1}' "${TAXA_FILE}" > "${workdir}/sample_lists/${target_group}.samples"

    for chr_num in {1..19}; do
        region=$(printf "VITVarB40-14_v2.0.hap1.chr%02d" "${chr_num}")
        grid_points=$(awk -v r="${region}" '$1 == r {print $2}' "${GRID_FILE}")
        [[ -n "${grid_points}" ]] || {
            echo "Skipping ${region}: no grid size found." >&2
            continue
        }

        vcf_in="${workdir}/vcf_split/${region}.${target_group}.vcf"
        bcftools view -r "${region}" \
            -S "${workdir}/sample_lists/${target_group}.samples" \
            -Ov -o "${vcf_in}" \
            "${RAW_VCF}"

        SweeD-P -name "${target_group}.${region}" \
            -input "${vcf_in}" \
            -grid "${grid_points}" \
            -folded

        mv "SweeD_Report.${target_group}.${region}" "${workdir}/sweed_out/"
        mv "SweeD_Info.${target_group}.${region}" "${workdir}/sweed_out/"
    done
done
