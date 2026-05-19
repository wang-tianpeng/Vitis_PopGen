#!/usr/bin/env bash
set -euo pipefail

GROUPS="${GROUPS:-arizonica berlandieri cinerea girdiana monticola mustangensis riparia rupestris}"

for target_group in ${GROUPS}; do
    workdir="sweed_simulations_${target_group}"
    [[ -d "${workdir}" ]] || {
        echo "Skipping missing directory: ${workdir}" >&2
        continue
    }

    output="${workdir}/max_clr_per_chromosome_${target_group}.txt"
    echo -e "Species\tChromosome\tMaxCLR" > "${output}"

    for chr_num in {1..19}; do
        region=$(printf "VITVarB40-14_v2.0.hap1.chr%02d" "${chr_num}")
        report_dir="${workdir}/sweed_simulations_${target_group}_${region}/sweed_out"
        [[ -d "${report_dir}" ]] || continue

        max_clr=$(awk 'NR > 1 && $0 !~ /^\\/\\// {if ($3 > max) max = $3} END {print max + 0}' \
            "${report_dir}"/SweeD_Report.sim_"${target_group}"_"${region}"_* 2>/dev/null)
        echo -e "${target_group}\t${region}\t${max_clr}" >> "${output}"
    done
done
