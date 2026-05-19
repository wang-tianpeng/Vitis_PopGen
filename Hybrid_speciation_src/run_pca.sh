#!/usr/bin/env bash
set -euo pipefail


plink --vcf doaniana_sample97.mac2.sample94.vcf.gz --aec --double-id --pca --out doaniana.sample94.mac


plink --vcf VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.sample97.vcf --keep doaniana.sample94.mac2.nosex --aec --double-id --pca --out doaniana.sample94.pruned
