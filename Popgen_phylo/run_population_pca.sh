#!/usr/bin/env bash
set -euo pipefail

## Tianpeng Wang
## 2025-05-14
## This script is used to run PCA & MDS analysis on the Vitis dataset using PLINK2/

plink2 --vcf VITVarB40-14_v2.0_hap1_filtered_final_pruned.vcf.gz --make-bed --double-id --aec --pca --out vitis639.pruned

plink --bfile vitis639.pruned --aec --double-id --mds-plot 3 --cluster --out vitis639.pruned.mds


grep -v -e "ssp.sylvestris" -e "ssp.vinifera" -e "asia" -e "outgroup" -e "\bNA\b" group_vitis639.pca_group_info.txt >group_vitis487.NorthAme.pca_group_info.txt

perl -alne 'print "$F[2]\t$F[2]"' group_vitis487.NorthAme.pca_group_info.txt > group_vitis487.NorthAme.pca_group_info.plink.id

plink2 --bfile vitis639.pruned --keep group_vitis487.NorthAme.pca_group_info.plink.id --aec --double-id --pca --out group_vitis487.NorthAme.pca
plink --bfile vitis639.pruned --keep group_vitis487.NorthAme.pca_group_info.plink.id --aec --double-id --mds-plot 3 --cluster --out group_vitis487.NorthAme.mds
