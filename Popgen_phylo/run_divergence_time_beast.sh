#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"

## 20250325


## 1. prepare the dataset
## 92 samples across 48 species.
raw=selections.txt
name=sample92spe48 
vcf=VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.vcf.gz

## 2. Data wrangling
cut -f 2 $raw |tail -n +2 > $name.id
perl -alne 'print "$F[1]\t$F[1]" ' $raw > $name.plink.id

## 3. Extract the 92 samples
bcftools view -S $name.id $vcf -Oz -o $name.vcf.gz &
tabix -p vcf $name.vcf.gz &

## vcf to fasta and nex files
vcf2phylip.py -i $name.vcf.gz -f -n -b --output-prefix $name


## After complex parameters adjustment by using BEAUti, the XML file is generated.
## run the BEAST2
${BASE_DIR}/beast -threads "${THREADS:-8}" Vitis_DivTime_v9_strict_iqtreeMLv2.xml

## run the treeannotator
${BASE_DIR}/treeannotator -burnin 10 -height mean Vitis_DivTime_v8.trees.txt Vitis_DivTime_v8.trees.annotated

${BASE_DIR}/treeannotator -burnin 10 -height mean Vitis_DivTime_v9.trees.txt Vitis_DivTime_v9.trees.annotated

