#!/bin/sh 
set -euo pipefail

file=VITVarB40-14_v2.0_hap1_filtered.4distes.vcf
#id_remove=nigra2_remove.id # if exists


prefix=${file/.vcf/}

## plink format transformation for Faststructure
plink --vcf $file --aec --set-missing-var-ids @_#  --double-id --out $prefix.faststr --make-bed
## plink formate transformation for admixture
plink --vcf $file --recode 12 --aec --set-missing-var-ids @_#  --double-id --out $prefix.admix

for i in `seq 2 18`
do
	echo "admixture --cv $prefix.admix.ped $i >>log${i}.txt"
done >cmd.admixture.$prefix.sh

for i in `seq 2 18`
do
	echo "structure.py -K $i --input=$prefix.faststr --output=$prefix.faststr.$i"
done >cmd.faststructure.$prefix.sh
