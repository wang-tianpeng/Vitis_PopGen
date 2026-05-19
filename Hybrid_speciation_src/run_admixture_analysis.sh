#!/bin/sh 
set -euo pipefail

file=doaniana_sample97.mac2.vcf.gz
prefix=$(basename "$file" .vcf.gz)
#id_remove=nigra2_remove.id # if exists
sample94=ind.sample94.id

## use bcftools to extract the sample of interest
bcftools view -S $sample94 $file -Oz -o $prefix.sample94.vcf.gz &
tabix -p vcf $prefix.sample94.vcf.gz &


## plink format transformation for Faststructure
plink --vcf $prefix.sample94.vcf.gz --aec --set-missing-var-ids @_#  --double-id --out $prefix.sample94.faststr --make-bed
## plink formate transformation for admixture
plink --vcf $prefix.sample94.vcf.gz --recode 12 --aec --set-missing-var-ids @_#  --double-id --out $prefix.sample94.admix

for i in `seq 2 6`
do
	echo "admixture --cv $prefix.admix.ped $i >>$prefix.log${i}.txt"
done >cmd.admixture.$prefix.sh
