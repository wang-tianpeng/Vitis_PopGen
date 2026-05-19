#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"



group=20250225_SamplesSelection.txt
name=doaniana_sample97

cut -f 2,4 $group > $name.groupinfo
cut -f 2 $group > $name.id
perl -alne 'print "$F[1]\t$F[1]" ' $group > $name.plink.id


vcf=VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz


bcftools view -S $name.id $vcf -Oz -o $name.vcf.gz &
tabix -p vcf $name.vcf.gz &
#wait

## 2.1 filter by mac=3 /// 2025-03-02
plink --aec --double-id --set-missing-var-ids @__# \
    --vcf $name.vcf.gz --make-bed --mac 2 --recode vcf-iid --out $name.mac2






nohup pixy --stats pi fst dxy \
    --vcf $name.mac2.vcf.gz \
    --populations $name.groupinfo \
    --window_size 20000 \
    --n_cores 5 \
    --bypass_invariant_check 'yes' \
    --output_prefix $name.mac2.pixy.output &

python ${BASE_DIR}/parseVCF.py -i $name.mac2.vcf.gz -o $name.mac2.geno.gz &

python ${BASE_DIR}/popgenWindows.py -g $name.geno.gz \
    --popsFile $name.groupinfo \
    -f phased \
    -w 20000 \
    -s 10000 \
    -T 10 \
    -p doaniana -p mustangensis -p acerifolia -p outgroup \
    -o $name.popgenWindows.output

while read -r sample group; do
    if [[ $group != "group" ]]; then
        echo "$sample    $sample" >> sample.${group}.plink.id 
    fi
done < $name.groupinfo

for group_file in sample.*.plink.id; do
    group_name=$(basename $group_file .id)
    plink --vcf $name.mac2.vcf.gz \
      --keep $group_file \
      --double-id \
      --set-missing-var-ids @:# \
      --aec \
      --het \
      --out ${group_name}
done
