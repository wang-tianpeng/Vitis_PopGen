#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"

# ${BASE_DIR}/10_pi_fst_dxy

rawgeno=VITVarB40-14_v2.0_hap1_filtered_final.geno.gz
group=taxa.sample451.group20

cut -f 2 $group |sort |uniq -c

popgenWindows.py -g $rawgeno -w 100000 -f phased -T 20 \
    --popsFile $group \
    -o sample351_pi_fst_dxy_100k.csv \
    -p acerifolia -p aestivalis -p amurensis -p arizonica -p berlandieri -p californica \
    -p cinerea -p girdiana -p monticola -p mustangensis -p labrusca -p riparia -p rupestris \
    -p sylvestris -p vinifera -p vulpina -p x_champinii -p x_champiniib -p x_doaniana


cut -f 2 $group |sort |uniq -c


pixy --stats pi fst dxy \
    --vcf VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz \
    --populations taxa.sample451.group20 \
    --window_size 10000 \
    --n_cores 20 \
    --output_folder pixy_output_10k \
    --output_prefix sample451_group20_10k
