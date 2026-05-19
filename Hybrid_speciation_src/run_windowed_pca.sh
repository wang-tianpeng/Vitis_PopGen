#!/bin/env bash
set -euo pipefail

# 2025-05-09
# Note: modified config.py file in modules and set VCF_PASS_FILTER = False

vcf=doaniana_sample94.pruned.mac2.vcf.gz
group=doaniana_sample94.groupinfo
chrom=doaniana_sample94.chrom_sizes.tsv

bcftools query -h $vcf | grep -v "^#" | perl -alne '$_=~/ID=(.*?),length=(\d+?)\>/' > $chrom

# VITVarB40-14_v2.0.hap1.chr01:1-24898126
winpca pca $vcf VITVarB40-14_v2.0.hap1.chr01:1-24898126 VITVarB40-14_v2.0.hap1.chr01
winpca chromplot VITVarB40-14_v2.0.hap1.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -m $group -g population -i 5
winpca chromplot VITVarB40-14_v2.0.hap1.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -p het -m $group -g population -i 5
winpca flip VITVarB40-14_v2.0.hap1.chr12 -r
winpca flip VITVarB40-14_v2.0.hap1.chr12 -w 11100000-26200000

COLORS='mustangensis:0A9396,acerifolia:BB3E03,riparia:DDCC77,doaniana:117733'
COLORS='mustangensis:0A9396,acerifolia:BB3E03,riparia:BB3E03,doaniana:DDCC77'
winpca flip VITVarB40-14_v2.0.hap1.chr12 -w 11100000-26200000

winpca chromplot VITVarB40-14_v2.0.hap1.chr12 VITVarB40-14_v2.0.hap1.chr01:1-26925356 -m $group -g population -i 5 -c $COLORS


COLORS='mustangensis:44AA99,acerrip:0077BB,doaniana:EE7733,insilicoF1:888888'

while IFS= read -r LINE ; do
  CHROM=$(echo "$LINE" | cut -f1)
  SIZE=$(echo "$LINE" | cut -f2)
  winpca pca doaniana_sample94.pruned.mac2.vcf.gz ${CHROM}:1-${SIZE} ${CHROM} -t 5 
  winpca chromplot ${CHROM} ${CHROM}:1-${SIZE} -m ../group.doaniana_sample94.f1.groupinfo -g population -i 5 -c $COLORS
  winpca chromplot ${CHROM} ${CHROM}:1-${SIZE} -m ../group.doaniana_sample94.f1.groupinfo -g population -i 5 -p het -c $COLORS
done < <(cat ../doaniana_sample94.chrom_sizes.tsv)

cut -f 1 doaniana_sample94.chrom_sizes.tsv


CHROMS='VITVarB40-14_v2.0.hap1.chr01,VITVarB40-14_v2.0.hap1.chr02,VITVarB40-14_v2.0.hap1.chr03,VITVarB40-14_v2.0.hap1.chr04,VITVarB40-14_v2.0.hap1.chr05,VITVarB40-14_v2.0.hap1.chr06,VITVarB40-14_v2.0.hap1.chr07,VITVarB40-14_v2.0.hap1.chr08,VITVarB40-14_v2.0.hap1.chr09,VITVarB40-14_v2.0.hap1.chr10,VITVarB40-14_v2.0.hap1.chr11,VITVarB40-14_v2.0.hap1.chr12,VITVarB40-14_v2.0.hap1.chr13,VITVarB40-14_v2.0.hap1.chr14,VITVarB40-14_v2.0.hap1.chr15,VITVarB40-14_v2.0.hap1.chr16,VITVarB40-14_v2.0.hap1.chr17,VITVarB40-14_v2.0.hap1.chr18,VITVarB40-14_v2.0.hap1.chr19'

winpca genomeplot ./ $CHROMS -m group.doaniana_sample94.f1.groupinfo -g population -i 20 -f HTML -c $COLORS
