#!/usr/bin/env bash
set -euo pipefail

# 2025-06-02


vcf=doaniana_sample94.pruned.mac2.vcf.gz

acerip=id.ind5.acerip
mus=id.ind5.mus

bcftools view -S ${acerip} ${vcf} -Ov -o acerip_samples.vcf
bcftools view -S ${mus} ${vcf} -Ov -o mus_samples.vcf

perl simuF1.pl --vcf_a acerip_samples.vcf --vcf_b mus_samples.vcf --out simu_acerip_musF1.vcf
vcf_index.sh simu_acerip_musF1.vcf

bcftools merge -O z -o simu_acerip_musF1_concat.vcf.gz doaniana_sample94.pruned.mac2.vcf.gz simu_acerip_musF1.vcf.gz    

group=group.doaniana_sample94.f1.groupinfo
### chr1
winpca pca simu_acerip_musF1_concat.vcf.gz VITVarB40-14_v2.0.hap1.chr01:1-24898126 simu_acerip_musF1_concat.chr01
winpca chromplot simu_acerip_musF1_concat.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -m $group -g population -i 5
winpca chromplot simu_acerip_musF1_concat.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -p het -m $group -g population -i 5
### chr2
winpca pca simu_acerip_musF1_concat.vcf.gz VITVarB40-14_v2.0.hap1.chr02:1-21551818 simu_acerip_musF1_concat.chr02
winpca chromplot simu_acerip_musF1_concat.chr02 VITVarB40-14_v2.0.hap1.chr02:1-21551818 -m $group -g population -i 5
winpca chromplot simu_acerip_musF1_concat.chr02 VITVarB40-14_v2.0.hap1.chr02:1-21551818 -p het -m $group -g population -i 5



COLORS='mustangensis:44AA99,acerrip:0077BB,doaniana:EE7733,insilicoF1:888888'

while IFS= read -r LINE ; do
  CHROM=$(echo "$LINE" | cut -f1)
  SIZE=$(echo "$LINE" | cut -f2)
  winpca chromplot ${CHROM} ${CHROM}:1-${SIZE} -m group.doaniana_sample94.f1.groupinfo -g population -i 5 -c $COLORS
  winpca chromplot ${CHROM} ${CHROM}:1-${SIZE} -m group.doaniana_sample94.f1.groupinfo -g population -i 5 -p het -c $COLORS
done < <(cat doaniana_sample94.chrom_sizes.tsv)

cut -f 1 doaniana_sample94.chrom_sizes.tsv


CHROMS='VITVarB40-14_v2.0.hap1.chr01,VITVarB40-14_v2.0.hap1.chr02,VITVarB40-14_v2.0.hap1.chr03,VITVarB40-14_v2.0.hap1.chr04,VITVarB40-14_v2.0.hap1.chr05,VITVarB40-14_v2.0.hap1.chr06,VITVarB40-14_v2.0.hap1.chr07,VITVarB40-14_v2.0.hap1.chr08,VITVarB40-14_v2.0.hap1.chr09,VITVarB40-14_v2.0.hap1.chr10,VITVarB40-14_v2.0.hap1.chr11,VITVarB40-14_v2.0.hap1.chr12,VITVarB40-14_v2.0.hap1.chr13,VITVarB40-14_v2.0.hap1.chr14,VITVarB40-14_v2.0.hap1.chr15,VITVarB40-14_v2.0.hap1.chr16,VITVarB40-14_v2.0.hap1.chr17,VITVarB40-14_v2.0.hap1.chr18,VITVarB40-14_v2.0.hap1.chr19'

winpca genomeplot ./ $CHROMS -m ../group.doaniana_sample94.f1.groupinfo -g population -i 20 -f HTML -c $COLORS

winpca flip VITVarB40-14_v2.0.hap1.chr05 -w 18300000-24900000
winpca flip VITVarB40-14_v2.0.hap1.chr12 -w 11100000-26100000


# 20250608
