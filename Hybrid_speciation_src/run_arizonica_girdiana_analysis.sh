#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"




#taxa=group.arizgir.sample187.filter.groupinfo # group.arizgir.sample212.mix.groupinfo
taxa=group.arizgir.sample134.filter.groupinfo
name=arizgir.sample134.filter

perl -alne 'print "$F[0]\t$F[0]"' $taxa > $taxa.plinkid


plink --vcf VITVarB40-14_v2.0_hap1_filtered_final_pruned.vcf.gz \
  --make-bed \
  --double-id \
  --aec \
  --keep-allele-order \
  --set-missing-var-ids @:# \
  --keep $taxa.plinkid \
  --recode vcf-iid \
  --out $name.pruned




plink2 --vcf $name.pruned.vcf --keep $taxa.plinkid --aec --double-id --pca --out $name.pca
  ### downstream analysis in local folder


vcf=arizgir.sample187.filter.pruned.vcf
group=group.arizgir.sample187.filter.groupinfo.winpca
chrom=doaniana_sample94.chrom_sizes.tsv

winpca pca $vcf VITVarB40-14_v2.0.hap1.chr01:1-24898126 arizgir.sample187.filter.pruned.chr01
winpca chromplot arizgir.sample187.filter.pruned.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -m $group -g population -i 5
winpca chromplot arizgir.sample187.filter.pruned.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -p het -m $group -g population -i 5

perl -alne 'if ($F[1] =~ /x.arizonica.girdiana/  ) {print "$F[0]\t$F[0]"}else {print} ' $group > $group.individuals
winpca chromplot arizgir.sample187.filter.pruned.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -m $group.individuals -g population -i 5
winpca chromplot arizgir.sample187.filter.pruned.chr01 VITVarB40-14_v2.0.hap1.chr01:1-24898126 -p het -m $group.individuals -g population -i 5

  ### loop through all chromosomes
while IFS= read -r LINE ; do
  CHROM=$(echo "$LINE" | cut -f1)
  SIZE=$(echo "$LINE" | cut -f2)
  winpca pca $vcf ${CHROM}:1-${SIZE} ${CHROM} -t 5 
  winpca chromplot ${CHROM} ${CHROM}:1-${SIZE} -m $group.individuals -g population -i 5 #-c $COLORS
  winpca chromplot ${CHROM} ${CHROM}:1-${SIZE} -m $group.individuals -g population -i 5 -p het #-c $COLORS
done < <(cat doaniana_sample94.chrom_sizes.tsv)



### 3.1 use the non-pruned vcf file

rawvcf=VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz
chrom=vitis.chrom.sizes

taxa=group.arizgir.sample187.filter.groupinfo
perl -alne 'print "$F[0]\t$F[0]"' $taxa > $taxa.plinkid
echo $taxa.plinkid

name=arizgir.sample187.filter

plink --vcf $rawvcf \
  --make-bed \
  --double-id \
  --aec \
  --keep-allele-order \
  --set-missing-var-ids @:# \
  --mac 3 \
  --keep $taxa.plinkid \
  --recode vcf-iid \
  --out $name.raw

### 3.2 convert to EIGENSTRAT format

plink --vcf $name.raw.vcf --make-bed --aec --double-id --keep-allele-order --set-missing-var-ids @:# --recode 12 --out $name.raw.dates
  ##### use converted ped and map file, and the suffix should be map/ped.
perl -alne '$new=$F[3]*0.00000261; print"$F[0]\t$F[1]\t$new\t$F[3]"' $name.raw.dates.map | perl -p -ne 's/VITVarB40-14_v2.0.hap1.chr0?//' | perl -pe 's/VITVarB40-14_v2.0.hap1.//' >1 ; mv 1 $name.raw.dates.map
awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t"; for(i=7;i<=NF;++i) printf $i"\t";ind++;printf "\n"}' $name.raw.dates.ped > 2; mv 2 $name.raw.dates.ped

file=$name.raw.dates

echo "genotypename:    ${file}.ped" > par.PED.EIGENSTRAT
echo "snpname:         ${file}.map" >> par.PED.EIGENSTRAT
echo "indivname:       ${file}.ped" >> par.PED.EIGENSTRAT
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT
echo "genotypeoutname: ${file}.geno" >> par.PED.EIGENSTRAT
echo "snpoutname:      ${file}.snp" >> par.PED.EIGENSTRAT
echo "indivoutname:    ${file}.ind" >> par.PED.EIGENSTRAT
echo "familynames:     NO" >> par.PED.EIGENSTRAT

convertf -p par.PED.EIGENSTRAT

perl -i.bak -pe 's/^\s+//' $name.raw.dates.ind 
  ## remove leading spaces


perl -e '
        open my $fh, "$ARGV[0]" or die "Cannot open $ARGV[0]: $!";
        while (<$fh>) {
            chomp;
            my @cols = split/\s+/;
            $group_info{$cols[0]} = $cols[1];
        }
        close $fh;

    open my $fh2, "$ARGV[1]" or die "Cannot open $ARGV[1]: $!";
    while (<$fh2>) {
        chomp;
        my @F = split/\s+/;

         if (exists $group_info{$F[0]}) {
        print "$F[0]\t$F[1]\t$group_info{$F[0]}\n";
    } else {
        print STDERR "Warning: ID $F[0] not found in $group. Using original 3rd column: $F[2]";
        print "$F[0]\t$F[1]\t$F[2]\n";
    }
    }

' $taxa $name.raw.dates.ind >$name.raw.dates.ind.newgroup

perl -alne 'if ($F[2]=~/x.arizonica.girdiana/) {print "$F[0]\t$F[1]\t$F[0]"} else {print} ' \
 $name.raw.dates.ind.newgroup >$name.raw.dates.ind.newgroup.ind
  ### admix individuals level.
grep -v -e "arizonica" -e "girdiana" $name.raw.dates.ind.newgroup.ind |perl -pe 's/U\t//' | perl -alne '$name="output_$F[0]"; print"arizonica\tgirdiana\t$F[0]\t$name"' >admix.ind

export LD_LIBRARY_PATH=${BASE_DIR}/lib:$LD_LIBRARY_PATH
export PATH=${BASE_DIR}/bin:$PATH

echo 'DIR: ./
indivname: DIR/arizgir.sample187.filter.raw.dates.ind.newgroup.ind
snpname: DIR/arizgir.sample187.filter.raw.dates.snp
genotypename:   DIR/arizgir.sample187.filter.raw.dates.geno
admixlist: admix.ind
binsize:        0.001 #0.00005
maxdis:         1.0
seed:           777
runmode:        1
chithresh:      0.0
mincount:       1
zdipcorrmode:   YES
qbin:           50 # 50
runfit:  YES 
afffit:  YES 
lovalfit:  0.2 #0.1' > par.dates.arizgir.sample187.filter.raw

nohup dates -p par.dates.arizgir.sample187.filter.raw >20250609.arizgir.DATES.log &

for file in output_*/*.jout; do   if [ -f "$file" ]; then     echo -n "$file";     cat "$file";   fi; done

  ### results show the majority of the individuals are 10~15 generations ago.
  ### Consider using dataset from simulated F1 data and introgression-free samples to test again.






taxa=group.arizgir.sample134.filter.groupinfo
name=arizgir.sample134.filter

vcf=arizgir.sample187.filter.raw.vcf

ariz=id.arizonica.sample8
gir=id.girdiana.sample8

bcftools view -S ${ariz} ${vcf} -Ov -o sample.arizonica.sample8.vcf
bcftools view -S ${gir} ${vcf} -Ov -o sample.girdiana.sample8.vcf


perl simuF1.pl --vcf_a sample.arizonica.sample8.vcf --vcf_b sample.girdiana.sample8.vcf --out sample_simu_arizgirF1.vcf
vcf_index.sh sample_simu_arizgirF1.vcf

bcftools merge -O z -o arizgir.sample140.withF1.raw.vcf.gz arizgir.sample132.filter.raw.vcf.gz sample_simu_arizgirF1.vcf.gz

f1vcf=arizgir.sample140.withF1.raw.vcf.gz
name=arizgir.sample140.withF1
plink2 --vcf $f1vcf --aec --double-id --pca --out $name


plink --vcf $f1vcf --make-bed --aec --double-id --keep-allele-order --set-missing-var-ids @:# --recode 12 --out $name.raw.dates
  ##### use converted ped and map file, and the suffix should be map/ped.
perl -alne '$new=$F[3]*0.00000261; print"$F[0]\t$F[1]\t$new\t$F[3]"' $name.raw.dates.map | perl -p -ne 's/VITVarB40-14_v2.0.hap1.chr0?//' | perl -pe 's/VITVarB40-14_v2.0.hap1.//' >1 ; mv 1 $name.raw.dates.map
awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t"; for(i=7;i<=NF;++i) printf $i"\t";ind++;printf "\n"}' $name.raw.dates.ped > 2; mv 2 $name.raw.dates.ped

file=$name.raw.dates

echo "genotypename:    ${file}.ped" > par.PED.EIGENSTRAT.sample140.withF1
echo "snpname:         ${file}.map" >> par.PED.EIGENSTRAT.sample140.withF1
echo "indivname:       ${file}.ped" >> par.PED.EIGENSTRAT.sample140.withF1
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT.sample140.withF1
echo "genotypeoutname: ${file}.geno" >> par.PED.EIGENSTRAT.sample140.withF1
echo "snpoutname:      ${file}.snp" >> par.PED.EIGENSTRAT.sample140.withF1
echo "indivoutname:    ${file}.ind" >> par.PED.EIGENSTRAT.sample140.withF1
echo "familynames:     NO" >> par.PED.EIGENSTRAT.sample140.withF1

convertf -p par.PED.EIGENSTRAT.sample140.withF1

perl -i.bak -pe 's/^\s+//' $name.raw.dates.ind 
  ## remove leading spaces


perl -e '
        open my $fh, "$ARGV[0]" or die "Cannot open $ARGV[0]: $!";
        while (<$fh>) {
            chomp;
            my @cols = split/\s+/;
            $group_info{$cols[0]} = $cols[1];
        }
        close $fh;

    open my $fh2, "$ARGV[1]" or die "Cannot open $ARGV[1]: $!";
    while (<$fh2>) {
        chomp;
        my @F = split/\s+/;

         if (exists $group_info{$F[0]}) {
        print "$F[0]\t$F[1]\t$group_info{$F[0]}\n";
    } else {
        print STDERR "Warning: ID $F[0] not found in $group. Using original 3rd column: $F[2]";
        print "$F[0]\t$F[1]\t$F[2]\n";
    }
    }

' $taxa $name.raw.dates.ind >$name.raw.dates.ind.newgroup

perl -alne 'if ($F[2]=~/x.arizonica.girdiana/) {print "$F[0]\t$F[1]\t$F[0]"} else {print} ' \
 $name.raw.dates.ind.newgroup | perl -alne 'if ($F[2]=~/insilicoF1/) {print "$F[0]\t$F[1]\t$F[0]"} else {print} ' >$name.raw.dates.ind.newgroup.ind
  ### admix individuals level.
grep -v -e "arizonica" -e "girdiana" $name.raw.dates.ind.newgroup.ind |perl -pe 's/U\t//' | perl -alne '$name="output_$F[0]"; print"arizonica\tgirdiana\t$F[0]\t$name"' >admix.ind.withF1

export LD_LIBRARY_PATH=${BASE_DIR}/lib:$LD_LIBRARY_PATH
export PATH=${BASE_DIR}/bin:$PATH

echo 'DIR: ./
indivname: DIR/arizgir.sample140.withF1.raw.dates.ind.newgroup.ind
snpname: DIR/arizgir.sample140.withF1.raw.dates.snp
genotypename:   DIR/arizgir.sample140.withF1.raw.dates.geno
admixlist: admix.ind.withF1
binsize:        0.001 #0.00005
maxdis:         1.0
seed:           777
runmode:        1
chithresh:      0.0
mincount:       1
zdipcorrmode:   YES
qbin:           50 # 50
runfit:  YES 
afffit:  YES 
lovalfit:  0.2 #0.1' > par.dates.arizgir.sample140.withF1

nohup dates -p par.dates.arizgir.sample140.withF1 >20250609.arizgir.DATES.log &

for file in output_*/*.jout; do   if [ -f "$file" ]; then     echo -n "$file";     cat "$file";   fi; done


bcftools view -T <(bcftools query -f '%CHROM\t%POS\n' VITVarB40-14_v2.0_hap1_filtered_final_pruned.vcf.gz) \
  -Oz -o arizgir.sample140.withF1.pruned.vcf.gz arizgir.sample140.withF1.raw.vcf.gz

bcftools index arizgir.sample140.withF1.pruned.vcf.gz
