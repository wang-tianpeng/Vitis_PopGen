#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"
# 2025-05-29




## ${BASE_DIR}/example

export LD_LIBRARY_PATH=${BASE_DIR}/lib:$LD_LIBRARY_PATH
export PATH=${BASE_DIR}/bin:$PATH

DIR: ./data
indivname: DIR/family_packed.ind
snpname: DIR/family_packed.snp
genotypename:   DIR/family_packed.geno
admixlist: admix
binsize:        0.001
maxdis:         1.0
seed:           77
runmode:        1
chithresh:      0.0
mincount:       1
zdipcorrmode:   YES
qbin:           10
runfit:  YES 
afffit:  YES 
lovalfit:  0.45 

dates -p par.dates >20250529_test.DATES.log

### after running the DATES, we set test for each individual, by modifying the admix.test, family_packed.ind.test

dates -p par.dates.test >20250529_test.DATES.test.log

## output NA00001.jout, generations ago, se.
## 60.114   9.254




## 1. Data wrangling


### Use 0.00000261 (cM/bp) as the recombination rate for the chromosomes use Morales-Cruz et al 2021

vcf=doaniana_sample94.pruned.mac2.vcf.gz
chrom=vitis.chrom.sizes
group=doaniana_sample94.4g.groupinfo
name=${vcf%.vcf.gz}

## 2. converting to eigensoft format
plink --vcf $vcf --make-bed --aec --double-id --keep-allele-order --set-missing-var-ids @:# --recode 12 --out $name.dates
### use converted ped and map file, and the suffix should be map/ped.
perl -alne '$new=$F[3]*0.00000261; print"$F[0]\t$F[1]\t$new\t$F[3]"' $name.dates.map | perl -p -ne 's/VITVarB40-14_v2.0.hap1.chr0?//' | perl -pe 's/VITVarB40-14_v2.0.hap1.//' >1 ; mv 1 $name.dates.map
awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t"; for(i=7;i<=NF;++i) printf $i"\t";ind++;printf "\n"}' $name.dates.ped > 2; mv 2 $name.dates.ped

file=$name.dates

echo "genotypename:    ${file}.ped" > par.PED.EIGENSTRAT
echo "snpname:         ${file}.map" >> par.PED.EIGENSTRAT
echo "indivname:       ${file}.ped" >> par.PED.EIGENSTRAT
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT
echo "genotypeoutname: ${file}.geno" >> par.PED.EIGENSTRAT
echo "snpoutname:      ${file}.snp" >> par.PED.EIGENSTRAT
echo "indivoutname:    ${file}.ind" >> par.PED.EIGENSTRAT
echo "familynames:     NO" >> par.PED.EIGENSTRAT

convertf -p par.PED.EIGENSTRAT
perl -i.bak -pe 's/^\s+//'  $name.ind

##### for the doaniana_sample94.pruned.mac2.dates.ind file, 
##### it contain 3 column (ind, u, control). 
##### we need to replace the 3rd control column with the group information, 
##### which is the 2nd column in the doaniana_sample94.4g.groupinfo file. 
##### the connection of two files is the first column in the doaniana_sample94.pruned.mac2.dates.ind file 
##### and the first column in the doaniana_sample94.4g.groupinfo file is ind.
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

' $group $name.dates.ind >$name.dates.ind.newgroup

echo $name

## 3. Do the DATES analysis

export LD_LIBRARY_PATH=${BASE_DIR}/lib:$LD_LIBRARY_PATH
export PATH=${BASE_DIR}/bin:$PATH

echo 'DIR: ./
indivname: DIR/doaniana_sample94.pruned.mac2.dates.ind.newgroup
snpname: DIR/doaniana_sample94.pruned.mac2.dates.snp
genotypename:   DIR/doaniana_sample94.pruned.mac2.dates.geno
admixlist: admix
binsize:        0.001 #0.00005
maxdis:         1.0
seed:           777
runmode:        1
chithresh:      0.0
mincount:       1
zdipcorrmode:   YES
qbin:           10 # 50
runfit:  YES 
afffit:  YES 
lovalfit:  0.45 #0.1' > par.dates.doaniana_sample94.pruned.mac2

dates -p par.dates.$name >20250529_doaniana.DATES.log

### after running the DATES, we set test for each individual, by modifying the admix.test
perl -alne 'if ($F[2]=~/doaniana/) {print "$F[0]\t$F[1]\t$F[0]"} else {print} ' doaniana_sample94.mac2.dates.ind.newgroup >doaniana_sample94.mac2.dates.ind.newgroup.ind

grep -v -e "acer" -e "mus" doaniana_sample94.pruned.mac2.dates.ind.group.ind |perl -pe 's/U\t//' | perl -alne '$name="output_$F[0]"; print"acer_rip\tmustangensis\t$F[0]\t$name"' >admix.ind

dates -p par.dates.$name >20250529_doaniana.DATES.log



for file in output_*/*.jout; do
  if [ -f "$file" ]; then
    echo -n "$file"
    cat "$file"
  fi
done > all.output.20250608.jout
