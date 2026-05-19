#!/usr/bin/env bash
set -euo pipefail
# date: 2025-04-09

## update: 2025-04-13

tree=202504_topology_21group_3x.newick
group=taxa.sample639.group20.n3

perl -alne 'print "$F[0]\t$F[0]"' $group > $group.plink.id
## for treemix clust file
perl -alne 'print "$F[0]\t$F[0]\t$F[1]"' taxa.sample639.group20.n3 >taxa.sample639.group20.n3.clust

## 1. filter and select samples from vcf
raw=VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.vcf.gz
prefix=$(basename "$raw" .vcf.gz)
clust=taxa.sample639.group20.n3.clust

plink --vcf $raw --keep $group.plink.id --aec --double-id --export vcf-iid --out $prefix.group20

plink --vcf $prefix.group20 --make-bed --out $prefix.group20 --allow-no-sex --allow-extra-chr
plink --bfile $prefix.group20 --freq --missing --within $clust --out $prefix.group20 --allow-no-sex --allow-extra-chr

file=$prefix.group20
gzip $file".frq.strat"

plink2treemix.py $file".frq.strat.gz" $file".treemix.frq.gz"

gunzip $file".treemix.frq.gz"
gunzip $file".frq.strat.gz"

awk 'BEGIN{print "scaffold_pos\tscaffold\tpos"}{split($2,pos,":");print $2"\t"pos[1]"\t"pos[2]}' $file".map" > $file".positions"
paste $file".positions" $file".treemix.frq" > $file".frequencies"


<<EOC
awk '{
    printf $0;
    for(i = 4; i <= NF; i++){
        split($i,values,",");
        if((values[1] + values[2]) > 0) 
            freq=values[1]/(values[1] + values[2]);
            else
                freq=0;
            printf freq"\t";
            }
	printf "\n";
    }' $file".frequencies">$file".frequencies2"
mv $file".frequencies2" $file".frequencies"

awk 'BEGIN{scaffold="";pos=0;newpos=0} {if($2==scaffold){newpos=pos+$3}else{scaffold=$2;pos=newpos} chpos=pos+$3; print $0, chpos}' "$file.frequencies2" > "$file.frequencies.newpos"


awk 'BEGIN{scaffold="";pos=0;newpos=0}{
if($2==scaffold){
newpos=pos+$3;
}else
{scaffold=$2;pos=newpos};chpos=pos+$3;print $0,chpos}'$file".frequencies2">$file".frequencies.newpos"
EOC

gzip $file".treemix.frq"


### TREEMIX ANALYSIS, for each run, use 5 replicates


for i in {1..5}; do for m in {1..20}; do echo "treemix -i $file.treemix.frq.gz -tr 202504_topology_21group_3x.newick -m $m -o res_v2/$file.treemix.frq.$i.$m -root Outgroup -bootstrap -k 10 1>res_v2/treemix_${i}_${m}.log 2>res_v2/treemix_${i}_${m}.err"; done; done >cmd.treemix_all.sh

ParaFly -c cmd.treemix_all.sh -CPU 15 -failed_cmds cmd.treemix_all.failed_cmds

zip -r output_treemix.zip VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.group20.treemix.frq.*
