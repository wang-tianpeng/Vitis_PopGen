#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"

## select groups for plotting

cut -f 2 taxa.sample639.group |sort |uniq -c |perl -alne 'print $F[1] if $F[0] > 1' >taxa.sample639.group30.select
perl filter_Dsuite.pl taxa.sample639.group30.select alldata_BBAA.txt >alldata_group30_BBAA.txt

## plot introgression
### plot D
ruby plot_d.rb alldata_group30_BBAA.txt plot_order.txt 0.3 test.d.03.svg &
ruby plot_d.rb alldata_group30_BBAA.txt plot_order.txt 0.5 test.d.05.svg &
ruby plot_d.rb alldata_group30_BBAA.txt plot_order.txt 0.7 test.d.07.svg &
ruby plot_d.rb alldata_group30_BBAA.txt plot_order.txt 0.9 test.d.09.svg &

### plot f4-ratio
ruby plot_f4ratio.rb alldata_group30_BBAA.txt plot_order.txt 0.3 test.f4.03.svg &
ruby plot_f4ratio.rb alldata_group30_BBAA.txt plot_order.txt 0.5 test.f4.05.svg &
ruby plot_f4ratio.rb alldata_group30_BBAA.txt plot_order.txt 0.7 test.f4.07.svg &
ruby plot_f4ratio.rb alldata_group30_BBAA.txt plot_order.txt 0.9 test.f4.09.svg &


# 20250409
### filter and select samples for n>=5

cut -f 2 taxa.sample639.group |sort |uniq -c| perl -alne 'print $F[1] if $F[0] >= 3' | perl -alne 'print unless $_=~/xxx/' >taxa.group.n3.id
grep -f taxa.group.n3.id taxa.sample639.group > taxa.sample639.group.n3

### use Dsuite to calculate D
vcf=VITVarB40-14_v2.0_hap1_filtered_final.vcf.gz
taxa=taxa.sample639.group20.n3
tree=202504_topology_21group_3x.newick.v3

perl -alne ' print "$F[0]\t$F[0]"' $taxa > $taxa.plink.id
plink --vcf $vcf --keep $taxa.plink.id \
    --mac 3 --aec --double-id --recode vcf-iid \
    --out VITVarB40-14_v2.0_hap1_filtered_final.sample639.group20.n3

Dsuite Dtrios -n output.Dsuite.sample639.group20.n3 -t $tree VITVarB40-14_v2.0_hap1_filtered_final.sample639.group20.n3.vcf $taxa
Dsuite Dtrios --ABBAclustering -n output.Dsuite.sample639.group20.n3.ABBAclustering -t $tree VITVarB40-14_v2.0_hap1_filtered_final.sample639.group20.n3.vcf $taxa

## for the v3 ABBA results, use the Fbranch to plot
tree=202504_topology_20group_3x.newick.v5
Dsuite Fbranch 202504_topology_20group_3x.newick.v5 taxa.sample639.group20_output.Dsuite.sample639.group20.n3.v3_tree.txt >taxa.sample639.group20_output.Dsuite.sample639.group20.n3.v3_tree.Fbranch.new

Fbranch=taxa.sample639.group20_output.Dsuite.sample639.group20.n3.v3_tree.Fbranch.new

python3 ${BASE_DIR}/dtools.py \
    -n plot.Fbranch.group20.v3.tree.v5g \
    $Fbranch $tree



## filter and select 21 samples representing 21 species 

cut -f 1 202507_sample21.spe21.group.id | seqkit grep -f - VITVarB40-14_v2.0_hap1_filtered.4distes.min10.fasta.outgroup.discard.fasta |seqkit seq -w 0 >sample21.spe21.4distes.fasta

vcf2phylip.py -i VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.vcf.gz -f -m 4
cut -f 1 202507_sample21.spe21.group.id | seqkit grep -f - VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.min4.fasta |seqkit seq -w 0 >sample21.spe21.thinned.fasta

## 20250713 - Run IQ-TREE on 4distes dataset
iqtree -s sample21.spe21.4distes.fasta -st DNA -B 1000 -nt AUTO -pre sample21.spe21.4distes

iqtree -s sample21.spe21.thinned.fasta -st DNA -B 1000 -nt AUTO -pre sample21.spe21.thinned

## convert the indiviuals name in sample21.spe21.4distes.treefile to the group name, with the group file 202507_sample21.spe21.group.id
awk '{print $1"\t"$2}' 202507_sample21.spe21.group.id > sample21.name_to_group.map

perl -e '
    open(MAP, "<sample21.name_to_group.map") or die "Cannot open mapping file: $!";
    while(<MAP>) {
        chomp;
        my ($ind, $group) = split(/\t/, $_);
        $map{$ind} = $group;
    }
    close(MAP);
    
    open(TREE, "<sample21.spe21.4distes.treefile") or die "Cannot open tree file: $!";
    my $tree = <TREE>;
    close(TREE);
    
    foreach my $ind (keys %map) {
        $tree =~ s/\b$ind\b/$map{$ind}/g;
    }
    
    open(OUT, ">sample21.spe21.4distes.group_names.treefile") or die "Cannot open output file: $!";
    print OUT $tree;
    close(OUT);
    
    print "Tree with group names saved to sample21.spe21.4distes.group_names.treefile\n";
' 

perl -e '
    open(MAP, "<sample21.name_to_group.map") or die "Cannot open mapping file: $!";
    while(<MAP>) {
        chomp;
        my ($ind, $group) = split(/\t/, $_);
        $map{$ind} = $group;
    }
    close(MAP);
    
    open(TREE, "<sample21.spe21.thinned.treefile") or die "Cannot open tree file: $!";
    my $tree = <TREE>;
    close(TREE);
    
    foreach my $ind (keys %map) {
        $tree =~ s/\b$ind\b/$map{$ind}/g;
    }
    
    open(OUT, ">sample21.spe21.thinned.group_names.treefile") or die "Cannot open output file: $!";
    print OUT $tree;
    close(OUT);
    
    print "Tree with group names saved to sample21.spe21.thinned.group_names.treefile\n";
'
### results not good. using one sample per species as representative is not good.
### considering mannually constructing a tree with 21 species




########## 2025-08-20

tree=202504_topology_20group_3x.newick.v5
taxa=taxa.sample451.group20
vcf=VITVarB40-14_v2.0_hap1_filtered_final.sample451.group20.vcf.gz

Dsuite Dtrios -n output.Dsuite.sample451.group20 -t $tree $vcf $taxa

Dsuite Fbranch -P $tree output.Dsuite.sample451.group20_tree.txt >output.Dsuite.sample451.group20_tree.Fbranch.new

python3 ${BASE_DIR}/dtools.py -n output.Dsuite.sample451.group20_tree.Fbranch.new.svg output.Dsuite.sample451.group20_tree.Fbranch.new $tree

Dsuite Dinvestigate -w 100,20 $vcf $taxa trios_riparia.txt



Rscript adjust_pvalues_dsuite.R taxa.sample639.group20_output.Dsuite.sample639.group20.n3.v3_tree.txt output.txt
