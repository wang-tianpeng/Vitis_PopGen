#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"
########## 2025-08-20

tree=202504_topology_20group_3x.newick.v5
taxa=taxa.sample451.group20
vcf=VITVarB40-14_v2.0_hap1_filtered_final.sample451.group20.vcf.gz

Dsuite Dtrios -n output.Dsuite.sample451.group20 -t $tree $vcf $taxa

Dsuite Fbranch -P $tree output.Dsuite.sample451.group20_tree.txt >output.Dsuite.sample451.group20_tree.Fbranch.new

python3 ${BASE_DIR}/dtools.py -n output.Dsuite.sample451.group20_tree.Fbranch.new.svg output.Dsuite.sample451.group20_tree.Fbranch.new $tree


Dsuite Dinvestigate -w 100,20 $vcf $taxa trios_riparia.txt

## output files:
# acerifolia_riparia_rupestris_localFstats__100_20.txt
# girdiana_riparia_aestivalis_localFstats__100_20.txt
# rupestris_riparia_arizonica_localFstats__100_20.txt

## we calculate the top 5% quantile fdM value for each of the three tests
for file in *__100_10.txt; do
    echo $file
    Rscript -e "data <- read.table('$file', header=TRUE); q <- quantile(data\$f_dM, 0.95, na.rm=TRUE); write(q, file=paste0('$file','.top5perc.txt'))"
    Rscript -e "data <- read.table('$file', header=TRUE); q <- quantile(data\$f_dM, 0.9, na.rm=TRUE); write(q, file=paste0('$file','.top10perc.txt'))"
done

for file in *__100_10.txt; do
    q=$(cat $file.top5perc.txt)
    echo "Filtering $file with threshold $q"
    awk -v threshold=$q '$6 >= threshold' $file > ${file%.txt}_top5perc_windows.txt
done

for file in *__100_10.txt; do
    q=$(cat $file.top10perc.txt)
    echo "Filtering $file with threshold $q"
    awk -v threshold=$q '$6 >= threshold' $file > ${file%.txt}_top10perc_windows.txt
done

for file in *_top5perc_windows.txt; do
    echo "Merging $file"
    tail -n +2 $file | bedtools merge -d 20000 -c 6,6 -o mean,count -i - > ${file%.txt}_merged.bed
    wc ${file%.txt}_merged.bed
    ## filter the merged windows with count > 1
    awk '$5 > 1' ${file%.txt}_merged.bed > ${file%.txt}_merged_filtered.bed
done


for file in *_top10perc_windows.txt; do
    echo "Merging $file"
    tail -n +2 $file | bedtools merge -d 20000 -c 6,6 -o mean,count -i - > ${file%.txt}_merged.bed
    wc ${file%.txt}_merged.bed
    ## filter the merged windows with count > 1
    awk '$5 > 1' ${file%.txt}_merged.bed > ${file%.txt}_merged_filtered.bed
done


## for the three tests of introgression to riparia, we have: 
# rupestris_riparia_arizonica_localFstats__100_20_top5perc_windows_merged_filtered.bed
# acerifolia_riparia_rupestris_localFstats__100_20_top5perc_windows_merged_filtered.bed
# girdiana_riparia_aestivalis_localFstats__100_20_top5perc_windows_merged_filtered.bed
cat rupestris_riparia_arizonica_localFstats__100_20_top5perc_windows_merged_filtered.bed acerifolia_riparia_rupestris_localFstats__100_20_top5perc_windows_merged_filtered.bed girdiana_riparia_aestivalis_localFstats__100_20_top5perc_windows_merged_filtered.bed | bedtools sort -i - | bedtools merge -d 20000 -c 4,5 -o mean,sum > riparia_introgression_merge_windows.bed
bedtools intersect -a rupestris_riparia_arizonica_localFstats__100_20_top5perc_windows_merged_filtered.bed -b acerifolia_riparia_rupestris_localFstats__100_20_top5perc_windows_merged_filtered.bed | bedtools intersect -a - -b girdiana_riparia_aestivalis_localFstats__100_20_top5perc_windows_merged_filtered.bed > riparia_introgression_intersect_windows.bed
