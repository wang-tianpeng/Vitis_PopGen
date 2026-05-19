#!/bin/bash
set -euo pipefail


if [ $# -ne 2 ]
 then
 echo "Please provide the following arguments: <vcf file> <clust file>"
 echo "The .clust file contains three columns: samplename\tsamplename\tgroup"
 exit 1
fi

clust=$2
file=${1%.gz}
file=${file%.vcf}

plink --vcf $file.vcf --aec --double-id --recode --out $file

awk -F"\t" '{
        split($2,chr,":")
	$1="1"
	$2="1:"chr[2]
        print $0
}' ${file}.map > better.map
mv better.map ${file}.map

plink --file $file --make-bed --out $file --allow-no-sex --allow-extra-chr
plink --bfile $file --freq --missing --within $2 --out $file --allow-no-sex --allow-extra-chr

gzip $file".frq.strat"

plink2treemix.py $file".frq.strat.gz" $file".treemix.frq.gz"

gunzip $file".treemix.frq.gz"
gunzip $file".frq.strat.gz"

awk 'BEGIN{print "scaffold_pos\tscaffold\tpos"}{split($2,pos,":");print $2"\t"pos[1]"\t"pos[2]}' $file".map" > $file".positions"
paste $file".positions" $file".treemix.frq" > $file".frequencies"

awk '{printf $0
	for(i = 4; i <= NF; i++){
		split($i,values,",")
		if((values[1]+values[2])>0) freq=values[1]/(values[1]+values[2])
		else freq=0
		printf freq"\t"
	}
	printf "\n"}' $file".frequencies" > $file".frequencies2"
mv $file".frequencies2" $file".frequencies"

awk 'BEGIN{scaffold="";pos=0;newpos=0}
	{if($2==scaffold){newpos=pos+$3}else{scaffold=$2;pos=newpos};chpos=pos+$3;print $0,chpos}' \
	$file".frequencies" > $file".frequencies.newpos"

gzip $file".treemix.frq"
