#!/bin/env bash
set -euo pipefail

## 20250318

## 1. prepare the dataset
vcf=doaniana_sample97.mac2.vcf.gz
taxa=group.species.select.txt
name=$(basename "$vcf" .vcf.gz) 


## 2. Data wrangling
vcf2phylip.py -i $vcf -f --output-prefix $name
head -1 $name.min4.phy
###### 97 1684639
###### make sure the taxa sample order is consistent with the vcf file.


## 3. run the Hyde analysis by run_hyde.py
nohup run_hyde_mp.py -i $name.min4.phy -m $taxa -o outgroup -t 4 -n 97 -s 1684639 -j 5 --prefix $name.hyde &

### 3.2 run the individual hyde test

individual_hyde.py -i $name.min4.phy -m $taxa -o outgroup -tr $name.hyde-out-filtered.txt -t 4 -n 97 -s 1684639 -j 5 --prefix $name.hyde.individual 


### 3.3 bootstrap resampling of individuals within hybrid populations.

bootstrap_hyde.py -i $name.min4.phy -m $taxa -o outgroup -tr $name.hyde-out-filtered.txt -t 4 -n 97 -s 1684639 --reps 100 --prefix $name.hyde.bootstrap100 
bootstrap_hyde_mp.py -i $name.min4.phy -m $taxa -o outgroup -tr $name.hyde-out-filtered.txt -t 4 -n 97 -s 1684639 -j 10 --reps 500 --prefix $name.hyde.bootstrap 
## plus
bootstrap_hyde.py -i VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.sample97.min10.phy \
    -m group.species.select.txt \
    -o outgroup \
    -tr doaniana_sample97.mac2.hyde-out-filtered.txt \
    -t 4 -n 97 -s 35567 --reps 100 \
    --prefix VITVarB40-14_v2.0_hap1_filtered_final_pruned_thinned.sample97.hyde.boot100
