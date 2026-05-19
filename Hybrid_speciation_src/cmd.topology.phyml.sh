#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"


file=doaniana_sample97.mac2.geno.gz
name=doaniana_sample97.mac2


for x in 50 100 200 500 1000
do
echo "Inferring trees with window snp $x"

python ./phyml_sliding_windows.py \
    --threads "${THREADS:-8}" \
    --windType sites -w $x --model GTR --optimise n \
    -g $file \
    --prefix $name.phyml.snp.$x \
    --tmp ${BASE_DIR}/phyml_tmp

done

## 2.2 run sliding window method based on coordinate

for x in 10 20 50 100 200
do
echo "Inferring trees with window snp $x"

python ./phyml_sliding_windows.py \
    --threads "${THREADS:-8}" \
    --windType coordinate -w ${x}000 -S 10000 --model GTR --optimise n \
    -g $file \
    --prefix $name.phyml.win${x}s10 \
    --tmp ${BASE_DIR}/phyml_tmp

done

###  The output files conatains too many NAs, after testing the mac2 filtering methods,
###  It seems phyml will automatically remove sequences if two or more are identical. 
###  Continue to debug by using the raxml method to generate trees.


for x in 50 100 200 500 1000
do
echo "Inferring trees with window snp $x"

python2 ./raxml_sliding_windows.py \
    --threads "${THREADS:-8}" \
    --windType sites -w $x --model GTR \
    -g $file \
    --prefix $name.raxml.snp.$x \
    --raxml ${BASE_DIR}/raxmlHPC \
    --log ${BASE_DIR}/raxml_tmp

done

###   again, the output files conatains too many NAs, which is problematic.



## dir in topology_manually/
# generate_windows.sh
#!/usr/bin/env bash

# generate_windows.sh

window_size=20000    # 20 kb
step_size=20000       # 20 kb step

while read chrom length _; do
  # We'll slide in increments of 10 kb from 1..length
  for start in $(seq 1 $step_size $length); do
    end=$((start + window_size - 1))
    if [ $end -gt $length ]; then
      end=$length
    fi
    echo -e "${chrom}\t${start}\t${end}"
    if [ $end -ge $length ]; then
      break
    fi
  done
done < genome.chr.txt > windows_20kb.bed

# step2.sh, step3.sh, step4.sh, step5.sh, step6.sh, step7.twisst.sh

### results meet some bugs. Try use the phased geno data to debug.




perl -pe 's/\//\|/g' doaniana_sample97.mac2.geno > doaniana_sample97.mac2.geno.phased

file=doaniana_sample97.mac2.geno.phased
name=doaniana_sample97.mac2.phased

for x in 50 100 200 500
do
echo "Inferring trees with window snp $x"

python ./phyml_sliding_windows.py \
    --threads "${THREADS:-8}" \
    --windType sites -w $x --model GTR --optimise n \
    -g $file \
    --prefix $name.phyml.snp.$x \
    --tmp ${BASE_DIR}/phyml_tmp \
    --log ${BASE_DIR}/phyml_tmp.log \
    --verbose

done
