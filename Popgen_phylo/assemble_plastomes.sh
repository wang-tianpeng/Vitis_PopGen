#!/bin/bash -l
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"


#getOrganelle
#trimmomatic
#mash

## software
PATH=${BASE_DIR}/bin:$PATH
PATH=${BASE_DIR}/bin/:$PATH
TRIMMOMATIC=${BASE_DIR}/trimmomatic-0.39.jar
PATH=${BASE_DIR}/bin:$PATH
source activate getorg

THREADS=4
ADAPTERSPE=${BASE_DIR}/TruSeq3-PE.fa
SEQLIST=${BASE_DIR}/redo2.txt
TEMP_DIR=${BASE_DIR}/temp
RESULTS=${BASE_DIR}/results

ID=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n 1 | cut -f1)
SAMPLE=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n 1 | cut -f2)
PREFIX=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n 1 | cut -f3)
FTP=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n 1 | cut -f4)
MD5=$(head -n "$SLURM_ARRAY_TASK_ID" "$SEQLIST" | tail -n 1 | cut -f5)

TEMP_DIR="$TEMP_DIR"/temp_"$ID"
mkdir -pv "$TEMP_DIR"
echo "$TEMP_DIR"
cd "$TEMP_DIR"

if [ "$FTP" = "NA" ]; then
    java -jar "$TRIMMOMATIC" PE -threads "$THREADS" "$PREFIX"-READ1-Sequences.txt.gz "$PREFIX"-READ2-Sequences.txt.gz "$ID"_1_trimmed_paired.fq.gz "$ID"_1_unpaired.fq.gz "$ID"_2_trimmed_paired.fq.gz "$ID"_2_unpaired.fq.gz ILLUMINACLIP:"$ADAPTERSPE":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:60

else
    INDEX=1
    for i in $(echo "$FTP" | tr ";" "\n")
    do
        echo "downloading" "$i"
        wget -O "$ID"_"$INDEX".fastq.gz "$i"
        INDEX=$((INDEX + 1))
    done

    # checksums
    INDEX=1
    for i in $(echo "$MD5" | tr ";" "\n")
    do
        echo "$i" "$ID"_"$INDEX".fastq.gz >> chk.md5
        INDEX=$((INDEX + 1))
    done

    if md5sum --status -c chk.md5; then
        echo "SUMS CHECKOUT"
    else
        echo "SUMS BAD"
        rm -r "$TEMP_DIR"
        exit 1
    fi

    # trim
    java -jar "$TRIMMOMATIC" PE -threads "$THREADS" "$ID"_1.fastq.gz "$ID"_2.fastq.gz "$ID"_1_trimmed_paired.fq.gz "$ID"_1_unpaired.fq.gz "$ID"_2_trimmed_paired.fq.gz "$ID"_2_unpaired.fq.gz ILLUMINACLIP:"$ADAPTERSPE":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:60 

fi


get_organelle_from_reads.py \
-1 "$ID"_1_trimmed_paired.fq.gz \
-2 "$ID"_2_trimmed_paired.fq.gz \
-s ${BASE_DIR}/V_vinifera_Pt.fasta \
--reduce-reads-for-coverage inf --max-reads inf \
-R 30 \
-w 95 \
--disentangle-time-limit 7200 \
--overwrite \
-t "$THREADS" \
-o "$RESULTS"/assemblies/"$ID" \
-F embplant_pt \
--which-spades ${BASE_DIR}/bin \
--prefix "$ID" \
--config-dir ${BASE_DIR}/getorg

# cleanup
rm -r "$TEMP_DIR"
