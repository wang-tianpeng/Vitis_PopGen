#!/usr/bin/env bash
set -euo pipefail
GENE_BED="VITVarB40-14_v2.0.gene.bed"



for FILE in *_sweep_regions_simu_quant92.tsv; do
    SPECIES=${FILE%_sweep_regions_simu_quant92.tsv}
    
    echo "Processing ${SPECIES}..."
    
    awk 'NR>1 {print $2"\t"int($3)"\t"int($4)}' "$FILE" | \
    bedtools intersect -a "$GENE_BED" -b - -wa -u > "${SPECIES}_sweep_genes.bed"
    
done



echo -e "Species\tNum_Sweeps\tTotal_Length_bp\tMin_Length_bp\tMax_Length_bp\tAvg_Length_bp" > sweeps_regions_summary.tsv

for FILE in *_sweep_regions_simu_quant92.tsv; do
    SPECIES=${FILE%_sweep_regions_simu_quant92.tsv}
    
    # $3 is Start_Pos, $4 is End_Pos.
    awk -v sp="$SPECIES" '
    NR>1 {
        len = int($4) - int($3);
        count++;
        total_len += len;
        
        if (count == 1) {
            min_len = len;
            max_len = len;
        } else {
            if (len < min_len) min_len = len;
            if (len > max_len) max_len = len;
        }
    }
    END {
        if (count > 0) 
            printf "%s\t%d\t%.0f\t%d\t%d\t%.2f\n", sp, count, total_len, min_len, max_len, total_len/count;
        else
            printf "%s\t0\t0\t0\t0\t0\n", sp;
    }
    ' "$FILE" >> sweeps_regions_summary.tsv
done



echo -e "Species\tNum_Sweep_Genes" > sweep_genes_summary.tsv
for GENE_FILE in *_sweep_genes.bed; do
    SPECIES=${GENE_FILE%_sweep_genes.bed}
    NUM_GENES=$(wc -l < "$GENE_FILE")
    echo -e "${SPECIES}\t${NUM_GENES}" >> sweep_genes_summary.tsv
done


## test how many sweeps have genes for each species, how many sweeps have no genes for each species
echo -e "Species\tNum_Sweeps_With_Genes\tNum_Sweeps_Without_Genes" > sweeps_genes_overlap_summary.tsv
GENE_BED="VITVarB40-14_v2.0.gene.bed"

for FILE in *_sweep_regions_simu_quant92.tsv; do
    SPECIES=${FILE%_sweep_regions_simu_quant92.tsv}
    
    awk 'NR>1 {print $2"\t"int($3)"\t"int($4)}' "$FILE" > "tmp_sweeps_${SPECIES}.bed"
    
    COUNT_WITH=$(bedtools intersect -a "tmp_sweeps_${SPECIES}.bed" -b "$GENE_BED" -u | wc -l)
    
    COUNT_WITHOUT=$(bedtools intersect -a "tmp_sweeps_${SPECIES}.bed" -b "$GENE_BED" -v | wc -l)
    
    echo -e "${SPECIES}\t${COUNT_WITH}\t${COUNT_WITHOUT}" >> sweeps_genes_overlap_summary.tsv
    
    rm "tmp_sweeps_${SPECIES}.bed"
done
