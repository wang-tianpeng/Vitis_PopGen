#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"
VCF="${VCF:-data/VITVarB40-14_v2.0_hap1_filtered_final.sample451.group20.vcf.gz}"
TAXA_FILE="${TAXA_FILE:-data/taxa.sample451.group20}"
PREFIX="${PREFIX:-vitis_451}"
export PREFIX TAXA_FILE

cd "${BASE_DIR}"
mkdir -p data eigenstrat results/f2 results/tables

plink --vcf "${VCF}" \
    --double-id --aec --recode 12 --set-missing-var-ids @:# \
    --keep-allele-order --out "data/${PREFIX}"

perl -alne '$new=$F[3] * 0.00000261; print "$F[0]\t$F[1]\t$new\t$F[3]"' "data/${PREFIX}.map" \
    | perl -pe 's/VITVarB40-14_v2.0.hap1.chr0?//' \
    | perl -pe 's/VITVarB40-14_v2.0.hap1.//' \
    > "data/${PREFIX}.map.tmp"
mv "data/${PREFIX}.map" "data/${PREFIX}.map.raw"
mv "data/${PREFIX}.map.tmp" "data/${PREFIX}.map"

awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t"; for(i=7;i<=NF;i++) printf $i"\t"; ind++; printf "\n"}' \
    "data/${PREFIX}.ped" > "data/${PREFIX}.ped.tmp"
mv "data/${PREFIX}.ped" "data/${PREFIX}.ped.raw"
mv "data/${PREFIX}.ped.tmp" "data/${PREFIX}.ped"

cat > par.convertf <<EOF
genotypename:    data/${PREFIX}.ped
snpname:         data/${PREFIX}.map
indivname:       data/${PREFIX}.ped
outputformat:    EIGENSTRAT
genotypeoutname: eigenstrat/${PREFIX}.geno
snpoutname:      eigenstrat/${PREFIX}.snp
indivoutname:    eigenstrat/${PREFIX}.ind
familynames:     YES
EOF

convertf -p par.convertf

paste <(cut -f1 "data/${PREFIX}.fam") \
      <(cut -f2 "${TAXA_FILE}") \
      <(yes U | head -n "$(wc -l < "data/${PREFIX}.fam")") \
      > "eigenstrat/${PREFIX}.ind"

Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
    library(admixtools)
    library(dplyr)
    library(readr)
})

prefix <- file.path("eigenstrat", Sys.getenv("PREFIX", "vitis_451"))
taxa_file <- Sys.getenv("TAXA_FILE", "data/taxa.sample451.group20")
f2_dir <- "results/f2"
dir.create(f2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

popinfo <- read_tsv(taxa_file, col_names = c("sample", "pop"), show_col_types = FALSE)
pops <- unique(popinfo$pop)

extract_f2(prefix, f2_dir, pops = pops)
f2_blocks <- f2_from_geno(
    pref = prefix,
    pops = pops,
    blgsize = 0.05,
    maxmiss = 0.1
)

saveRDS(f2_blocks, "results/f2/f2_blocks.rds")

f3_results <- f3(f2_blocks)
write_tsv(arrange(f3_results, z), "results/tables/f3_scan.tsv")
RSCRIPT
