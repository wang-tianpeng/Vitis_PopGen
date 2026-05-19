# Population Genetics and Phylogenomics

This folder contains scripts for population-genetic summaries, phylogenetic data preparation, plastome assembly, and demographic inference.

## Contents

- `assemble_plastomes.sh`: plastome assembly workflow from paired-end reads.
- `calculate_pi_fst_dxy.sh`: computes nucleotide diversity, FST, and Dxy-style summaries.
- `cmd.phylo.sh`: phylogenetic tree-building command wrapper.
- `run_divergence_time_beast.sh`: BEAST divergence-time workflow.
- `run_population_pca.sh`: PCA workflow for population structure.
- `run_structure_admixture.sh`: converts VCFs and runs structure/ADMIXTURE-style analyses.
- `run_mushi_vitis.py` and `run_mushi_vitis_Ver4.py`: demographic inference with mushi from folded SFS data.
- `vcf2phylip.py`: converts VCF alignments to PHYLIP/FASTA/NEXUS formats.
- `WGSpipeline_Step1_V1.pl` and `WGSpipeline_Step2_multisamples.pl`: whole-genome resequencing mapping and multisample variant-calling helper scripts.
