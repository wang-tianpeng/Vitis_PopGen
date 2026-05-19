# Hybrid Speciation Workflows

This folder contains scripts for testing hybrid-speciation signals in Vitis populations.

## Contents

- `calculate_heterozygosity_pi_fst.sh`: calculates heterozygosity, nucleotide diversity, and FST-style population summaries.
- `run_admixture_analysis.sh`: prepares PLINK inputs and runs ancestry clustering workflows.
- `run_arizonica_girdiana_analysis.sh`: population-specific analysis workflow for the arizonica-girdiana comparison.
- `run_dates_analysis.sh`: prepares and runs DATES-style admixture timing analyses.
- `run_hyde_analysis.sh`: runs HyDe tests for hybridization.
- `run_pca.sh`: runs genome-wide PCA.
- `run_windowed_pca.sh`: performs windowed PCA analyses.
- `run_triangulaR.sh`: entry point for triangulaR-based tests.
- `cmd.topology.phyml.sh`: phylogenetic topology testing commands.
- `simulate_f1_hybrids.sh` and `simuF1.pl`: simulate F1 hybrid genotypes from parental VCFs.
- `11.2_phy_pt.slurm`: Slurm job template retained for cluster execution.
