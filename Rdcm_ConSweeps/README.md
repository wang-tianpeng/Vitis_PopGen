# Convergent Sweeps and rdmc Analyses

This folder contains scripts for testing overlap between introgression windows, shared selective sweeps, and rdmc model-comparison signals.

## Contents

- `extract_neutral_allele_frequencies.sh`: extracts neutral allele-frequency panels from VCFs.
- `run_neutral_frequency_extraction.sh`: wrapper for neutral-frequency extraction.
- `run_overlap_permutation.sh`: performs permutation tests for overlap between introgression and sweep intervals.
- `plot_permutation_results.R`: plots permutation-test distributions and observed overlap values.
- `run_rdmc_workflow.sh`: end-to-end rdmc workflow for a population pair and shared-sweep BED file.
- `run_rdmc_windows.sh`: window-level rdmc workflow.
- `generate_rdmc_slurm_jobs.sh`: generates Slurm jobs for rdmc population-pair runs.
- `plot_rdmc_manual_profiles.R`: manual plotting of rdmc profiles.
- `tidy_plot_rdmc_results.R`: tidies rdmc outputs and creates summary plots.

## Usage Notes

Most scripts expect BED intervals, a filtered VCF, population sample files, and a genetic map. Adjust population names, file paths, SNP thresholds, and Slurm settings at the top of each script.
