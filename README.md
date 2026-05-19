# Vitis Population Genomics Workflows

Code repository for paper: **The dynamics of introgression and parallel adaptation across North American Vitis species (https://doi.org/10.64898/2026.02.23.707523)**

This directory contains analysis scripts for Vitis population genomics, hybrid speciation, introgression, selection scans, convergent sweeps, phylogenomics, and species distribution modeling.

## Subdirectories

- `Hybrid_speciation_src/`: hybrid-speciation analyses, including PCA, ADMIXTURE, HyDe, DATES, topology tests, and F1 simulations.
- `Popgen_phylo/`: population-genetic summary statistics, phylogenetic conversion, plastome assembly, demographic inference, and divergence-time workflows.
- `Pop_Introgression/`: introgression tests using window statistics, TreeMix, ADMIXTOOLS, and VCF-to-TreeMix conversion.
- `Rdcm_ConSweeps/`: workflows for shared sweeps, neutral allele frequencies, rdmc analyses, and permutation plots.
- `selection_sweeps/`: SweeD and msprime workflows for selection scans and simulation-based thresholds.
- `species_sdm/`: species distribution modeling, map generation, overlap statistics, and genotype-environment comparison plots.

## Reproducibility Notes

Scripts assume standard command-line genomics tools are available on `PATH` where needed, including `bcftools`, `plink`, `vcftools`, `samtools`, `Rscript`, and analysis-specific tools such as SweeD, BEAST, ADMIXTOOLS, TreeMix, HyDe, and triangulaR. Most scripts define configurable variables near the top; adjust paths, thread counts, Slurm resources, and input filenames for your computing environment before running.
