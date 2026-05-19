# Species Distribution Modeling

This folder contains R scripts for species distribution modeling, map generation, overlap calculations, and comparisons between ecological overlap and genomic statistics.

## Contents

- `sdm_model_functions.R`: reusable functions for SDM data preparation, model fitting, and prediction.
- `run_sdm_models.R`: main SDM model-fitting workflow.
- `create_sdm_maps.R`: creates distribution maps from model outputs.
- `create_present_lgm_ig_maps.R`: creates present, Last Glacial Maximum, and interglacial map products.
- `plot_doaniana_overlap_maps.R`: creates focused maps for doaniana-related comparisons.
- `plot_general_sdm_maps.R`: creates general SDM and overlap maps.
- `compare_sdm_overlap_metrics.R`: compares geographic overlap metrics.
- `compare_d_statistics_with_sdm.R`: compares D-statistics with SDM-derived overlap and distance metrics.
- `plot_all_sdm_genomic_comparisons.R`: combines SDM and genomic comparison plots.
- `run_sdm_statistical_comparisons.R`: runs broader SDM statistical comparisons.
- `plot_individual_sdm_statistics.R`: plots individual-level SDM/genomic summary relationships.

These scripts require R geospatial and plotting packages such as `terra`, `biomod2`, `dplyr`, `ggplot2`, and related dependencies. Keep raster inputs, occurrence tables, and output directories relative to the project root where possible.
