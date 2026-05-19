#!/usr/bin/env Rscript

# Usage: Rscript step4_PLOT_rdmc_generic.R <analysis_directory> <pop1> <pop2>
# Example: Rscript step4_PLOT_rdmc_generic.R rdmc_analysis_advanced.rupestris_vs_riparia rupestris riparia


library(tidyverse)

args_cli <- commandArgs(trailingOnly = TRUE)

if (length(args_cli) < 3) {
  stop("Usage: Rscript step4_PLOT_rdmc_generic.R <analysis_directory> <pop1> <pop2>\nExample: Rscript step4_PLOT_rdmc_generic.R rdmc_analysis_advanced.rupestris_vs_riparia rupestris riparia")
}

analysis_dir <- args_cli[1]
pop1 <- args_cli[2]
pop2 <- args_cli[3]

message(paste("Processing Directory:", analysis_dir))
message(paste("Population 1:", pop1))
message(paste("Population 2:", pop2))

gmap_file <- "vitis_sample451_phased.rdmc.map"
# Note: Using the filename from your snippet. Ensure this file exists.
neutral_file <- "neutral_freqs_group8.v2.filter.tsv" 
taxa_file <- "taxa.sample451.group8.v2"

results_dir <- file.path(analysis_dir, "results")
sweep_dir <- file.path(analysis_dir, "sweep_freq_files")

if (!dir.exists(results_dir)) {
  stop(paste("Results directory not found:", results_dir))
}

message("Loading result files...")
result_files <- list.files(results_dir, pattern = "^result_", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files found in ", results_dir)
}

all_results_raw <- result_files %>%
    map_dfr(read_tsv, show_col_types = FALSE) %>% 
    mutate(locus = str_remove(locus, "sweep_VITVarB40-14_v2.0.hap1."))

message("Processing raw data...")
all_results_raw_data <- all_results_raw %>%
    filter(!is.na(cle)) %>% 
    filter(is.finite(cle))  %>% 
    group_by(selected_sites,locus, model) %>%
    slice_max(cle, n = 1, with_ties = FALSE)  %>% 
    ungroup()

raw_data_outfile <- file.path(analysis_dir, paste0("all_results_raw_data_", pop1, "_", pop2, ".tsv"))
write_tsv(all_results_raw_data, raw_data_outfile)
message("Saved raw data summary to: ", raw_data_outfile)

message("Calculating Delta CLE values...")

plot_data <- all_results_raw_data %>%
    filter(!is.na(cle)) %>% 
    filter(is.finite(cle))  %>% 
    group_by(locus, model) %>%
    # slice_max(cle, n = 1, with_ties = FALSE) %>%
     summarise(max_cle = max(cle, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = model, values_from = max_cle) %>%
    left_join(
        all_results_raw %>% distinct(locus, neut_cle),
        by = "locus"
    ) %>%
    rowwise() %>%
    mutate(
        independent = if("independent" %in% names(.)) independent else NA_real_,
        migration = if("migration" %in% names(.)) migration else NA_real_,
        standing = if("standing" %in% names(.)) standing else NA_real_,
        
        cle_values = list(c(independent, migration, standing)),
        model_names = list(c("independent", "migration", "standing")),
        
        cle_values_clean = list(replace(cle_values, is.na(cle_values), -Inf)),
        sorted_indices = list(order(unlist(cle_values_clean), decreasing = TRUE)),
        
        best_model = model_names[sorted_indices[1]],
        cle_best = cle_values[sorted_indices[1]],
        cle_second_best = cle_values[sorted_indices[2]],
        
        delta_cle_vs_second = cle_best - cle_second_best,
        delta_cle_vs_neutral = cle_best - neut_cle
    ) %>%
    ungroup() %>%
    arrange(best_model, desc(delta_cle_vs_neutral)) %>%
    mutate(plot_index = row_number())

plot_data_outfile <- file.path(analysis_dir, paste0("rdmc_delta_cle_plot_data_", pop1, "_", pop2, ".tsv"))
write_tsv(plot_data, plot_data_outfile)

rds_outfile <- file.path(analysis_dir, paste0("all_results_raw_data_", pop1, "_", pop2, ".rds"))
saveRDS(list(raw = all_results_raw_data, plot = plot_data), rds_outfile)
message("Saved plot data to: ", plot_data_outfile)

message("Generating plot...")

my_color <- c("independent" = "#f1a340", "migration" = "#30688e", "standing" = "#5c9b7c")

final_plot <- ggplot(plot_data, aes(x = plot_index)) +
    geom_point(aes(y = delta_cle_vs_neutral, color = best_model)) +
    scale_color_manual(values = my_color) +
    
    scale_y_log10(
        breaks = 10^(-1:3), 
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    
    facet_wrap(~best_model, scales = "free_x") +
    
    labs(
        title = paste("RDMC Analysis:", pop1, "vs", pop2),
        x = "Index",
        y = expression(Delta[CLE])
    ) +
    theme_bw() +
    theme(
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", face = "bold")
    )

plot_outfile <- file.path(analysis_dir, paste0("rdmc_delta_cle_plot_", pop1, "_", pop2, ".png"))
ggsave(plot_outfile, final_plot, width = 8, height = 4, dpi = 300)

message("Step complete.")



message("Generating locus-specific CLE profile plots...")

profile_plot_outfile <- file.path(analysis_dir, paste0("rdmc_cle_profiles_", pop1, "_", pop2, ".pdf"))

unique_loci <- unique(all_results_raw_data$locus)

plots_per_page <- 4

pdf(profile_plot_outfile, width = 12, height = 8)

for(i in seq(1, length(unique_loci), by = plots_per_page)) {
    
    batch_loci <- unique_loci[i : min(i + plots_per_page - 1, length(unique_loci))]
    
    batch_data <- all_results_raw_data %>% 
        filter(locus %in% batch_loci)
    
    p <- ggplot(batch_data, aes(x = selected_sites, y = cle, color = model, group = model)) +
        geom_line(alpha = 0.7) +
        geom_point(size = 2) +
        scale_color_manual(values = my_color) +
        facet_wrap(~locus, scales = "free", ncol = 2) + 
        labs(
            title = paste("CLE Profiles:", pop1, "vs", pop2),
            subtitle = paste("Loci", i, "to", min(i + plots_per_page - 1, length(unique_loci))),
            x = "Selected Site Position", 
            y = "Composite Likelihood Estimate (CLE)"
        ) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = "#f0f0f0"),
            strip.text = element_text(face = "bold", size = 10),
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    
    print(p)
}

dev.off()

message("Step complete.")

