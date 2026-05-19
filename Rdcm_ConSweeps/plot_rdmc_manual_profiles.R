
# .berlandieri_vs_cinerea
# .girdiana_vs_arizonica
# .rupestris_vs_riparia
# .mustangensis_vs_monticola


library(tidyverse)
args=list()
args$gmap_file="vitis_sample451_phased.rdmc.map"
args$sweep_dir="rdmc_analysis_advanced.rupestris_vs_riparia/sweep_freq_files/"
args$neutral_file="neutral_freqs_group8.v2.filter.tsv"
args$taxa_file="taxa.sample451.group8.v2"
args$output_dir="rdmc_analysis_advanced.rupestris_vs_riparia/results/"
args$pop1="rupestris"
args$pop2="riparia"



all_results_raw <- list.files(args$output_dir, pattern = "^result_", full.names = TRUE) %>%
    map_dfr(read_tsv)  %>% 
    mutate(locus = str_remove(locus, "sweep_VITVarB40-14_v2.0.hap1."))


#     map_dfr(read_tsv)  %>% 
#     mutate(locus = str_remove(locus, "sweep_VITVarB40-14_v2.0.hap1."))




all_results_raw_data <- all_results_raw %>%
    filter(!is.na(cle)) %>% 
    filter(is.finite(cle))  %>% 
    group_by(locus, model) %>%
    slice_max(cle, n = 1, with_ties = FALSE)  %>% ungroup()

write_tsv(all_results_raw_data, file.path(dirname(args$output_dir), paste0("all_results_raw_data_",args$pop1,".tsv") ))



plot_data <- all_results_raw %>%
    filter(!is.na(cle)) %>% 
    filter(is.finite(cle))  %>% 
    group_by(locus, model) %>%
    # slice_max(cle, n = 1, with_ties = FALSE)
    summarise(max_cle = max(cle, na.rm = T)) %>%
    ungroup() %>%
    pivot_wider(names_from = model, values_from = max_cle) %>%
    left_join(
        all_results_raw %>% distinct(locus, neut_cle),
        by = "locus"
    ) %>%
    # Now, for each locus (row), find the best and second best model
    rowwise() %>%
    mutate(
        cle_values = list(c(independent, migration, standing)),
        model_names = list(c("independent", "migration", "standing")),
        
        sorted_indices = list(order(cle_values, decreasing = TRUE)),
        
        best_model = model_names[sorted_indices[1]],
        cle_best = cle_values[sorted_indices[1]],
        cle_second_best = cle_values[sorted_indices[2]],
        
        delta_cle_vs_second = cle_best - cle_second_best,
        delta_cle_vs_neutral = cle_best - neut_cle
    ) %>%
    ungroup() %>%
    # select(locus, best_model, delta_cle_vs_second, delta_cle_vs_neutral) %>%
    arrange(best_model, desc(delta_cle_vs_neutral)) %>%
    mutate(plot_index = row_number())

write_tsv(plot_data, file.path(dirname(args$output_dir), paste0("rdmc_delta_cle_plot_data_",args$pop1,"_",args$pop2,".tsv") ))
saveRDS(c(all_results_raw_data,plot_data), file.path(dirname(args$output_dir), paste0("all_results_raw_data_",args$pop1,"_",args$pop2,".rds") ))

######### DECIDE TO PLOT IN LOCAL R




# library(ggplot2)


plot_data

my_color <- c("independent" = "#f1a340", "migration" = "#30688e", "standing" = "#5c9b7c")

final_plot <- ggplot(plot_data, aes(x = plot_index)) +
    # geom_segment(
    #     aes(y = delta_cle_vs_second,  yend = delta_cle_vs_neutral, xend = plot_index, color = best_model),
    # ) +
    geom_point(aes(y = delta_cle_vs_neutral, color = best_model)) +
    scale_color_manual(values = my_color) +
    
    scale_y_log10(
        breaks = 10^(-1:3), 
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    
    facet_wrap(~best_model, scales = "free_x") +
    
    labs(
        x = "Index",
        y = expression(Delta[CLE]) # Using expression for the delta symbol
    ) +
    theme_bw() +
    theme(
        legend.position = "none", # Hide the legend as facets serve the purpose
        strip.background = element_rect(fill = "white"), # Facet title background
        strip.text = element_text(color = "black", face = "bold") # Facet title text
    )
ggsave(final_plot,filename="test.png")
plot_file <- file.path(dirname(args$output_dir), "rdmc_delta_cle_plot.png")
ggsave(plot_file, final_plot, width = 8, height = 4, dpi = 300)

message("Step complete.")

# # The original summary file generation can be simplified now
#     select(locus, best_model, delta_cle_vs_second, delta_cle_vs_neutral)

# write_tsv(best_models_summary, summary_file)
