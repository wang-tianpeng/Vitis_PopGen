library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

aba_baba <- read_excel("inputs/samples_SDM_INFO.xlsx", sheet = 4)
areas_dist <- read_csv("email_FollowUp/areas_dist_new.csv")
areas_group <- read_csv("email_FollowUp/aba_baba_overlaps_group.csv")

aba_baba_clean <- as_tibble(aba_baba[, c("P2", "P3", "D statistic")]) %>% 
  filter( !( str_detect(P2, "x_") | str_detect(P3, "x_") ))

aba_baba_clean$P2 <- sub("x_", "", aba_baba_clean$P2) %>% paste("Vitis.", ., sep = "")
aba_baba_clean$P3 <- sub("x_", "", aba_baba_clean$P3) %>% paste("Vitis.", ., sep = "")

aba_baba_clean$dist_cent <- NA
aba_baba_clean$Area_Overlap <- NA

for(i in 1:nrow(aba_baba_clean)){
  pair_key <- paste(sort(c(aba_baba_clean$P2[i], aba_baba_clean$P3[i])), collapse = "_")
  
  match_idx <- which(apply(areas_dist[, c("sp1", "sp2")], 1, function(x) paste(sort(x), collapse = "_")) == pair_key)
  
  if(length(match_idx) > 0){
    aba_baba_clean$dist_cent[i] <- areas_dist$dist_cent[match_idx]
    aba_baba_clean$Area_Overlap[i] <- areas_dist$Area_Overlap[match_idx]
  }
}


#  filter(! ( str_detect(sp1, "doaniana") | str_detect(sp2, "doaniana") ) )  %>% 
#  filter(! ( str_detect(sp1, "champinii") | str_detect(sp2, "champinii") ) ) %>% 
#  select(1,2,)

plot_data <- na.omit(aba_baba_clean)

plot_correlation <- function(data, x_var, y_var, x_label, y_label) {
  
  formula <- as.formula(paste0("`", y_var, "` ~ `", x_var, "`"))
  fit <- lm(formula, data = data) %>% summary()
  
  r_val <- cor(data[[x_var]], data[[y_var]], use = "complete.obs")
  p_val <- fit$coefficients[2, 4]
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001\nR = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nR = ", round(r_val, 2))
  }
  
  p <- ggplot(data, aes_string(x = x_var, y = paste0("`", y_var, "`"))) +
    geom_point(alpha = 0.6, color = "#2c7bb6", size = 3) +
    geom_smooth(method = "lm", color = "#E76254", fill = "#E76254", alpha = 0.2) + #loess
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.1, vjust = 1.2, size = 5, fontface = "italic") +
    labs(x = x_label, y = y_label) +
    ylim(0, 0.6) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  return(p)
}

plot_a <- plot_correlation(plot_data, "Area_Overlap", "D statistic", 

plot_b <- plot_correlation(plot_data, "dist_cent", "D statistic", 
                           "Species Centroid Distance (km)", "D statistic")

final_plot <- ggarrange(plot_a, plot_b, nrow = 2, labels = c("A", "B"))

print(final_plot)

ggsave("./D_statistic_correlation_NoHyb.pdf", final_plot, width = 6, height = 8)
ggsave("./D_statistic_correlation_NoHyb.png", final_plot, width = 6, height = 8, dpi = 300)

ggsave("./D_statistic_FIG2D_NoHyb.pdf", plot_b, width = 6, height = 5)
ggsave("./D_statistic_FIG2D_NoHyb.png", plot_b, width = 6, height = 5, dpi = 300)

message("Plots saved to output/D_statistic_correlation.pdf and .png")


## new log regression model


plot_log_correlation <- function(data, x_var, y_var, x_label, y_label) {
  
  
  data_log <- data[data[[x_var]] > 0, ]
  
  formula_log <- as.formula(paste0("`", y_var, "` ~ log(`", x_var, "`)"))
  fit <- lm(formula_log, data = data_log) %>% summary()
  
  p_val <- fit$coefficients[2, 4]
  r_squared <- fit$r.squared
  slope_sign <- sign(fit$coefficients[2, 1])
  r_val <- slope_sign * sqrt(r_squared)
  
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001","\nr = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nr = ", round(r_val, 2))
  }
  
  p <- ggplot(data_log, aes_string(x = x_var, y = paste0("`", y_var, "`"))) +
    geom_point(alpha = 0.6, color = "#2c7bb6", size = 3) +
    
    geom_smooth(method = "lm", formula = y ~ log(x), 
                color = "#E76254", fill = "#E76254", alpha = 0.2) +
    
    annotate("text", x = Inf, y = Inf, label = stats_text, 
         hjust = 1.2, vjust = 1.2, size = 5, fontface = "italic") +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  return(p)
}

plot_log_dist <- plot_log_correlation(plot_data, "dist_cent", "D statistic", 
                                      "Centroid Distance (km)", "D statistic")

print(plot_log_dist)

ggsave("output/D_statistic_Centroid_LogFit.pdf", plot_log_dist, width = 6, height = 5)
ggsave("output/D_statistic_Centroid_LogFit.png", plot_log_dist, width = 6, height = 5)


### FIG 2C: Boxplot of D statistic by Overlap Category (No Hybrids) ###
areas_group_clean <- areas_group %>% 
  filter(! ( str_detect(P2, "doaniana") | str_detect(P3, "doaniana") ) )  %>% 
  filter(! ( str_detect(P2, "champinii") | str_detect(P3, "champinii") ) ) %>% 
  mutate(cat = str_replace(cat, "no_Overlap", "No_Overlap")) %>% 
  mutate(cat = str_replace(cat, "Other", "Eurasian_taxa")) %>%
  mutate(cat = factor(cat, levels = c("Overlap", "No_Overlap", "Eurasian_taxa")))


  plot_boxplot_group <- function(data, x_var, y_var, y_label) {
    
    my_colors <- c("Overlap" = "#E76254", "No_Overlap" = "#486A8C", "Eurasian_taxa" = "gray50")
    
    p <- ggplot(data, aes_string(x = x_var, y = paste0("`", y_var, "`"), fill = x_var, color = x_var)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(position = position_jitter(0.1), size = 2, alpha = 0.6) +
    
    stat_compare_means(method = "anova", label.y = 0.55, size = 5, color = "black") + 
    

    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    ylim(0, 0.6) +
    labs(x = "Species Niche Overlap Category", y = y_label) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    )
    return(p)
  }

    # stat_compare_means(comparisons = my_comparisons, method = "t.test", 
    # # stat_compare_means(label.y.npc = "top", size = 4, color = "black") + 
    
plot_box <- plot_boxplot_group(areas_group_clean, "cat", "D statistic", "D statistic")
ggsave("output/D_statistic_Group_Boxplot_noHyb.pdf", plot_box, width = 6, height = 5)
ggsave("output/D_statistic_Group_Boxplot_noHyb.png", plot_box, width = 6, height = 5)




# ==============================================================================
# ==============================================================================

anova_model <- aov(`D statistic` ~ cat, data = areas_group_clean)
print(summary(anova_model))
cat("\n=== ANOVA Results ===\n")


tukey_result <- TukeyHSD(anova_model)

cat("\n=== Tukey HSD Post-hoc Test Results ===\n")
print(tukey_result)

tukey_df <- as.data.frame(tukey_result$cat)

cat("\n=== Specific Contrasts involving 'Overlap' ===\n")
print(tukey_df[grep("Overlap", rownames(tukey_df)), ])
