library(tidyverse)
library(ggpubr)

df_ind_data <- read_csv("email_FollowUp/df_coords_clean.csv")

df_ind_data_fil <- df_ind_data %>% filter( !str_detect(Species_label, "x ") )

y_var <- "admix_1minusmax"
y_label <- "Admixture (1-Q)"

print(names(df_ind_data))

plot_geo_corr <- function(y_var, y_label) {
  
  x_val <- as.numeric(df_ind_data_fil$geo_cent)
  y_val <- as.numeric(df_ind_data_fil[[y_var]])
  
  valid_idx <- !is.na(x_val) & !is.na(y_val)
  x_val <- x_val[valid_idx]
  y_val <- y_val[valid_idx]
  plot_df <- data.frame(x = x_val, y = y_val)
  
  fit <- lm(y ~ x, data = plot_df) %>% summary()
  p_val <- fit$coefficients[2, 4]
  r_val <- cor(x_val, y_val)
  
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001\nR = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nR = ", round(r_val, 2))
  }
  
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(alpha = 0.5, color = "#2c7bb6", size = 2) +
    geom_smooth(method = "lm", color = "#E76254", fill = "#E76254", alpha = 0.2) +
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.2, vjust = 1.2, size = 5, fontface = "italic") +
    labs(x = "Sample Relative Distance to Centroid", y = y_label) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  return(p)
}

p1 <- plot_geo_corr("maxprop", "Ancestry proportion (Q)")
p2 <- plot_geo_corr("admix_1minusmax", "Admixture (1 - Q)")
p3 <- plot_geo_corr("H_raw", "Shannon Entropy (H)")
p4 <- plot_geo_corr("effN", "Effective Ancestry (EffN)")
p5 <- plot_geo_corr("simpson_norm", "Normalized Simpson Index")

final_plot <- ggarrange(p1, p2, p3, p4, p5, ncol = 3, nrow = 2, labels = "AUTO")

print(final_plot)

ggsave("output/Individual_Geo_Correlations_new.pdf", final_plot, width = 18, height = 10)
ggsave("output/Individual_Geo_Correlations_new.png", final_plot, width = 18, height = 10, dpi = 300)

ggsave("output/Individual_Geo_Correlations_new_p2_FIG2E.pdf", p2, width = 6, height = 5)
ggsave("output/Individual_Geo_Correlations_new_p2_FIG2E.png", p2, width = 6, height = 5, dpi = 300)


message("Correlation plots saved to output/Individual_Geo_Correlations.pdf")




library(tidyverse)
library(ggpubr)
library(ggExtra)
library(cowplot)

df_ind_data <- read_csv("email_FollowUp/df_coords_clean.csv")

df_ind_data_fil <- df_ind_data %>% filter( !str_detect(Species_label, "x ") )

y_var <- "admix_1minusmax"
y_label <- "Admixture (1-Q)"
plot_geo_corr_marginal <- function(y_var, y_label) {
  
  x_val <- as.numeric(df_ind_data_fil$geo_cent)
  y_val <- as.numeric(df_ind_data_fil[[y_var]])
  
  valid_idx <- !is.na(x_val) & !is.na(y_val)
  plot_df <- data.frame(x = x_val[valid_idx], y = y_val[valid_idx])
  
  fit <- lm(y ~ x, data = plot_df) %>% summary()
  p_val <- fit$coefficients[2, 4]
  r_val <- cor(plot_df$x, plot_df$y)
  
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001\nr = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nr = ", round(r_val, 2))
  }
  
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(alpha = 0.5, color = "#2c7bb6", size = 2) +
    geom_smooth(method = "lm", color = "#d7191c", fill = "#d7191c", alpha = 0.2) +
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.1, vjust = 1.2, size = 4, fontface = "italic") +
    labs(x = "Relative Distance to Centroid", y = y_label) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    )
  
  p_marg <- ggMarginal(p, type = "density", fill = "gray80", color = "gray50", size = 4)
  
  return(p_marg)
}

p1 <- plot_geo_corr_marginal("maxprop", "Max Admixture Prop")
p2 <- plot_geo_corr_marginal("admix_1minusmax", "Admixture (1-Q)")
p3 <- plot_geo_corr_marginal("H_raw", "Shannon Entropy")
p4 <- plot_geo_corr_marginal("effN", "Effective Ancestry")
p5 <- plot_geo_corr_marginal("simpson_norm", "Norm. Simpson Index")

final_plot <- plot_grid(p1, p2, p3, p4, p5, ncol = 3, labels = "AUTO")

print(final_plot)

ggsave("output/Individual_Geo_Correlations_Marginal.pdf", final_plot, width = 15, height = 10)
ggsave("output/Individual_Geo_Correlations_Marginal.png", final_plot, width = 15, height = 10, dpi = 300)

message("Marginal plots saved.")








library(ggpmisc)

plot_geo_corr_facet_all <- function(y_var, y_label) {
  
  plot_data <- df_ind_data_fil
  plot_data$x <- as.numeric(plot_data$geo_cent)
  plot_data[[y_var]] <- as.numeric(plot_data[[y_var]])
  
  p <- ggplot(plot_data, aes_string(x = "x", y = y_var)) +
    
    geom_smooth(method = "lm", color = "black", fill = "gray50", alpha = 0.2) +
    geom_point(alpha = 0.8, aes(color = Species), size = 2) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top",
             size = 3.5, p.accuracy = 0.001, r.accuracy = 0.01,
             label.sep = "\n") +
    
    facet_wrap(~Species, scales = "free") +
    
    labs(x = "Sample Relative Distance to Centroid", y = y_label, 
         title = paste(y_label)) +
    
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey40"),
      legend.position = "none",
      strip.text = element_text(face = "bold")
  
  return(p)
}

pdf("output/Individual_Geo_Correlations_Facet_Species.pdf", width = 16, height = 12)

print(plot_geo_corr_facet_all("maxprop", "Ancestry proportion (Q)"))

print(plot_geo_corr_facet_all("admix_1minusmax", "Admixture (1 - Q)"))

print(plot_geo_corr_facet_all("H_raw", "Shannon Entropy (H)"))

print(plot_geo_corr_facet_all("effN", "Effective Ancestry (EffN)"))

print(plot_geo_corr_facet_all("simpson_norm", "Normalized Simpson Index"))

dev.off()


p_facet_2 <- plot_geo_corr_facet_all("admix_1minusmax", "Admixture (1 - Q)")
ggsave("output/Individual_Geo_Correlations_Facet_Species_Admix1minusmax.pdf", 
       p_facet_2, width = 16, height = 12)
ggsave("output/Individual_Geo_Correlations_Facet_Species_Admix1minusmax.png", 
       p_facet_2, width = 16, height = 12)
