library(tidyverse)
library(readxl)
library(ggpubr)

# ==============================================================================
# ==============================================================================

admixture <- read.csv("centroids/df_overlaps.csv")
names(admixture)[14] <- "overlap" # Geographic Overlap
names(admixture)[c(5,7,9)] <- c("Shannon_H", "Simpson_Index", "Eff_N")

aba_baba <- read_excel("inputs/samples_SDM_INFO.xlsx", sheet = 4)
areas_dist <- read.csv("centroids/areas_dist.csv")

aba_baba_clean <- as.data.frame(aba_baba[, c("P2", "P3", "D statistic")])
aba_baba_clean$P2 <- sub("x_", "", aba_baba_clean$P2) %>% paste("Vitis.", ., sep = "")
aba_baba_clean$P3 <- sub("x_", "", aba_baba_clean$P3) %>% paste("Vitis.", ., sep = "")
aba_baba_clean$dist_cent <- NA

for(i in 1:nrow(aba_baba_clean)){
  pair_key <- paste(sort(c(aba_baba_clean$P2[i], aba_baba_clean$P3[i])), collapse = "_")
  match_idx <- which(apply(areas_dist[, c("sp1", "sp2")], 1, function(x) paste(sort(x), collapse = "_")) == pair_key)
  if(length(match_idx) > 0){
    aba_baba_clean$dist_cent[i] <- areas_dist$dist_cent[match_idx]
  }
}
d_stat_data <- na.omit(aba_baba_clean)

# ==============================================================================
# ==============================================================================
p4 <- plot_corr(d_stat_data, "dist_cent", "D statistic", 
                "Centroid Distance (km)", "D statistic")
data = d_stat_data
x_var = "dist_cent"
y_var = "D statistic"
x_lab = "Centroid Distance (km)"
y_lab = "D statistic"

plot_corr <- function(data, x_var, y_var, x_lab, y_lab) {
  
  formula <- as.formula(paste0("`", y_var, "` ~ `", x_var, "`"))
  fit <- lm(formula, data = data) %>% summary()
  
  p_val <- fit$coefficients[2, 4]
  r_val <- cor(data[[x_var]], data[[y_var]], use = "complete.obs")
  
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001\nR = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nR = ", round(r_val, 2))
  }
  
  p <- ggplot(data, aes_string(x = x_var, y = paste0("`", y_var, "`"))) +
    geom_point(alpha = 0.6, color = "#528fad", size = 2.5) +
    geom_smooth(method = "lm", color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.3, vjust = 1.2, size = 4.5, fontface = "italic") +
    labs(x = x_lab, y = y_lab) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  return(p)
}

# ==============================================================================
# ==============================================================================

p1 <- plot_corr(admixture, "overlap", "Shannon_H", 

p2 <- plot_corr(admixture, "overlap", "Simpson_Index", 

p3 <- plot_corr(admixture, "overlap", "Eff_N", 

p4 <- plot_corr(d_stat_data, "dist_cent", "D statistic", 
                "Centroid Distance (km)", "D statistic")

# ==============================================================================
# ==============================================================================

final_composite <- ggarrange(p4, p1, p2, p3,
                             ncol = 2, nrow = 2 
                             ) # labels = c("A", "B", "C", "D")

print(final_composite)

ggsave("output/Final_Composite_Correlation.pdf", final_composite, width = 10, height = 8)
ggsave("output/Final_Composite_Correlation.png", final_composite, width = 10, height = 8, dpi = 300)

message("Composite plot saved to output/Final_Composite_Correlation.pdf")
