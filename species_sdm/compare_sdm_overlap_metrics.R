library(tidyverse)


overlaps <- read.csv("centroids/areas_dist.csv")
# flextable::flextable(overlaps[c(1:5,15:20,35:40),-1])




load("centroids/cambio.R")
# flextable::flextable(output[c(2:4,14:16,35:40),])


admixture <- read.csv("centroids/df_overlaps.csv")
admixture[sample(1:nrow(admixture),10,F),c("X.1","D","I","sobrelape")]

names(admixture)[14] <- "overlap"

names(admixture)[c(3,5,7,9)]
# "maxprop"     "H_raw"       "simpson_raw" "effN"

names(admixture)[c(3,5,7,9)] <- c( "Max_Prop","Shannon_H","Simpson_Index","Eff_N")

variables <- names(admixture)[c(3,5,7,9)]

pdf("./Comparison_statistics.pdf",width = 21,height = 5)

point_col <- "#528fad" 
line_col <- "#e74c3c"

for(i in variables){
  
  fit <- lm(admixture[,i]~admixture$D) %>% summary() 
  p_val <- fit$coefficients[2,4]
  r_val <- cor(admixture[,i], admixture$D, use = "complete.obs")
  
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001\nr = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nr = ", round(r_val, 2))
  }
  
  a <- ggplot(admixture, mapping = aes_string("D", i)) +
    geom_point(alpha = 0.6, color = point_col, size = 2) +
    geom_smooth(method = "lm", color = line_col, fill = line_col, alpha = 0.2) +
    labs(x = "D (Suitability Overlap)", y = i) + 
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.1, vjust = 1.2, size = 5, fontface = "italic") + 
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  fit <- lm(admixture[,i]~admixture$I) %>% summary() 
  p_val <- fit$coefficients[2,4]
  r_val <- cor(admixture[,i], admixture$I, use = "complete.obs")
  
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001\nr = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nr = ", round(r_val, 2))
  }
  
  b <- ggplot(admixture, mapping = aes_string("I", i)) +
    geom_point(alpha = 0.6, color = point_col, size = 2) +
    geom_smooth(method = "lm", color = line_col, fill = line_col, alpha = 0.2) +
    labs(x = "I (Suitability Overlap)", y = i) + 
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.1, vjust = 1.2, size = 5, fontface = "italic") + 
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  fit <- lm(admixture[,i]~admixture$overlap) %>% summary() 
  p_val <- fit$coefficients[2,4]
  r_val <- cor(admixture[,i], admixture$overlap, use = "complete.obs")
  
  if(p_val < 0.001){
    stats_text <- paste0("p < 0.001\nr = ", round(r_val, 2))
  } else {
    stats_text <- paste0("p = ", round(p_val, 3), "\nr = ", round(r_val, 2))
  }
  
  c <- ggplot(admixture, mapping = aes_string("overlap", i)) +
    geom_point(alpha = 0.6, color = point_col, size = 2) +
    geom_smooth(method = "lm", color = line_col, fill = line_col, alpha = 0.2) +
    labs(x = "Geographic Overlap", y = i) + 
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.1, vjust = 1.2, size = 5, fontface = "italic") + 
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  print(ggpubr::ggarrange(a, b, c, ncol = 3))
}

dev.off()
