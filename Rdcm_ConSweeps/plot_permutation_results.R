#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(readr)

  make_option(c("-r", "--results_file"), type="character", help="Permutation results file (one value per line)"),
  make_option(c("-o", "--observed_value"), type="numeric", help="Observed overlap value"),
  make_option(c("-p", "--output_plot"), type="character", help="Output plot filename"),
  make_option(c("-t", "--title"), type="character", default="Permutation Test", help="Plot title")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

perm_values <- read_table(opt$results_file, col_names = "overlap", show_col_types = FALSE)
obs_val <- opt$observed_value

n_perms <- nrow(perm_values)
n_greater <- sum(perm_values$overlap >= obs_val)
p_value <- (n_greater + 1) / (n_perms + 1)

p <- ggplot(perm_values, aes(x = overlap)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = obs_val, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = opt$title,
    subtitle = paste0("P-value: ", format(p_value, digits = 4), " (N=", n_perms, ")"),
    x = "Overlapping Base Pairs (bp)",
    y = "Frequency"
  ) +
  theme_bw() +
  annotate("text", x = obs_val, y = 0, label = "Observed", vjust = -1, color = "red", angle = 90, hjust = -0.1)

ggsave(opt$output_plot, plot = p, width = 6, height = 4)
cat(paste0("P-value: ", p_value, "\n"))
