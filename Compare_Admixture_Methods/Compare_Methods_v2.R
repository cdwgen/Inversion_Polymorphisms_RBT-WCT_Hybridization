# ==== Load Libraries ====
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(purrr)
library(stringr)
library(Metrics)
library(patchwork) 

# ==== File paths ====
admix_file <- "ADMIXTURE.tsv"
norm_file  <- "diagnostics_normalized.tsv"
elai_file  <- "ELAI.tsv"

# ==== Read data ====
admix <- read_tsv(admix_file, col_types = cols())
norm  <- read_tsv(norm_file, col_types = cols())
elai  <- read_tsv(elai_file, col_types = cols())

# ==== Dynamic Data Wrangling ====
admix_long <- admix %>%
  rename_with(~ str_remove(., "_admix"), ends_with("_admix")) %>%
  pivot_longer(cols = -Sample_ID, names_to = "Ancestry", values_to = "ADMIXTURE")

norm_long <- norm %>%
  rename_with(~ str_remove(., "norm_p"), starts_with("norm_p")) %>%
  pivot_longer(cols = -Sample_ID, names_to = "Ancestry", values_to = "Diagnostics")

elai_long <- elai %>%
  pivot_longer(cols = -Sample_ID, names_to = "Ancestry", values_to = "ELAI")

# ==== Merge into one long-format table ====
combined <- norm_long %>%
  left_join(admix_long, by = c("Sample_ID", "Ancestry")) %>%
  left_join(elai_long, by = c("Sample_ID", "Ancestry"))

write_tsv(combined, "ancestry_estimate_comparison.tsv")

# ==== Calculate Statistics ====
methods <- c("Diagnostics", "ADMIXTURE", "ELAI")

# Summary Stats (FIXED: Diagnostics is now the baseline)
summary_stats <- combined %>%
  group_by(Ancestry) %>%
  summarise(
    across(all_of(methods),
           list(
             cor = ~ cor(.x, Diagnostics, use = "complete.obs"),
             rmse = ~ rmse(Diagnostics, .x)
           ),
           .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )
write_tsv(summary_stats, "ancestry_summary_statistics.tsv")

# Pairwise Stats Calculation
pairwise_stats <- function(df, ancestry) {
  df_sub <- df %>% filter(Ancestry == ancestry) %>% select(all_of(methods))
  combs <- combn(methods, 2, simplify = FALSE)
  
  map_dfr(combs, ~ {
    tibble(
      Ancestry = ancestry,
      Method1 = .x[1],
      Method2 = .x[2],
      Correlation = cor(df_sub[[.x[1]]], df_sub[[.x[2]]], use = "complete.obs"),
      RMSE = rmse(df_sub[[.x[1]]], df_sub[[.x[2]]])
    )
  })
}

ancestries <- unique(combined$Ancestry)
pairwise_stats_df <- map_dfr(ancestries, ~ pairwise_stats(combined, .x)) %>%
  mutate(MethodPair = paste(Method1, "vs", Method2)) %>%
  mutate(MethodPair = factor(MethodPair))

write_tsv(pairwise_stats_df, "ancestry_pairwise_stats.tsv")

# ==== Visualization Function ====
create_plot <- function(data, y_var, y_label) {
  
  p <- ggplot(data, aes(x = Ancestry, y = .data[[y_var]], fill = MethodPair)) +
    # STYLE UPDATE: 
    # 1. Removed "color = black" (no heavy border)
    # 2. Increased width to 0.8
    # 3. Dodge 0.9 ensures a small gap between bars in a group
    geom_col(position = position_dodge(width = 0.9), width = 0.8) + 
    
    scale_fill_brewer(palette = "Dark2") +
    
    # STYLE UPDATE: switched to 'sans' (Arial)
    theme_classic(base_size = 12, base_family = "sans") + 
    
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      
      # TEXT UPDATES: Pure black, size 11 for readability
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      
      # FIX: Changed face from "bold" to "plain"
      axis.title = element_text(face = "plain", size = 13),
      
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.title = element_blank(), 
      axis.title.x = element_blank() 
    ) +
    labs(y = y_label)
  
  # Specific Y-axis handling
  if(y_var == "Correlation") {
    # Zoom in for correlation
    p <- p + coord_cartesian(ylim = c(0.95, 1.00))
  } else {
    # Anchor to 0 for RMSE (removes floating gap)
    p <- p + scale_y_continuous(expand = c(0, 0))
  }
  
  return(p)
}

# ==== Generate Individual Plots ====
plot_a <- create_plot(pairwise_stats_df, "Correlation", "Pearson Correlation (r)")

plot_b <- create_plot(pairwise_stats_df, "RMSE", "RMSE") +
  theme(axis.title.x = element_text()) + 
  labs(x = "Ancestry Component")

# ==== Combine with Patchwork ====
combined_plot <- (plot_a / plot_b) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'a', tag_suffix = '.') & 
  theme(
    legend.position = 'top',
    # Match tags to the rest of the plot (Sans, Black, Size 12)
    plot.tag = element_text(face = 'bold', size = 12, family = "sans", color = "black"),
    text = element_text(family = "sans", color = "black")
  )

# ==== Save Final Figure ====
ggsave("Figure_Ancestry_Comparison.pdf", combined_plot, width = 7, height = 9, dpi = 300)
ggsave("Figure_Ancestry_Comparison.png", combined_plot, width = 7, height = 9, dpi = 600)

print(combined_plot)
message("Figure saved with consistent formatting (Sans-serif, Plain titles).")