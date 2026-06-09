#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

# ==== 1. LOAD & FILTER GLOBAL ANCESTRY DATA ====
ancestry <- read_tsv("ancestry_table.tsv")
whitelist <- read_lines("whitelist_v2.txt")

# Filter to whitelist AND exclude YCT ancestry
yct_threshold <- 0.001 

ancestry_filtered <- ancestry %>% 
  filter(Sample_ID %in% whitelist) %>%
  filter(pYCT <= yct_threshold) 

# ==== 2. LOAD & MERGE RAW DIAGNOSTIC DATA ====
diagnostic_data <- read_tsv("counts_output/RBT_diagnostic_metrics.tsv")

# Define ground-truth F1s based on inversion data
true_f1s <- c("GOC1_03_18", "LSM_14_11", "ParmenterCr_20_007", 
              "W_14_015", "KRBY_22_06", "KRBY_22_02", 
              "W_14_008", "SML_20_157")

plot_data <- left_join(ancestry_filtered, diagnostic_data, by = "Sample_ID") %>%
  filter(!is.na(Heterozygosity)) %>%
  # Categorization logic prioritizing biological ground-truth
  mutate(
    Hybrid_Class = case_when(
      # Biological Override: Force known F1s first
      Sample_ID %in% true_f1s ~ "F1",
      
      # Unadmixed: Strict 99% threshold
      Raw_Ancestry <= 0.01 ~ "Unadmixed WCT",
      Raw_Ancestry >= 0.99 ~ "Unadmixed RBT",
      
      # F1s: Centered hybrid index, high heterozygosity for any remaining un-flagged F1s
      Raw_Ancestry >= 0.40 & Raw_Ancestry <= 0.60 & Heterozygosity > 0.80 ~ "F1",
      
      # F2s: Centered hybrid index, intermediate heterozygosity
      Raw_Ancestry >= 0.40 & Raw_Ancestry <= 0.60 & Heterozygosity >= 0.40 & Heterozygosity <= 0.60 ~ "F2",
      
      # Backcrosses to WCT (Bounded along the left diagonal with +/- 0.1 tolerance)
      (Heterozygosity - 2 * Raw_Ancestry) >= -0.1 & (Heterozygosity - 2 * Raw_Ancestry) <= 0.1 ~ "Backcross (to WCT)",
      
      # Backcrosses to RBT (Bounded along the right diagonal with +/- 0.1 tolerance)
      (Heterozygosity - 2 * (1 - Raw_Ancestry)) >= -0.1 & (Heterozygosity - 2 * (1 - Raw_Ancestry)) <= 0.1 ~ "Backcross (to RBT)",
      
      # Anything that escapes these strict bounds is a later-generation introgressant
      TRUE ~ "Complex Admixture"
    ),
    # UPDATED: Levels now match the order in inversion/fusion plots
    Hybrid_Class = factor(Hybrid_Class, levels = c(
      "Unadmixed RBT",
      "Unadmixed WCT",
      "Unadmixed YCT",
      "F1",
      "F2",
      "Complex Admixture",
      "Backcross (to RBT)",
      "Backcross (to WCT)",
      "Backcross (to YCT)"
    ))
  )

# ==== 3. UNIFIED THEME & COLOR PALETTE ====
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right"
    )
}

# Master Manuscript Palette (Okabe-Ito) - Ordered to match factor levels
class_colors <- c(
  "Unadmixed RBT"      = "#D55E00", # Vermilion
  "Unadmixed WCT"      = "#0072B2", # Blue
  "Unadmixed YCT"      = "#009E73", # Teal
  "F1"                 = "#F0E442", # Yellow
  "F2"                 = "#CC79A7", # Pink/Rose
  "Complex Admixture"  = "#999999", # Gray
  "Backcross (to RBT)" = "#E69F00", # Gold
  "Backcross (to WCT)" = "#56B4E9", # Sky Blue
  "Backcross (to YCT)" = "#000000"  # Black
)

# ==== 4. DEFINE THE TRIANGLE BOUNDARIES ====
triangle_bounds <- data.frame(
  x = c(0, 0.5, 1, 0),
  y = c(0, 1, 0, 0)
)

# ==== 5. BUILD THE STATIC PLOT ====
p_static <- ggplot() +
  geom_path(data = triangle_bounds, aes(x = x, y = y), color = "black", linewidth = 1) +
  geom_point(data = plot_data, 
             aes(x = Raw_Ancestry, y = Heterozygosity, fill = Hybrid_Class), 
             shape = 21, color = "black", size = 3.5, stroke = 0.5, alpha = 0.85) +
  # drop = TRUE here will hide the YCT classes from the legend since they are filtered out
  scale_fill_manual(values = class_colors, name = "Inferred Hybrid Class", drop = TRUE) +
  labs(
    x = "Hybrid Index (Proportion RBT Ancestry)",
    y = "Interspecific Heterozygosity"
  ) +
  theme_pub()

ggsave("hybrid_triangle_categorized.png", plot = p_static, width = 8, height = 6, dpi = 300)
ggsave("hybrid_triangle_categorized.pdf", plot = p_static, width = 8, height = 6, device = "pdf", dpi = 300)

# ==== 6. BUILD THE INTERACTIVE PLOT ====
p_interactive <- ggplot() +
  geom_path(data = triangle_bounds, aes(x = x, y = y), color = "black", linewidth = 1) +
  geom_point(data = plot_data, 
             aes(x = Raw_Ancestry, y = Heterozygosity, text = Sample_ID, color = Hybrid_Class), 
             size = 3, alpha = 0.8) +
  scale_color_manual(values = class_colors, name = "Inferred Hybrid Class", drop = TRUE) +
  labs(
    x = "Hybrid Index (Proportion RBT Ancestry)",
    y = "Interspecific Heterozygosity"
  ) +
  theme_pub()

interactive_triangle <- ggplotly(p_interactive, tooltip = c("text", "x", "y"))
saveWidget(interactive_triangle, "hybrid_triangle_plot_Final.html")

# ==== 7. EXPORT HYBRID CLASSES ====
plot_data %>%
  select(Sample_ID, Hybrid_Class) %>%
  write_tsv("triangle_hybrid_classifications.tsv")
