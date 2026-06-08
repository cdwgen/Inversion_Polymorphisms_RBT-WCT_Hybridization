library(tidyverse)
library(patchwork)

# 1. Load your data
# Note: Ensure the file name matches your actual inversion counts file
df <- read_tsv("inversion_counts.txt", show_col_types = FALSE) %>%
  rename(Het_Inversions = `# Heterokaryotypic Inversions`)

# 2. LOAD AND APPLY WHITELIST
# Updated to match the filename you provided
whitelist <- read_lines("whitelist_triangle_plot.txt")

df <- df %>%
  # Updated "Sample" to "Sample_ID" to match your file header
  filter(Sample_ID %in% whitelist)

# 3. Calculate "Minor Parent Ancestry"
df <- df %>%
  mutate(
    Minor_Ancestry = pmin(RBT_ELAI, WCT_ELAI),
    # Create bins for the boxplot panel
    Ancestry_Bin = cut(Minor_Ancestry, 
                       breaks = c(-Inf, 0.05, 0.15, 0.35, Inf), 
                       labels = c("<5%", "5-15%", "15-35%", ">35%"))
  )

# === THEME ===
theme_molecol <- function() {
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = "topleft"
    )
}

# === PANEL A: Continuous Jitter Plot + Trendline ===
pA <- ggplot(df, aes(x = Minor_Ancestry, y = Het_Inversions)) +
  geom_jitter(width = 0.01, height = 0.15, alpha = 0.4, color = "#4575b4", size = 2) +
  geom_smooth(method = "loess", color = "#d73027", fill = "gray80", linewidth = 1) +
  scale_y_continuous(breaks = 0:6, limits = c(-0.2, 6.2)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Minor Parent Ancestry",
    y = "Number of Heterokaryotypic Inversions",
    tag = "a."
  ) +
  theme_molecol()

# === PANEL B: Boxplots by Generation/Bin ===
pB <- ggplot(df, aes(x = Ancestry_Bin, y = Het_Inversions)) +
  geom_boxplot(fill = "gray90", color = "black", outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, height = 0.15, alpha = 0.3, color = "#4575b4", size = 1.5) +
  scale_y_continuous(breaks = 0:6, limits = c(-0.2, 6.2)) +
  labs(
    x = "Admixture Proportion",
    y = "Number of Heterokaryotypic Inversions",
    tag = "b."
  ) +
  theme_molecol() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# === COMBINE AND SAVE ===
combined_plot <- pA + pB + plot_layout(widths = c(1, 1))

ggsave("Figure_S10_Inversion_Decay.pdf", combined_plot, width = 10, height = 5, dpi = 300)
ggsave("Figure_S10_Inversion_Decay.png", combined_plot, width = 10, height = 5, dpi = 300)

print(combined_plot)